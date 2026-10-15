"""Vet the BLS signals of one TESS sector with LEO-Vetter (Kunimoto et al. 2025).

The script reads the signals that search_qlp_sector.py found, and keeps those with a BLS SNR of at least --min-snr
and at least --min-transits transits. For each kept signal it:
1. loads the light curve and detrends it as the search did;
2. runs the LEO-Vetter flux-level tests;
3. checks the results against the LEO-Vetter default thresholds, in two groups:
   - false alarm (FA) tests: is the signal a real astrophysical event, not noise or a systematic?
   - false positive (FP) tests: is the event a transiting planet, not an eclipsing binary?

A signal that passes both groups is a planet candidate.
LEO-Vetter also has pixel-level tests; this script does not run them.

Run from the repository root, for example for Sector 98:
    .venv/bin/python analysis/vet_qlp_sector_signals.py 98

--workers sets the number of processes. The script writes each row as soon as it has it. It keeps the rows of an
earlier run without an error, and vets only the other signals. --restart vets all signals again.

Output: results/qlp_sector_<nnnn>_leo_vetting/vetting.csv, one row per vetted signal, with the search values,
the main LEO-Vetter metrics, passes_false_alarm_tests, passes_false_positive_tests and the names of the failed tests.
"""

import argparse
import csv
import os
import warnings
from multiprocessing import Pool

import numpy as np
from leo_vetter import thresholds
from leo_vetter.main import TCELightCurve
from leo_vetter.stellar import quadratic_ldc

from build_qlp_sector_star_list import REPOSITORY_ROOT, sector_directory
from download_qlp_sector_light_curves import light_curve_path
from incremental_csv_output import IncrementalCsvWriter, earlier_rows_to_keep
from pipeline_status import PipelineStatus
from search_qlp_sector import (
    detrend, load_light_curve, output_directory as search_output_directory, robust_standard_deviation,
)

SOLAR_DENSITY_GRAMS_PER_CUBIC_CENTIMETRE = 1.41
# The star list has no parameter errors. These fractional errors stand in for them.
FRACTIONAL_RADIUS_ERROR = 0.1
FRACTIONAL_MASS_ERROR = 0.1
TEFF_ERROR_KELVIN = 150.0
# Used when the TIC has no value.
DEFAULT_STELLAR_PARAMETERS = {"radius": 1.0, "mass": 1.0, "teff": 5800.0, "logg": 4.44}

FALSE_ALARM_TESTS = [
    thresholds.weak, thresholds.invalid_transits, thresholds.bad_shape, thresholds.non_unique,
    thresholds.chases, thresholds.dmm, thresholds.single_event, thresholds.bad_fit, thresholds.sinusoidal,
    thresholds.unphysical_duration, thresholds.asymmetric, thresholds.chi, thresholds.data_gapped,
]
FALSE_POSITIVE_TESTS = [thresholds.odd_even, thresholds.vshaped, thresholds.large, thresholds.secondary]
# transit_aRs is a/R* from the transit fit. aRs is a/R* from the period and the TIC stellar mass and radius.
REPORTED_METRICS = [
    "MES", "SHP", "CHI", "Rp", "transit_b", "transit_RpRs", "transit_aRs", "aRs", "trap_dep", "sec_dep", "pri_dep",
]

OUTPUT_COLUMNS = (
    ["tic_id", "signal_number", "period_days", "epoch_btjd", "duration_hours", "depth_ppm", "snr", "sde",
     "transit_count", "tmag", "used_default_stellar_parameters"]
    + REPORTED_METRICS
    + ["passes_false_alarm_tests", "passes_false_positive_tests", "failed_tests", "error"]
)


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_leo_vetting"


def star_parameters(star_row):
    """Return the LEO-Vetter star dict, and whether any TIC value was missing."""
    values, used_default = {}, False
    for name, default in DEFAULT_STELLAR_PARAMETERS.items():
        text = star_row.get(name, "")
        if text in ("", "nan") or not np.isfinite(float(text)):
            values[name], used_default = default, True
        else:
            values[name] = float(text)
    u1, u2 = quadratic_ldc(values["teff"], values["logg"])
    return {
        "rad": values["radius"],
        "mass": values["mass"],
        "Teff": values["teff"],
        "e_rad": FRACTIONAL_RADIUS_ERROR * values["radius"],
        "e_mass": FRACTIONAL_MASS_ERROR * values["mass"],
        "e_Teff": TEFF_ERROR_KELVIN,
        "rho": SOLAR_DENSITY_GRAMS_PER_CUBIC_CENTIMETRE * values["mass"] / values["radius"] ** 3,
        "u1": float(u1),
        "u2": float(u2),
    }, used_default


def failed_test_names(metrics, tests):
    return [test.__name__ for test in tests if test(metrics, thresholds._default_thresholds)[0]]


def vet_signal(arguments):
    sector, signal, star_row, path = arguments
    row = {column: signal.get(column, "") for column in OUTPUT_COLUMNS}
    row["tmag"] = star_row.get("tmag", "")
    try:
        time, raw_flux = load_light_curve(path)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            detrended_time, detrended_flux = detrend(time, raw_flux)
            raw_flux = np.interp(detrended_time, time, raw_flux)
            flux_error = np.full_like(detrended_flux, robust_standard_deviation(detrended_flux))
            star, row["used_default_stellar_parameters"] = star_parameters(star_row)
            light_curve = TCELightCurve(
                signal["tic_id"], detrended_time, raw_flux, detrended_flux, flux_error,
                float(signal["period_days"]), float(signal["epoch_btjd"]), float(signal["duration_hours"]) / 24,
                planetno=int(signal["signal_number"]),
            )
            light_curve.compute_flux_metrics(star, verbose=False)
        metrics = light_curve.metrics
        for name in REPORTED_METRICS:
            value = metrics.get(name, "")
            row[name] = round(float(value), 5) if isinstance(value, (int, float, np.floating)) else value
        failed_false_alarm = failed_test_names(metrics, FALSE_ALARM_TESTS)
        failed_false_positive = failed_test_names(metrics, FALSE_POSITIVE_TESTS)
        row["passes_false_alarm_tests"] = not failed_false_alarm
        row["passes_false_positive_tests"] = not failed_false_positive
        row["failed_tests"] = " ".join(failed_false_alarm + failed_false_positive)
    except Exception as error:
        row["error"] = f"{type(error).__name__}: {error}"
    return row


def selected_signals(sector, minimum_snr, minimum_transits):
    for chunk_path in sorted(search_output_directory(sector).glob("chunk_*[0-9].csv")):
        for signal in csv.DictReader(open(chunk_path)):
            if signal["signal_number"] == "0":
                continue
            if float(signal["snr"]) >= minimum_snr and int(signal["transit_count"]) >= minimum_transits:
                yield signal


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    parser.add_argument("--min-snr", type=float, default=7.0)
    parser.add_argument("--min-transits", type=int, default=3)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 2))
    parser.add_argument("--restart", action="store_true",
                        help="vet all signals again, not only those without a row in vetting.csv")
    arguments = parser.parse_args()
    sector = arguments.sector

    output_directory(sector).mkdir(parents=True, exist_ok=True)
    output_path = output_directory(sector) / "vetting.csv"
    star_rows = {int(row["tic_id"]): row for row in csv.DictReader(open(sector_directory(sector) / "star_list.csv"))}
    signals = list(selected_signals(sector, arguments.min_snr, arguments.min_transits))
    previous_rows = earlier_rows_to_keep(output_path, arguments.restart,
                                         {(signal["tic_id"], signal["signal_number"]) for signal in signals},
                                         key_columns=("tic_id", "signal_number"))
    vetted = {(row["tic_id"], row["signal_number"]) for row in previous_rows}
    tasks = [
        (sector, signal, star_rows[int(signal["tic_id"])], light_curve_path(sector, int(signal["tic_id"])))
        for signal in signals if (signal["tic_id"], signal["signal_number"]) not in vetted
    ]
    print(f"Sector {sector}: {len(signals)} signals with SNR >= {arguments.min_snr} "
          f"and >= {arguments.min_transits} transits, {len(previous_rows)} vetted in an earlier run", flush=True)

    status = PipelineStatus(sector, 4, "leo_vetting", "vet_qlp_sector_signals.py", len(signals), "signals")
    with Pool(arguments.workers) as pool, IncrementalCsvWriter(output_path, OUTPUT_COLUMNS, previous_rows) as writer, \
            status:
        status.set_items_done(len(previous_rows))
        for count, row in enumerate(pool.imap_unordered(vet_signal, tasks, chunksize=2), start=len(previous_rows) + 1):
            writer.write_row(row)
            status.set_items_done(count)
            if count % 500 == 0:
                print(f"  {count} of {len(signals)} signals vetted", flush=True)
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
