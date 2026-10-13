"""Compute the false positive probability of each candidate of one sector that passed the all-sector verification.

The script runs TRICERATOPS (Giacalone et al. 2021) on each candidate with the verdict "passed" in
all_sector_verification.csv. It uses the light curves of all the sectors of that verification, folded on the fitted
period and epoch, and the fitted depth. TRICERATOPS takes one exposure time, so the script bins each sector to the
longest exposure time of the sectors. TRICERATOPS fits the folded transit with each false positive scenario: an eclipsing binary on the
target, on an unresolved companion, on a nearby star or on a background star. It weighs the fits by the stars in
the TIC and Gaia, and gives:
- FPP, the false positive probability: the probability that the transit is not a planet on the target;
- NFPP, the nearby false positive probability: the probability that the transit is on a nearby star.

Giacalone et al. (2021) call a candidate with FPP < 0.015 and NFPP < 0.001 validated, and a candidate with
FPP < 0.5 and NFPP < 0.001 a likely planet. Validation of a TESS candidate also needs high-resolution imaging, which
this script does not have, so it reports the second test only.

For a candidate with pixel verdict "on target", the difference image rules out the stars more than
UNRESOLVABLE_SEPARATION_PIXELS from the target, so the script removes their nearby-star scenarios. Giacalone et al.
(2021) do the same for stars that ground-based photometry rules out.

TRICERATOPS needs the photometric aperture of each sector. The QLP files do not give it, so the script uses the pixels
within APERTURE_RADIUS_PIXELS of the target in each sector. TRICERATOPS takes one dilution setting for all sectors.
The QLP SAP flux is not corrected for the light of other stars, so the script tells TRICERATOPS that the depth is
not corrected.
TRICERATOPS draws random samples, so the script runs it RUNS_PER_CANDIDATE times and gives the mean and the
standard deviation.

Run from the repository root, with the TRICERATOPS environment, after verify_qlp_sector_candidates_over_all_sectors.py:
    .venv_triceratops/bin/python analysis/validate_qlp_sector_candidates_with_triceratops.py 98

Output: results/qlp_sector_<nnnn>_triceratops_validation/triceratops_validation.csv, one row per candidate.
The script writes each row as soon as it has it. It keeps the rows of an earlier run without an error, and validates
only the other candidates. --restart validates all candidates again.
The script also appends the result of each candidate to its verification record.
"""

import argparse
import csv
import warnings

import numpy as np

from build_qlp_sector_star_list import REPOSITORY_ROOT
from candidate_verification_records import append_verification_record_entry
from check_qlp_sector_candidate_pixels import PIXEL_SIZE_ARCSEC, UNRESOLVABLE_SEPARATION_PIXELS
from incremental_csv_output import IncrementalCsvWriter, earlier_rows_to_keep
from search_qlp_sector import robust_standard_deviation
from verify_qlp_sector_candidates_over_all_sectors import (
    load_all_sectors, output_directory as all_sector_verification_output_directory, phase_time,
)

APERTURE_RADIUS_PIXELS = 2.5
# The folded light curve keeps the cadences within this many transit durations of mid-transit.
FOLDED_WINDOW_DURATIONS = 3.0
SIMULATIONS_PER_SCENARIO = 1_000_000
RUNS_PER_CANDIDATE = 5
LIKELY_PLANET_FPP = 0.5
VALIDATED_FPP = 0.015
NEARBY_FPP_LIMIT = 0.001

OUTPUT_COLUMNS = [
    "tic_id", "period_days", "tmag", "pixel_verdict", "sectors_used", "fit_period_days", "fit_epoch_btjd",
    "fit_depth_ppm", "nearby_scenarios_excluded",
    "fpp_mean", "fpp_standard_deviation", "nfpp_mean", "nfpp_standard_deviation", "most_likely_scenario",
    "triceratops_verdict", "error",
]


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_triceratops_validation"


def folded_light_curve(sector, candidate):
    """Return the time from mid-transit, the detrended flux and the flux error of all sectors, and the exposure time.

    TRICERATOPS takes one exposure time, so the script bins each sector to the longest exposure time of the sectors.
    """
    period, epoch = float(candidate["fit_period_days"]), float(candidate["fit_epoch_btjd"])
    duration = float(candidate["fit_duration_hours"]) / 24
    sectors = load_all_sectors(sector, candidate, period, epoch, duration)
    exposure = max(sector_exposure for _, _, sector_exposure, _, _ in sectors)
    times, fluxes = [], []
    for _, _, _, time, flux in sectors:
        bin_numbers = np.floor((time - time.min()) / (exposure / 86400)).astype(int)
        counts = np.bincount(bin_numbers)
        filled = counts > 0
        times.append(np.bincount(bin_numbers, time)[filled] / counts[filled])
        fluxes.append(np.bincount(bin_numbers, flux)[filled] / counts[filled])
    time, flux = np.concatenate(times), np.concatenate(fluxes)
    time_from_mid_transit = phase_time(time, period, epoch)
    near = np.abs(time_from_mid_transit) < FOLDED_WINDOW_DURATIONS * duration
    out_of_transit = np.abs(time_from_mid_transit) > duration
    return (time_from_mid_transit[near], flux[near], robust_standard_deviation(flux[out_of_transit]),
            exposure / 86400)


def circular_aperture(target_pixel):
    centre = np.round(target_pixel)
    offsets = np.arange(-np.ceil(APERTURE_RADIUS_PIXELS), np.ceil(APERTURE_RADIUS_PIXELS) + 1)
    return np.array([[centre[0] + column, centre[1] + row] for column in offsets for row in offsets
                     if np.hypot(column, row) <= APERTURE_RADIUS_PIXELS])


def validate(sector, candidate):
    import triceratops.triceratops as triceratops

    sectors = np.array([int(value) for value in candidate["sectors_used"].split()])
    target = triceratops.target(ID=int(candidate["tic_id"]), sectors=sectors)
    target.calc_depths(float(candidate["fit_depth_ppm"]) * 1e-6,
                       all_ap_pixels=[circular_aperture(coordinates[0]) for coordinates in target.pix_coords],
                       dilution_corrected=False)
    # The pixel check rules out the stars more than UNRESOLVABLE_SEPARATION_PIXELS from an "on target" dip.
    # TRICERATOPS skips the nearby-star scenarios of a star with a tdepth of 0.
    nearby_scenarios_excluded = candidate["pixel_verdict"] == "on target"
    if nearby_scenarios_excluded:
        ruled_out = target.stars["sep (arcsec)"] > UNRESOLVABLE_SEPARATION_PIXELS * PIXEL_SIZE_ARCSEC
        target.stars.loc[ruled_out, "tdepth"] = 0
    time, flux, flux_error, exposure_days = folded_light_curve(sector, candidate)
    fpps, nfpps, scenarios = [], [], []
    for _ in range(RUNS_PER_CANDIDATE):
        target.calc_probs(time=time, flux_0=flux, flux_err_0=flux_error, P_orb=float(candidate["fit_period_days"]),
                          N=SIMULATIONS_PER_SCENARIO, exptime=exposure_days, dilution_corrected=False,
                          verbose=0)
        fpps.append(target.FPP)
        nfpps.append(target.NFPP)
        probabilities = target.probs
        best = probabilities.loc[probabilities["prob"].idxmax()]
        scenarios.append(f"{best['scenario']} on TIC {best['ID']}")
    fpp, nfpp = float(np.mean(fpps)), float(np.mean(nfpps))
    if fpp < VALIDATED_FPP and nfpp < NEARBY_FPP_LIMIT:
        verdict = "validated without imaging"
    elif fpp < LIKELY_PLANET_FPP and nfpp < NEARBY_FPP_LIMIT:
        verdict = "likely planet"
    else:
        verdict = "not validated"
    return {
        "nearby_scenarios_excluded": nearby_scenarios_excluded,
        "fpp_mean": f"{fpp:.4g}", "fpp_standard_deviation": f"{np.std(fpps):.2g}",
        "nfpp_mean": f"{nfpp:.4g}", "nfpp_standard_deviation": f"{np.std(nfpps):.2g}",
        "most_likely_scenario": max(set(scenarios), key=scenarios.count),
        "triceratops_verdict": verdict,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    parser.add_argument("--restart", action="store_true",
                        help="validate all candidates again, not only those without a row in triceratops_validation.csv")
    arguments = parser.parse_args()
    sector = arguments.sector

    directory = output_directory(sector)
    directory.mkdir(parents=True, exist_ok=True)
    output_path = directory / "triceratops_validation.csv"
    earlier_rows = earlier_rows_to_keep(output_path, arguments.restart)
    validated = {(row["tic_id"], row["period_days"]) for row in earlier_rows}
    verifications = list(csv.DictReader(
        open(all_sector_verification_output_directory(sector) / "all_sector_verification.csv")))
    passed = [row for row in verifications if row["verdict"] == "passed"]
    candidates = [row for row in passed if (row["tic_id"], row["period_days"]) not in validated]
    print(f"Sector {sector}: {len(passed)} candidates passed the all-sector verification, "
          f"{len(earlier_rows)} validated in an earlier run", flush=True)

    with IncrementalCsvWriter(output_path, OUTPUT_COLUMNS, earlier_rows) as writer:
        for count, candidate in enumerate(candidates, start=1):
            row = {column: candidate.get(column, "") for column in OUTPUT_COLUMNS}
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    row.update(validate(sector, candidate))
            except Exception as error:
                row["error"] = f"{type(error).__name__}: {error}"
            result = row["error"] or (
                f"sectors {candidate['sectors_used']}; FPP {row['fpp_mean']} +/- {row['fpp_standard_deviation']}; "
                f"NFPP {row['nfpp_mean']} +/- {row['nfpp_standard_deviation']}; most likely scenario "
                f"{row['most_likely_scenario']}; nearby scenarios excluded: {row['nearby_scenarios_excluded']}; "
                f"{RUNS_PER_CANDIDATE} runs"
            )
            append_verification_record_entry(sector, candidate["tic_id"], candidate["period_days"], "TRICERATOPS",
                                             row["triceratops_verdict"] or "error", result)
            writer.write_row(row)
            print(f"  {count} of {len(candidates)}: TIC {row['tic_id']} {row['triceratops_verdict'] or 'error'}: "
                  f"{result}", flush=True)
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
