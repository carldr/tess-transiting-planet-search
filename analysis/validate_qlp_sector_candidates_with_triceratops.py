"""Compute the false positive probability of each candidate of one sector that the pixel check puts on the target.

The script runs TRICERATOPS (Giacalone et al. 2021) on each candidate with pixel verdict "on target" or
"unresolved". TRICERATOPS fits the folded transit with each false positive scenario: an eclipsing binary on the
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

TRICERATOPS needs the photometric aperture. The QLP files do not give it, so the script uses the pixels within
APERTURE_RADIUS_PIXELS of the target. The QLP SAP flux is not corrected for the light of other stars.
TRICERATOPS draws random samples, so the script runs it RUNS_PER_CANDIDATE times and gives the mean and the
standard deviation.

Run from the repository root, with the TRICERATOPS environment, after check_qlp_sector_candidate_pixels.py:
    .venv_triceratops/bin/python analysis/validate_qlp_sector_candidates_with_triceratops.py 98

Output: results/qlp_sector_<nnnn>_triceratops_validation/triceratops_validation.csv, one row per candidate.
"""

import argparse
import csv
import warnings

import numpy as np

from build_qlp_sector_star_list import REPOSITORY_ROOT
from check_qlp_sector_candidate_pixels import (
    PIXEL_SIZE_ARCSEC, UNRESOLVABLE_SEPARATION_PIXELS, output_directory as pixel_checks_output_directory,
)
from download_qlp_sector_light_curves import light_curve_path
from search_qlp_sector import detrend, load_light_curve, robust_standard_deviation

PIXEL_VERDICTS_TO_VALIDATE = ("on target", "unresolved")
APERTURE_RADIUS_PIXELS = 2.5
# The folded light curve keeps the cadences within this many transit durations of mid-transit.
FOLDED_WINDOW_DURATIONS = 3.0
QLP_EXPOSURE_DAYS = 200 / 86400
SIMULATIONS_PER_SCENARIO = 1_000_000
RUNS_PER_CANDIDATE = 5
LIKELY_PLANET_FPP = 0.5
VALIDATED_FPP = 0.015
NEARBY_FPP_LIMIT = 0.001

OUTPUT_COLUMNS = [
    "tic_id", "period_days", "depth_ppm", "tmag", "check_verdict", "pixel_verdict", "nearby_scenarios_excluded",
    "fpp_mean", "fpp_standard_deviation", "nfpp_mean", "nfpp_standard_deviation", "most_likely_scenario",
    "triceratops_verdict", "error",
]


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_triceratops_validation"


def folded_light_curve(sector, candidate):
    """Return the time from mid-transit, the detrended flux, and the flux error, near the transits."""
    time, flux = load_light_curve(light_curve_path(sector, int(candidate["tic_id"])))
    time, flux = detrend(time, flux)
    period, epoch = float(candidate["period_days"]), float(candidate["epoch_btjd"])
    duration = float(candidate["duration_hours"]) / 24
    phase_time = (time - epoch + 0.5 * period) % period - 0.5 * period
    near = np.abs(phase_time) < FOLDED_WINDOW_DURATIONS * duration
    out_of_transit = ~(np.abs(phase_time) < duration)
    return phase_time[near], flux[near], robust_standard_deviation(flux[out_of_transit])


def circular_aperture(target_pixel):
    centre = np.round(target_pixel)
    offsets = np.arange(-np.ceil(APERTURE_RADIUS_PIXELS), np.ceil(APERTURE_RADIUS_PIXELS) + 1)
    return np.array([[centre[0] + column, centre[1] + row] for column in offsets for row in offsets
                     if np.hypot(column, row) <= APERTURE_RADIUS_PIXELS])


def validate(sector, candidate):
    import triceratops.triceratops as triceratops

    target = triceratops.target(ID=int(candidate["tic_id"]), sectors=np.array([sector]))
    target.calc_depths(float(candidate["depth_ppm"]) * 1e-6,
                       all_ap_pixels=[circular_aperture(target.pix_coords[0][0])], dilution_corrected=False)
    # The pixel check rules out the stars more than UNRESOLVABLE_SEPARATION_PIXELS from an "on target" dip.
    # TRICERATOPS skips the nearby-star scenarios of a star with a tdepth of 0.
    nearby_scenarios_excluded = candidate["pixel_verdict"] == "on target"
    if nearby_scenarios_excluded:
        ruled_out = target.stars["sep (arcsec)"] > UNRESOLVABLE_SEPARATION_PIXELS * PIXEL_SIZE_ARCSEC
        target.stars.loc[ruled_out, "tdepth"] = 0
    time, flux, flux_error = folded_light_curve(sector, candidate)
    fpps, nfpps, scenarios = [], [], []
    for _ in range(RUNS_PER_CANDIDATE):
        target.calc_probs(time=time, flux_0=flux, flux_err_0=flux_error, P_orb=float(candidate["period_days"]),
                          N=SIMULATIONS_PER_SCENARIO, exptime=QLP_EXPOSURE_DAYS, dilution_corrected=False,
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
    sector = parser.parse_args().sector

    pixel_checks = list(csv.DictReader(open(pixel_checks_output_directory(sector) / "pixel_checks.csv")))
    candidates = [row for row in pixel_checks if row["pixel_verdict"] in PIXEL_VERDICTS_TO_VALIDATE]
    print(f"Sector {sector}: {len(candidates)} candidates on the target or unresolved", flush=True)

    rows = []
    for count, candidate in enumerate(candidates, start=1):
        row = {column: candidate.get(column, "") for column in OUTPUT_COLUMNS}
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                row.update(validate(sector, candidate))
        except Exception as error:
            row["error"] = f"{type(error).__name__}: {error}"
        rows.append(row)
        print(f"  {count} of {len(candidates)}: TIC {row['tic_id']} FPP {row['fpp_mean']} NFPP {row['nfpp_mean']} "
              f"{row['triceratops_verdict'] or row['error']}", flush=True)

    directory = output_directory(sector)
    directory.mkdir(parents=True, exist_ok=True)
    output_path = directory / "triceratops_validation.csv"
    with open(output_path, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
