"""Find which star causes the dip of each surviving candidate of one sector, from TESS full-frame image pixels.

TESS pixels are 21 arcsec wide, so the light of nearby stars falls on the target pixels. A dip on a nearby star
can look like a transit on the target. For each candidate that check_qlp_sector_new_candidates.py did not refute,
the script checks each sector in which the star falls on the science pixels and TESScut has full-frame images.
In each sector it:
1. downloads a 21x21-pixel TESScut cutout of the full-frame images;
2. makes a difference image: for each transit, the mean image just before and after the transit minus the mean
   image in the transit, averaged over the transits. The pixels that lose light in the transit are bright in it;
3. gets the Gaia DR3 sources within SOURCE_SEARCH_RADIUS_PIXELS pixels of the target that are bright enough to cause
   the dip, and computes their pixel positions from the cutout WCS;
4. fits the difference image with the TESS pixel response function (PRF) at the position of each of those stars,
   and also at a free position;
5. compares the chi-square of the fits. The star with the lowest chi-square is the most likely source.
   The chi-square values are scaled so that the best fit has a reduced chi-square of at least 1.

The fit cannot separate the target from a star within UNRESOLVABLE_SEPARATION_PIXELS, so "target" below includes
those stars. The pixel verdict is:
- "on target": the target fits better than every other star by DELTA_CHI_SQUARE_DECISION or more;
- "on neighbour": a nearby star fits better than the target by DELTA_CHI_SQUARE_DECISION or more;
- "unresolved": the target and a nearby star fit about equally well;
- "too faint": the difference image has an SNR below MINIMUM_DIFFERENCE_IMAGE_SNR at the best star.

The verdict of the candidate combines the verdicts of the sectors. It is "on neighbour" when any sector puts the dip
on a neighbour. Otherwise it is the first of "on target", "unresolved" and "too faint" that a sector gives.

Run from the repository root, after check_qlp_sector_new_candidates.py:
    .venv/bin/python analysis/check_qlp_sector_candidate_pixels.py 98

The script saves the cutouts in data/candidate_tesscut_cutouts/<TIC ID>/.
Output, in results/qlp_sector_<nnnn>_candidate_pixel_checks/:
- pixel_checks.csv: one row per candidate, with the combined verdict;
- pixel_checks_by_sector.csv: one row per candidate and sector;
- tic<TIC ID>_<period>_s<sector>.png: the out-of-transit image and the difference image, with the stars marked.
The script writes the rows of each candidate as soon as it has them. It keeps the candidates of an earlier run with a
pixel verdict, and checks only the other candidates. --restart checks all candidates again.
The script also appends the result of each candidate to its verification record.
"""

import argparse
import csv
import warnings

import astropy.units as units
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.gaia import Gaia
from astroquery.mast import Tesscut
from scipy.optimize import minimize
from tess_stars2px import tess_stars2px_function_entry
from transitDiffImage.tessprfmodel import SimpleTessPRF

from build_qlp_sector_star_list import REPOSITORY_ROOT, sector_directory
from candidate_verification_records import append_verification_record_entry, verification_stamp
from check_qlp_sector_new_candidates import output_directory as candidate_checks_output_directory
from incremental_csv_output import IncrementalCsvWriter, complete_csv_rows, earlier_rows_to_keep

CUTOUT_DIRECTORY = REPOSITORY_ROOT / "data/candidate_tesscut_cutouts"
CUTOUT_SIZE_PIXELS = 21
PIXEL_SIZE_ARCSEC = 21.0
# The fit uses the pixels within this many pixels of the target, and the stars inside them.
FIT_HALF_WIDTH_PIXELS = 5
SOURCE_SEARCH_RADIUS_PIXELS = FIT_HALF_WIDTH_PIXELS * np.sqrt(2)
# In-transit cadences are within IN_TRANSIT_FRACTION of the half duration from mid-transit. Out-of-transit cadences
# are between 1 and OUT_OF_TRANSIT_OUTER_DURATIONS durations from mid-transit.
IN_TRANSIT_FRACTION = 0.75
OUT_OF_TRANSIT_OUTER_DURATIONS = 2.0
MINIMUM_IN_TRANSIT_CADENCES = 3
MINIMUM_OUT_OF_TRANSIT_CADENCES = 6
MINIMUM_TRANSITS = 2
MINIMUM_DIFFERENCE_IMAGE_SNR = 5.0
# A chi-square difference of 9 is 3 sigma for one parameter.
DELTA_CHI_SQUARE_DECISION = 9.0
# The WCS and the PRF disagree at the target by up to about half a pixel, so the fit cannot separate stars closer
# than this.
UNRESOLVABLE_SEPARATION_PIXELS = 0.5

OUTPUT_COLUMNS = [
    "tic_id", "period_days", "epoch_btjd", "duration_hours", "depth_ppm", "tmag", "check_verdict", "sector", "transits_used",
    "target_position_offset_pixels", "stars_tested", "stars_within_half_pixel_of_target", "difference_image_snr", "best_source", "best_source_is_target",
    "target_chi_square_above_best", "best_other_star_chi_square_above_target", "free_fit_offset_from_target_arcsec",
    "pixel_verdict", "error",
]
CANDIDATE_OUTPUT_COLUMNS = [
    "tic_id", "period_days", "epoch_btjd", "duration_hours", "depth_ppm", "tmag", "check_verdict", "sectors_checked",
    "pixel_verdicts_by_sector", "pixel_verdict", "verification_date_utc", "code_version",
]


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_candidate_pixel_checks"


def cutout_path(tic_id, ra, dec, sector):
    directory = CUTOUT_DIRECTORY / str(tic_id)
    existing = sorted(directory.glob(f"tess-s{sector:04d}-*_astrocut.fits"))
    if existing:
        return existing[0]
    directory.mkdir(parents=True, exist_ok=True)
    manifest = Tesscut.download_cutouts(
        coordinates=SkyCoord(ra, dec, unit="deg"), size=CUTOUT_SIZE_PIXELS, sector=sector, path=str(directory)
    )
    return directory / manifest["Local Path"][0].split("/")[-1]


def load_cutout(path):
    """Return time, flux cube, WCS, camera, CCD and detector corner column and row of a TESScut cutout."""
    with fits.open(path) as hdul:
        table = hdul[1].data
        time, flux, quality = np.array(table["TIME"]), np.array(table["FLUX"]), np.array(table["QUALITY"])
        keep = (quality == 0) & np.isfinite(time) & np.all(np.isfinite(flux), axis=(1, 2))
        return (time[keep], flux[keep], WCS(hdul[2].header), hdul[0].header["CAMERA"], hdul[0].header["CCD"],
                hdul[1].header["1CRV4P"], hdul[1].header["2CRV4P"])


def difference_image(time, flux, period, epoch, duration):
    """Return the mean difference image, its per-pixel error, the out-of-transit image and the number of transits."""
    phase = (time - epoch + 0.5 * period) % period - 0.5 * period
    transit_numbers = np.round((time - epoch) / period)
    in_transit = np.abs(phase) < IN_TRANSIT_FRACTION * 0.5 * duration
    out_of_transit = (np.abs(phase) > duration) & (np.abs(phase) < OUT_OF_TRANSIT_OUTER_DURATIONS * duration)
    differences = []
    for number in np.unique(transit_numbers[in_transit]):
        this_in = in_transit & (transit_numbers == number)
        this_out = out_of_transit & (transit_numbers == number)
        if this_in.sum() >= MINIMUM_IN_TRANSIT_CADENCES and this_out.sum() >= MINIMUM_OUT_OF_TRANSIT_CADENCES:
            # Subtract the median of each transit's difference image. Scattered light changes the background of
            # the whole cutout from transit to transit, and would add to the scatter of every pixel.
            transit_difference = flux[this_out].mean(axis=0) - flux[this_in].mean(axis=0)
            differences.append(transit_difference - np.median(transit_difference))
    if len(differences) < MINIMUM_TRANSITS:
        raise ValueError(f"only {len(differences)} usable transits in the cutout")
    differences = np.array(differences)
    error = differences.std(axis=0, ddof=1) / np.sqrt(len(differences))
    return differences.mean(axis=0), error, np.median(flux[out_of_transit], axis=0), len(differences)


class PixelResponseModel:
    """The TESS PRF on the cutout pixel grid, placed by cutout pixel coordinates (0 at the first pixel centre)."""

    def __init__(self, sector, camera, ccd, corner_column, corner_row):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.prf = SimpleTessPRF(
                shape=(CUTOUT_SIZE_PIXELS, CUTOUT_SIZE_PIXELS), sector=sector, camera=camera, ccd=ccd,
                column=corner_column, row=corner_row,
            )
        # Find the shift from cutout pixel coordinates to PRF coordinates, from the centroid of a centred model.
        centre = CUTOUT_SIZE_PIXELS // 2
        self.column_shift, self.row_shift = 0.0, 0.0
        image = self.image(centre, centre)
        rows, columns = np.indices(image.shape)
        self.column_shift = centre - (image * columns).sum() / image.sum()
        self.row_shift = centre - (image * rows).sum() / image.sum()

    def image(self, column, row):
        image = self.prf.evaluate(self.prf.col_coord[0] + column + self.column_shift,
                                  self.prf.row_coord[0] + row + self.row_shift)
        return image / image.sum()


def fit_source(data, error, model_image):
    """Fit data = amplitude * model_image + offset by weighted least squares; return chi-square, amplitude, its error."""
    weights = 1 / error**2
    design = np.vstack([model_image.ravel(), np.ones(model_image.size)]).T
    weighted_design = design * np.sqrt(weights.ravel())[:, None]
    weighted_data = data.ravel() * np.sqrt(weights.ravel())
    coefficients, *_ = np.linalg.lstsq(weighted_design, weighted_data, rcond=None)
    covariance = np.linalg.inv(weighted_design.T @ weighted_design)
    chi_square = float(np.sum((weighted_data - weighted_design @ coefficients) ** 2))
    return chi_square, coefficients[0], np.sqrt(covariance[0, 0])


def gaia_sources(ra, dec, radius_arcsec):
    query = f"""
        SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius_arcsec / 3600}))
        AND phot_g_mean_mag IS NOT NULL
    """
    return Gaia.launch_job(query).get_results()


def stars_able_to_cause_dip(candidate, ra, dec):
    """Return the target and the Gaia sources bright enough to cause the dip: [(name, ra, dec, G)], target first."""
    sources = gaia_sources(ra, dec, SOURCE_SEARCH_RADIUS_PIXELS * PIXEL_SIZE_ARCSEC)
    separations = SkyCoord(sources["ra"], sources["dec"], unit="deg").separation(SkyCoord(ra, dec, unit="deg"))
    target_source_id = candidate.get("gaia_source_id", "")
    if target_source_id and int(target_source_id) in set(int(value) for value in sources["source_id"]):
        target_index = [int(value) for value in sources["source_id"]].index(int(target_source_id))
    else:
        target_index = int(np.argmin(separations))
    target_g = float(sources["phot_g_mean_mag"][target_index])
    depth = float(candidate["depth_ppm"]) * 1e-6
    stars = [(f"TIC {candidate['tic_id']}", float(sources["ra"][target_index]), float(sources["dec"][target_index]),
              target_g)]
    for index, source in enumerate(sources):
        if index == target_index:
            continue
        # A neighbour causes the dip only when its eclipse can remove depth times the target flux.
        if 10 ** (-0.4 * (float(source["phot_g_mean_mag"]) - target_g)) >= depth:
            stars.append((f"Gaia DR3 {source['source_id']}", float(source["ra"]), float(source["dec"]),
                          float(source["phot_g_mean_mag"])))
    return stars


def check_pixels(candidate, star_row, sector, figure_directory):
    ra, dec = float(star_row["ra"]), float(star_row["dec"])
    period, epoch = float(candidate["period_days"]), float(candidate["epoch_btjd"])
    duration = float(candidate["duration_hours"]) / 24
    time, flux, wcs, camera, ccd, corner_column, corner_row = load_cutout(
        cutout_path(candidate["tic_id"], ra, dec, sector)
    )
    difference, error, out_of_transit_image, transits_used = difference_image(time, flux, period, epoch, duration)
    model = PixelResponseModel(sector, camera, ccd, corner_column, corner_row)
    stars = stars_able_to_cause_dip(candidate, ra, dec)
    positions = [wcs.world_to_pixel(SkyCoord(star_ra, star_dec, unit="deg")) for _, star_ra, star_dec, _ in stars]
    target_column, target_row = (float(value) for value in positions[0])

    # Fit only the pixels near the target, and only the stars inside them.
    centre_column, centre_row = int(round(target_column)), int(round(target_row))
    window = (slice(max(centre_row - FIT_HALF_WIDTH_PIXELS, 0), centre_row + FIT_HALF_WIDTH_PIXELS + 1),
              slice(max(centre_column - FIT_HALF_WIDTH_PIXELS, 0), centre_column + FIT_HALF_WIDTH_PIXELS + 1))
    data, data_error = difference[window], error[window]
    inside = [index for index, (column, row) in enumerate(positions)
              if abs(column - centre_column) <= FIT_HALF_WIDTH_PIXELS
              and abs(row - centre_row) <= FIT_HALF_WIDTH_PIXELS]

    # Measure the offset between the WCS and the PRF at the target, from a free fit to the out-of-transit image.
    def out_of_transit_chi_square(position):
        return fit_source(out_of_transit_image[window], np.sqrt(np.abs(out_of_transit_image[window])) + 1,
                          model.image(*position)[window])[0]
    target_fit = minimize(out_of_transit_chi_square, [target_column, target_row], method="Nelder-Mead")
    target_position_offset = float(np.hypot(*(target_fit.x - [target_column, target_row])))

    fits_by_star = {index: fit_source(data, data_error, model.image(*positions[index])[window]) for index in inside}
    free_fit = minimize(lambda position: fit_source(data, data_error, model.image(*position)[window])[0],
                        [target_column, target_row], method="Nelder-Mead")
    # Only a star whose fit loses light in the transit can cause the dip.
    dimming = [index for index in fits_by_star if fits_by_star[index][1] > 0] or list(fits_by_star)
    best_index = min(dimming, key=lambda index: fits_by_star[index][0])
    degrees_of_freedom = data.size - 2
    scale = max(1.0, fits_by_star[best_index][0] / degrees_of_freedom)
    best_chi_square, best_amplitude, best_amplitude_error = fits_by_star[best_index]
    snr = best_amplitude / best_amplitude_error / np.sqrt(scale)
    # The fit cannot separate the target from a star closer than UNRESOLVABLE_SEPARATION_PIXELS, so it treats them
    # as one group.
    target_group = [index for index in inside
                    if np.hypot(positions[index][0] - target_column, positions[index][1] - target_row)
                    < UNRESOLVABLE_SEPARATION_PIXELS]
    target_group_chi_square = min(fits_by_star[index][0] for index in target_group)
    target_above_best = (target_group_chi_square - best_chi_square) / scale
    others = [fits_by_star[index][0] for index in inside if index not in target_group]
    best_other_above_target = (min(others) - target_group_chi_square) / scale if others else np.inf

    # When no star fits a dip, the amplitude and the SNR are negative, and the verdict is too faint.
    if snr < MINIMUM_DIFFERENCE_IMAGE_SNR:
        pixel_verdict = "too faint"
    elif best_index not in target_group and target_above_best >= DELTA_CHI_SQUARE_DECISION:
        pixel_verdict = "on neighbour"
    elif best_index in target_group and best_other_above_target >= DELTA_CHI_SQUARE_DECISION:
        pixel_verdict = "on target"
    else:
        pixel_verdict = "unresolved"

    draw_figure(figure_directory / f"tic{candidate['tic_id']}_{period:.5f}_s{sector:04d}.png", out_of_transit_image,
                difference, stars, positions, inside, best_index, free_fit.x, candidate, sector, pixel_verdict)
    return {
        "sector": sector,
        "transits_used": transits_used,
        "target_position_offset_pixels": round(target_position_offset, 2),
        "stars_tested": len(inside),
        "stars_within_half_pixel_of_target": len(target_group) - 1,
        "difference_image_snr": round(float(snr), 1),
        "best_source": stars[best_index][0],
        "best_source_is_target": best_index == 0,
        "target_chi_square_above_best": round(float(target_above_best), 1),
        "best_other_star_chi_square_above_target": round(float(best_other_above_target), 1),
        "free_fit_offset_from_target_arcsec": round(
            float(np.hypot(*(free_fit.x - [target_column, target_row]))) * PIXEL_SIZE_ARCSEC, 1
        ),
        "pixel_verdict": pixel_verdict,
    }


def draw_figure(path, out_of_transit_image, difference, stars, positions, inside, best_index, free_position,
                candidate, sector, pixel_verdict):
    figure, axes = plt.subplots(1, 2, figsize=(11, 5))
    for axis, image, title in ((axes[0], out_of_transit_image, "Out of transit"),
                               (axes[1], difference, "Out of transit minus in transit")):
        shown = axis.imshow(image, origin="lower", cmap="viridis")
        figure.colorbar(shown, ax=axis, fraction=0.046)
        for index in inside:
            column, row = positions[index]
            marker = "*" if index == 0 else "o"
            colour = "red" if index == best_index else "white"
            axis.scatter(column, row, marker=marker, s=160 if index == 0 else 50, facecolors="none",
                         edgecolors=colour, linewidths=1.5)
        axis.scatter(*free_position, marker="x", color="orange", s=80)
        axis.set_title(title)
    figure.suptitle(f"TIC {candidate['tic_id']}, P = {float(candidate['period_days']):.5f} d, Sector {sector}: "
                    f"{pixel_verdict}. "
                    f"Star: target. Red: best star. Cross: free fit.")
    figure.tight_layout()
    figure.savefig(path, dpi=90)
    plt.close(figure)


def sectors_with_full_frame_images(tic_id, ra, dec):
    """Return the sectors in which the star falls on the science pixels and TESScut has full-frame images."""
    on_science_pixels = set(int(value) for value in tess_stars2px_function_entry(int(tic_id), ra, dec)[3])
    available = set(int(value) for value in Tesscut.get_sectors(coordinates=SkyCoord(ra, dec, unit="deg"))["sector"])
    return sorted(on_science_pixels & available)


def combined_pixel_verdict(sector_verdicts):
    """Combine the verdicts of the sectors. A dip on a neighbour in any sector outweighs the other sectors."""
    for verdict in ("on neighbour", "on target", "unresolved", "too faint"):
        if verdict in sector_verdicts:
            return verdict
    return ""


def append_verification_record(sector, candidate_row):
    """Append the row of a candidate, as written to pixel_checks.csv, to its verification record.

    A candidate with no pixel verdict gets no record row: the next run checks that candidate again.
    """
    if not candidate_row["pixel_verdict"]:
        return
    append_verification_record_entry(sector, candidate_row["tic_id"], candidate_row["period_days"], "pixel check",
                                     candidate_row["pixel_verdict"], candidate_row["pixel_verdicts_by_sector"],
                                     candidate_row)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    parser.add_argument("--restart", action="store_true",
                        help="check all candidates again, not only those without a row in pixel_checks.csv")
    arguments = parser.parse_args()
    sector = arguments.sector

    directory = output_directory(sector)
    directory.mkdir(parents=True, exist_ok=True)
    sector_output_path, candidate_output_path = directory / "pixel_checks_by_sector.csv", directory / "pixel_checks.csv"
    # A candidate row with no pixel verdict has an error in each sector, so the script checks that candidate again.
    earlier_candidate_rows = [row for row in earlier_rows_to_keep(candidate_output_path, arguments.restart)
                              if row["pixel_verdict"]]
    checked = {(row["tic_id"], row["period_days"]) for row in earlier_candidate_rows}
    # The script keeps the sector rows of the kept candidates only. A stop between the two files can leave sector rows
    # of a candidate with no candidate row.
    earlier_sector_rows = [row for row in complete_csv_rows(sector_output_path)
                           if (row["tic_id"], row["period_days"]) in checked]
    checks = list(csv.DictReader(open(candidate_checks_output_directory(sector) / "candidate_checks.csv")))
    not_refuted = [row for row in checks if row["verdict"] != "refuted"]
    candidates = [row for row in not_refuted if (row["tic_id"], row["period_days"]) not in checked]
    star_rows = {row["tic_id"]: row for row in csv.DictReader(open(sector_directory(sector) / "star_list.csv"))}
    print(f"Sector {sector}: {len(not_refuted)} candidates not refuted, {len(checked)} checked in an earlier run",
          flush=True)

    # A stop between a CSV row and its record row leaves the record row missing. This appends it.
    for candidate_row in earlier_candidate_rows:
        append_verification_record(sector, candidate_row)

    sector_writer = IncrementalCsvWriter(sector_output_path, OUTPUT_COLUMNS, earlier_sector_rows)
    candidate_writer = IncrementalCsvWriter(candidate_output_path, CANDIDATE_OUTPUT_COLUMNS, earlier_candidate_rows)
    with sector_writer, candidate_writer:
        for count, candidate in enumerate(candidates, start=1):
            star_row = star_rows[candidate["tic_id"]]
            try:
                sectors = sectors_with_full_frame_images(candidate["tic_id"], float(star_row["ra"]),
                                                         float(star_row["dec"]))
            except Exception as error:
                sectors = [sector]
                print(f"  TIC {candidate['tic_id']}: sector list failed ({error}); Sector {sector} only", flush=True)
            rows = []
            for pixel_sector in sectors:
                row = {column: candidate.get(column, "") for column in OUTPUT_COLUMNS}
                row["check_verdict"] = candidate["verdict"]
                row["sector"] = pixel_sector
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        row.update(check_pixels(candidate, star_row, pixel_sector, directory))
                except Exception as error:
                    row["error"] = f"{type(error).__name__}: {error}"
                rows.append(row)
            summary = "; ".join(
                f"S{row['sector']} {row.get('pixel_verdict') or 'error: ' + row['error']}"
                + (f" (SNR {row['difference_image_snr']}, best star {row['best_source']})"
                   if row.get("pixel_verdict") else "")
                for row in rows
            )
            candidate_row = {column: candidate.get(column, "") for column in CANDIDATE_OUTPUT_COLUMNS}
            candidate_row["check_verdict"] = candidate["verdict"]
            candidate_row["sectors_checked"] = " ".join(str(row["sector"]) for row in rows)
            candidate_row["pixel_verdicts_by_sector"] = summary
            candidate_row["pixel_verdict"] = combined_pixel_verdict([row.get("pixel_verdict", "") for row in rows])
            candidate_row.update(verification_stamp())
            for row in rows:
                sector_writer.write_row(row)
            candidate_row = candidate_writer.write_row(candidate_row)
            append_verification_record(sector, candidate_row)
            print(f"  {count} of {len(candidates)}: TIC {candidate['tic_id']} {candidate_row['pixel_verdict']}: "
                  f"{summary}", flush=True)
    print(f"Wrote {sector_output_path.relative_to(REPOSITORY_ROOT)} and "
          f"{candidate_output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
