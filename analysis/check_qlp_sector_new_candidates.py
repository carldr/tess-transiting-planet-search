"""Check each new candidate of one sector against independent data, to support or refute it.

The script reads new_candidates.csv from the cross-match, and runs these checks on each candidate:
1. Other TESS sectors: it gets the star's SPOC, TESS-SPOC or QLP light curve of each other sector.
   - It searches each sector alone near the candidate period. A sector refutes the candidate when the transit should
     have an SNR of at least ABSENT_EXPECTED_SNR there, but the search does not detect it.
   - It searches all sectors together, with this sector, in a narrow window around the candidate period. A real
     transit stays in phase over the years between sectors, so its SNR grows. The script searches the same light
     curves at null periods away from the candidate period. The joint search detects the transit when its SNR is
     JOINT_SEARCH_SIGNIFICANCE_FACTOR times the highest null SNR.
   A single-sector search alone cannot support a candidate: near a wrong period, it reaches SNR 7 in 12% of sectors.
2. Period pile-up: a period at which the sector search finds many more signals than at nearby periods is a
   systematic of the sector, and refutes the candidate.
3. Ephemeris match: it looks for a signal at the same period on another star within EPHEMERIS_MATCH_RADIUS_ARCMIN.
   It compares the signals that the sector search found, the TOIs and the TESS EBs. When the other star's dip
   removes more flux than the candidate's dip, the other star is the likely source, and it refutes the candidate.
4. Gaia DR3 target: a RUWE above 1.4 or a non-single-star solution means an unresolved companion.
   A Gaia eclipsing binary classification refutes the candidate.
5. Gaia DR3 neighbours within NEIGHBOUR_RADIUS_ARCSEC: a neighbour that could cause the dip with an eclipse depth of
   at most 100% is a possible source. Such a neighbour with a Gaia eclipsing binary classification refutes the
   candidate.
6. VSX: it lists the known variable stars within NEIGHBOUR_RADIUS_ARCSEC. An eclipsing or ellipsoidal binary
   refutes the candidate.
7. Companion: it reads the LEO-Vetter transit fit. A companion radius above COMPANION_MAXIMUM_RADIUS_EARTH_RADII,
   with complete TIC stellar parameters, refutes the candidate, unless the transit is grazing: a grazing fit leaves
   the radius almost unconstrained. A grazing transit, or an a/R* that does not fit the TIC star, is a caution.

The verdict is "refuted" when any check refutes the candidate, "supported" when the joint search detects the transit
and no check refutes it, and "unconfirmed" otherwise.

Run from the repository root, after cross_match_qlp_sector_signals_with_known_catalogues.py:
    .venv/bin/python analysis/check_qlp_sector_new_candidates.py 98

The script saves the other-sector light curves in data/candidate_other_sector_light_curves/.
Output: results/qlp_sector_<nnnn>_new_candidate_checks/candidate_checks.csv, one row per candidate.
The script writes each row as soon as it has it. It keeps the rows of an earlier run without an error, and checks
only the other candidates. --restart checks all candidates again.
The script also appends the result of each candidate to its verification record.
"""

import argparse
import csv
import warnings

import astropy.units as units
import lightkurve
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.timeseries import BoxLeastSquares
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astroquery.vizier import Vizier

from build_qlp_sector_star_list import REPOSITORY_ROOT, sector_directory
from candidate_verification_records import append_verification_record_entry, verification_stamp
from cross_match_qlp_sector_signals_with_known_catalogues import (
    CATALOGUE_DIRECTORY, output_directory as cross_match_output_directory,
)
from download_qlp_sector_light_curves import light_curve_path
from incremental_csv_output import IncrementalCsvWriter, earlier_rows_to_keep
from pipeline_status import PipelineStatus
from search_qlp_sector import (
    DURATIONS_DAYS, FREQUENCY_FACTOR, MAXIMUM_PERIOD_DAYS, MINIMUM_PERIOD_DAYS, detrend, load_light_curve,
    output_directory as search_output_directory, robust_standard_deviation, select_usable_cadences,
)

OTHER_SECTOR_LIGHT_CURVE_DIRECTORY = REPOSITORY_ROOT / "data/candidate_other_sector_light_curves"
# Light curve products in order of preference, and the flux column of each.
FLUX_COLUMN_BY_AUTHOR = {"SPOC": "pdcsap_flux", "TESS-SPOC": "pdcsap_flux", "QLP": "sap_flux"}
OTHER_SECTOR_PERIOD_FRACTION = 0.01
DETECTION_SNR = 7.0
ABSENT_EXPECTED_SNR = 10.0
JOINT_SEARCH_BIN_DAYS = 20 / 1440
JOINT_SEARCH_PERIOD_FRACTION = 0.003
# Null periods: the same light curves searched at these multiples of the candidate period, away from its harmonics.
JOINT_SEARCH_NULL_PERIOD_FACTORS = (0.71, 0.77, 0.83, 0.89, 0.94, 1.06, 1.13, 1.21, 1.29, 1.41)
JOINT_SEARCH_SIGNIFICANCE_FACTOR = 1.5

# A period at which the sector search finds this many times more signals than at nearby periods is a systematic.
PERIOD_PILE_UP_RATIO = 5.0
PERIOD_PILE_UP_WINDOW_FRACTION = 0.001
PERIOD_PILE_UP_REFERENCE_FRACTION = 0.2

EPHEMERIS_MATCH_RADIUS_ARCMIN = 10.0
EPHEMERIS_PERIOD_RATIOS = (0.5, 1.0, 2.0)
EPHEMERIS_PERIOD_TOLERANCE = 0.001
SEARCH_SIGNAL_MINIMUM_SNR = 7.0

NEIGHBOUR_RADIUS_ARCSEC = 60.0
RUWE_BINARY_THRESHOLD = 1.4
# VSX types of eclipsing and ellipsoidal binaries.
VSX_BINARY_TYPE_PREFIXES = ("E", "ELL")

# 2 Jupiter radii. A larger companion is a star: the largest known planets are about 2 Jupiter radii.
COMPANION_MAXIMUM_RADIUS_EARTH_RADII = 22.4
# a/R* from the transit fit, divided by a/R* from the TIC star, outside 1/2 to 2 means a stellar density wrong by a
# factor of 8 or more. The transit then does not fit the TIC star, or the orbit is very eccentric.
SCALED_SEMI_MAJOR_AXIS_RATIO_LIMIT = 2.0

OUTPUT_COLUMNS = [
    "tic_id", "period_days", "epoch_btjd", "duration_hours", "depth_ppm", "MES", "tmag",
    "other_sectors_searched", "other_sectors_absent", "other_sector_results",
    "joint_search_snr", "joint_search_period_days", "joint_search_highest_null_snr", "joint_search_detected",
    "period_pile_up_ratio",
    "ephemeris_match_likely_sources", "ephemeris_match_contaminated_stars",
    "gaia_source_id", "gaia_ruwe", "gaia_non_single_star", "gaia_eclipsing_binary",
    "gaia_neighbours", "gaia_neighbours_able_to_cause_dip", "gaia_neighbour_eclipsing_binaries",
    "vsx_variables", "vsx_binaries",
    "companion_radius_earth_radii", "tic_stellar_parameters_complete", "transit_grazing",
    "fitted_to_expected_scaled_semi_major_axis",
    "refuted_by", "cautions", "verdict", "error", "verification_date_utc", "code_version",
]
VERIFICATION_RECORD_STEP = "independent data checks"


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_new_candidate_checks"


# 1. Other TESS sectors.

def other_sector_products(tic_id, sector):
    """Return [(sector, author, search result row)] for each other sector, with the preferred product of each."""
    search = lightkurve.search_lightcurve(f"TIC {tic_id}", author=list(FLUX_COLUMN_BY_AUTHOR))
    search = search[search.table["target_name"] == str(tic_id)]
    best_product_by_sector = {}
    for index, row in enumerate(search.table):
        product_sector = int(row["mission"].split()[-1])
        rank = list(FLUX_COLUMN_BY_AUTHOR).index(row["author"])
        if product_sector != sector and rank < best_product_by_sector.get(product_sector, (99, None))[0]:
            best_product_by_sector[product_sector] = (rank, index)
    return [(product_sector, list(FLUX_COLUMN_BY_AUTHOR)[rank], search[index])
            for product_sector, (rank, index) in sorted(best_product_by_sector.items())]


def load_other_sector_light_curve(product, author):
    """Return time and normalised flux of the usable cadences. Download the file again if the saved copy is broken."""
    arguments = {"download_dir": str(OTHER_SECTOR_LIGHT_CURVE_DIRECTORY), "quality_bitmask": "none",
                 "flux_column": FLUX_COLUMN_BY_AUTHOR[author]}
    try:
        light_curve = product.download(**arguments)
    except Exception:
        light_curve = product.download(cache=False, **arguments)
    return select_usable_cadences(
        light_curve.time.value, np.ma.filled(light_curve.flux.value, np.nan), np.asarray(light_curve.quality.value),
    )


def search_one_sector(time, flux, candidate):
    """Search one detrended sector near the candidate period. Return the BLS SNR, and the SNR the transit should have."""
    period, depth = float(candidate["period_days"]), float(candidate["depth_ppm"]) * 1e-6
    duration = float(candidate["duration_hours"]) / 24
    noise = robust_standard_deviation(flux)
    expected_snr = depth / noise * np.sqrt(len(time) * duration / period)
    model = BoxLeastSquares(time, flux, np.full_like(flux, noise))
    result = model.autopower(
        DURATIONS_DAYS[DURATIONS_DAYS < period / 2],
        minimum_period=period * (1 - OTHER_SECTOR_PERIOD_FRACTION),
        maximum_period=period * (1 + OTHER_SECTOR_PERIOD_FRACTION),
        frequency_factor=FREQUENCY_FACTOR,
        objective="likelihood",
    )
    return float(result.depth_snr[int(np.argmax(result.power))]), float(expected_snr)


def binned(time, flux):
    """Return the mean time and flux in bins of JOINT_SEARCH_BIN_DAYS."""
    bin_numbers = np.floor((time - time.min()) / JOINT_SEARCH_BIN_DAYS).astype(int)
    counts = np.bincount(bin_numbers)
    filled = counts > 0
    return (np.bincount(bin_numbers, time)[filled] / counts[filled],
            np.bincount(bin_numbers, flux)[filled] / counts[filled])


def joint_search(time, flux, period, duration):
    """Search all sectors together, within JOINT_SEARCH_PERIOD_FRACTION of period. Return the best SNR and period."""
    noise = robust_standard_deviation(flux)
    model = BoxLeastSquares(time, flux, np.full_like(flux, noise))
    # The frequency step keeps the phase drift over the whole baseline below 0.3 transit durations.
    frequency_step = 0.3 * duration / (time.max() - time.min()) ** 2
    frequencies = np.arange(1 / (period * (1 + JOINT_SEARCH_PERIOD_FRACTION)),
                            1 / (period * (1 - JOINT_SEARCH_PERIOD_FRACTION)), frequency_step)
    result = model.power(1 / frequencies, duration, objective="likelihood")
    best = int(np.argmax(result.power))
    return float(result.depth_snr[best]), float(result.period[best])


def check_other_sectors(candidate, sector):
    """Search each other sector alone, then all sectors together.

    A sector refutes the candidate when the transit should have an SNR of at least ABSENT_EXPECTED_SNR there, but the
    sector search finds less than DETECTION_SNR. The joint search supports the candidate when its SNR at the candidate
    period is JOINT_SEARCH_SIGNIFICANCE_FACTOR times the highest SNR at the null periods of the same light curves.
    """
    tic_id = int(candidate["tic_id"])
    period, duration = float(candidate["period_days"]), float(candidate["duration_hours"]) / 24
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        binned_sectors = [binned(*detrend(*load_light_curve(light_curve_path(sector, tic_id))))]
        searched, absent, results = [], [], []
        for product_sector, author, product in other_sector_products(tic_id, sector):
            try:
                time, flux = detrend(*load_other_sector_light_curve(product, author))
                snr, expected_snr = search_one_sector(time, flux, candidate)
            except Exception as error:
                results.append(f"S{product_sector} {author}: {type(error).__name__}: {error}")
                continue
            searched.append(product_sector)
            binned_sectors.append(binned(time, flux))
            if snr < DETECTION_SNR and expected_snr >= ABSENT_EXPECTED_SNR:
                absent.append(product_sector)
            results.append(f"S{product_sector} {author}: SNR {snr:.1f}, expected SNR {expected_snr:.1f}")
        row = {
            "other_sectors_searched": " ".join(map(str, searched)),
            "other_sectors_absent": " ".join(map(str, absent)),
            "other_sector_results": "; ".join(results),
        }
        if searched:
            time = np.concatenate([times for times, _ in binned_sectors])
            flux = np.concatenate([fluxes for _, fluxes in binned_sectors])
            joint_snr, joint_period = joint_search(time, flux, period, duration)
            null_snrs = [joint_search(time, flux, period * factor, duration)[0]
                         for factor in JOINT_SEARCH_NULL_PERIOD_FACTORS]
            row.update({
                "joint_search_snr": round(joint_snr, 1),
                "joint_search_period_days": round(joint_period, 6),
                "joint_search_highest_null_snr": round(max(null_snrs), 1),
                "joint_search_detected": joint_snr >= JOINT_SEARCH_SIGNIFICANCE_FACTOR * max(null_snrs),
            })
    return row


# 2. Period pile-up.

def check_period_pile_up(candidate, search_signals):
    """Return how many times more sector signals are at the candidate period than at nearby periods."""
    period, periods = float(candidate["period_days"]), search_signals["periods"]
    local_count = np.sum(np.abs(periods / period - 1) < PERIOD_PILE_UP_WINDOW_FRACTION)
    reference_low = max(period * (1 - PERIOD_PILE_UP_REFERENCE_FRACTION), MINIMUM_PERIOD_DAYS)
    reference_high = min(period * (1 + PERIOD_PILE_UP_REFERENCE_FRACTION), MAXIMUM_PERIOD_DAYS)
    reference_count = np.sum((periods >= reference_low) & (periods <= reference_high))
    local_density = local_count / (2 * PERIOD_PILE_UP_WINDOW_FRACTION * period)
    reference_density = reference_count / (reference_high - reference_low)
    return {"period_pile_up_ratio": round(float(local_density / reference_density), 1)}


# 3. Ephemeris match.

def periods_match_within_tolerance(period, other_periods):
    return np.any([np.abs(other_periods / period / ratio - 1) < EPHEMERIS_PERIOD_TOLERANCE
                   for ratio in EPHEMERIS_PERIOD_RATIOS], axis=0)


def flux_deficit(depth_ppm, tmag):
    """Return the flux that a dip removes, in units of the flux of a star of TESS magnitude 0."""
    return np.asarray(depth_ppm, dtype=float) * 1e-6 * 10 ** (-0.4 * np.asarray(tmag, dtype=float))


def load_ephemeris_sources(sector):
    """Return the signals to match against, as a list of dicts of arrays, and {TIC ID: (ra, dec, tmag)}.

    Each dict has label, tic_ids, coordinates, periods, epochs (None for catalogues) and flux_deficits
    (infinite for EBs, whose eclipses are deep).
    """
    stars = {}
    for row in csv.DictReader(open(sector_directory(sector) / "star_list.csv")):
        stars[int(row["tic_id"])] = (float(row["ra"]), float(row["dec"]), float(row["tmag"]))
    signals = []
    for chunk_path in sorted(search_output_directory(sector).glob("chunk_*[0-9].csv")):
        for row in csv.DictReader(open(chunk_path)):
            if row["signal_number"] != "0" and float(row["snr"]) >= SEARCH_SIGNAL_MINIMUM_SNR:
                signals.append((int(row["tic_id"]), float(row["period_days"]), float(row["epoch_btjd"]),
                                float(row["depth_ppm"])))
    tic_ids = np.array([signal[0] for signal in signals])
    star_values = np.array([stars[tic_id] for tic_id in tic_ids])
    search_signals = {
        "label": f"S{sector} search", "tic_ids": tic_ids,
        "coordinates": SkyCoord(star_values[:, 0], star_values[:, 1], unit="deg"),
        "periods": np.array([signal[1] for signal in signals]),
        "epochs": np.array([signal[2] for signal in signals]),
        "flux_deficits": flux_deficit([signal[3] for signal in signals], star_values[:, 2]),
    }

    tois = [row for row in csv.DictReader(open(CATALOGUE_DIRECTORY / "toi.csv"))
            if row["Period (days)"] and row["Depth (ppm)"] and row["TESS Mag"]]
    toi_signals = {
        "label": "TOI", "tic_ids": np.array([int(row["TIC ID"]) for row in tois]),
        "coordinates": SkyCoord([row["RA"] for row in tois], [row["Dec"] for row in tois], unit=("hourangle", "deg")),
        "periods": np.array([float(row["Period (days)"]) for row in tois]), "epochs": None,
        "flux_deficits": flux_deficit([row["Depth (ppm)"] for row in tois], [row["TESS Mag"] for row in tois]),
    }

    ebs = [row for row in csv.DictReader(open(CATALOGUE_DIRECTORY / "tess_eb.csv")) if row["period_days"]]
    eb_signals = {
        "label": "TESS EB", "tic_ids": np.array([int(row["tic_id"]) for row in ebs]),
        "coordinates": SkyCoord([float(row["ra"]) for row in ebs], [float(row["dec"]) for row in ebs], unit="deg"),
        "periods": np.array([float(row["period_days"]) for row in ebs]), "epochs": None,
        "flux_deficits": np.full(len(ebs), np.inf),
    }
    return [search_signals, toi_signals, eb_signals], stars


def check_ephemeris_match(candidate, ephemeris_sources, stars):
    """Find signals at the candidate period on other stars nearby.

    The star whose dip removes more flux is the likely source. A nearby star that is the likely source refutes the
    candidate. A nearby star that the candidate's star contaminates is only listed.
    """
    tic_id = int(candidate["tic_id"])
    period, epoch = float(candidate["period_days"]), float(candidate["epoch_btjd"])
    duration = float(candidate["duration_hours"]) / 24
    ra, dec, tmag = stars[tic_id]
    position = SkyCoord(ra, dec, unit="deg")
    candidate_flux_deficit = flux_deficit(candidate["depth_ppm"], tmag)
    sources, contaminated = [], []
    for signals in ephemeris_sources:
        separations = position.separation(signals["coordinates"]).arcmin
        selected = (separations < EPHEMERIS_MATCH_RADIUS_ARCMIN) & (signals["tic_ids"] != tic_id)
        selected &= periods_match_within_tolerance(period, signals["periods"])
        if signals["epochs"] is not None:
            phase_offset = (signals["epochs"] - epoch + period / 2) % period - period / 2
            selected &= np.abs(phase_offset) < duration
        for index in np.flatnonzero(selected):
            description = (f"{signals['label']} TIC {signals['tic_ids'][index]} at {separations[index]:.1f}' "
                           f"P {signals['periods'][index]:.5f} d")
            if signals["flux_deficits"][index] >= candidate_flux_deficit:
                sources.append(description)
            else:
                contaminated.append(description)
    return {"ephemeris_match_likely_sources": "; ".join(sources),
            "ephemeris_match_contaminated_stars": "; ".join(contaminated)}


# 4 and 5. Gaia DR3 target and neighbours.

def gaia_dr2_source_ids_by_tic_id(tic_ids):
    """Return {TIC ID: Gaia DR2 source ID}. The TIC gives Gaia DR2 IDs, and some differ from the DR3 IDs."""
    # MAST refuses a query with no IDs. A run with every candidate done has no IDs.
    if not tic_ids:
        return {}
    table =Catalogs.query_criteria(catalog="Tic", ID=[int(tic_id) for tic_id in tic_ids])
    return {int(row["ID"]): int(row["GAIA"]) for row in table if str(row["GAIA"]).isdigit()}


def gaia_dr3_source_id(gaia_dr2_source_id):
    query = f"""
        SELECT dr3_source_id FROM gaiadr3.dr2_neighbourhood
        WHERE dr2_source_id = {gaia_dr2_source_id} ORDER BY angular_distance
    """
    matches = Gaia.launch_job(query).get_results()
    if len(matches) == 0:
        raise ValueError(f"Gaia DR2 {gaia_dr2_source_id} has no Gaia DR3 counterpart")
    return int(matches["dr3_source_id"][0])


def check_gaia(candidate, gaia_dr2_source_id, stars):
    ra, dec, _ = stars[int(candidate["tic_id"])]
    gaia_source_id = gaia_dr3_source_id(gaia_dr2_source_id)
    query = f"""
        SELECT g.source_id, g.phot_g_mean_mag, g.ruwe, g.non_single_star, eb.source_id AS eclipsing_binary_source_id,
               DISTANCE(POINT(g.ra, g.dec), POINT({ra}, {dec})) * 3600 AS separation_arcsec
        FROM gaiadr3.gaia_source AS g
        LEFT JOIN gaiadr3.vari_eclipsing_binary AS eb ON eb.source_id = g.source_id
        WHERE 1 = CONTAINS(POINT(g.ra, g.dec), CIRCLE({ra}, {dec}, {NEIGHBOUR_RADIUS_ARCSEC / 3600}))
    """
    sources = Gaia.launch_job(query).get_results()
    is_target = sources["source_id"] == gaia_source_id
    if not is_target.any():
        raise ValueError(f"Gaia DR3 {gaia_source_id} is not within {NEIGHBOUR_RADIUS_ARCSEC}\" of the TIC position")
    target = sources[is_target][0]
    neighbours = sources[~is_target]
    depth = float(candidate["depth_ppm"]) * 1e-6
    able, eclipsing = [], []
    for neighbour in neighbours:
        if np.ma.is_masked(neighbour["phot_g_mean_mag"]):
            continue
        flux_ratio = 10 ** (-0.4 * (float(target["phot_g_mean_mag"]) - float(neighbour["phot_g_mean_mag"])))
        required_depth = depth * flux_ratio
        description = (f"Gaia DR3 {neighbour['source_id']} at {float(neighbour['separation_arcsec']):.0f}\" "
                       f"G {float(neighbour['phot_g_mean_mag']):.1f} needs {required_depth:.0%}")
        if required_depth <= 1:
            able.append(description)
            if not np.ma.is_masked(neighbour["eclipsing_binary_source_id"]):
                eclipsing.append(description)
    ruwe = float(target["ruwe"]) if not np.ma.is_masked(target["ruwe"]) else np.nan
    return {
        "gaia_source_id": gaia_source_id,
        "gaia_ruwe": round(ruwe, 2),
        "gaia_non_single_star": int(target["non_single_star"]),
        "gaia_eclipsing_binary": not np.ma.is_masked(target["eclipsing_binary_source_id"]),
        "gaia_neighbours": len(neighbours),
        "gaia_neighbours_able_to_cause_dip": "; ".join(able),
        "gaia_neighbour_eclipsing_binaries": "; ".join(eclipsing),
    }


# 6. VSX.

def check_vsx(candidate, stars):
    ra, dec, _ = stars[int(candidate["tic_id"])]
    position = SkyCoord(ra, dec, unit="deg")
    tables = Vizier(columns=["Name", "Type", "Period", "RAJ2000", "DEJ2000"]).query_region(
        position, radius=NEIGHBOUR_RADIUS_ARCSEC * units.arcsec, catalog="B/vsx/vsx",
    )
    variables, binaries = [], []
    if tables:
        for row in tables[0]:
            separation = position.separation(SkyCoord(row["RAJ2000"], row["DEJ2000"], unit="deg")).arcsec
            if separation > NEIGHBOUR_RADIUS_ARCSEC:
                continue
            period = "" if np.ma.is_masked(row["Period"]) else f" P {float(row['Period']):.5f} d"
            description = f"{row['Name']} {row['Type']}{period} at {separation:.0f}\""
            variables.append(description)
            if str(row["Type"]).startswith(VSX_BINARY_TYPE_PREFIXES):
                binaries.append(description)
    return {"vsx_variables": "; ".join(variables), "vsx_binaries": "; ".join(binaries)}


def to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return np.nan


def check_companion(candidate):
    """Check that the transit fit gives a planet-sized companion that fits the TIC star."""
    radius = to_float(candidate.get("Rp"))
    impact_parameter = to_float(candidate.get("transit_b"))
    radius_ratio = to_float(candidate.get("transit_RpRs"))
    fitted_scaled_semi_major_axis = to_float(candidate.get("transit_aRs"))
    expected_scaled_semi_major_axis = to_float(candidate.get("aRs"))
    return {
        "companion_radius_earth_radii": radius,
        "tic_stellar_parameters_complete": candidate.get("used_default_stellar_parameters") == "False",
        "transit_grazing": bool(impact_parameter + radius_ratio > 1),
        "fitted_to_expected_scaled_semi_major_axis": round(
            fitted_scaled_semi_major_axis / expected_scaled_semi_major_axis, 2
        ),
    }


def verdict(row):
    refuted_by, cautions = [], []
    if row.get("other_sectors_absent"):
        refuted_by.append("absent in other sectors")
    if isinstance(row.get("period_pile_up_ratio"), float) and row["period_pile_up_ratio"] >= PERIOD_PILE_UP_RATIO:
        refuted_by.append("period shared by many stars in the sector")
    if row.get("ephemeris_match_likely_sources"):
        refuted_by.append("nearby star is the likely source")
    if row.get("gaia_eclipsing_binary") is True:
        refuted_by.append("Gaia eclipsing binary")
    if row.get("gaia_neighbour_eclipsing_binaries"):
        refuted_by.append("Gaia eclipsing binary neighbour")
    if row.get("vsx_binaries"):
        refuted_by.append("VSX binary")
    # The LEO-Vetter fit gives no radius error, and a grazing transit leaves the radius almost unconstrained, so a
    # large radius of a grazing transit is a caution.
    if row.get("tic_stellar_parameters_complete") is True and (
        to_float(row.get("companion_radius_earth_radii")) > COMPANION_MAXIMUM_RADIUS_EARTH_RADII
    ):
        if row.get("transit_grazing") is True:
            cautions.append("companion radius above 2 Jupiter radii, from a grazing fit")
        else:
            refuted_by.append("companion too large for a planet")
    if row.get("tic_stellar_parameters_complete") is False:
        cautions.append("TIC stellar parameters missing")
    if row.get("transit_grazing") is True:
        cautions.append("grazing transit")
    ratio = to_float(row.get("fitted_to_expected_scaled_semi_major_axis"))
    if np.isfinite(ratio) and not 1 / SCALED_SEMI_MAJOR_AXIS_RATIO_LIMIT <= ratio <= SCALED_SEMI_MAJOR_AXIS_RATIO_LIMIT:
        cautions.append("duration does not fit the star")
    if isinstance(row.get("gaia_ruwe"), float) and row["gaia_ruwe"] > RUWE_BINARY_THRESHOLD:
        cautions.append("RUWE above 1.4")
    if row.get("gaia_non_single_star", 0):
        cautions.append("Gaia non-single star")
    if row.get("gaia_neighbours_able_to_cause_dip"):
        cautions.append("neighbour able to cause dip")
    if not row.get("other_sectors_searched"):
        cautions.append("no other sector")
    row["refuted_by"] = "; ".join(refuted_by)
    row["cautions"] = "; ".join(cautions)
    if refuted_by:
        row["verdict"] = "refuted"
    elif row.get("joint_search_detected") is True:
        row["verdict"] = "supported"
    else:
        row["verdict"] = "unconfirmed"


def check_candidate(candidate, sector, ephemeris_sources, stars, gaia_dr2_source_id):
    row = {column: candidate.get(column, "") for column in OUTPUT_COLUMNS}
    errors = []
    checks = [
        lambda: check_other_sectors(candidate, sector),
        lambda: check_period_pile_up(candidate, ephemeris_sources[0]),
        lambda: check_ephemeris_match(candidate, ephemeris_sources, stars),
        lambda: check_gaia(candidate, gaia_dr2_source_id, stars),
        lambda: check_vsx(candidate, stars),
        lambda: check_companion(candidate),
    ]
    for check in checks:
        try:
            row.update(check())
        except Exception as error:
            errors.append(f"{type(error).__name__}: {error}")
    row["error"] = "; ".join(errors)
    verdict(row)
    return row


def append_verification_record(sector, row):
    """Append the result row of a candidate, as written to candidate_checks.csv, to its verification record.

    A row with an error gets no record row: the next run checks that candidate again.
    """
    if row["error"]:
        return
    append_verification_record_entry(
        sector, row["tic_id"], row["period_days"], VERIFICATION_RECORD_STEP, row["verdict"],
        f"other sectors: {row['other_sector_results'] or 'none'}; joint search SNR {row['joint_search_snr']} "
        f"against highest null SNR {row['joint_search_highest_null_snr']}; Gaia RUWE {row['gaia_ruwe']}; "
        f"companion radius {row['companion_radius_earth_radii']} Earth radii; "
        f"refuted by: {row['refuted_by'] or 'nothing'}; cautions: {row['cautions'] or 'none'}",
        row,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    parser.add_argument("--restart", action="store_true",
                        help="check all candidates again, not only those without a row in candidate_checks.csv")
    parser.add_argument("--tic-id", dest="tic_ids", action="append", default=[],
                        help="check only the candidates of this TIC ID; give it more than once for more TIC IDs")
    arguments = parser.parse_args()
    sector = arguments.sector

    directory = output_directory(sector)
    directory.mkdir(parents=True, exist_ok=True)
    output_path = directory / "candidate_checks.csv"
    all_candidates = list(csv.DictReader(open(cross_match_output_directory(sector) / "new_candidates.csv")))
    earlier_rows = earlier_rows_to_keep(
        output_path, arguments.restart, {(row["tic_id"], row["period_days"]) for row in all_candidates},
        tic_ids=arguments.tic_ids)
    checked = {(row["tic_id"], row["period_days"]) for row in earlier_rows}
    candidates = [candidate for candidate in all_candidates
                  if (candidate["tic_id"], candidate["period_days"]) not in checked
                  and (not arguments.tic_ids or candidate["tic_id"] in arguments.tic_ids)]
    print(f"Sector {sector}: {len(all_candidates)} new candidates, {len(earlier_rows)} checked in an earlier run",
          flush=True)
    ephemeris_sources, stars = load_ephemeris_sources(sector)
    gaia_dr2_source_ids = gaia_dr2_source_ids_by_tic_id({candidate["tic_id"] for candidate in candidates})

    # A stop between a CSV row and its record row leaves the record row missing. This appends it.
    for row in earlier_rows:
        append_verification_record(sector, row)

    # The candidates run one at a time: lightkurve redirects stdout during a download, which breaks other threads.
    rows = list(earlier_rows)
    status = PipelineStatus(sector, 7, "independent_checks", "check_qlp_sector_new_candidates.py", len(candidates))
    with IncrementalCsvWriter(output_path, OUTPUT_COLUMNS, earlier_rows) as writer, status:
        for count, candidate in enumerate(candidates, start=1):
            status.start_candidate(candidate["tic_id"], candidate["period_days"])
            row = check_candidate(candidate, sector, ephemeris_sources, stars,
                                  gaia_dr2_source_ids.get(int(candidate["tic_id"])))
            row.update(verification_stamp())
            row = writer.write_row(row)
            append_verification_record(sector, row)
            status.finish_candidate()
            rows.append(row)
            print(f"  {count} of {len(candidates)}: TIC {row['tic_id']} {row['verdict']}"
                  f"{' (' + row['refuted_by'] + ')' if row['refuted_by'] else ''}", flush=True)

    counts = {name: sum(row["verdict"] == name for row in rows) for name in ("supported", "unconfirmed", "refuted")}
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}: {counts}")


if __name__ == "__main__":
    main()
