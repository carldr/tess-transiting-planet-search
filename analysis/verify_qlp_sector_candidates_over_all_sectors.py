"""Verify each candidate of one sector that the pixel check puts on the target, with the data of all TESS sectors.

The script takes each candidate with pixel verdict "on target" or "unresolved", and:
1. loads the light curve of the sector and of each other sector with a SPOC, TESS-SPOC or QLP light curve. It
   detrends each light curve with the transits masked, so the detrending does not make the transits shallower;
2. fits one transit model to all sectors with batman and emcee. Each sector model is integrated over the exposure
   time of that sector. The fit gives the period, the epoch, the depth, Rp/R*, the impact parameter, a/R* and the
   duration, each with its 1-sigma error. The planet radius uses the TIC stellar radius and its error;
3. compares the stellar density from a/R* with the TIC stellar density. An eclipsing binary on another star often
   gives a transit shape that does not fit the target;
4. measures the depth of the odd and the even transits with the fitted shape. Different depths mean an eclipsing
   binary at twice the period;
5. measures the depth of each sector with the fitted shape. Different depths can mean a contaminating star that the
   apertures of the sectors include differently;
6. measures the depth at phase 0.5. It compares it with the largest secondary eclipse a planet can give: full
   reflection plus the day side thermal emission without heat redistribution, at 800 nm;
7. reads the Gaia DR3 radial velocity statistics. The star is RV-variable when rv_chisq_pvalue <= 0.01,
   rv_renormalised_gof > 4 and rv_nb_transits >= 10 (Katz et al. 2023). A planet cannot cause a variation that
   Gaia can detect, so this means a stellar companion or a pulsating star;
8. lists the LAMOST DR7 low-resolution spectra within LAMOST_RADIUS_ARCSEC, with their class and radial velocity.

The verdict is "refuted" when the odd and even depths differ by ODD_EVEN_SIGMA_LIMIT sigma or more, when the
secondary eclipse is deeper than a planet can give by SECONDARY_SIGMA_LIMIT sigma, or when the companion radius is
above COMPANION_MAXIMUM_RADIUS_EARTH_RADII. Otherwise it is "passed". The cautions column lists weaker warnings.

Run from the repository root, after check_qlp_sector_candidate_pixels.py:
    .venv/bin/python analysis/verify_qlp_sector_candidates_over_all_sectors.py 98

Output, in results/qlp_sector_<nnnn>_all_sector_verification/:
- all_sector_verification.csv: one row per candidate;
- tic<TIC ID>_<period>.png: the folded transit of each sector with the model, and the folded light curve at phase 0.5.
The script also appends the result of each candidate to its verification record.
"""

import argparse
import csv
import warnings

import astropy.units as units
import batman
import emcee
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import G, R_earth, R_sun, M_sun, c, h, k_B
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astroquery.vizier import Vizier
from scipy import stats
from scipy.optimize import minimize
from wotan import flatten

from build_qlp_sector_star_list import REPOSITORY_ROOT, sector_directory
from candidate_verification_records import append_verification_record_entry
from check_qlp_sector_candidate_pixels import output_directory as pixel_checks_output_directory
from check_qlp_sector_new_candidates import (
    COMPANION_MAXIMUM_RADIUS_EARTH_RADII, SCALED_SEMI_MAJOR_AXIS_RATIO_LIMIT, load_other_sector_light_curve,
    other_sector_products, output_directory as candidate_checks_output_directory,
)
from download_qlp_sector_light_curves import light_curve_path
from search_qlp_sector import (
    BIWEIGHT_WINDOW_LENGTH_DAYS, load_light_curve, robust_standard_deviation,
)

PIXEL_VERDICTS_TO_VERIFY = ("on target", "unresolved")
# The detrending masks the cadences within this many transit durations of mid-transit, and uses a window of at
# least this many transit durations.
DETRENDING_MASK_DURATIONS = 1.0
DETRENDING_WINDOW_DURATIONS = 3.0
# The fit uses the cadences within this many transit durations of mid-transit.
FIT_WINDOW_DURATIONS = 2.5
# The red noise factor compares the scatter of bin means with white noise, for bins of these fractions of the
# transit duration.
RED_NOISE_BIN_FRACTIONS = (0.25, 0.5, 0.75, 1.0)
# The exposure time of the full-frame images, by the first sector of each exposure time.
FULL_FRAME_IMAGE_EXPOSURE_SECONDS_BY_FIRST_SECTOR = {1: 1800, 27: 600, 56: 200}
# batman integrates the model over the exposure in steps of at most this many seconds.
SUPERSAMPLE_STEP_SECONDS = 120
MCMC_WALKERS = 40
MCMC_STEPS = 4000
MCMC_BURN_IN_STEPS = 1500
# The depth at mid-transit comes from the model of this many random samples of the chain.
DEPTH_SAMPLES = 2000
ODD_EVEN_SIGMA_LIMIT = 3.0
SECONDARY_SIGMA_LIMIT = 3.0
SECTOR_DEPTH_CONSISTENCY_P_VALUE = 0.001
SECONDARY_ECLIPSE_WAVELENGTH_METRES = 800e-9
GAIA_RV_VARIABLE_P_VALUE = 0.01
GAIA_RV_VARIABLE_RENORMALISED_GOF = 4.0
GAIA_RV_MINIMUM_TRANSITS = 10
LAMOST_RADIUS_ARCSEC = 3.0

OUTPUT_COLUMNS = [
    "tic_id", "period_days", "depth_ppm", "tmag", "pixel_verdict", "sectors_used", "exposure_seconds_by_sector",
    "fit_period_days", "fit_period_error_days", "fit_epoch_btjd", "fit_epoch_error_days",
    "fit_depth_ppm", "fit_depth_error_ppm", "fit_radius_ratio", "fit_radius_ratio_error",
    "fit_impact_parameter", "fit_impact_parameter_error", "fit_scaled_semi_major_axis",
    "fit_scaled_semi_major_axis_error", "fit_duration_hours", "fit_duration_error_hours",
    "planet_radius_earth_radii", "planet_radius_error_earth_radii", "stellar_radius_solar_radii",
    "stellar_density_fit_to_tic_ratio", "odd_depth_ppm", "even_depth_ppm", "odd_even_difference_sigma",
    "depth_ppm_by_sector", "sector_depth_consistency_p_value", "secondary_depth_ppm", "secondary_depth_error_ppm",
    "secondary_maximum_planet_depth_ppm", "gaia_source_id", "gaia_rv_nb_transits", "gaia_rv_chisq_pvalue",
    "gaia_rv_renormalised_gof", "gaia_rv_amplitude_robust_km_s", "gaia_vbroad_km_s", "gaia_rv_variable",
    "lamost_spectra", "refuted_by", "cautions", "verdict", "error",
]


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_all_sector_verification"


def full_frame_image_exposure_seconds(sector):
    return [seconds for first, seconds in FULL_FRAME_IMAGE_EXPOSURE_SECONDS_BY_FIRST_SECTOR.items() if sector >= first][-1]


def phase_time(time, period, epoch):
    """Return the time from the nearest mid-transit."""
    return (time - epoch + 0.5 * period) % period - 0.5 * period


def detrend_with_transits_masked(time, flux, period, epoch, duration):
    """Divide the flux by a biweight trend fitted to the cadences out of transit."""
    out_of_transit = np.abs(phase_time(time, period, epoch)) > DETRENDING_MASK_DURATIONS * duration
    window = max(BIWEIGHT_WINDOW_LENGTH_DAYS, DETRENDING_WINDOW_DURATIONS * duration)
    _, trend = flatten(time[out_of_transit], flux[out_of_transit], method="biweight", window_length=window,
                       edge_cutoff=0.0, return_trend=True)
    good = np.isfinite(trend)
    trend = np.interp(time, time[out_of_transit][good], trend[good])
    detrended = flux / trend
    keep = np.isfinite(detrended)
    time, detrended = time[keep], detrended[keep]
    keep = detrended < 1 + 5 * robust_standard_deviation(detrended)
    return time[keep], detrended[keep]


def red_noise_factor(time, flux, duration):
    """Return the factor by which correlated noise raises the error of a mean over a transit (Winn et al. 2008).

    The factor is the mean ratio of the measured scatter of bin means to the white-noise scatter, for bins of
    RED_NOISE_BIN_FRACTIONS of the transit duration, and at least 1.
    """
    white = np.std(flux)
    ratios = []
    for fraction in RED_NOISE_BIN_FRACTIONS:
        bin_numbers = np.floor((time - time.min()) / (fraction * duration)).astype(int)
        counts = np.bincount(bin_numbers)
        full = counts >= 0.5 * np.median(counts[counts > 0])
        if full.sum() < 10:
            continue
        means = np.bincount(bin_numbers, flux)[full] / counts[full]
        ratios.append(np.std(means) / (white / np.sqrt(np.mean(counts[full]))))
    return max(1.0, float(np.mean(ratios))) if ratios else 1.0


def load_all_sectors(sector, candidate, period, epoch, duration):
    """Return [(sector, author, exposure seconds, time, detrended flux)] for this sector and each other sector."""
    tic_id = int(candidate["tic_id"])
    sectors = [(sector, "QLP", full_frame_image_exposure_seconds(sector),
                *load_light_curve(light_curve_path(sector, tic_id)))]
    for product_sector, author, product in other_sector_products(tic_id, sector):
        sectors.append((product_sector, author, float(product.exptime[0].to_value(units.s)),
                        *load_other_sector_light_curve(product, author)))
    return [(number, author, exposure, *detrend_with_transits_masked(time, flux, period, epoch, duration))
            for number, author, exposure, time, flux in sectors]


class TransitModel:
    """A batman transit model of all sectors, with the parameters of TransitModel.PARAMETER_NAMES."""

    PARAMETER_NAMES = ["epoch", "period", "radius_ratio", "impact_parameter", "log10_scaled_semi_major_axis",
                       "limb_darkening_q1", "limb_darkening_q2"]

    def __init__(self, sectors):
        self.sectors = sectors
        self.models = []
        for _, _, exposure, time, _ in sectors:
            parameters = self.batman_parameters([0, 1, 0.1, 0.3, 1, 0.4, 0.4])
            supersample = max(1, int(np.ceil(exposure / SUPERSAMPLE_STEP_SECONDS)))
            self.models.append(batman.TransitModel(parameters, time, supersample_factor=supersample,
                                                   exp_time=exposure / 86400))

    @staticmethod
    def batman_parameters(theta):
        epoch, period, radius_ratio, impact_parameter, log10_a, q1, q2 = theta
        parameters = batman.TransitParams()
        parameters.t0, parameters.per, parameters.rp = epoch, period, radius_ratio
        parameters.a = 10 ** log10_a
        parameters.inc = np.degrees(np.arccos(min(impact_parameter / parameters.a, 1.0)))
        parameters.ecc, parameters.w = 0.0, 90.0
        # Kipping (2013): q1 and q2 in [0, 1] give every physical quadratic limb darkening law.
        parameters.u = [2 * np.sqrt(q1) * q2, np.sqrt(q1) * (1 - 2 * q2)]
        parameters.limb_dark = "quadratic"
        return parameters

    def light_curves(self, theta):
        parameters = self.batman_parameters(theta)
        return [model.light_curve(parameters) for model in self.models]


def log_probability(theta, model, errors, period_guess, epoch_guess):
    epoch, period, radius_ratio, impact_parameter, log10_a, q1, q2 = theta
    if not (abs(period - period_guess) < 0.01 * period_guess and abs(epoch - epoch_guess) < 0.2
            and 0 < radius_ratio < 0.5 and 0 <= impact_parameter < 1 + radius_ratio and 0.2 < log10_a < 2.5
            and 0 < q1 < 1 and 0 < q2 < 1):
        return -np.inf
    return -0.5 * sum(np.sum(((flux - curve) / error) ** 2) for (_, _, _, _, flux), curve, error
                      in zip(model.sectors, model.light_curves(theta), errors))


def observed_depths(samples):
    """Return the depth at mid-transit of the model of each sample, without integration over the exposure."""
    depths = []
    for theta in samples:
        parameters = TransitModel.batman_parameters(theta)
        depths.append(1 - batman.TransitModel(parameters, np.array([theta[0]])).light_curve(parameters)[0])
    return np.array(depths)


def transit_duration_days(period, radius_ratio, impact_parameter, scaled_semi_major_axis):
    inclination = np.arccos(np.clip(impact_parameter / scaled_semi_major_axis, 0, 1))
    argument = np.sqrt(np.clip((1 + radius_ratio) ** 2 - impact_parameter ** 2, 0, None)) / (
        scaled_semi_major_axis * np.sin(inclination))
    return period / np.pi * np.arcsin(np.clip(argument, 0, 1))


def scaled_depths(groups, sectors, curves, errors):
    """Fit flux - 1 = scale * (curve - 1) for each group of cadences. Return {group: (depth scale, its error)}."""
    numerators, denominators = {}, {}
    for (_, _, _, _, flux), curve, error, group in zip(sectors, curves, errors, groups):
        for name in np.unique(group):
            chosen = group == name
            shape = curve[chosen] - 1
            numerators[name] = numerators.get(name, 0) + np.sum((flux[chosen] - 1) * shape) / error ** 2
            denominators[name] = denominators.get(name, 0) + np.sum(shape ** 2) / error ** 2
    return {name: (numerators[name] / denominators[name], 1 / np.sqrt(denominators[name]))
            for name in numerators if denominators[name] > 0}


def secondary_eclipse_depth(sectors, errors, period, epoch, duration):
    """Return the depth and its error of a box of the transit duration at phase 0.5."""
    inside, outside = [], []
    for (_, _, _, time, flux), error in zip(sectors, errors):
        offset = np.abs(phase_time(time, period, epoch + 0.5 * period))
        inside.append((flux[offset < 0.5 * duration], error))
        outside.append((flux[(offset > duration) & (offset < 2.5 * duration)], error))

    def weighted_mean(parts):
        weights = np.concatenate([np.full(len(flux), error ** -2) for flux, error in parts])
        values = np.concatenate([flux for flux, _ in parts])
        return np.sum(weights * values) / np.sum(weights), 1 / np.sqrt(np.sum(weights))

    (inside_mean, inside_error), (outside_mean, outside_error) = weighted_mean(inside), weighted_mean(outside)
    return outside_mean - inside_mean, np.hypot(inside_error, outside_error)


def planck(wavelength, temperature):
    return 1 / (np.exp((h * c / (wavelength * units.m * k_B * temperature * units.K)).decompose().value) - 1)


def maximum_planet_secondary_depth(radius_ratio, scaled_semi_major_axis, stellar_temperature):
    """Return the depth of full reflection plus the day side thermal emission without heat redistribution."""
    reflection = (radius_ratio / scaled_semi_major_axis) ** 2
    day_side_temperature = stellar_temperature * np.sqrt(1 / scaled_semi_major_axis) * (2 / 3) ** 0.25
    thermal = radius_ratio ** 2 * planck(SECONDARY_ECLIPSE_WAVELENGTH_METRES, day_side_temperature) / planck(
        SECONDARY_ECLIPSE_WAVELENGTH_METRES, stellar_temperature)
    return reflection + thermal


def tic_star(tic_id):
    row = Catalogs.query_criteria(catalog="TIC", ID=int(tic_id))[0]
    return {name: float(row[name]) if row[name] is not np.ma.masked else np.nan
            for name in ("rad", "e_rad", "mass", "Teff")}


def gaia_radial_velocity_statistics(gaia_source_id):
    if not gaia_source_id:
        return {}
    row = Gaia.launch_job(
        "SELECT rv_nb_transits, rv_chisq_pvalue, rv_renormalised_gof, rv_amplitude_robust, vbroad "
        f"FROM gaiadr3.gaia_source WHERE source_id = {int(gaia_source_id)}"
    ).get_results()[0]

    def value(name):
        return "" if row[name] is np.ma.masked else float(row[name])

    statistics = {
        "gaia_rv_nb_transits": value("rv_nb_transits"), "gaia_rv_chisq_pvalue": value("rv_chisq_pvalue"),
        "gaia_rv_renormalised_gof": value("rv_renormalised_gof"),
        "gaia_rv_amplitude_robust_km_s": value("rv_amplitude_robust"), "gaia_vbroad_km_s": value("vbroad"),
    }
    statistics["gaia_rv_variable"] = (
        all(statistics[name] != "" for name in ("gaia_rv_nb_transits", "gaia_rv_chisq_pvalue", "gaia_rv_renormalised_gof"))
        and statistics["gaia_rv_chisq_pvalue"] <= GAIA_RV_VARIABLE_P_VALUE
        and statistics["gaia_rv_renormalised_gof"] > GAIA_RV_VARIABLE_RENORMALISED_GOF
        and statistics["gaia_rv_nb_transits"] >= GAIA_RV_MINIMUM_TRANSITS
    )
    return statistics


def lamost_spectra(ra, dec):
    tables = Vizier(columns=["**", "_r"], row_limit=-1).query_region(
        SkyCoord(ra, dec, unit="deg"), radius=LAMOST_RADIUS_ARCSEC * units.arcsec, catalog="V/156/dr7lrs")
    if not tables:
        return ""
    return "; ".join(
        f"{row['Obs.Date']} {row['Class']} {row['subClass']}, RV {row['z'] * c.to_value('km/s'):.0f} "
        f"+/- {row['e_z'] * c.to_value('km/s'):.1f} km/s, SNR g {row['snrg']:.0f}"
        for row in tables[0]
    )


def verify(sector, candidate, star_row, gaia_source_id, figure_directory):
    period_guess = float(candidate["joint_search_period_days"] or candidate["period_days"])
    duration = float(candidate["duration_hours"]) / 24
    all_sectors = load_all_sectors(sector, candidate, period_guess, float(candidate["epoch_btjd"]), duration)
    # The reference epoch is the transit nearest the middle of the data, so the epoch and the period do not
    # correlate.
    all_times = np.concatenate([time for _, _, _, time, _ in all_sectors])
    middle_transit = np.round((np.median(all_times) - float(candidate["epoch_btjd"])) / period_guess)
    epoch_guess = float(candidate["epoch_btjd"]) + middle_transit * period_guess
    all_sectors = [sector_data for sector_data in all_sectors
                   if np.any(np.abs(phase_time(sector_data[3], period_guess, epoch_guess)) < FIT_WINDOW_DURATIONS * duration)]
    all_sector_errors = [
        robust_standard_deviation(out_of_transit_flux) * red_noise_factor(time[out_of_transit], out_of_transit_flux, duration)
        for _, _, _, time, flux in all_sectors
        for out_of_transit in [np.abs(phase_time(time, period_guess, epoch_guess)) > duration]
        for out_of_transit_flux in [flux[out_of_transit]]
    ]
    sectors = [(number, author, exposure, time[near], flux[near])
               for number, author, exposure, time, flux in all_sectors
               for near in [np.abs(phase_time(time, period_guess, epoch_guess)) < FIT_WINDOW_DURATIONS * duration]]
    errors = all_sector_errors
    star = tic_star(candidate["tic_id"])
    density = star["mass"] * M_sun / (4 / 3 * np.pi * (star["rad"] * R_sun) ** 3)
    a_guess = ((G * density * (period_guess * units.day) ** 2 / (3 * np.pi)) ** (1 / 3)).decompose().value
    model = TransitModel(sectors)
    start = [epoch_guess, period_guess, np.sqrt(float(candidate["depth_ppm"]) * 1e-6), 0.3,
             np.log10(a_guess) if np.isfinite(a_guess) else 1.0, 0.4, 0.4]
    arguments = (model, errors, period_guess, epoch_guess)
    best = minimize(lambda theta: -log_probability(theta, *arguments), start, method="Nelder-Mead",
                    options={"maxiter": 20000, "maxfev": 20000, "xatol": 1e-7, "fatol": 1e-3}).x
    rng = np.random.default_rng(int(candidate["tic_id"]))
    scales = np.array([1e-4, 1e-6, 1e-3, 1e-2, 1e-3, 1e-2, 1e-2])
    walkers = best + scales * rng.standard_normal((MCMC_WALKERS, len(best)))
    walkers[:, 3] = np.abs(walkers[:, 3])
    walkers[:, 5:] = np.clip(walkers[:, 5:], 0.01, 0.99)
    sampler = emcee.EnsembleSampler(MCMC_WALKERS, len(best), log_probability, args=arguments)
    sampler.run_mcmc(walkers, MCMC_STEPS)
    chain = sampler.get_chain(discard=MCMC_BURN_IN_STEPS, flat=True)
    epoch, period, radius_ratio, impact_parameter = chain[:, 0], chain[:, 1], chain[:, 2], chain[:, 3]
    scaled_semi_major_axis = 10 ** chain[:, 4]
    durations = transit_duration_days(period, radius_ratio, impact_parameter, scaled_semi_major_axis)
    depth = observed_depths(chain[rng.choice(len(chain), DEPTH_SAMPLES, replace=False)])
    stellar_radius = star["rad"] + (star["e_rad"] if np.isfinite(star["e_rad"]) else 0) * rng.standard_normal(len(chain))
    planet_radius = radius_ratio * stellar_radius * (R_sun / R_earth).decompose().value
    fit_density = 3 * np.pi * scaled_semi_major_axis ** 3 / (G * (period * units.day) ** 2)

    def median_and_error(samples):
        low, middle, high = np.percentile(samples, [15.87, 50, 84.13])
        return middle, 0.5 * (high - low)

    median = np.median(chain, axis=0)
    curves = model.light_curves(median)
    transit_numbers = [np.round((time - median[0]) / median[1]).astype(int) for _, _, _, time, _ in sectors]
    odd_even = scaled_depths([np.where(numbers % 2 == 0, "even", "odd") for numbers in transit_numbers],
                             sectors, curves, errors)
    by_sector = scaled_depths([np.full(len(time), str(number)) for number, _, _, time, _ in sectors],
                              sectors, curves, errors)
    model_depth = float(np.median(depth))
    if {"odd", "even"} <= set(odd_even):
        (odd, odd_error), (even, even_error) = odd_even["odd"], odd_even["even"]
        odd_even_sigma = abs(odd - even) / np.hypot(odd_error, even_error)
    else:
        odd = even = odd_even_sigma = np.nan
    sector_scales = np.array([value for value, _ in by_sector.values()])
    sector_errors = np.array([error for _, error in by_sector.values()])
    weighted_mean = np.sum(sector_scales / sector_errors ** 2) / np.sum(sector_errors ** -2)
    consistency_chi_square = np.sum(((sector_scales - weighted_mean) / sector_errors) ** 2)
    consistency_p_value = float(stats.chi2.sf(consistency_chi_square, len(by_sector) - 1)) if len(by_sector) > 1 else np.nan
    duration_median = float(np.median(durations))
    secondary, secondary_error = secondary_eclipse_depth(all_sectors, all_sector_errors, median[1], median[0],
                                                         duration_median)
    secondary_maximum = maximum_planet_secondary_depth(median[2], 10 ** median[4], star["Teff"])
    density_ratio = float(np.median(fit_density) / density) if np.isfinite(density.value) else np.nan
    radius_median, radius_error = median_and_error(planet_radius)

    row = {
        "sectors_used": " ".join(str(number) for number, *_ in sectors),
        "exposure_seconds_by_sector": " ".join(f"S{number}:{exposure:.0f}" for number, _, exposure, *_ in sectors),
        "stellar_radius_solar_radii": star["rad"],
        "stellar_density_fit_to_tic_ratio": f"{density_ratio:.3g}",
        "odd_depth_ppm": f"{odd * model_depth * 1e6:.0f}", "even_depth_ppm": f"{even * model_depth * 1e6:.0f}",
        "odd_even_difference_sigma": f"{odd_even_sigma:.2f}",
        "depth_ppm_by_sector": "; ".join(f"S{name}: {value * model_depth * 1e6:.0f} +/- {error * model_depth * 1e6:.0f}"
                                         for name, (value, error) in by_sector.items()),
        "sector_depth_consistency_p_value": f"{consistency_p_value:.3g}",
        "secondary_depth_ppm": f"{secondary * 1e6:.0f}", "secondary_depth_error_ppm": f"{secondary_error * 1e6:.0f}",
        "secondary_maximum_planet_depth_ppm": f"{secondary_maximum * 1e6:.0f}",
        "planet_radius_earth_radii": f"{radius_median:.2f}", "planet_radius_error_earth_radii": f"{radius_error:.2f}",
    }
    for name, samples, factor, digits in (
        ("fit_period", period, 1, 7), ("fit_epoch", epoch, 1, 5), ("fit_depth", depth, 1e6, 0),
        ("fit_radius_ratio", radius_ratio, 1, 4), ("fit_impact_parameter", impact_parameter, 1, 3),
        ("fit_scaled_semi_major_axis", scaled_semi_major_axis, 1, 2), ("fit_duration", durations, 24, 3),
    ):
        value, error = median_and_error(samples * factor)
        unit = {"fit_period": "_days", "fit_epoch": "_btjd", "fit_depth": "_ppm", "fit_duration": "_hours"}.get(name, "")
        error_unit = {"fit_period": "_days", "fit_epoch": "_days", "fit_depth": "_ppm", "fit_duration": "_hours"}.get(name, "")
        row[f"{name}{unit}"] = f"{value:.{digits}f}"
        row[f"{name}_error{error_unit}"] = f"{error:.{max(digits, 2) + 1}f}" if digits else f"{error:.0f}"

    row.update(gaia_radial_velocity_statistics(gaia_source_id))
    row["lamost_spectra"] = lamost_spectra(float(star_row["ra"]), float(star_row["dec"]))

    refuted_by, cautions = [], []
    if odd_even_sigma >= ODD_EVEN_SIGMA_LIMIT:
        refuted_by.append("odd and even depths differ")
    if secondary - SECONDARY_SIGMA_LIMIT * secondary_error > secondary_maximum:
        refuted_by.append("secondary eclipse too deep for a planet")
    elif secondary >= SECONDARY_SIGMA_LIMIT * secondary_error:
        cautions.append("secondary eclipse detected")
    if radius_median > COMPANION_MAXIMUM_RADIUS_EARTH_RADII:
        refuted_by.append("companion too large for a planet")
    if np.median(impact_parameter + radius_ratio) > 1:
        cautions.append("grazing transit")
    if not SCALED_SEMI_MAJOR_AXIS_RATIO_LIMIT ** -3 <= density_ratio <= SCALED_SEMI_MAJOR_AXIS_RATIO_LIMIT ** 3:
        cautions.append("stellar density does not fit the star")
    if consistency_p_value < SECTOR_DEPTH_CONSISTENCY_P_VALUE:
        cautions.append("depths differ between sectors")
    if row.get("gaia_rv_variable") is True:
        cautions.append("Gaia radial velocity variable")
    row["refuted_by"] = "; ".join(refuted_by)
    row["cautions"] = "; ".join(cautions)
    row["verdict"] = "refuted" if refuted_by else "passed"

    draw_figure(figure_directory / f"tic{candidate['tic_id']}_{float(candidate['period_days']):.5f}.png",
                sectors, all_sectors, model, median, duration_median, candidate, row["verdict"])
    return row


def draw_figure(path, sectors, all_sectors, model, theta, duration, candidate, verdict):
    figure, axes = plt.subplots(1, len(sectors) + 1, figsize=(4 * (len(sectors) + 1), 4), squeeze=False)
    for axis, (number, author, exposure, time, flux), curve in zip(axes[0], sectors, model.light_curves(theta)):
        hours = phase_time(time, theta[1], theta[0]) * 24
        order = np.argsort(hours)
        hours, flux_ppm, curve_ppm = hours[order], (flux[order] - 1) * 1e6, (curve[order] - 1) * 1e6
        axis.plot(hours, flux_ppm, ".", color="0.6", markersize=2)
        bin_index = np.digitize(hours, np.linspace(hours.min(), hours.max(), 40))
        binned = [(hours[bin_index == value].mean(), flux_ppm[bin_index == value].mean())
                  for value in np.unique(bin_index) if (bin_index == value).sum() > 2]
        axis.plot(*np.array(binned).T, "o", color="black", markersize=4)
        axis.plot(hours, curve_ppm, color="red")
        axis.set_title(f"Sector {number} {author}, {exposure:.0f} s")
        axis.set_xlabel("hours from mid-transit")
    axes[0][0].set_ylabel("flux - 1 (ppm)")
    axis = axes[0][-1]
    axis.set_title("Phase 0.5 (secondary eclipse)")
    axis.set_xlabel("hours from phase 0.5")
    for number, author, exposure, time, flux in all_sectors:
        hours = phase_time(time, theta[1], theta[0] + 0.5 * theta[1]) * 24
        near = np.abs(hours) < FIT_WINDOW_DURATIONS * duration * 24
        axis.plot(hours[near], (flux[near] - 1) * 1e6, ".", markersize=2, label=f"S{number}")
    axis.axvspan(-12 * duration, 12 * duration, color="orange", alpha=0.2)
    axis.set_ylim(axes[0][0].get_ylim())
    axis.legend()
    figure.suptitle(f"TIC {candidate['tic_id']}, P = {theta[1]:.6f} d: {verdict}")
    figure.tight_layout()
    figure.savefig(path, dpi=90)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    sector = parser.parse_args().sector

    pixel_checks = list(csv.DictReader(open(pixel_checks_output_directory(sector) / "pixel_checks.csv")))
    checks = {(row["tic_id"], row["period_days"]): row
              for row in csv.DictReader(open(candidate_checks_output_directory(sector) / "candidate_checks.csv"))}
    star_rows = {row["tic_id"]: row for row in csv.DictReader(open(sector_directory(sector) / "star_list.csv"))}
    candidates = [row for row in pixel_checks if row["pixel_verdict"] in PIXEL_VERDICTS_TO_VERIFY]
    directory = output_directory(sector)
    directory.mkdir(parents=True, exist_ok=True)
    print(f"Sector {sector}: {len(candidates)} candidates on the target or unresolved", flush=True)

    rows = []
    for count, pixel_check in enumerate(candidates, start=1):
        candidate = {**checks[(pixel_check["tic_id"], pixel_check["period_days"])], **pixel_check}
        row = {column: candidate.get(column, "") for column in OUTPUT_COLUMNS}
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                row.update(verify(sector, candidate, star_rows[candidate["tic_id"]], candidate["gaia_source_id"],
                                  directory))
        except Exception as error:
            row["error"] = f"{type(error).__name__}: {error}"
        rows.append(row)
        result = row["error"] or (
            f"sectors {row['sectors_used']}; P {row['fit_period_days']} +/- {row['fit_period_error_days']} d; "
            f"T0 {row['fit_epoch_btjd']} +/- {row['fit_epoch_error_days']} BTJD; depth {row['fit_depth_ppm']} "
            f"+/- {row['fit_depth_error_ppm']} ppm; duration {row['fit_duration_hours']} h; b {row['fit_impact_parameter']}; "
            f"Rp {row['planet_radius_earth_radii']} +/- {row['planet_radius_error_earth_radii']} Earth radii; "
            f"density ratio {row['stellar_density_fit_to_tic_ratio']}; odd-even {row['odd_even_difference_sigma']} sigma; "
            f"sector depths {row['depth_ppm_by_sector']} (p {row['sector_depth_consistency_p_value']}); "
            f"secondary {row['secondary_depth_ppm']} +/- {row['secondary_depth_error_ppm']} ppm "
            f"(planet maximum {row['secondary_maximum_planet_depth_ppm']}); Gaia RV variable {row.get('gaia_rv_variable', '')}; "
            f"LAMOST: {row['lamost_spectra'] or 'none'}; refuted by: {row['refuted_by'] or 'nothing'}; "
            f"cautions: {row['cautions'] or 'none'}"
        )
        append_verification_record_entry(sector, candidate["tic_id"], candidate["period_days"],
                                         "all-sector verification", row["verdict"] or "error", result)
        print(f"  {count} of {len(candidates)}: TIC {row['tic_id']} {row['verdict'] or row['error']}: {result}",
              flush=True)

    output_path = directory / "all_sector_verification.csv"
    with open(output_path, "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
