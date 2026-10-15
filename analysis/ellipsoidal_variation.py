"""Measure the ellipsoidal variation of a star at the period of a transiting companion.

The tides of a companion of mass M2 stretch the star, which gives a flux variation at half the orbital period, with
its minima at the conjunctions: A = ELLIPSOIDAL_COEFFICIENT (M2 / M1) (R1 / a)^3 sin^2 i. A planet gives a few ppm
at most. A low-mass star in a close orbit gives hundreds of ppm, so a detected variation is a warning of a stellar
companion.

The fit takes one sector of flux that is not detrended. It removes the transits and the phase-0.5 window, and fits
one polynomial of degree ELLIPSOIDAL_TREND_DEGREE per continuous segment, plus -cos(4 pi phase) (the ellipsoidal
variation), -cos(2 pi phase) (reflection) and sin(2 pi phase) (Doppler beaming). The error is the standard deviation
of the same fit at ELLIPSOIDAL_NULL_TRIALS wrong periods, which includes correlated noise.

calibrate_ellipsoidal_variation_test_on_known_planets.py runs this fit on known planets of hot stars.
"""

import numpy as np

ELLIPSOIDAL_COEFFICIENT = 1.1
ELLIPSOIDAL_TREND_DEGREE = 2
ELLIPSOIDAL_NULL_TRIALS = 200
# The wrong periods are random factors in this range times the period, each with a random epoch.
ELLIPSOIDAL_NULL_PERIOD_FACTORS = (0.6, 1.6)
# A gap longer than this starts a new continuous segment, with its own polynomial.
ELLIPSOIDAL_SEGMENT_GAP_DAYS = 1.0
ELLIPSOIDAL_MINIMUM_POINTS = 1000


def ellipsoidal_design_matrix(time, period, epoch):
    phase = ((time - epoch) / period) % 1
    columns = [-np.cos(4 * np.pi * phase), -np.cos(2 * np.pi * phase), np.sin(2 * np.pi * phase)]
    breaks = np.concatenate([[0], np.where(np.diff(time) > ELLIPSOIDAL_SEGMENT_GAP_DAYS)[0] + 1, [len(time)]])
    for start, end in zip(breaks[:-1], breaks[1:]):
        middle = np.mean(time[start:end])
        for power in range(ELLIPSOIDAL_TREND_DEGREE + 1):
            column = np.zeros(len(time))
            column[start:end] = (time[start:end] - middle) ** power
            columns.append(column)
    return np.array(columns).T


def fitted_ellipsoidal_amplitude(time, flux, period, epoch):
    return np.linalg.lstsq(ellipsoidal_design_matrix(time, period, epoch), flux - 1, rcond=None)[0][0]


def ellipsoidal_amplitude(time, flux, period, epoch, duration, rng):
    """Return the ellipsoidal amplitude and its error, as fractions of the flux, of one sector of normalised flux.

    Return None when fewer than ELLIPSOIDAL_MINIMUM_POINTS cadences are out of the transits and the phase-0.5 window.
    """
    transit_time = (time - epoch + 0.5 * period) % period - 0.5 * period
    secondary_time = (time - epoch) % period - 0.5 * period
    keep = (np.abs(transit_time) > duration) & (np.abs(secondary_time) > duration)
    time, flux = time[keep], flux[keep]
    if len(time) < ELLIPSOIDAL_MINIMUM_POINTS:
        return None
    amplitude = fitted_ellipsoidal_amplitude(time, flux, period, epoch)
    null = [fitted_ellipsoidal_amplitude(time, flux, period * factor, epoch + rng.uniform(0, period))
            for factor in rng.uniform(*ELLIPSOIDAL_NULL_PERIOD_FACTORS, ELLIPSOIDAL_NULL_TRIALS)]
    return float(amplitude), float(np.std(null))


def combined_ellipsoidal_amplitude(amplitudes_and_errors):
    """Return the weighted mean of the amplitudes of several sectors, and its error."""
    amplitudes, errors = np.array(amplitudes_and_errors).T
    weights = errors ** -2
    return float(np.sum(amplitudes * weights) / np.sum(weights)), float(np.sum(weights) ** -0.5)


def companion_mass_from_ellipsoidal_amplitude(amplitude, stellar_mass, scaled_semi_major_axis, impact_parameter):
    """Return the companion mass, in the unit of stellar_mass, that gives this ellipsoidal amplitude."""
    sin_squared_inclination = 1 - (impact_parameter / scaled_semi_major_axis) ** 2
    return amplitude * stellar_mass * scaled_semi_major_axis ** 3 / (ELLIPSOIDAL_COEFFICIENT * sin_squared_inclination)
