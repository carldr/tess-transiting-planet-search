"""Load and detrend the archival TESS light curves of TIC 4206066.

Each sector uses the product that Rabtsevich used:
- Sector 6: TESS-SPOC full-frame-image light curve, 30 min, PDCSAP flux.
- Sector 32: SPOC target light curve, 2 min, PDCSAP flux.
- Sector 98: QLP full-frame-image light curve, 200 s, SAP flux.

Every cadence with a non-zero QUALITY flag is removed.
A wotan biweight filter removes the stellar and instrumental trends.
"""

from pathlib import Path

import numpy as np
from astropy.io import fits
from wotan import flatten

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
LIGHT_CURVE_DIRECTORY = REPOSITORY_ROOT / "data/tic4206066_mast"

LIGHT_CURVE_FILE_AND_FLUX_COLUMN_BY_SECTOR = {
    6: ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX"),
    32: ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX"),
    98: ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SAP_FLUX"),
}

# Three times the 1.87 h transit duration that Rabtsevich reports is 0.23 d.
# A 0.5 d window keeps the transit shape and removes slower variability.
BIWEIGHT_WINDOW_LENGTH_DAYS = 0.5


def load_sector_light_curve(sector):
    """Return time (BTJD), normalised flux and flux error, with flagged and NaN cadences removed."""
    file_name, flux_column = LIGHT_CURVE_FILE_AND_FLUX_COLUMN_BY_SECTOR[sector]
    data = fits.open(LIGHT_CURVE_DIRECTORY / file_name)[1].data
    time = np.asarray(data["TIME"], dtype=float)
    flux = np.asarray(data[flux_column], dtype=float)
    if flux_column == "SAP_FLUX":
        # QLP gives no SAP flux error. The scatter of the detrended flux replaces it later.
        flux_error = np.full_like(flux, np.nan)
    else:
        flux_error = np.asarray(data[flux_column + "_ERR"], dtype=float)
    keep = (data["QUALITY"] == 0) & np.isfinite(time) & np.isfinite(flux)
    time, flux, flux_error = time[keep], flux[keep], flux_error[keep]
    median_flux = np.median(flux)
    return time, flux / median_flux, flux_error / median_flux


def detrend_sector_light_curve(sector):
    """Return time, detrended flux and flux error for one sector.

    The function removes points more than 5 sigma above the trend, and keeps low points because transits are low.
    Where the product has no flux error, the error is the robust scatter of the detrended flux.
    """
    time, flux, flux_error = load_sector_light_curve(sector)
    detrended_flux = flatten(
        time, flux, method="biweight", window_length=BIWEIGHT_WINDOW_LENGTH_DAYS, edge_cutoff=0.0
    )
    keep = np.isfinite(detrended_flux)
    time, detrended_flux, flux_error = time[keep], detrended_flux[keep], flux_error[keep]
    robust_scatter = 1.4826 * np.median(np.abs(detrended_flux - np.median(detrended_flux)))
    keep = detrended_flux < 1 + 5 * robust_scatter
    time, detrended_flux, flux_error = time[keep], detrended_flux[keep], flux_error[keep]
    if not np.all(np.isfinite(flux_error)):
        flux_error = np.full_like(detrended_flux, robust_scatter)
    return time, detrended_flux, flux_error
