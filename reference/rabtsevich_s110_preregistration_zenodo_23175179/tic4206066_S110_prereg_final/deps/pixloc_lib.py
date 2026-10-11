"""Fresh pixel-level localisation code for TIC 4206066 (re-derivation items 10 and 11).

Written 2026-09-30 from scratch. Uses only numpy, scipy, astropy and TESS_PRF 0.1.3.
No project analysis code is imported.

Conventions
-----------
* Pixel coordinates are 0-based, integer = pixel centre (astropy WCS origin=0 and TESS_PRF.locate agree).
* Images are arrays [row, col] = [y, x].
* Sky offsets are (dE, dN) in arcsec relative to the Gaia DR3 position of TIC 4206066 propagated
  with its proper motion to the mid-time of each sector.
* A difference image is  mean(out-of-transit) - mean(in-transit), so a dimming source is POSITIVE.
"""
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage
import PRF

# default locations of the archival inputs, next to the deposit folder (s110_localise.py sets PRF_DIR and FILES[32])
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
PRF_DIR = os.path.join(ROOT, "tess_prf/start_s0004/")
FILES = {
    32: os.path.join(ROOT, "archival_inputs/mast/tess2020324010417-s0032-0000000004206066-0200-s_tp.fits"),
    98: os.path.join(ROOT, "archival_inputs/tesscut/tess-s0098-1-2_80.859830_-1.326960_21x21_astrocut.fits"),
    6: os.path.join(ROOT, "archival_inputs/tesscut/tess-s0006-1-1_80.859830_-1.326960_21x21_astrocut.fits"),
}
QLP_S98 = os.path.join(ROOT, "archival_inputs/mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")

# Signal A ephemeris of the archival analysis (replaced at run time by s110_localise.py)
P = 3.182785293657263
T0 = 3998.7367930333035
T14 = 1.874 / 24.0  # days

# Gaia DR3 3220388198192519424, epoch J2016.0 (VizieR I/355/gaiadr3)
GAIA_RA, GAIA_DEC = 80.85982665567, -1.32696401034
PMRA, PMDEC = -3.406, 70.806  # mas/yr
# The TIC lists RA 80.859842, Dec -1.32728: that is the TIC (epoch J2000.0) position, i.e. the Gaia
# position moved back 16 yr by the proper motion (1.13 arcsec south). Both are the same star.

HARD_BITS = 1610619071  # QLP default mask: bits 1,2,4,8,16,32,128,2048,4096,2^29,2^30
BOX = 11  # side of the fitting box in pixels


def target_coords_at(btjd):
    """Gaia DR3 position propagated to a BTJD epoch (J2016.0 = JD 2457389.0)."""
    dt_yr = (np.asarray(btjd) + 2457000.0 - 2457389.0) / 365.25
    dec = GAIA_DEC + PMDEC / 3.6e6 * dt_yr
    ra = GAIA_RA + PMRA / 3.6e6 * dt_yr / np.cos(np.deg2rad(GAIA_DEC))
    return ra, dec


class SectorData:
    """Pixel cube of one sector, cut to an 11x11 box around the target."""

    def __init__(self, sector, mask_bits=HARD_BITS, qlp_bits=0, verbose=True):
        self.sector = sector
        h = fits.open(FILES[sector])
        d = h[1].data
        self.time = np.array(d["TIME"], float)
        self.quality = np.array(d["QUALITY"], int)
        self.cadenceno = np.array(d["CADENCENO"], int)
        flux = np.array(d["FLUX"], np.float32)
        ferr = np.array(d["FLUX_ERR"], np.float32)
        self.pc1 = np.array(d["POS_CORR1"], float) if "POS_CORR1" in d.names else None
        self.pc2 = np.array(d["POS_CORR2"], float) if "POS_CORR2" in d.names else None
        self.cadence_days = float(np.nanmedian(np.diff(self.time)))
        self.wcs_full = WCS(h[2].header)
        h1 = h[1].header
        self.col0 = int(h1["1CRV4P"])  # full-frame column of pixel x=0
        self.row0 = int(h1["2CRV4P"])
        self.camera, self.ccd = int(h[0].header["CAMERA"]), int(h[0].header["CCD"])
        self.aperture_hdu = np.array(h[2].data)
        h.close()

        # extra QLP flags (S98 only): the TESScut CADENCENO column is all zero, so match by TIME
        # (QLP TIME differs from the cutout TIME by -7 to +12 s; nearest neighbour within 60 s)
        self.qlp_mask_note = "none"
        qlp_flags = np.zeros_like(self.quality)
        if qlp_bits and sector == 98:
            q = fits.open(QLP_S98)[1].data
            qt = np.array(q["TIME"], float)
            qq = np.array(q["QUALITY"], int)
            order = np.argsort(qt)
            j = np.searchsorted(qt[order], self.time)
            j = np.clip(j, 1, len(qt) - 1)
            left, right = order[j - 1], order[j]
            nearest = np.where(np.abs(qt[left] - self.time) < np.abs(qt[right] - self.time), left, right)
            okm = np.abs(qt[nearest] - self.time) * 86400 < 60
            qlp_flags[okm] = qq[nearest[okm]]
            self.qlp_mask_note = f"QLP bits {qlp_bits} matched by TIME ({okm.sum()} of {len(okm)} cadences matched)"
        self.qlp_flags = qlp_flags

        # target pixel position (full cutout coordinates) at mid-sector epoch
        self.tmid = float(np.nanmean(self.time))
        self.ra_t, self.dec_t = target_coords_at(self.tmid)
        xt, yt = self.wcs_full.all_world2pix([[self.ra_t, self.dec_t]], 0)[0]
        ny, nx = flux.shape[1:]
        half = BOX // 2
        xc, yc = int(round(xt)), int(round(yt))
        x0, y0 = max(0, min(xc - half, nx - BOX)), max(0, min(yc - half, ny - BOX))
        self.x0, self.y0 = x0, y0
        self.xt, self.yt = xt - x0, yt - y0  # in box coordinates
        self.full_shape = (ny, nx)

        # background (S98/S6 FFI cutouts are not background subtracted; SPOC TPF S32 is)
        self.bkg = np.zeros(len(self.time))
        if sector != 32:
            mean_img = np.nanmedian(flux, axis=0)
            faint = mean_img <= np.nanpercentile(mean_img, 40)
            self.bkg = np.nanmedian(flux[:, faint], axis=1)
            flux = flux - self.bkg[:, None, None]
            self.bkg_note = f"per-cadence median of the {faint.sum()} faintest pixels (lowest 40%) of the {ny}x{nx} cutout"
        else:
            self.bkg_note = "SPOC FLUX is already background subtracted (FLUX_BKG column)"

        self.flux = np.array(flux[:, y0:y0 + BOX, x0:x0 + BOX], np.float64)
        self.ferr = np.array(ferr[:, y0:y0 + BOX, x0:x0 + BOX], np.float64)
        good = np.isfinite(self.time) & ((self.quality & mask_bits) == 0) & ((qlp_flags & qlp_bits) == 0)
        good &= np.all(np.isfinite(self.flux), axis=(1, 2)) & np.all(np.isfinite(self.ferr), axis=(1, 2))
        if self.pc1 is not None:
            good &= np.isfinite(self.pc1) & np.isfinite(self.pc2)
        self.good = good

        # PRF at the detector location of the box centre
        self.prf = PRF.TESS_PRF(self.camera, self.ccd, sector, self.col0 + x0 + half, self.row0 + y0 + half,
                                localdatadir=PRF_DIR)
        # photometric aperture for the depth measurement
        if sector == 32:
            ap = (self.aperture_hdu & 2) > 0
            self.ap = ap[y0:y0 + BOX, x0:x0 + BOX]
            self.ap_note = "SPOC optimal aperture (bit 2 of the APERTURE extension), 18 pixels"
        else:
            self.ap = np.zeros((BOX, BOX), bool)
            xi, yi = int(round(self.xt)), int(round(self.yt))
            self.ap[yi - 1:yi + 2, xi - 1:xi + 2] = True
            self.ap_note = "3x3 pixels centred on the target pixel"
        if verbose:
            print(f"S{sector}: {good.sum()}/{len(good)} good cadences, cadence {self.cadence_days*86400:.0f} s, "
                  f"target at box pixel ({self.xt:.3f},{self.yt:.3f}), box origin ({x0},{y0}), "
                  f"detector col/row {self.col0 + x0 + half}/{self.row0 + y0 + half}; {self.qlp_mask_note}")

    # ---- coordinate helpers -------------------------------------------------------------
    def pix_to_offset(self, x, y):
        """Box pixel -> (dE, dN) arcsec from the propagated target position."""
        ra, dec = self.wcs_full.all_pix2world(np.c_[np.atleast_1d(x) + self.x0, np.atleast_1d(y) + self.y0], 0).T
        dE = (ra - self.ra_t) * np.cos(np.deg2rad(self.dec_t)) * 3600.0
        dN = (dec - self.dec_t) * 3600.0
        return np.squeeze(dE), np.squeeze(dN)

    def offset_to_pix(self, dE, dN):
        ra = self.ra_t + np.atleast_1d(dE) / 3600.0 / np.cos(np.deg2rad(self.dec_t))
        dec = self.dec_t + np.atleast_1d(dN) / 3600.0
        x, y = self.wcs_full.all_world2pix(np.c_[ra, dec], 0).T
        return np.squeeze(x - self.x0), np.squeeze(y - self.y0)

    def offset_to_pixshift(self, dE, dN):
        """Pixel displacement that corresponds to a sky offset (dE, dN) at the target."""
        x, y = self.offset_to_pix(dE, dN)
        return float(x - self.xt), float(y - self.yt)

    # ---- transit bookkeeping -------------------------------------------------------------
    def transit_times(self, phase_shift=0.0):
        t = self.time[self.good]
        n = np.arange(np.floor((t.min() - T0) / P) - 1, np.ceil((t.max() - T0) / P) + 2)
        tc = T0 + (n + phase_shift) * P
        return tc[(tc > t.min() - 2 * T14) & (tc < t.max() + 2 * T14)]

    def windows(self, tc, W=1.0, gap=1.5):
        """Boolean masks (over all cadences) for in-transit and the two out-of-transit windows."""
        dt = self.time - tc
        w_in = self.good & (np.abs(dt) <= T14 / 2)
        w_pre = self.good & (dt < -gap * T14) & (dt >= -(gap + W) * T14)
        w_post = self.good & (dt > gap * T14) & (dt <= (gap + W) * T14)
        return w_in, w_pre, w_post

    def usable_transits(self, tc_list, W=1.0, gap=1.5, min_cov=0.8, min_side=0.3):
        out = []
        n_in_exp = T14 / self.cadence_days
        n_out_exp = W * T14 / self.cadence_days
        for tc in tc_list:
            w_in, w_pre, w_post = self.windows(tc, W, gap)
            cov_in = w_in.sum() / n_in_exp
            cov_pre, cov_post = w_pre.sum() / n_out_exp, w_post.sum() / n_out_exp
            # in-transit coverage >= min_cov; out-of-transit: both sides together >= 0.6 of expected,
            # each side >= min_side (so the baseline brackets the transit)
            ok = (cov_in >= min_cov) and ((cov_pre + cov_post) / 2 >= 0.6) and \
                 (cov_pre >= min_side) and (cov_post >= min_side)
            out.append(dict(tc=float(tc), cov_in=float(cov_in), cov_pre=float(cov_pre), cov_post=float(cov_post),
                            n_in=int(w_in.sum()), n_out=int(w_pre.sum() + w_post.sum()), ok=bool(ok)))
        return out

    # ---- pointing proxies and decorrelation ---------------------------------------------
    def pointing_proxy(self, kind, flux=None, half=2):
        """Two per-cadence pointing regressors.
        kind='poscorr'  : SPOC POS_CORR1/2 (S32 TPF; all zero in TESScut files).
        kind='centroid' : flux-weighted centroid of the target in a (2*half+1)^2 box (proxy for pointing;
                          in S32 it correlates with POS_CORR at r=0.993 (x) and 0.9998 (y))."""
        if kind == "poscorr":
            if self.pc1 is None or not np.any(self.pc1 != 0):
                raise ValueError(f"S{self.sector}: POS_CORR not available (all zero or missing)")
            return self.pc1.copy(), self.pc2.copy()
        F = self.flux if flux is None else flux
        xi, yi = int(round(self.xt)), int(round(self.yt))
        sl = (slice(yi - half, yi + half + 1), slice(xi - half, xi + half + 1))
        yy, xx = np.mgrid[sl]
        sub = F[:, sl[0], sl[1]]
        tot = sub.sum((1, 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            cx = (sub * xx).sum((1, 2)) / tot
            cy = (sub * yy).sum((1, 2)) / tot
        return cx, cy

    def decorrelated_flux(self, exclude_tc, kind="poscorr", W=1.0, gap=1.5, flux=None):
        """Per pixel linear fit  f = a + b*r1 + c*r2  on out-of-transit cadences (all good cadences
        farther than (gap+W)*T14 from every transit in exclude_tc); subtract b*r1 + c*r2 everywhere.
        r1, r2 = pointing_proxy(kind). Returns a new cube."""
        F = self.flux if flux is None else flux
        r1, r2 = self.pointing_proxy(kind, flux=F)
        oot = self.good & np.isfinite(r1) & np.isfinite(r2)
        for tc in exclude_tc:
            oot &= np.abs(self.time - tc) > (gap + W) * T14
        m1, m2 = r1[oot].mean(), r2[oot].mean()
        X = np.c_[np.ones(oot.sum()), r1[oot] - m1, r2[oot] - m2]
        Y = F[oot].reshape(oot.sum(), -1)
        coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
        model = np.nan_to_num((r1 - m1))[:, None] * coef[1][None, :] + np.nan_to_num((r2 - m2))[:, None] * coef[2][None, :]
        self.decorr_coef = coef
        return F - model.reshape(F.shape)


# ---- difference imaging ---------------------------------------------------------------------
def stack_difference(sd, tc_list, W=1.0, gap=1.5, flux=None, min_cov=0.8):
    """Per-transit difference images (out - in) and their weighted stack.
    Weight w_k = 1/(1/n_in + 1/n_out) (effective number of cadences).
    Returns dict with 'diff' (stack), 'sig' (per-pixel error from FLUX_ERR), 'per_transit', 'mean_out', 'transits'."""
    F = sd.flux if flux is None else flux
    info = sd.usable_transits(tc_list, W, gap, min_cov)
    diffs, sigs, ws, outs, used = [], [], [], [], []
    for tr in info:
        if not tr["ok"]:
            continue
        w_in, w_pre, w_post = sd.windows(tr["tc"], W, gap)
        w_out = w_pre | w_post
        m_in, m_out = F[w_in].mean(0), F[w_out].mean(0)
        v = (sd.ferr[w_in] ** 2).sum(0) / w_in.sum() ** 2 + (sd.ferr[w_out] ** 2).sum(0) / w_out.sum() ** 2
        w = 1.0 / (1.0 / w_in.sum() + 1.0 / w_out.sum())
        diffs.append(m_out - m_in); sigs.append(v); ws.append(w); outs.append(m_out); used.append(tr)
    if not diffs:
        return None
    diffs, sigs, ws, outs = np.array(diffs), np.array(sigs), np.array(ws), np.array(outs)
    stack = (ws[:, None, None] * diffs).sum(0) / ws.sum()
    sig = np.sqrt((ws[:, None, None] ** 2 * sigs).sum(0)) / ws.sum()
    mean_out = (ws[:, None, None] * outs).sum(0) / ws.sum()
    return dict(diff=stack, sig=sig, per_transit=diffs, weights=ws, mean_out=mean_out, transits=used,
                n_transits=len(used))


def box_depth(sd, tc_list, W=1.0, gap=1.5, flux=None, ap=None, min_cov=0.8):
    """Box depth of the aperture light curve: 1 - mean_in/mean_out per transit, weighted mean."""
    F = sd.flux if flux is None else flux
    ap = sd.ap if ap is None else ap
    lc = F[:, ap].sum(1)
    info = sd.usable_transits(tc_list, W, gap, min_cov)
    d, w = [], []
    for tr in info:
        if not tr["ok"]:
            continue
        w_in, w_pre, w_post = sd.windows(tr["tc"], W, gap)
        fo = lc[w_pre | w_post].mean()
        d.append(1 - lc[w_in].mean() / fo)
        w.append(1.0 / (1.0 / w_in.sum() + 1.0 / (w_pre.sum() + w_post.sum())))
    d, w = np.array(d), np.array(w)
    depth = (w * d).sum() / w.sum()
    err = np.sqrt((w ** 2 * (d - depth) ** 2).sum()) / w.sum() if len(d) > 1 else np.nan
    return dict(depth=float(depth), err=float(err), per_transit=d.tolist(), n=len(d))


# ---- PRF fitting ---------------------------------------------------------------------------
def _linear_amp(model, D, wgt):
    """Weighted LS for D = A*model + c. Returns A, c, chi2."""
    m, d, w = model.ravel(), D.ravel(), wgt.ravel()
    Sw, Sm, Sd = w.sum(), (w * m).sum(), (w * d).sum()
    Smm, Smd = (w * m * m).sum(), (w * m * d).sum()
    det = Sw * Smm - Sm * Sm
    if det <= 0:
        return 0.0, 0.0, np.inf
    A = (Sw * Smd - Sm * Sd) / det
    c = (Smm * Sd - Sm * Smd) / det
    r = d - A * m - c
    return A, c, float((w * r * r).sum())


def fit_prf(sd, D, sig, x_init=None, y_init=None, search=2.5, coarse=0.1, fine=0.01, fine_half=0.15):
    """Fit A*PRF(x,y)+c to image D with per-pixel errors sig by a profile-likelihood grid search:
    coarse grid (step `coarse` px within +-search px of the initial position) then a fine grid
    (step `fine` px within +-fine_half px of the coarse minimum). Amplitude and constant are solved
    linearly at each grid point. Returns dict with x, y, A, c, chi2, chi2_at_target, A_at_target, snr_at_target."""
    x_init = sd.xt if x_init is None else x_init
    y_init = sd.yt if y_init is None else y_init
    wgt = 1.0 / np.maximum(sig, 1e-12) ** 2
    shape = D.shape

    def chi2_at(x, y):
        m = sd.prf.locate(x, y, shape)
        return _linear_amp(m, D, wgt)

    best = (np.inf, None)
    xs = np.arange(x_init - search, x_init + search + 1e-9, coarse)
    ys = np.arange(y_init - search, y_init + search + 1e-9, coarse)
    xs = xs[(xs > -0.5) & (xs < shape[1] - 0.5)]
    ys = ys[(ys > -0.5) & (ys < shape[0] - 0.5)]
    for x in xs:
        for y in ys:
            A, c, chi2 = chi2_at(x, y)
            if chi2 < best[0] and A > 0:
                best = (chi2, (x, y, A, c))
    if best[1] is None:  # no positive-amplitude minimum: fall back to unconstrained
        for x in xs:
            for y in ys:
                A, c, chi2 = chi2_at(x, y)
                if chi2 < best[0]:
                    best = (chi2, (x, y, A, c))
    xb, yb = best[1][0], best[1][1]
    for x in np.arange(xb - fine_half, xb + fine_half + 1e-9, fine):
        for y in np.arange(yb - fine_half, yb + fine_half + 1e-9, fine):
            A, c, chi2 = chi2_at(x, y)
            if chi2 < best[0]:
                best = (chi2, (x, y, A, c))
    chi2, (x, y, A, c) = best
    # amplitude and its error at the fixed target position
    m_t = sd.prf.locate(sd.xt, sd.yt, shape)
    A_t, c_t, chi2_t = _linear_amp(m_t, D, wgt)
    mm, w = m_t.ravel(), wgt.ravel()
    Sw, Sm, Smm = w.sum(), (w * mm).sum(), (w * mm * mm).sum()
    varA = Sw / (Sw * Smm - Sm * Sm)
    return dict(x=float(x), y=float(y), A=float(A), c=float(c), chi2=float(chi2), chi2_at_target=float(chi2_t),
                A_at_target=float(A_t), snr_at_target=float(A_t / np.sqrt(varA)), dchi2_target=float(chi2_t - chi2))


def fit_star(sd, image):
    """PRF fit of the mean stellar image (position check of the WCS)."""
    sig = np.sqrt(np.maximum(image, 0) + np.median(sd.ferr[sd.good], axis=0) ** 2)
    return fit_prf(sd, image, sig, search=1.5)


# ---- centroid method ----------------------------------------------------------------------------
def centroid_method(sd, tc_list, W=1.0, gap=1.5, flux=None, half=2, depth=None, min_cov=0.8):
    """Flux-weighted centroid in a (2*half+1)^2 box around the target pixel, in minus out of transit,
    converted to an implied source position: s = c_out - dc*(1-delta)/delta, delta = box-flux depth."""
    F = sd.flux if flux is None else flux
    xi, yi = int(round(sd.xt)), int(round(sd.yt))
    sl = (slice(yi - half, yi + half + 1), slice(xi - half, xi + half + 1))
    yy, xx = np.mgrid[sl]
    sub = F[:, sl[0], sl[1]]
    tot = sub.sum((1, 2))
    cx = (sub * xx).sum((1, 2)) / tot
    cy = (sub * yy).sum((1, 2)) / tot
    info = sd.usable_transits(tc_list, W, gap, min_cov)
    dcx, dcy, cox, coy, dd, ws = [], [], [], [], [], []
    for tr in info:
        if not tr["ok"]:
            continue
        w_in, w_pre, w_post = sd.windows(tr["tc"], W, gap)
        w_out = w_pre | w_post
        dcx.append(cx[w_in].mean() - cx[w_out].mean()); dcy.append(cy[w_in].mean() - cy[w_out].mean())
        cox.append(cx[w_out].mean()); coy.append(cy[w_out].mean())
        dd.append(1 - tot[w_in].mean() / tot[w_out].mean())
        ws.append(1.0 / (1.0 / w_in.sum() + 1.0 / w_out.sum()))
    ws = np.array(ws)
    def wm(a): return float((ws * np.array(a)).sum() / ws.sum())
    dcx, dcy, cox, coy, delta = wm(dcx), wm(dcy), wm(cox), wm(coy), wm(dd)
    if depth is None:
        depth = delta
    sx = cox - dcx * (1 - depth) / depth
    sy = coy - dcy * (1 - depth) / depth
    return dict(dcx=dcx, dcy=dcy, cx_out=cox, cy_out=coy, depth_box=delta, x=float(sx), y=float(sy))


# ---- injections --------------------------------------------------------------------------------
def shifted_template(image, dx, dy):
    """Shift an image by (dx, dy) pixels (cubic spline, edges padded with nearest value)."""
    return ndimage.shift(image, (dy, dx), order=3, mode="nearest")


def inject(sd, flux, template, depth, tc_list):
    """Subtract depth*template from every good cadence within +-T14/2 of each tc (box transit)."""
    F = flux.copy()
    for tc in tc_list:
        w = sd.good & (np.abs(sd.time - tc) <= T14 / 2)
        F[w] -= depth * template[None, :, :]
    return F


# ---- statistics ---------------------------------------------------------------------------------
def cov_stats(pts):
    pts = np.asarray(pts, float)
    mu = pts.mean(0)
    C = np.cov(pts.T)
    return mu, C


def mahalanobis(x, mu, C):
    d = np.asarray(x, float) - np.asarray(mu, float)
    return float(np.sqrt(d @ np.linalg.solve(C, d)))


def chi2_2dof_p(d):
    return float(np.exp(-0.5 * d * d))


def gauss_loglike(x, mu, C):
    d = np.asarray(x, float) - np.asarray(mu, float)
    return float(-0.5 * d @ np.linalg.solve(C, d) - 0.5 * np.log(np.linalg.det(2 * np.pi * C)))


def combine(xs, Cs):
    """Inverse-covariance weighted combination of independent 2-D estimates."""
    Ci = [np.linalg.inv(C) for C in Cs]
    Cc = np.linalg.inv(sum(Ci))
    xc = Cc @ sum(ci @ np.asarray(x) for ci, x in zip(Ci, xs))
    return xc, Cc


def ellipse_points(mu, C, nsig=1.0, n=200):
    vals, vecs = np.linalg.eigh(C)
    t = np.linspace(0, 2 * np.pi, n)
    circ = np.c_[np.cos(t), np.sin(t)] * np.sqrt(np.maximum(vals, 0)) * nsig
    return (vecs @ circ.T).T + np.asarray(mu)
