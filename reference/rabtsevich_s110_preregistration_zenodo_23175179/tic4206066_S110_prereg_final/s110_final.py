"""TIC 4206066, TESS Sector 110 test: the analysis of the FINAL pre-registration (written 2026-10-02, before any
Sector 110 data exist). See PREREG_S110_FINAL.md for the rules in words; this script is the rules in code.

It is run once on the Sector 110 light curve, unchanged, with the primary settings, before any plot or other look at
any Sector 110 product. The variants listed in the document are run afterwards and reported next to it.

Per signal (A and B are tested separately; there is no combined verdict):
  * predicted mid-times: ephemeris_final.json (A: final re-derivation; B: unchanged from v1/v2). When the light curve
    covers the frozen Sector 110 epochs, every frozen epoch is listed, observed or not;
  * EDGE ZONES: for every detrending (real or fake ephemeris) the series actually passed to wotan (cleaned light curve
    with that run's transit masks removed) is split into continuous segments at gaps longer than BREAK_TOL (0.3 d, the
    break tolerance lcio.detrend passes to wotan), i.e. exactly the segments wotan fits; in each segment the first and
    last W/2 are edge zones (W = detrending window: 0.75 d for the primary, so W/2 = 0.375 d = 9 h), because there the
    centred filter has less than half a window of data on one side (wotan: edge cut-off is at most W/2, "as this fills
    the window completely"). A transit whose mask joins a data gap into a break lies between segments: not tested;
  * a predicted transit is COUNTED if (i) >= 80 % of its T14 window has cadences after cleaning, (ii) it has >= 2
    cadences on each side between 0.75 and 3 T14 from the predicted mid-time (baseline rule "T3" of the re-derivation,
    item 03), and (iii) its whole T14 window lies inside one segment and outside its edge zones; (iv) in Sector 110
    only the frozen epochs can be counted. A transit that fails any of these is NOT
    TESTED: listed with its status and depth, never in the statistic, and never counted as a miss;
  * single-transit box depth = 1 - mean(detrended flux inside |t - tc| < T14/2); statistic d = mean over the n counted
    transits (no period or time scan);
  * NOISE MODEL (null pool): the predicted windows of both signals are cut (+-1.5 T14); then for each fake phase offset
    phi in [0.15, 0.85] P (steps of T14/2) the light curve is detrended with the fake ephemeris masked exactly as the
    real one is (+-1.5 T14), and the windows at the fake times that pass the same counting rule (with the segments and
    edge zones of that fake detrending) form the pool. So the noise windows are processed exactly like the real windows. In addition the whole set
    of windows at each fake phase gives a full-set statistic d_phi; the fraction of fake phases with d_phi >= d is
    reported (supplementary, not part of the outcome: the phases resolve only a few per cent);
  * det = 99.9th percentile of the pool's n-window mean (bootstrap); lo = 0.1th percentile of noise + N(D, D_err);
    hi = 99.9th percentile of noise + N(D98, D98_err). (The bootstrap of a finite pool of overlapping windows
    underestimates the tails: with 99.7 / 0.3 the realised false-alarm level measured on archival fake ephemerides was
    0.8 % (A) and 1.2 % (B) with the earlier segment definition, audit 2. 99.9 / 0.1 is the stricter level adopted;
    the realised levels per product and per number of transits are in calib/false_alarm_*.json and section 4.)
    P(detect | real) = fraction of noise + N(D, D_err) above det;
    medians over N_SEEDS bootstrap seeds, spread reported. D, D98: calib/expected_depths.json (same statistic, same
    counting rule, same noise model, archival sectors; D from S6 + S32, D98 from the discovery sector S98);
  * timing (only if d > det): common shift dt in -30..+30 min (0.5 min steps) maximising the mean depth of the counted
    transits; if several shifts share the maximum, the middle one of them (ordered by shift) is taken; sigma_meas = SD
    of the recovered shift over N_TIM injections of a box of depth d into n random pool windows (true shift 0);
    sigma_pred = mean of the frozen per-epoch prediction sigmas of the counted epochs; consistent if
    |dt| <= 3 sqrt(sigma_pred^2 + sigma_meas^2). Individual transit shifts are reported; they do not enter the outcome.
Outcomes:
  O1 RECOVERED      d > det  and lo <= d <= hi  and timing consistent
  O2 INCONSISTENT   d > det  and (d < lo or d > hi or timing inconsistent)
  O3 EXCLUDED       d <= det and d < lo and P(detect | real) >= 0.95
  O4 INSUFFICIENT   everything else; the JSON gives one fixed reason string (O4_REASONS)
Registered runs (no --dry-run) refuse to start unless the input header is Sector 110 / TIC 4206066 and the installed
package versions equal requirements.txt; a refused run writes nothing and is not a run. --dry-run (archival files only)
lifts the sector check, never the environment check, and marks the JSON as a dry run.

Usage (exact command lines for Sector 110: PREREG_S110_FINAL.md section 9):
  python s110_final.py --lc FILE --product spoc|qlp [--flux COL] [--window 0.75] [--drop-desat-hours 0]
                       [--count-rule final|v2] [--t0-shift-frac 0] [--cut-true] --out OUT.json
  python s110_final.py --calibrate --out calib/expected_depths.json     (archival expected depths; run once, frozen)
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import sys
from importlib import metadata

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "deps"))
import lcio  # noqa: E402  (frozen audit loader/detrender; same code as the v2 package copy, one comment and one docstring made neutral on 2026-10-05)

EPHF = json.load(open(os.path.join(HERE, "ephemeris_final.json")))
EPH = {k: dict(P=float(EPHF[k]["P"]), T0=float(EPHF[k]["T0"]), T14=float(EPHF[k]["T14"])) for k in ("A", "B")}
SIGMA_PRED = {k: {int(e["epoch"]): float(e["sigma_pred_min"]) for e in EPHF[k]["epochs"]} for k in ("A", "B")}
FROZEN_EPOCHS = {k: [(int(e["epoch"]), float(e["tc_btjd"])) for e in EPHF[k]["epochs"]] for k in ("A", "B")}
_C = np.array(EPHF["A"]["cov_T0_P"])

N_BOOT, N_SEEDS, SEED0 = 20000, 20, 110210
DET_PCT, LO_PCT, HI_PCT = 99.9, 0.1, 99.9           # see docstring (audit 2)
N_TIM, SEED_TIM = 300, 220110
PHI_MIN, PHI_MAX, PHI_STEP_T14 = 0.15, 0.85, 0.5   # fake phases, step in units of T14
CUT_FAC = 1.5                      # real windows cut from the light curve before the null pool is built
MIN_POOL = 100                     # fewer pool windows -> O4
SENS_MIN = 0.95
COV_MIN = 0.8                      # in-transit coverage of a counted transit
BASE_IN, BASE_OUT, BASE_MIN = 0.75, 3.0, 2   # baseline rule T3: >= 2 cadences per side in 0.75..3 T14
BREAK_TOL = 0.3                    # segment split at gaps > 0.3 d (= lcio.detrend break_tol, passed to wotan)
SHIFTS_MIN = np.arange(-30, 30.01, 0.5)
SEG_HALF_MIN = 30.5                # pool segments keep |t - tc| < T14/2 + 30.5 min (enough for the shift scan)
FROZEN_QLP_MASK = lcio.FROZEN_QLP_MASK   # 1610619071
SPOC_MASK = lcio.DEFAULT_BITMASK         # 17087
DESAT = 32                               # SPOC quality bit "desaturation event" (momentum dump)
PRIMARY_ARGS = dict(product="spoc", flux="PDCSAP_FLUX", window=0.75, drop_desat_hours=0.0, count_rule="final",
                    t0_shift_frac=0.0, cut_true=False)
FALLBACK1_ARGS = dict(PRIMARY_ARGS, product="qlp", flux="QLP_DEFAULT")   # fallback 1 (QLP); fallback 2 = PRIMARY_ARGS
TARGET_TIC, TARGET_SECTOR = 4206066, 110
import re  # noqa: E402
FILE_PATTERNS = {   # registered input products for Sector 110 (file names as distributed by MAST)
    "spoc_2min": re.compile(r"^tess\d+-s0110-0000000004206066-\d+-s_lc\.fits$"),
    "tessspoc": re.compile(r"^hlsp_tess-spoc_tess_phot_0000000004206066-s0110_tess_v\d+_lc\.fits$"),
    "qlp": re.compile(r"^hlsp_qlp_tess_ffi_s0110-0000000004206066_tess_v\d+_llc\.fits$"),
}
O4_REASONS = ("fewer than 10 usable cadences", "no registered QLP flux column", "no counted transit",
              "noise pool has fewer than 100 windows", "d not above det and not below lo", "sensitivity below 0.95",
              "d not above det and not below lo; sensitivity below 0.95")

# archival products for --calibrate (same files and masks as v1/v2)
MAST = os.path.join(HERE, "..", "archival_inputs", "mast")
ARCHIVAL = {
    "S6": ("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "spoc", "PDCSAP_FLUX"),
    "S32": ("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "spoc", "PDCSAP_FLUX"),
    "S98": ("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "qlp", "SYS_RM_FLUX"),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def versions():
    v = dict(python=platform.python_version(), machine=platform.machine(), system=platform.system())
    for pkg in ("numpy", "astropy", "scipy", "wotan", "numba", "llvmlite", "TESS_PRF"):
        try:
            v[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            v[pkg] = "not installed"
    return v


def pinned_versions():
    """Versions in requirements.txt (python from its first line)."""
    req, py = {}, None
    for line in open(os.path.join(HERE, "requirements.txt")):
        line = line.strip()
        if line.startswith("#") and "python " in line:
            py = line.split("python ")[1].split(")")[0].strip()
        elif line.startswith("# machine:") or line.startswith("# system:"):
            k, v = line[2:].split(":", 1)
            req[k.strip()] = v.strip()
        elif "==" in line:
            k, v = line.split("==")
            req[k.strip()] = v.strip()
    if py:
        req["python"] = py
    return req


def environment_mismatch():
    have = versions()
    return {k: (v, have.get(k)) for k, v in pinned_versions().items() if have.get(k) != v}


def registered_settings(a):
    """True if the arguments are the registered primary (SPOC) or fallback-1 (QLP) settings."""
    return any(all(vars(a).get(q) == v for q, v in reg.items()) for reg in (PRIMARY_ARGS, FALLBACK1_ARGS))


def provenance(script_path, input_paths=()):
    p = dict(run_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
             script=os.path.basename(script_path), script_sha256=sha256(script_path),
             deps_sha256={f: sha256(os.path.join(HERE, "deps", f)) for f in sorted(os.listdir(os.path.join(HERE, "deps")))
                          if f.endswith(".py")},
             ephemeris_sha256=sha256(os.path.join(HERE, "ephemeris_final.json")), versions=versions())
    cal = os.path.join(HERE, "calib", "expected_depths.json")
    if os.path.exists(cal):
        p["expected_depths_sha256"] = sha256(cal)
    p["inputs"] = {os.path.basename(x): sha256(x) for x in input_paths}
    return p


QLP_FLUX_ORDER = ("SYS_RM_FLUX", "DET_FLUX")   # QLP primary column, then its pre-declared replacement


class NoQLPColumn(Exception):
    """The QLP file has neither registered flux column (O4 for both signals)."""


def load(path, product, flux):
    h = fits.open(path)
    d = h[1].data
    if flux == "QLP_DEFAULT":
        have = [c for c in QLP_FLUX_ORDER if c in d.columns.names]
        if not have:
            raise NoQLPColumn("no QLP flux column of %s in %s" % (QLP_FLUX_ORDER, os.path.basename(path)))
        flux = have[0]
    t = np.asarray(d["TIME"], float)
    f = np.asarray(d[flux], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    desat_t = t[np.isfinite(t) & ((q & DESAT) != 0)]
    good = (q & (SPOC_MASK if product == "spoc" else FROZEN_QLP_MASK)) == 0
    m = np.isfinite(t) & np.isfinite(f) & good & (f > 0)
    t, f = t[m], f[m] / np.nanmedian(f[m])
    return dict(t=t, f=f, q=q[m], cad=float(np.median(np.diff(t))) if len(t) > 1 else np.nan,
                sector=int(h[0].header.get("SECTOR", -1)), tic=h[0].header.get("TICID", None),
                file=os.path.basename(path), flux_column=flux), desat_t


def segments(t):
    """Continuous segments of a time series: (starts, ends), split at gaps > BREAK_TOL."""
    i = np.flatnonzero(np.diff(t) > BREAK_TOL)
    return t[np.r_[0, i + 1]], t[np.r_[i, len(t) - 1]]


def wotan_segments(t, eph, mask_fac=1.5):
    """The segments wotan fits inside lcio.detrend: the series with both signals' +-mask_fac*T14 windows removed
    (lcio.detrend passes t[~msk] to wotan with break_tolerance BREAK_TOL)."""
    msk = lcio.in_transit(t, "A", mask_fac, eph) | lcio.in_transit(t, "B", mask_fac, eph)
    return segments(t[~msk])


def in_clear_zone(tc, T14, segs, window):
    """True if the whole T14 window lies in one segment, at least window/2 from both of its ends."""
    starts, ends = segs
    lo, hi = tc - T14 / 2, tc + T14 / 2
    return bool(np.any((lo >= starts + window / 2) & (hi <= ends - window / 2)))


def window_ok(t, tc, T14, cad, rule, segs, window):
    """(coverage, n_left, n_right, counted, status) for a window centred on tc, with the counting rule."""
    dt = t - tc
    cov = (np.abs(dt) < T14 / 2).sum() * cad / T14
    left = int(((dt < -BASE_IN * T14) & (dt > -BASE_OUT * T14)).sum())
    right = int(((dt > BASE_IN * T14) & (dt < BASE_OUT * T14)).sum())
    if cov < COV_MIN:
        status = "not tested: coverage < 80 %"
    elif rule == "v2":
        status = "counted"
    elif left < BASE_MIN or right < BASE_MIN:
        status = "not tested: no baseline on one side (T3)"
    elif not in_clear_zone(tc, T14, segs, window):
        status = "not tested: in an edge zone (within W/2 of a segment end) or between segments"
    else:
        status = "counted"
    return float(cov), left, right, status == "counted", status


def window_depth(t, fd, tc, T14):
    inn = np.abs(t - tc) < T14 / 2
    return (1 - fd[inn].mean()) * 1e6 if inn.any() else np.nan


def predicted(t, e):
    n0, n1 = np.ceil((t.min() - e["T0"]) / e["P"]), np.floor((t.max() - e["T0"]) / e["P"])
    n = np.arange(n0, n1 + 1).astype(int)
    return n, e["T0"] + e["P"] * n


def best_shift(t, fd, tcs, T14):
    vals = np.array([np.nanmean([window_depth(t, fd, tc + s / 1440, T14) for tc in tcs]) for s in SHIFTS_MIN])
    top = np.flatnonzero(vals == np.nanmax(vals))
    return float(SHIFTS_MIN[top[len(top) // 2]])


def clean(a, lc, desat_t, true_eph):
    """Quality cut is done in load(); here the dry-run cut and the desaturation variant."""
    if a.cut_true:   # dry runs only: remove the real transits so that a shifted ephemeris is a clean null
        keep = ~(lcio.in_transit(lc["t"], "A", 1.5, true_eph) | lcio.in_transit(lc["t"], "B", 1.5, true_eph))
        lc = {**lc, "t": lc["t"][keep], "f": lc["f"][keep], "q": lc["q"][keep]}
    if a.drop_desat_hours > 0 and len(desat_t):
        near = np.min(np.abs(lc["t"][:, None] - desat_t[None, :]), axis=1) < a.drop_desat_hours / 24
        lc = {**lc, "t": lc["t"][~near], "f": lc["f"][~near], "q": lc["q"][~near]}
    return lc


def events(t, fd, cad, k, e, rule, list_frozen, segs, window):
    rows, seen = [], set()
    epochs, tcs = predicted(t, e)
    frozen = {n for n, _ in FROZEN_EPOCHS[k]}
    for n, tc in zip(epochs, tcs):
        cov, left, right, counted, status = window_ok(t, tc, e["T14"], cad, rule, segs, window)
        if list_frozen and int(n) not in frozen:
            counted, status = False, "not tested: not a frozen Sector 110 epoch"
        dep = window_depth(t, fd, tc, e["T14"])
        rows.append(dict(epoch=int(n), tc=float(tc), coverage=cov, baseline_left=left, baseline_right=right,
                         depth_ppm=None if not np.isfinite(dep) else float(dep), counted=counted, status=status))
        seen.add(int(n))
    if list_frozen:
        for n, tc in FROZEN_EPOCHS[k]:
            if n not in seen:
                rows.append(dict(epoch=n, tc=tc, coverage=0.0, baseline_left=0, baseline_right=0, depth_ppm=None,
                                 counted=False, status="not tested: mid-time outside the light curve"))
        rows.sort(key=lambda r: r["tc"])
    return rows


def null_pool(lc, k, eph, window, rule):
    """Pool of noise windows processed exactly like the real ones (see module docstring).
    Returns window depths, segments for the timing injections, and the full-set statistic per fake phase."""
    t = lc["t"]
    cut = lcio.in_transit(t, "A", CUT_FAC, eph) | lcio.in_transit(t, "B", CUT_FAC, eph)
    lcc = {**lc, "t": t[~cut], "f": lc["f"][~cut], "q": lc["q"][~cut]}
    e = eph[k]
    step = PHI_STEP_T14 * e["T14"] / e["P"]
    depths, segs, sets = [], [], []
    for phi in np.arange(PHI_MIN, PHI_MAX + 1e-9, step):
        fake = {kk: dict(v) for kk, v in eph.items()}
        fake[k]["T0"] = e["T0"] + phi * e["P"]
        d = lcio.detrend(lcc, window=window, eph=fake)
        tt, ff = d["t"], d["fd"]
        segs_fake = wotan_segments(lcc["t"], fake)   # the segments this detrending actually used
        here = []
        for r in events(tt, ff, lc["cad"], k, fake[k], rule, False, segs_fake, window):
            if r["counted"] and r["depth_ppm"] is not None:
                s = np.abs(tt - r["tc"]) < e["T14"] / 2 + SEG_HALF_MIN / 1440
                depths.append(r["depth_ppm"])
                here.append(r["depth_ppm"])
                segs.append((r["tc"], tt[s].copy(), ff[s].copy()))
        if here:
            sets.append(dict(phi=float(phi), n=len(here), d=float(np.mean(here))))
    return np.asarray(depths), segs, sets


def thresholds(rng, W, n, D, De, D98, De98):
    null = W[rng.integers(0, len(W), (N_BOOT, n))].mean(1)
    det = float(np.percentile(null, DET_PCT))
    sig = null + rng.normal(D, De, N_BOOT)
    sig98 = null + rng.normal(D98, De98, N_BOOT)
    return dict(det=det, lo=float(np.percentile(sig, LO_PCT)), hi=float(np.percentile(sig98, HI_PCT)),
                sens=float(np.mean(sig > det)), null_sd=float(np.std(null)))


def sigma_meas_inj(rng, segs, T14, n, depth_ppm):
    """SD of the recovered common shift for a box of depth_ppm injected into n random pool windows (true shift 0)."""
    rec = []
    for _ in range(N_TIM):
        pick = rng.choice(len(segs), n, replace=False)
        vals = np.zeros(len(SHIFTS_MIN))
        for j, s in enumerate(SHIFTS_MIN):
            dd = []
            for i in pick:
                tc, ts, fs = segs[i]
                fi = fs - (depth_ppm / 1e6) * (np.abs(ts - tc) < T14 / 2)
                dd.append(window_depth(ts, fi, tc + s / 1440, T14))
            vals[j] = np.nanmean(dd)
        top = np.flatnonzero(vals == np.nanmax(vals))
        rec.append(SHIFTS_MIN[top[len(top) // 2]])
    rec = np.asarray(rec)
    return float(np.std(rec, ddof=1)), float(np.mean(rec))


def outcome(stat, det, lo, hi, sens, timing):
    """(outcome label, fixed reason string or None)."""
    if stat > det:
        if lo <= stat <= hi and timing["consistent"]:
            return "O1_RECOVERED", None
        why = []
        if stat < lo:
            why.append("shallower than predicted")
        if stat > hi:
            why.append("deeper than predicted")
        if not timing["consistent"]:
            why.append("timing off")
        return "O2_INCONSISTENT (" + ", ".join(why) + ")", None
    if stat < lo and sens >= SENS_MIN:
        return "O3_EXCLUDED", None
    r1, r2 = stat >= lo, sens < SENS_MIN
    return "O4_INSUFFICIENT", (O4_REASONS[6] if (r1 and r2) else (O4_REASONS[4] if r1 else O4_REASONS[5]))


def run(a):
    is_primary = registered_settings(a)
    mism = environment_mismatch()
    if mism:
        raise SystemExit("run refused before reading any data: environment differs from requirements.txt %s" % mism)
    hdr = fits.getheader(a.lc, 0)
    if not a.dry_run and (int(hdr.get("SECTOR", -1)) != TARGET_SECTOR or hdr.get("TICID") not in (TARGET_TIC, str(TARGET_TIC))):
        raise SystemExit("registered run refused before reading any data: input is not Sector 110 / TIC 4206066 "
                         "(use --dry-run for archival files)")
    if a.dry_run and int(hdr.get("SECTOR", -1)) == TARGET_SECTOR:
        raise SystemExit("--dry-run is for archival files only")
    ptype = [k for k, pat in FILE_PATTERNS.items() if pat.match(os.path.basename(a.lc))]
    if not a.dry_run:
        allowed = ("spoc_2min", "tessspoc") if a.product == "spoc" else ("qlp",)
        if not ptype or ptype[0] not in allowed:
            raise SystemExit("registered run refused before reading any data: file name is not a registered Sector 110 "
                             "product for --product %s" % a.product)
    eph = {k: dict(v) for k, v in EPH.items()}
    for k in eph:
        eph[k]["T0"] += a.t0_shift_frac * eph[k]["P"]
    exp = json.load(open(os.path.join(HERE, "calib", "expected_depths.json")))["expected"]
    out = dict(script="s110_final.py (final pre-registration, 2026-10-02)", args=vars(a), dry_run=a.dry_run,
               registered_settings=is_primary, product_type=(ptype[0] if ptype else None),
               provenance=provenance(os.path.abspath(__file__), [a.lc]),
               sector=int(hdr.get("SECTOR", -1)), tic=hdr.get("TICID", None), ephemeris_used=eph, expected=exp,
               n_seeds=N_SEEDS, results={})
    try:
        lc, desat_t = load(a.lc, a.product, a.flux)
    except NoQLPColumn as err:
        out["results"] = {k: dict(n_counted=0, outcome="O4_INSUFFICIENT", reason=O4_REASONS[1], detail=str(err))
                          for k in ("A", "B")}
        json.dump(out, open(a.out, "w"), indent=1)
        print("A O4_INSUFFICIENT; B O4_INSUFFICIENT (%s)" % O4_REASONS[1])
        return
    out.update(file=lc["file"], flux_column=lc["flux_column"])
    eph = {k: dict(v) for k, v in EPH.items()}
    for k in eph:
        eph[k]["T0"] += a.t0_shift_frac * eph[k]["P"]
    lc = clean(a, lc, desat_t, EPH)
    msk_all = lcio.in_transit(lc["t"], "A", 1.5, eph) | lcio.in_transit(lc["t"], "B", 1.5, eph)
    if len(lc["t"]) < 10 or (~msk_all).sum() < 10:
        out["results"] = {k: dict(n_counted=0, outcome="O4_INSUFFICIENT", reason=O4_REASONS[0]) for k in ("A", "B")}
        json.dump(out, open(a.out, "w"), indent=1)
        print("A O4_INSUFFICIENT; B O4_INSUFFICIENT (%s)" % O4_REASONS[0])
        return
    segs_lc = wotan_segments(lc["t"], eph)       # the segments the real detrending actually uses
    d = lcio.detrend(lc, window=a.window, eph=eph)
    t, fd, cad = d["t"], d["fd"], lc["cad"]
    out.update(cadence_s=cad * 86400, n_cadences=int(len(t)))
    for k in ("A", "B"):
        e = eph[k]
        T14 = e["T14"]
        list_frozen = a.t0_shift_frac == 0 and any(t.min() <= tc <= t.max() for _, tc in FROZEN_EPOCHS[k])
        rows = events(t, fd, cad, k, e, a.count_rule, list_frozen, segs_lc, a.window)
        cnt = [r for r in rows if r["counted"]]
        n = len(cnt)
        res = dict(predicted=rows, n_counted=n, n_not_tested=int(len(rows) - n),
                   edge_zone_window_d=a.window / 2, segments=[list(map(float, x)) for x in zip(*segs_lc)])
        if n == 0:
            res.update(outcome="O4_INSUFFICIENT", reason=O4_REASONS[2])
            out["results"][k] = res
            continue
        stat = float(np.mean([r["depth_ppm"] for r in cnt]))
        W, segs, sets = null_pool(lc, k, eph, a.window, a.count_rule)
        res.update(statistic_ppm=stat, n_pool_windows=int(len(W)))
        if len(W) < MIN_POOL:
            res.update(outcome="O4_INSUFFICIENT", reason=O4_REASONS[3], n_pool_windows=int(len(W)))
            out["results"][k] = res
            continue
        D, De = exp[k]["primary_S6_S32"]
        D98, De98 = exp[k]["upper_S98"]
        per_seed = [thresholds(np.random.default_rng(SEED0 + i), W, n, D, De, D98, De98) for i in range(N_SEEDS)]
        med = {q: float(np.median([p[q] for p in per_seed])) for q in ("det", "lo", "hi", "sens", "null_sd")}
        spread = {q: [float(np.min([p[q] for p in per_seed])), float(np.max([p[q] for p in per_seed]))]
                  for q in ("det", "lo", "hi", "sens")}
        det, lo, hi, sens = med["det"], med["lo"], med["hi"], med["sens"]
        timing = None
        if stat > det:
            tcs = [r["tc"] for r in cnt]
            dt_hat = best_shift(t, fd, tcs, T14)
            s_meas, inj_bias = sigma_meas_inj(np.random.default_rng(SEED_TIM), segs, T14, n, stat)
            sp = [SIGMA_PRED[k].get(r["epoch"]) for r in cnt]
            if all(x is not None for x in sp):
                s_pred, s_src = float(np.mean(sp)), "frozen per-epoch table"
            elif k == "A":   # archival dry runs only: epochs outside the frozen Sector 110 table
                ne = np.mean([r["epoch"] for r in cnt])
                s_pred = float(1440 * np.sqrt(_C[0, 0] + 2 * ne * _C[0, 1] + ne * ne * _C[1, 1]))
                s_src = "covariance at the mean counted epoch (archival dry run)"
            else:
                s_pred, s_src = float(np.mean(list(SIGMA_PRED["B"].values()))), "mean of the frozen B table (archival dry run)"
            s_tot = float(np.hypot(s_pred, s_meas))
            indiv = [dict(epoch=r["epoch"], tc=r["tc"], shift_min=best_shift(t, fd, [r["tc"]], T14)) for r in cnt]
            timing = dict(dt_min=dt_hat, sigma_meas_min=s_meas, injection_mean_shift_min=inj_bias,
                          sigma_pred_min=s_pred, sigma_pred_source=s_src, sigma_total_min=s_tot,
                          tolerance_min=3 * s_tot, consistent=bool(abs(dt_hat) <= 3 * s_tot), individual=indiv)
        oc, reason = outcome(stat, det, lo, hi, sens, timing)
        res.update(detect_above_ppm=det, consistent_band_ppm=[lo, hi], expected_primary_ppm=[D, De],
                   expected_S98_ppm=[D98, De98], P_detect_if_real=sens, sensitivity_ok=bool(sens >= SENS_MIN),
                   null_mean_sd_ppm=med["null_sd"],
                   null_window_sigma_mad_ppm=float(1.4826 * np.median(np.abs(W - np.median(W)))),
                   full_set_phase_null=dict(n_phases=len(sets), n_phases_d_ge_observed=int(sum(x["d"] >= stat for x in sets)),
                                            frac_d_ge_observed=float(np.mean([x["d"] >= stat for x in sets])) if sets else None,
                                            n_phases_d_gt_det=int(sum(x["d"] > det for x in sets)),
                                            note="supplementary; fake phases in steps of T14/2 overlap, so they resolve only a few per cent"),
                   seed_spread=spread, timing=timing, outcome=oc, reason=reason)
        out["results"][k] = res
    json.dump(out, open(a.out, "w"), indent=1)
    for k, r in out["results"].items():
        if "detect_above_ppm" not in r:
            print(k, "n=%d" % r["n_counted"], r["outcome"], "(%s)" % r.get("reason", ""))
            continue
        tm = r.get("timing")
        fp = r.get("full_set_phase_null") or {}
        print(k, "n=%d not_tested=%d pool=%s phases d>=obs %s/%s d>det %s" % (r["n_counted"], r["n_not_tested"], r.get("n_pool_windows"),
              fp.get("n_phases_d_ge_observed"), fp.get("n_phases"), fp.get("n_phases_d_gt_det")), r.get("outcome"),
              "d=%.0f det>%.0f band=%s sens=%.2f" % (r.get("statistic_ppm", np.nan), r.get("detect_above_ppm", np.nan),
                                                    [round(x) for x in r.get("consistent_band_ppm", [np.nan, np.nan])],
                                                    r.get("P_detect_if_real", np.nan)),
              "det spread %s" % ([round(x) for x in r["seed_spread"]["det"]] if "seed_spread" in r else None),
              "timing=%s" % (None if not tm else "%+.1f min, tol %.1f (meas %.1f, pred %.1f)" % (
                  tm["dt_min"], tm["tolerance_min"], tm["sigma_meas_min"], tm["sigma_pred_min"])))


def calibrate(a):
    """Expected depths with the test statistic itself (counted transits, event mean) on the archival sectors.
    Error of a sector value = SD of the n-window mean of that sector's null pool (median over seeds)."""
    files = [os.path.join(MAST, fn) for fn, _, _ in ARCHIVAL.values()]
    out = dict(script="s110_final.py --calibrate", provenance=provenance(os.path.abspath(__file__), files),
               ephemeris=EPH, per_sector={}, expected={})
    for s, (fn, product, flux) in ARCHIVAL.items():
        lc, desat_t = load(os.path.join(MAST, fn), product, flux)
        lc = clean(a, lc, desat_t, EPH)
        segs_lc = wotan_segments(lc["t"], EPH)
        d = lcio.detrend(lc, window=a.window, eph=EPH)
        t, fd, cad = d["t"], d["fd"], lc["cad"]
        out["per_sector"][s] = dict(file=fn, product=product, flux=flux)
        for k in ("A", "B"):
            rows = events(t, fd, cad, k, EPH[k], "final", False, segs_lc, a.window)
            cnt = [r for r in rows if r["counted"]]
            n = len(cnt)
            dep = float(np.mean([r["depth_ppm"] for r in cnt])) if n else None
            W, _, _ = null_pool(lc, k, EPH, a.window, "final")
            sds = []
            for i in range(N_SEEDS if n else 0):
                rng = np.random.default_rng(SEED0 + i)
                sds.append(float(np.std(W[rng.integers(0, len(W), (N_BOOT, n))].mean(1))))
            out["per_sector"][s][k] = dict(depth_ppm=dep, err_ppm=float(np.median(sds)) if n else None, n_counted=n,
                                           n_not_tested=int(len(rows) - n), n_pool_windows=int(len(W)),
                                           events=rows)
            print(s, k, "n=%d not_tested=%d pool=%d depth %s +- %s ppm" % (n, out["per_sector"][s][k]["n_not_tested"], len(W),
                  None if dep is None else round(dep, 1), None if not n else round(float(np.median(sds)), 1)))
    for k in ("A", "B"):
        use = [s for s in ("S6", "S32") if out["per_sector"][s][k]["n_counted"] > 0]   # a sector with no counted transit drops out
        dd = np.array([out["per_sector"][s][k]["depth_ppm"] for s in use])
        ee = np.array([out["per_sector"][s][k]["err_ppm"] for s in use])
        w = 1 / ee ** 2
        out["expected"][k] = dict(primary_S6_S32=[float(np.sum(w * dd) / w.sum()), float(1 / np.sqrt(w.sum()))],
                                  primary_sectors_used=use,
                                  upper_S98=[out["per_sector"]["S98"][k]["depth_ppm"], out["per_sector"]["S98"][k]["err_ppm"]])
        print(k, "expected (S6+S32) %.1f +- %.1f; S98 %.1f +- %.1f" % (*out["expected"][k]["primary_S6_S32"],
                                                                        *out["expected"][k]["upper_S98"]))
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--lc")
    p.add_argument("--product", choices=["spoc", "qlp"])
    p.add_argument("--flux", default=None)
    p.add_argument("--window", type=float, default=0.75)
    p.add_argument("--drop-desat-hours", type=float, default=0.0)
    p.add_argument("--count-rule", choices=["final", "v2"], default="final")
    p.add_argument("--t0-shift-frac", type=float, default=0.0)
    p.add_argument("--cut-true", action="store_true")
    p.add_argument("--calibrate", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="archival files only: lifts the Sector 110 / TIC check")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    if a.calibrate:
        if environment_mismatch():
            raise SystemExit("calibration refused: environment differs from requirements.txt %s" % environment_mismatch())
        calibrate(a)
    else:
        if not a.lc or not a.product:
            p.error("--lc and --product are required unless --calibrate")
        if a.flux is None:
            a.flux = "PDCSAP_FLUX" if a.product == "spoc" else "QLP_DEFAULT"
        run(a)
