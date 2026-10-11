"""Modes A, B (all three archival products), A_S32, B_S32 (SPOC 2-min only, added after audit 3) and A_S32_a3, B_S32_a3
(the SPOC 2-min trial set run by audit 3; added after audit 4 so that every trial set that was run is included).
Re-measure the realised false-alarm and lower-tail levels with the FINAL code (audit-2 trial script, adapted):
trial ephemerides with arbitrary periods on archival sectors with the true A/B transits removed; pool exactly as
s110_final.null_pool (segments of each fake detrending); leave-out trials: each fake phase i is a trial, thresholds
from the pool without phases |j - i| <= 2. Same jobs and seeds as audit 2 (rng 777). Since audit 6 the "final" thresholds are the median of N_SEEDS seeds, computed by
s110_final.thresholds itself, as in the S110 analysis (before audit 6 one seed SEED0 was used here)."""
import json, os, sys, time
import numpy as np
PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the deposit folder
sys.path.insert(0, PKG); sys.path.insert(0, os.path.join(PKG, "deps"))
sys.argv_save = list(sys.argv)
import s110_final as S
import lcio
EXP = json.load(open(os.path.join(PKG, "calib", "expected_depths.json")))["expected"]

def load_cut(sec):
    fn, product, flux = S.ARCHIVAL[sec]
    lc, _ = S.load(os.path.join(S.MAST, fn), product, flux)
    keep = ~(lcio.in_transit(lc["t"], "A", 1.5, S.EPH) | lcio.in_transit(lc["t"], "B", 1.5, S.EPH))
    return {**lc, "t": lc["t"][keep], "f": lc["f"][keep], "q": lc["q"][keep]}

def pool_with_phase(lc, k, eph, window):
    t = lc["t"]
    cut = lcio.in_transit(t, "A", S.CUT_FAC, eph) | lcio.in_transit(t, "B", S.CUT_FAC, eph)
    lcc = {**lc, "t": t[~cut], "f": lc["f"][~cut], "q": lc["q"][~cut]}
    e = eph[k]; step = S.PHI_STEP_T14 * e["T14"] / e["P"]
    dep, ph = [], []
    for i, phi in enumerate(np.arange(S.PHI_MIN, S.PHI_MAX + 1e-9, step)):
        fake = {kk: dict(v) for kk, v in eph.items()}; fake[k]["T0"] = e["T0"] + phi * e["P"]
        d = lcio.detrend(lcc, window=window, eph=fake)
        segs = S.wotan_segments(lcc["t"], fake)
        for r in S.events(d["t"], d["fd"], lc["cad"], k, fake[k], "final", False, segs, window):
            if r["counted"] and r["depth_ppm"] is not None:
                dep.append(r["depth_ppm"]); ph.append(i)
    return np.array(dep), np.array(ph)

def one(args):
    sec, k, P, T0frac, window = args
    lc = load_cut(sec)
    eph = {kk: dict(v) for kk, v in S.EPH.items()}
    eph[k]["P"] = P; eph[k]["T0"] = lc["t"].min() + T0frac * P
    W, ph = pool_with_phase(lc, k, eph, window)
    D, De = EXP[k]["primary_S6_S32"]
    D98, De98 = EXP[k]["upper_S98"]
    rng_p = np.random.default_rng(int(P * 1e6) % 2**31)
    out = []
    for i in np.unique(ph):
        own = W[ph == i]; Wl = W[np.abs(ph - i) > 2]; n = len(own)
        if len(Wl) < S.MIN_POOL:
            continue
        row = dict(i=int(i), n=n, d=float(own.mean()))
        # "final": exactly the S110 procedure (s110_final.thresholds, median over N_SEEDS seeds SEED0 + s; audit 6).
        # "old": the pre-audit-2 percentiles 99.7 / 0.3 with the single seed SEED0, kept only as a historical comparison.
        per_seed = [S.thresholds(np.random.default_rng(S.SEED0 + s), Wl, n, D, De, D98, De98) for s in range(S.N_SEEDS)]
        med = {q: float(np.median([p[q] for p in per_seed])) for q in ("det", "lo", "sens")}
        r = np.random.default_rng(S.SEED0)
        null = Wl[r.integers(0, len(Wl), (S.N_BOOT, n))].mean(1)
        sig = null + r.normal(D, De, S.N_BOOT)
        old_det = np.percentile(null, 99.7)
        cases = (("final", med["det"], med["lo"], med["sens"]),
                 ("old", old_det, np.percentile(sig, 0.3), float(np.mean(sig > old_det))))
        for tag, det, lo, sens in cases:
            dd = own.mean() + rng_p.normal(D, De, 200)
            row[tag] = dict(fa=bool(own.mean() > det), below_lo=float(np.mean(dd < lo)), sens=sens, det=float(det),
                            false_O3=float(np.mean((dd < lo) & (dd <= det))) if sens >= S.SENS_MIN else None)
        out.append(row)
    return dict(sec=sec, k=k, P=P, pool=len(W), trials=out)

if __name__ == "__main__":
    from multiprocessing import Pool
    mode, outp, nproc = sys.argv_save[1], sys.argv_save[2], int(sys.argv_save[3])
    rng = np.random.default_rng(777); jobs = []
    if mode == "A":
        for sec in ("S32", "S98", "S6"):
            for P in np.linspace(2.35, 4.45, 10):
                jobs.append((sec, "A", float(P * (1 + rng.uniform(-0.01, 0.01))), float(rng.uniform()), 0.75))
    elif mode == "B":
        for sec in ("S32", "S98"):
            for P in (7.3, 8.9, 10.2, 12.6, 14.1):
                jobs.append((sec, "B", float(P), float(rng.uniform()), 0.75))
    elif mode == "A_S32":   # added after audit 3: more trials on the SPOC 2-min product (the S110 primary product)
        rng = np.random.default_rng(2026)
        for P in np.linspace(2.35, 4.45, 40):
            jobs.append(("S32", "A", float(P * (1 + rng.uniform(-0.01, 0.01))), float(rng.uniform()), 0.75))
    elif mode == "B_S32":
        rng = np.random.default_rng(2027)
        for P in np.linspace(7.0, 14.5, 24):
            jobs.append(("S32", "B", float(P * (1 + rng.uniform(-0.01, 0.01))), float(rng.uniform()), 0.75))
    elif mode == "A_S32_a3":   # the trial set run by audit 3 (rng 31337), same one() function
        rng = np.random.default_rng(31337)
        for P in np.linspace(2.30, 4.50, 40):
            jobs.append(("S32", "A", float(P * (1 + rng.uniform(-0.01, 0.01))), float(rng.uniform()), 0.75))
    elif mode == "B_S32_a3":
        rng = np.random.default_rng(31337)
        for P in np.linspace(6.5, 14.5, 24):
            jobs.append(("S32", "B", float(P * (1 + rng.uniform(-0.01, 0.01))), float(rng.uniform()), 0.75))
    t0 = time.time(); R = []
    with Pool(nproc) as p:
        for r in p.imap_unordered(one, jobs):
            R.append(r); json.dump(R, open(outp, "w"))
            print(r["sec"], round(r["P"], 3), len(r["trials"]), round(time.time() - t0), flush=True)
    rng_b = np.random.default_rng(1)
    for tag in ("old", "final"):
        fa = np.array([t[tag]["fa"] for r in R for t in r["trials"]], float)
        bl = np.array([t[tag]["below_lo"] for r in R for t in r["trials"]])
        o3 = [t[tag]["false_O3"] for r in R for t in r["trials"] if t[tag]["false_O3"] is not None]
        # cluster bootstrap over trial ephemerides
        cl = [np.array([t[tag]["fa"] for t in r["trials"]], float) for r in R if r["trials"]]
        bs = []
        for _ in range(4000):
            pick = rng_b.integers(0, len(cl), len(cl)); x = np.concatenate([cl[j] for j in pick]); bs.append(x.mean())
        print(mode, tag, "trials %d  false alarm %.2f%% (95%% CI %.2f-%.2f)  below lo %.2f%%  false O3 %.2f%% (%d)" % (
            len(fa), 100 * fa.mean(), 100 * np.percentile(bs, 2.5), 100 * np.percentile(bs, 97.5), 100 * bl.mean(),
            100 * np.mean(o3) if o3 else np.nan, len(o3)))
