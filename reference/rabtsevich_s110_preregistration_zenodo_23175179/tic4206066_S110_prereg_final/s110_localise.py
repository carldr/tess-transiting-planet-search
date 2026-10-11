"""TIC 4206066, Sector 110: pixel localisation, FINAL pre-registration (written 2026-10-02, before any Sector 110 data
exist). See PREREG_S110_FINAL.md section 6.

Method = the headline ("main") recipe of the final re-derivation, items 10 and 11, run through the frozen library
deps/pixloc_lib.py (the code of the re-derivation's library; only its comments and default file paths were made
neutral on 2026-10-05, sha256 in SHA256SUMS). The library is not otherwise edited; this wrapper only
  * registers the target pixel file under the library's SPOC-TPF branch (pixloc_lib treats sector 32 as "SPOC 2-min
    TPF": FLUX already background-subtracted, SPOC optimal aperture, POS_CORR present; the Sector 110 product is the
    same kind of file). The PRF model depends on the sector only through "sector < 4", so the label does not change
    it; camera and CCD are read from the file;
  * sets the library's ephemeris to the frozen final ephemeris (ephemeris_final.json);
  * takes the transits from the PRIMARY photometric run (counted transits of s110_final.py); the library's own
    usability rule (>= 80 % in transit, out-of-transit windows 1.5..2.5 T14 on both sides) is applied on top, and any
    transit it drops is listed. As in item 10, the POS_CORR decorrelation is fitted away from the localised signal's
    own predicted transits;
  * adds a calibration cloud of 3000 on-target injections (seed 4004, same recipe) for the empirical p (audit 5);
  * repeats the item 11 injections (300 on target, 300 at +5" E; real-image template; seeds 1001 / 2002; phase
    U(0.12, 0.88); POS_CORR decorrelation refitted per trial; depth = the pixel-aperture box depth of the real
    transits, as in item 11) with one pre-declared change: a fake transit is dropped from its trial if the other
    signal has an event within 2.5 T14 + T14_other / 2 of it (item 11 B-contamination check, "any window").
Statistics: offset x = (dE, dN) of the PRF fit to the stacked difference image; d = Mahalanobis distance of x from the
target position (0, 0) with the covariance of the 300-trial on-target cloud (as in item 10). The label uses an
EMPIRICAL p (revised after audit 5, which showed that the Gaussian tail exp(-d^2/2) is not calibrated: 12 of the 300
on-target dry-run trials had exp(-d^2/2) < 0.0027): a further N_CAL = 3000 on-target injections (seed 4004, same
recipe) give the calibration distances d_i, measured in the same metric, and p = (1 + #{d_i >= d}) / (N_CAL + 1)
(split-conformal / Monte Carlo p-value: Davison & Hinkley 1997; North et al. 2002, AJHG 71, 439; Bates et al. 2023,
Ann. Stat. 51, 149). exp(-d^2/2) and the p from the cloud mean are reported for information only.
L(5" E)/L(target) = ratio of the Gaussian densities of the two clouds at x (a ratio, not a probability).
Labels (signal A, Sector 110 alone): L1 compatible p >= 0.05; L2 tension 0.0027 <= p < 0.05; L3 off-target p < 0.0027.
"Not possible" if: the signal's photometric outcome is not O1/O2; POS_CORR1/2 missing, all zero, or non-finite on more
than 50 % of cadences; no transit passes the pixel rule; fewer than 90 % of the trials succeed in any cloud (270 of
300, 2700 of 3000); or d is not finite. Combined S32 + S98 + S110: inverse-covariance mean with
calib/localisation_archival.json, descriptive only.

A labelled run refuses to start unless: the photometry JSON comes from the primary settings, its sha256 equals
--primary-sha256 (the publicly posted hash of the primary output), its input is the SPOC 2-min light curve of Sector 110
(file name tessNNNN-s0110-0000000004206066-NNNN-s_lc.fits), the TPF is the matching SPOC 2-min target pixel file,
both are TIC 4206066, the PRF folder cam1_ccd2 matches PRF_SHA256SUMS, --n-trials 300, --n-cal 3000, --ephem final,
--events photometry, clearing on, and the installed package versions equal requirements.txt (as for the photometry).
The environment is checked first, before any file is read; all other checks come before any rule that writes a
result (revised after audit 5). Under a fallback primary product, localisation is not possible. --dry-run is for archival files only (refused on Sector 110 files); it lifts the
sector checks (never the environment check) and suppresses every label.
Equal stacks (labelled run and the dry run that mimics it): every fake stack has exactly as many transits as the real
stack (n_real). A candidate fake must (i) not lie within 2.5 T14 + T14_other / 2 of the other signal, (ii) pass the
pixel usability rule, (iii) have its T14 window inside one segment of the primary photometric run and at least W/2
from its ends (the photometric edge-zone rule, segments read from the photometry JSON), and (iv) lie more than 3 T14
from every fake already chosen in that trial. For each trial a phase is drawn and its candidates are taken; while
fewer than n_real are collected, further phases are drawn and their candidates added (no phase is rejected); then
n_real of the collected candidates are chosen at random. At most MAX_DRAWS phases per trial, after which the trial
fails. Each cloud records the number of distinct T14-wide phase slots of its fakes; fewer than MIN_SLOTS in either
cloud (including the calibration cloud), or fewer than 90 % successful trials, makes the localisation "not possible".
The TPF must be on camera 1, CCD 2 (the PRF files of PRF_SHA256SUMS); otherwise a labelled run that has passed every
check above gives "not possible".

Usage:
  python s110_localise.py --tpf S110_TPF.fits --photometry s110_primary.json --primary-sha256 HASH --prf-dir DIR --out loc_A.json
  python s110_localise.py --tpf S110_TPF.fits --photometry s110_primary.json --primary-sha256 HASH --prf-dir DIR --signal B --out loc_B.json
"""
import argparse
import json
import os
import re
import sys
import time

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "deps"))
import pixloc_lib as L  # noqa: E402  (frozen, unedited)
from s110_final import provenance, sha256, EPH, PRIMARY_ARGS, environment_mismatch  # noqa: E402

N_TRIALS, SEED_T, SEED_E5 = 300, 1001, 2002
N_CAL, SEED_CAL = 3000, 4004   # on-target calibration cloud of the empirical p (audit 5)
MAX_DRAWS = 100            # phase draws per trial before the trial fails (equal stacks)
MIN_SLOTS = 15             # distinct T14-wide phase slots per cloud
SEP_T14 = 3.0              # minimum separation of two fakes of one trial, in T14
MIN_OK_FRAC = 0.9          # >= 90 % successful trials per cloud (270 of 300, 2700 of 3000)
POSCORR_MIN_FRAC = 0.5     # POS_CORR must be finite and non-zero on at least half of the cadences
TARGET_TIC, TARGET_SECTOR = 4206066, 110
LC_PATTERN = re.compile(r"^tess\d+-s0110-0000000004206066-\d+-s_lc\.fits$")
TP_PATTERN = re.compile(r"^tess\d+-s0110-0000000004206066-\d+-s_tp\.fits$")


def other_events(sd, other):
    """Mid-times of the other signal with at least one good cadence within T14_other/2 (as in the archival analysis)."""
    t = sd.time[sd.good]
    m = np.arange(np.floor((t.min() - other["T0"]) / other["P"]) - 1, np.ceil((t.max() - other["T0"]) / other["P"]) + 2)
    tb = other["T0"] + m * other["P"]
    return np.array([x for x in tb if np.any(np.abs(t - x) < other["T14"] / 2)])


def in_clear_zone(tc, segs, half_w):
    lo, hi = tc - L.T14 / 2, tc + L.T14 / 2
    return any(lo >= a + half_w and hi <= b - half_w for a, b in segs)


def run_trials(sd, flux_raw, tmpl_fn, depth, dE, dN, n_trials, seed, tc_real, t_other, T14_other, clear, n_real=None,
               segs=None, half_w=None):
    """n_real None: item-11 recipe (fidelity check). n_real given: equal stacks (module docstring)."""
    rng = np.random.default_rng(seed)
    dx, dy = sd.offset_to_pixshift(dE, dN)
    tmpl = tmpl_fn(dx, dy)
    trials, n_dropped, n_failed, t0 = [], 0, 0, time.time()
    for _ in range(n_trials):
        if n_real is None:
            phi = rng.uniform(0.12, 0.88)
            tc_fake = sd.transit_times(phase_shift=phi)
            if clear and t_other.size:
                bad = np.any(np.abs(t_other[:, None] - tc_fake[None, :]) < 2.5 * L.T14 + T14_other / 2, axis=0)
                n_dropped += int(bad.sum())
                tc_fake = tc_fake[~bad]
            phis = [phi]
        else:
            chosen, phis = [], []
            while len(chosen) < n_real and len(phis) < MAX_DRAWS:
                phi = rng.uniform(0.12, 0.88)
                phis.append(phi)
                cand = sd.transit_times(phase_shift=phi)
                if clear and t_other.size:
                    bad = np.any(np.abs(t_other[:, None] - cand[None, :]) < 2.5 * L.T14 + T14_other / 2, axis=0)
                    n_dropped += int(bad.sum())
                    cand = cand[~bad]
                for u in sd.usable_transits(cand):
                    tc = u["tc"]
                    if u["ok"] and in_clear_zone(tc, segs, half_w) and all(abs(tc - c) > SEP_T14 * L.T14 for c in chosen):
                        chosen.append(tc)
            if len(chosen) < n_real:
                n_failed += 1
                continue
            chosen = np.array(chosen)
            tc_fake = np.sort(chosen[rng.choice(len(chosen), n_real, replace=False)])
            phi = phis[0]
        F = L.inject(sd, flux_raw, tmpl, depth, tc_fake)
        F = sd.decorrelated_flux(np.r_[tc_real, tc_fake], kind="poscorr", flux=F)
        st = L.stack_difference(sd, tc_fake, flux=F)
        if st is None:
            continue
        r = L.fit_prf(sd, st["diff"], st["sig"])
        e, n = sd.pix_to_offset(r["x"], r["y"])
        slots = sorted({int(np.floor((((tc - L.T0) / L.P) % 1.0) / (L.T14 / L.P))) for tc in tc_fake})
        trials.append(dict(phi=float(phi), n_phases=len(phis), n_events=st["n_transits"], dE=float(e), dN=float(n),
                           slots=slots))
    all_slots = sorted({s_ for t in trials for s_ in t["slots"]})
    return dict(true_dE=float(dE), true_dN=float(dN), n_requested=n_trials, n_trials=len(trials), n_real=n_real,
                n_failed_draw_cap=n_failed, n_fake_dropped_other_signal=n_dropped, distinct_phase_slots=len(all_slots),
                runtime_s=time.time() - t0, trials=trials)


def empirical_p(d, d_cal):
    """(1 + number of calibration distances >= d) / (N + 1); nan if d is not finite."""
    if not np.isfinite(d):
        return float("nan")
    return float((1 + int(np.sum(np.asarray(d_cal) >= d))) / (len(d_cal) + 1))


def hyp_stats(x, ct, c5, ccal):
    mt, Ct = L.cov_stats(ct)
    m5, C5 = L.cov_stats(c5)
    d = L.mahalanobis(x, (0, 0), Ct)
    d_cal = [L.mahalanobis(c, (0, 0), Ct) for c in ccal]
    lr = np.exp(L.gauss_loglike(x, m5, C5) - L.gauss_loglike(x, mt, Ct))
    return dict(x=list(map(float, x)), inj_mean_target=mt.tolist(), inj_cov_target=Ct.tolist(),
                inj_sd_target=np.sqrt(np.diag(Ct)).tolist(), inj_mean_E5=m5.tolist(), inj_cov_E5=C5.tolist(),
                d_raw=d, p_emp=empirical_p(d, d_cal), n_cal=len(d_cal), n_cal_at_least_d=int(np.sum(np.asarray(d_cal) >= d)),
                p_gauss_info_only=L.chi2_2dof_p(d), d_bc_info_only=L.mahalanobis(x, mt, Ct),
                p_bc_info_only=L.chi2_2dof_p(L.mahalanobis(x, mt, Ct)), LR_E5_over_target=float(lr))


def label(p):
    if not np.isfinite(p):
        return "NOT_POSSIBLE"
    return "L1_COMPATIBLE" if p >= 0.05 else ("L2_TENSION" if p >= 0.0027 else "L3_OFF_TARGET")


def finish(out, a, msg):
    out["result"] = msg
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print(msg)


def main(a):
    sig, oth = a.signal, ("B" if a.signal == "A" else "A")
    mism = environment_mismatch()
    if mism:
        raise SystemExit("run refused before reading any data: environment differs from requirements.txt %s" % mism)
    hdr0 = fits.getheader(a.tpf, 0)
    out = dict(script="s110_localise.py (final pre-registration, 2026-10-02, revised after audit 5)", args=vars(a),
               dry_run=a.dry_run,
               provenance=provenance(os.path.abspath(__file__), [a.tpf] + ([a.photometry] if a.photometry else [])),
               file_sector=int(hdr0.get("SECTOR", -1)), file_tic=hdr0.get("TICID", None),
               camera_ccd=[int(hdr0.get("CAMERA", -1)), int(hdr0.get("CCD", -1))])
    if a.dry_run and out["file_sector"] == TARGET_SECTOR:
        raise SystemExit("--dry-run is for archival files only; Sector 110 is localised only by the registered labelled run")
    # PRF folder: the registered one (camera 1, CCD 2) for a labelled run, the file's own one for a dry run
    prf_sub = os.path.join(a.prf_dir, "cam%d_ccd%d" % (tuple(out["camera_ccd"]) if a.dry_run else (1, 2)))
    if not os.path.isdir(prf_sub):
        raise SystemExit("PRF directory %s missing (download MAST prf_fitsfiles/start_s0004/%s/)" % (prf_sub, os.path.basename(prf_sub)))
    out["provenance"]["prf_files_sha256"] = {f: sha256(os.path.join(prf_sub, f)) for f in sorted(os.listdir(prf_sub))
                                             if f.endswith(".fits")}
    phot = json.load(open(a.photometry)) if a.photometry else None
    if not a.dry_run:
        # every check of origin and settings comes before any rule that writes a result (audit 5)
        problems = []
        prf_ref = {}
        for line in open(os.path.join(HERE, "PRF_SHA256SUMS")):
            if line.strip() and not line.startswith("#"):
                h, f = line.split()
                prf_ref[os.path.basename(f)] = h
        have = out["provenance"]["prf_files_sha256"]
        if set(have) != set(prf_ref) or any(prf_ref[f] != h for f, h in have.items()):
            problems.append("PRF files are not exactly the set listed in PRF_SHA256SUMS")
        if not TP_PATTERN.match(os.path.basename(a.tpf)):
            problems.append("TPF is not the SPOC 2-min Sector 110 target pixel file of TIC 4206066")
        if a.primary_sha256 is None or (a.photometry and sha256(a.photometry) != a.primary_sha256):
            problems.append("photometry JSON sha256 does not equal --primary-sha256 (the posted primary hash)")
        if phot is None:
            problems.append("no photometry JSON")
        else:
            if not all(phot["args"].get(q) == v for q, v in PRIMARY_ARGS.items()):
                problems.append("photometry JSON is not from the primary settings")
            if not LC_PATTERN.match(str(phot.get("file", ""))):
                problems.append("photometry input is not the SPOC 2-min Sector 110 light curve (fallback primary: localisation not possible)")
            if phot.get("sector") != TARGET_SECTOR or phot.get("tic") not in (TARGET_TIC, str(TARGET_TIC)):
                problems.append("photometry JSON is not Sector 110 / TIC 4206066")
        if out["file_sector"] != TARGET_SECTOR or out["file_tic"] not in (TARGET_TIC, str(TARGET_TIC)):
            problems.append("TPF is not Sector 110 / TIC 4206066")
        if (a.n_trials != N_TRIALS or a.n_cal != N_CAL or a.ephem != "final" or a.events != "photometry"
                or a.no_clear):
            problems.append("non-registered switches (--n-trials/--n-cal/--ephem/--events/--no-clear)")
        if problems:
            raise SystemExit("labelled run refused: " + "; ".join(problems) + " (use --dry-run for validation)")
        if out["camera_ccd"] != [1, 2]:
            return finish(out, a, "localisation not possible: the TPF is not on camera 1, CCD 2 (no registered PRF files)")
    if phot is not None:
        oc = phot["results"][sig].get("outcome", "")
        out["photometric_outcome"] = oc
        if not (oc.startswith("O1") or oc.startswith("O2")):
            return finish(out, a, "localisation not possible: signal %s not detected (photometric outcome %s)" % (sig, oc))
    # POS_CORR check before the library is called
    d1 = fits.getdata(a.tpf, 1)
    ok_pc = False
    if "POS_CORR1" in d1.columns.names and "POS_CORR2" in d1.columns.names:
        pc1, pc2 = np.asarray(d1["POS_CORR1"], float), np.asarray(d1["POS_CORR2"], float)
        usable = np.isfinite(pc1) & np.isfinite(pc2) & ((pc1 != 0) | (pc2 != 0))
        out["poscorr_usable_fraction"] = float(usable.mean())
        ok_pc = usable.mean() >= POSCORR_MIN_FRAC
    if not ok_pc:
        return finish(out, a, "localisation not possible: POS_CORR missing, zero or non-finite on > 50 % of cadences")
    # library configuration (no edits to the library)
    L.PRF_DIR = a.prf_dir
    if a.ephem == "final":
        L.P, L.T0, L.T14 = EPH[sig]["P"], EPH[sig]["T0"], EPH[sig]["T14"]
    elif sig != "A":
        raise SystemExit("--ephem items10 exists only for signal A")
    other = EPH[oth]
    L.FILES[32] = os.path.abspath(a.tpf)          # SPOC 2-min TPF branch of pixloc_lib
    sd = L.SectorData(32)
    out["ephemeris_used"] = dict(P=L.P, T0=L.T0, T14=L.T14)
    tc_all = sd.transit_times()
    if a.events == "photometry":
        if phot is None:
            raise SystemExit("--events photometry needs --photometry")
        tc_ev = np.array([r["tc"] for r in phot["results"][sig]["predicted"] if r["counted"]])
    else:
        tc_ev = tc_all
    flux = sd.decorrelated_flux(tc_all, kind="poscorr")
    usable = sd.usable_transits(tc_ev)
    out["events"] = dict(given=tc_ev.tolist(), library_rule=usable,
                         dropped_by_library_rule=[u["tc"] for u in usable if not u["ok"]])
    st = L.stack_difference(sd, tc_ev, flux=flux)
    if st is None:
        return finish(out, a, "localisation not possible: no transit passes the pixel usability rule")
    depth = L.box_depth(sd, tc_ev, flux=flux)
    r = L.fit_prf(sd, st["diff"], st["sig"])
    dE, dN = sd.pix_to_offset(r["x"], r["y"])
    star = L.fit_star(sd, st["mean_out"])
    out["real"] = dict(aperture_box_depth=depth, n_transits=st["n_transits"], transits=[t["tc"] for t in st["transits"]],
                       prf_fit=dict(r, dE=float(dE), dN=float(dN)), star_fit=star)
    base = L.stack_difference(sd, tc_all, flux=sd.flux)["mean_out"] - star["c"]

    def tmpl_fn(dx, dy):
        return L.shifted_template(base, dx, dy) if (dx or dy) else base.copy()

    t_oth = other_events(sd, other)
    clear = not a.no_clear
    sets = {}
    n_real = st["n_transits"] if a.events == "photometry" else None   # equal stacks (labelled path and its dry run)
    segs = phot["results"][sig]["segments"] if n_real is not None else None
    half_w = phot["results"][sig]["edge_zone_window_d"] if n_real is not None else None
    out["injections_n_real"] = n_real
    clouds = dict(target=(0.0, SEED_T, a.n_trials), E5=(5.0, SEED_E5, a.n_trials), target_cal=(0.0, SEED_CAL, a.n_cal))
    for name, (e5, seed, n_tr) in clouds.items():
        sets[name] = run_trials(sd, sd.flux, tmpl_fn, depth["depth"], e5, 0.0, n_tr, seed, tc_all,
                                t_oth, other["T14"], clear, n_real, segs, half_w)
        print(name, "done", round(sets[name]["runtime_s"]), "s", flush=True)
    out["injections"] = dict(other_signal_events=t_oth.tolist(), clear_other_signal=clear, sets=sets)
    for name, (_, _, n_tr) in clouds.items():
        n_min = int(np.ceil(MIN_OK_FRAC * n_tr))
        if sets[name]["n_trials"] < n_min:
            return finish(out, a, "localisation not possible: fewer than %d successful trials in the %s cloud" % (n_min, name))
    if n_real is not None and min(sets[k]["distinct_phase_slots"] for k in clouds) < MIN_SLOTS:
        return finish(out, a, "localisation not possible: fewer than %d distinct phase slots in a cloud" % MIN_SLOTS)
    ct, c5, ccal = (np.array([[t["dE"], t["dN"]] for t in sets[k]["trials"]]) for k in ("target", "E5", "target_cal"))
    hs = hyp_stats(np.array([dE, dN]), ct, c5, ccal)
    out["stats"] = hs
    lab = label(hs["p_emp"])
    if a.dry_run:
        out["label_if_this_were_sector_110"] = lab + " (dry run, not a label)"
    elif sig == "A":
        out["label"] = lab
    if sig == "A" and a.ephem == "final" and np.isfinite(hs["d_raw"]):
        arch = json.load(open(os.path.join(HERE, "calib", "localisation_archival.json")))["inputs"]
        xs = [np.array(arch[k]["x"]) for k in ("S32_poscorr", "S98_hard")] + [np.array([dE, dN])]
        Cs = [np.array(arch[k]["inj_cov_target"]) for k in ("S32_poscorr", "S98_hard")] + [np.array(hs["inj_cov_target"])]
        xc, Cc = L.combine(xs, Cs)
        dc = L.mahalanobis(xc, (0, 0), Cc)
        out["combined_S32_S98_S110_descriptive"] = dict(
            x=xc.tolist(), cov=Cc.tolist(), sd=np.sqrt(np.diag(Cc)).tolist(), d=dc, p=L.chi2_2dof_p(dc),
            LR_product_of_sectors=float(arch["S32_poscorr"]["LR_E5_over_target"] * arch["S98_hard"]["LR_E5_over_target"]
                                        * hs["LR_E5_over_target"]),
            note="descriptive; never changes the Sector 110 label" + ("; DRY RUN: meaningless if the input is S32" if a.dry_run else ""))
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    print("offset (%+.2f, %+.2f)\" n=%d d=%.2f p_emp=%.4g LR=%.2f %s" % (dE, dN, st["n_transits"], hs["d_raw"], hs["p_emp"],
                                                                hs["LR_E5_over_target"],
                                                                out.get("label", out.get("label_if_this_were_sector_110", ""))))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--tpf", required=True)
    p.add_argument("--photometry", default=None)
    p.add_argument("--signal", choices=["A", "B"], default="A")
    p.add_argument("--prf-dir", required=True, help="directory with MAST prf_fitsfiles/start_s0004/camC_ccdD/ subfolders")
    p.add_argument("--ephem", choices=["final", "items10"], default="final")
    p.add_argument("--events", choices=["photometry", "own"], default="photometry")
    p.add_argument("--n-trials", type=int, default=N_TRIALS)
    p.add_argument("--n-cal", type=int, default=N_CAL, help="on-target calibration trials of the empirical p")
    p.add_argument("--no-clear", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="archival validation: lifts the checks, suppresses labels")
    p.add_argument("--primary-sha256", default=None, help="publicly posted sha256 of the primary photometry JSON")
    p.add_argument("--out", required=True)
    main(p.parse_args())
