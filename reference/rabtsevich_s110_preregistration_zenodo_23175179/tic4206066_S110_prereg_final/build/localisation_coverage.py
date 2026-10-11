"""Coverage check of the empirical localisation p (added after audit 5), on the Sector 32 dry run of signal A.

The dry-run localisation is run once more, unchanged except for the seed of the calibration cloud (5005 instead of the
registered 4004); the fit cloud (seed 1001, which defines the metric) is the same. Every trial of the new 3000-trial
on-target cloud is then treated as an observation of a signal that truly sits on the target, and its empirical p is
computed against the 3000 calibration distances of the registered seed (dryrun/loc_s32_final_pipeline.json), and the
other way round. For a calibrated p about 5 % of on-target signals fall below 0.05 and about 0.27 % below 0.0027.
The Gaussian tail exp(-d^2/2) of the earlier rule is counted on the same trials for comparison.
All trials use the same Sector 32 pixels and overlapping windows, so these fractions check the method on this data
set; they are not independent binomial counts and not a measured rate for Sector 110.

Usage (from the deposit root, as in build/run_all.sh):
  python build/localisation_coverage.py --tpf S32_TPF --prf-dir DIR --photometry dryrun/s32_real.json \
      --reference dryrun/loc_s32_final_pipeline.json --out-loc dryrun/loc_s32_coverage_seed5005.json \
      --out calib/localisation_coverage.json
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "deps"))
import pixloc_lib as L  # noqa: E402
import s110_localise as SL  # noqa: E402
from s110_final import sha256  # noqa: E402

SEED_CHECK = 5005
LEVELS = (0.05, 0.0027)


def dists(cloud, C):
    return np.array([L.mahalanobis([t["dE"], t["dN"]], (0, 0), C) for t in cloud["trials"]])


def counts(p):
    p = np.asarray(p)
    return {"n": int(p.size), **{"below_%g" % q: int(np.sum(p < q)) for q in LEVELS},
            **{"percent_below_%g" % q: float(100 * np.mean(p < q)) for q in LEVELS}}


def main(a):
    SL.SEED_CAL = SEED_CHECK
    ns = argparse.Namespace(tpf=a.tpf, photometry=a.photometry, signal="A", prf_dir=a.prf_dir, ephem="final",
                            events="photometry", n_trials=SL.N_TRIALS, n_cal=SL.N_CAL, no_clear=False, dry_run=True,
                            primary_sha256=None, out=a.out_loc, coverage_check_seed_cal=SEED_CHECK)
    SL.main(ns)
    ref, new = json.load(open(a.reference)), json.load(open(a.out_loc))
    rs, ns_ = ref["injections"]["sets"], new["injections"]["sets"]
    if rs["target"]["trials"] != ns_["target"]["trials"]:
        raise SystemExit("fit clouds differ: the check run is not the registered dry run with another calibration seed")
    C = np.array(ref["stats"]["inj_cov_target"])
    d_ref, d_new = dists(rs["target_cal"], C), dists(ns_["target_cal"], C)
    p_new = [SL.empirical_p(d, d_ref) for d in d_new]
    p_ref = [SL.empirical_p(d, d_new) for d in d_ref]
    out = dict(script="build/localisation_coverage.py", script_sha256=sha256(os.path.abspath(__file__)),
               reference=os.path.basename(a.reference), reference_sha256=sha256(a.reference),
               check_run=os.path.basename(a.out_loc), check_run_sha256=sha256(a.out_loc),
               seeds=dict(fit=SL.SEED_T, calibration=4004, check=SEED_CHECK),
               empirical_p_check_cloud_vs_registered_calibration=counts(p_new),
               empirical_p_registered_calibration_vs_check_cloud=counts(p_ref),
               gaussian_tail_check_cloud=counts(np.exp(-0.5 * d_new ** 2)),
               gaussian_tail_registered_calibration=counts(np.exp(-0.5 * d_ref ** 2)),
               note=("on-target trials only: every count below a threshold is a false 'tension' or 'off-target' for a "
                     "signal that is on the target; trials share pixels and overlapping windows (see docstring)"))
    json.dump(out, open(a.out, "w"), indent=1)
    for k in ("empirical_p_check_cloud_vs_registered_calibration", "empirical_p_registered_calibration_vs_check_cloud",
              "gaussian_tail_check_cloud", "gaussian_tail_registered_calibration"):
        c = out[k]
        print("%-52s n=%d  <0.05: %d (%.2f %%)  <0.0027: %d (%.2f %%)" % (
            k, c["n"], c["below_0.05"], c["percent_below_0.05"], c["below_0.0027"], c["percent_below_0.0027"]))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--tpf", required=True)
    p.add_argument("--prf-dir", required=True)
    p.add_argument("--photometry", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--out-loc", required=True)
    p.add_argument("--out", required=True)
    main(p.parse_args())
