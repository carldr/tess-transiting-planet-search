"""Builds ephemeris_final.json from the final re-derivation outputs (2026-10-01, decision D13). Run once; the JSON is
then frozen. Sources (read-only):
  A: ../final_rederivation/03_signalA/out/linear_ephemeris_calibrated.json  ('adopted': WLS of the 22 item-03
     transit times with injection-calibrated errors) and 03_signalA/out/mcmc_summary.json (T14_h median);
  B: ../final_rederivation/14_signalB/out/signalB.json 'ephemeris_B_used' (= v1/v2 joint-fit ephemeris, unchanged)
     and the per-epoch B sigma of 05_ephemeris/out/s110_predicted_transits.csv (joint_fit_MCMC, 'sigma_min_corr').
"""
import csv, json, os
import numpy as np
from astropy.time import Time

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "final_rederivation")   # the author's re-derivation folder, not in this deposit
lin = json.load(open(os.path.join(R, "03_signalA/out/linear_ephemeris_calibrated.json")))["adopted"]
mc = json.load(open(os.path.join(R, "03_signalA/out/mcmc_summary.json")))
eb = json.load(open(os.path.join(R, "14_signalB/out/signalB.json")))["ephemeris_B_used"]
rows = list(csv.DictReader(open(os.path.join(R, "05_ephemeris/out/s110_predicted_transits.csv"))))

def utc(btjd):
    return Time(btjd + 2457000.0, format="jd", scale="tdb").utc.iso

A = dict(P=lin["P"], T0=lin["T0"], T14=mc["T14_h"]["median"] / 24.0, cov_T0_P=lin["cov"],
         source="item 03 transit times (T3 rule, B modelled jointly, D13), WLS with injection-calibrated errors "
                "(03_signalA/out/linear_ephemeris_calibrated.json 'adopted'); T14 = MCMC median (mcmc_summary.json)",
         sigma_pred_rule="per-epoch sigma from the full (T0, P) covariance, no extra factor (errors are injection-calibrated)")
C = np.array(lin["cov"])
A["epochs"] = []
for n in range(109, 118):
    tc = A["T0"] + n * A["P"]
    s = 1440 * np.sqrt(C[0, 0] + 2 * n * C[0, 1] + n * n * C[1, 1])
    A["epochs"].append(dict(epoch=n, tc_btjd=tc, tc_utc_calendar=utc(tc), sigma_pred_min=float(s)))
B = dict(P=float(eb["P"]), T0=float(eb["T0"]), T14=float(eb["T14"]),
         source="joint-fit ephemeris of v1/v2, unchanged (14_signalB/out/signalB.json 'ephemeris_B_used')",
         sigma_pred_rule="1.5 x the joint-fit MCMC per-epoch sigma, as in v2 (B timing errors are not injection-calibrated)")
B["epochs"] = []
for r in rows:
    if r["signal"] == "B" and r["error_model"] == "joint_fit_MCMC":
        n = int(r["epoch"]); tc = B["T0"] + n * B["P"]
        B["epochs"].append(dict(epoch=n, tc_btjd=tc, tc_utc_calendar=utc(tc), sigma_pred_min=1.5 * float(r["sigma_min_corr"])))
out = dict(note="Frozen ephemeris of the final S110 pre-registration. BTJD = BJD_TDB - 2457000. UTC strings are calendar "
                "labels of the barycentric times, not arrival times at the spacecraft.", A=A, B=B)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ephemeris_final.json"), "w"), indent=1)
for k in ("A", "B"):
    for e in out[k]["epochs"]:
        print(k, e["epoch"], round(e["tc_btjd"], 5), e["tc_utc_calendar"], round(e["sigma_pred_min"], 2))
