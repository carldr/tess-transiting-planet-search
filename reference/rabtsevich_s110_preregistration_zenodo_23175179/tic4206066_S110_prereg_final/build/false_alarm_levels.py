"""Builds calib/false_alarm_levels.json from every false-alarm trial file in calib/ (outputs of
build/false_alarm_check.py, all modes that were run): the realised false-detection level ("d > det" under the noise
model, final percentiles) per product, for A-like and B-like trials, overall and for each exact number of counted
transits n. Product "spoc_2min" = Sector 32 trials (SPOC 2-min, the S110 primary product); "qlp" = Sector 98 trials
(fallback 1); "tessspoc_30min" = Sector 6 trials. Every level is given as a point estimate with a one-sided 95 %
Clopper-Pearson upper bound; the overall level also with a 95 % interval from a bootstrap over trial ephemerides.
No grouping of n is used (an earlier grouping was chosen after seeing per-n counts and was dropped after audit 4)."""
import glob
import json
import os

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CAL = os.path.join(HERE, "..", "calib")
PROD = {"S32": "spoc_2min", "S98": "qlp", "S6": "tessspoc_30min"}
files = sorted(f for f in glob.glob(os.path.join(CAL, "false_alarm_*.json")) if not f.endswith("levels.json"))
rows = []   # (product, signal, ephemeris id, n, exceed)
for f in files:
    for j, r in enumerate(json.load(open(f))):
        for t in r["trials"]:
            rows.append((PROD[r["sec"]], r["k"], "%s:%d" % (os.path.basename(f), j), int(t["n"]), bool(t["final"]["fa"])))
rng = np.random.default_rng(5)


def level(sel, cluster_ci=False):
    k, n = int(sum(r[4] for r in sel)), len(sel)
    out = dict(trials=n, exceed=k, rate_percent=100.0 * k / n if n else None,
               upper95_percent=float(100 * stats.beta.ppf(0.95, k + 1, n - k)) if n else None)
    if cluster_ci:
        ids = sorted({r[2] for r in sel})
        cl = [np.array([r[4] for r in sel if r[2] == i], float) for i in ids]
        bs = [np.concatenate([cl[m] for m in rng.integers(0, len(cl), len(cl))]).mean() for _ in range(4000)]
        out["ci95_percent"] = [float(100 * np.percentile(bs, 2.5)), float(100 * np.percentile(bs, 97.5))]
        out["ephemerides"] = len(ids)
    return out


res = dict(note=__doc__, files=[os.path.basename(f) for f in files])
for prod in sorted({r[0] for r in rows}):
    res[prod] = {}
    for sig in ("A", "B"):
        sel = [r for r in rows if r[0] == prod and r[1] == sig]
        if sel:
            res[prod][sig] = dict(overall=level(sel, True),
                                  by_n={str(n): level([r for r in sel if r[3] == n]) for n in sorted({r[3] for r in sel})})
json.dump(res, open(os.path.join(CAL, "false_alarm_levels.json"), "w"), indent=1)
for prod in res:
    if prod in ("note", "files"):
        continue
    for sig, v in res[prod].items():
        o = v["overall"]
        print(prod, sig, "overall %.2f %% (%d/%d), CI %s, upper %.2f %%" % (o["rate_percent"], o["exceed"], o["trials"],
              [round(x, 2) for x in o["ci95_percent"]], o["upper95_percent"]))
        print("   ", {n: "%.2f%% (%d/%d) <%.2f%%" % (x["rate_percent"], x["exceed"], x["trials"], x["upper95_percent"])
                      for n, x in v["by_n"].items()})
