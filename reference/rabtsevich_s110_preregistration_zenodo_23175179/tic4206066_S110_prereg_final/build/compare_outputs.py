"""Compares a rebuilt copy of this deposit (build/reproduce.sh) with the deposited outputs (added 2026-10-03).

Compared: expected depths; the `results` block of every photometric dry run; the `stats`, `real`, `events`, `result`
and every injection trial list of every localisation dry run; the coverage counts; the fixed sentences of the summary
dry runs; the text of every refusal log and of build/calibrate.log. Not compared: run timestamps, the date printed
in one refusal message, and the hashes of output files, which change with every run.
Usage: python build/compare_outputs.py DEPOSIT REBUILT_COPY     (exit status 0 only if everything matches)
"""
import glob
import json
import os
import re
import sys

dep, new = sys.argv[1], sys.argv[2]
bad, n = [], 0


def check(name, a, b):
    global n
    n += 1
    if a != b:
        bad.append(name)
    print("%-48s %s" % (name, "same" if a == b else "DIFFERENT"))


def load(root, rel):
    return json.load(open(os.path.join(root, rel)))


e0, e1 = load(dep, "calib/expected_depths.json"), load(new, "calib/expected_depths.json")
check("calib/expected_depths.json", (e0["expected"], e0["per_sector"]), (e1["expected"], e1["per_sector"]))
c0, c1 = load(dep, "calib/localisation_coverage.json"), load(new, "calib/localisation_coverage.json")
keys = [k for k in c0 if k.startswith(("empirical_p", "gaussian_tail"))]
check("calib/localisation_coverage.json", [c0[k] for k in keys], [c1[k] for k in keys])
for f in sorted(glob.glob(os.path.join(dep, "dryrun", "*.json"))):
    rel = os.path.relpath(f, dep)
    if rel.endswith("s32_sap_other_version.json"):
        continue   # fixture with a script hash set to zeros on purpose (refusal check)
    a, b = load(dep, rel), load(new, rel)
    if "sentences" in a:
        check(rel, a["sentences"], b["sentences"])
    elif "results" in a:
        check(rel, a["results"], b["results"])
    else:
        part = lambda j: [j.get(k) for k in ("stats", "real", "events", "result")] + [
            v["trials"] for _, v in sorted(j.get("injections", {}).get("sets", {}).items())]
        check(rel, part(a), part(b))
for f in sorted(glob.glob(os.path.join(dep, "dryrun", "*refused.log")) + [os.path.join(dep, "build", "calibrate.log")]):
    rel = os.path.relpath(f, dep)
    date = lambda x: re.sub(r"today \d{4}-\d{2}-\d{2}", "today <date>", x)   # the refusal prints the run date
    check(rel, date(open(f).read()), date(open(os.path.join(new, rel)).read()))
print("\n%d of %d comparisons identical" % (n - len(bad), n) + ("" if not bad else "; different: " + ", ".join(bad)))
sys.exit(1 if bad else 0)
