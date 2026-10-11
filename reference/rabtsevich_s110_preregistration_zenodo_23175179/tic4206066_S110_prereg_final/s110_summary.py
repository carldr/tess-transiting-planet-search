"""TIC 4206066, Sector 110: the fixed sentences of PREREG_S110_FINAL.md section 9b, produced mechanically from the
output JSON files, so that the description of every outcome is fixed before the data (written 2026-10-02, revised
after audits 4, 5 and 6).

Checks before any sentence is written (the script refuses otherwise):
  * the primary JSON is a registered Sector 110 run (or, with --dry-run, an archival dry run);
  * every variant JSON is exactly the registered variant of that primary: same input file (sha256) except (d), same
    arguments except the variant's own switch; keys a to f only; (c) and (d) only for a SPOC primary;
  * every applicable variant is given, or declared not available with --not-available KEY=FILE, where FILE is the saved,
    non-empty error output of the attempted registered command; its sha256 and text are copied into the output
    (audit 5); every result (primary, variants, localisation) must carry the script, ephemeris, expected-depth and
    helper-module hashes of this package (audit 6);
  * variant (d), SPOC primary only (audit 6): --d-status pending while no QLP Sector 110 release is public (the saved
    MAST search may be given with --not-available d=FILE; the qualifier stays provisional; refused on or after
    2027-02-28); --d-status final when (d) was run (--variant d=FILE), or with --not-available d=FILE and
    --d-final-reason crashed (a public release was attempted and failed) or not-public-on-deadline (accepted only on or
    after 2027-02-28, by the UTC date of the summary run); for a QLP primary --d-status is final automatically;
  * a localisation JSON is required for each signal whose outcome is O1 or O2 when the primary product is the SPOC 2-min
    light curve; it must be for that signal, of the same kind (registered or dry run) and computed from this primary
    JSON (its sha256 in the localisation provenance), and its pixel file must be of the primary's sector and of
    TIC 4206066 (audit 5). For a fallback primary (QLP or TESS-SPOC) localisation is not
    possible (section 6) and no localisation file is accepted.
Robustness (section 3): qualifying variants are (a) to (e) with P(detect | real) >= 0.95; none qualifying -> "robustness
not testable"; O1 robust if no qualifying variant returns O3; O3 robust if no qualifying variant returns O1 or O2;
every qualifier except "not robust" is provisional while variant (d) is pending. Only these sentences may describe the
outcomes; combined values, variants, reprocessings, later products and individual transits appear only in tables.

Usage:
  python s110_summary.py --primary s110_primary.json [--variant a=FILE ...] [--not-available c=ERRORFILE ...]
                         [--loc-A FILE] [--loc-B FILE] [--d-status pending|final] [--d-final-reason R]
                         [--dry-run [--today YYYY-MM-DD]] --out s110_summary.json
"""
import argparse
import datetime
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET_TIC = 4206066
LEVELS_FILE = os.path.join(HERE, "calib", "false_alarm_levels.json")
FA = json.load(open(LEVELS_FILE))   # frozen realised levels (section 4)
SENS_MIN = 0.95
D_DEADLINE = "2027-02-28"   # section 3: the qualifier of a SPOC primary becomes final on this date
SWITCH = {"a": dict(window=0.5), "b": dict(window=1.0), "c": dict(drop_desat_hours=1.0), "e": dict(flux="SAP_FLUX"),
          "f": dict(count_rule="v2"), "d": dict(product="qlp", flux="QLP_DEFAULT")}
ARGKEYS = ("product", "flux", "window", "drop_desat_hours", "count_rule", "t0_shift_frac", "cut_true")


def sha(f):
    return hashlib.sha256(open(f, "rb").read()).hexdigest()


def short(oc):
    return oc.split("_")[0] if oc else "O4"


def fail(msg):
    raise SystemExit("summary refused: " + msg)


VERSION_KEYS = ("script_sha256", "ephemeris_sha256", "expected_depths_sha256", "deps_sha256")
DEPS = ("lcio.py", "pixloc_lib.py")


def package_versions(script):
    """The hashes a result must carry to belong to this package (audit 6): analysis script, ephemeris, expected
    depths and both helper modules."""
    return dict(script_sha256=sha(os.path.join(HERE, script)), ephemeris_sha256=sha(os.path.join(HERE, "ephemeris_final.json")),
                expected_depths_sha256=sha(os.path.join(HERE, "calib", "expected_depths.json")),
                deps_sha256={f: sha(os.path.join(HERE, "deps", f)) for f in DEPS})


def check_versions(js, what, script):
    pv, want = js.get("provenance", {}), package_versions(script)
    bad = [q for q in VERSION_KEYS if pv.get(q) != want[q]]
    if pv.get("script") != script or bad:
        fail("%s was not produced by this package's %s (differs: %s)" % (what, script, ", ".join(bad) or "script name"))


def input_sha(js):
    v = list(js["provenance"]["inputs"].values())
    return v[0] if v else None


def check_variant(k, v, prim):
    if v.get("dry_run") != prim.get("dry_run") or v.get("sector") != prim.get("sector"):
        fail("variant %s is not of the same kind / sector as the primary" % k)
    want = {q: prim["args"][q] for q in ARGKEYS}
    want.update(SWITCH[k])
    got = {q: v["args"].get(q) for q in ARGKEYS}
    if got != want:
        fail("variant %s arguments %s differ from the registered %s" % (k, got, want))
    if k != "d" and input_sha(v) != input_sha(prim):
        fail("variant %s was not run on the primary input file" % k)
    check_versions(v, "variant %s" % k, "s110_final.py")   # (d) may differ only in its input file (audit 6)


def fa_clause(sig, n, prim):
    ptype = prim.get("product_type") or ("spoc_2min" if prim["args"]["product"] == "spoc" else "qlp")
    key = "spoc_2min" if ptype in ("spoc_2min", "tessspoc") else "qlp"
    tab = FA[key][sig]
    o, x = tab["overall"], tab["by_n"].get(str(n))
    s = ("measured false-detection level for this product about %.2f %% overall (%d of %d trials; binomial 95 %% upper "
         "bound %.2f %%, approximate because the trials share data)" % (
             o["rate_percent"], o["exceed"], o["trials"], o["upper95_percent"]))
    if x:
        s += " and %.2f %% for exactly %d counted transits (%d of %d; binomial bound %.2f %%)" % (
            x["rate_percent"], n, x["exceed"], x["trials"], x["upper95_percent"])
    else:
        s += "; no trials with exactly %d counted transits" % n
    if ptype == "tessspoc":
        s += "; measured on SPOC 2-min data, no trials exist for the TESS-SPOC product"
    return s


def robustness(sig, oc, variants, d_pending):
    kind = short(oc)
    if kind not in ("O1", "O3"):
        return None, []
    qual = {k: v["results"][sig] for k, v in variants.items() if k in ("a", "b", "c", "d", "e")
            and v["results"][sig].get("P_detect_if_real", 0) >= SENS_MIN}
    against = sorted(k for k, r in qual.items() if (kind == "O1" and short(r["outcome"]) == "O3") or
                     (kind == "O3" and short(r["outcome"]) in ("O1", "O2")))
    if against:
        return "not robust", against
    q = "robust" if qual else "robustness not testable"
    return (q + " (provisional until 2027-02-28)" if d_pending else q), []


def p_text(p, n):
    s = "empirical p = %.4g from %d on-target injections, " % (p, n)
    if p >= 0.05:
        return s + "not below 0.05"
    if p >= 0.0027:
        return s + "below 0.05 and not below 0.0027"
    return s + "below 0.0027"


FALLBACK_LOC = {"result": "fallback product: localisation not possible"}


def loc_clause_A(loc):
    if "label" not in loc and "label_if_this_were_sector_110" not in loc:
        return "its position could not be measured (%s)." % loc.get("result", "no label written")
    lab = loc.get("label") or loc["label_if_this_were_sector_110"].split()[0]   # the second key exists only in dry runs
    p, n = loc["stats"]["p_emp"], loc["stats"]["n_cal"]
    if lab == "L1_COMPATIBLE":
        return "its pixel position is compatible with TIC 4206066 (%s)." % p_text(p, n)
    if lab == "L2_TENSION":
        return "its pixel position is in tension with TIC 4206066 (%s)." % p_text(p, n)
    if lab == "L3_OFF_TARGET":
        return ("its pixel position is off target (%s): the signal most likely comes from another source, and a planet "
                "on TIC 4206066 is disfavoured." % p_text(p, n))
    return "its position could not be measured (label %s)." % lab


def loc_clause_B(loc):
    if "stats" not in loc:
        return "B's localisation was not possible (%s)." % loc.get("result", "no result written")
    return "B's pixel offset is reported in the table without a label and is not interpreted in this test."


def sentence(sig, res, rob, against, loc, prim, where):
    oc = res.get("outcome", "O4_INSUFFICIENT")
    kind = short(oc)
    if kind == "O1":
        s = "Signal %s was recovered in %s at the pre-registered times and depth (%d counted transits; %s)" % (
            sig, where, res["n_counted"], fa_clause(sig, res["n_counted"], prim))
    elif kind == "O2":
        why = oc[oc.find("(") + 1:oc.rfind(")")] if "(" in oc else "condition not recorded"
        s = ("A dip was found at the predicted times of signal %s, but it is not the predicted signal (%s), as "
             "pre-registered; this contradicts the archival model of signal %s" % (sig, why, sig))
    elif kind == "O3":
        s = "The predicted signal %s is excluded by %s (pre-registered sensitivity %.3f >= 0.95)" % (
            sig, where, res.get("P_detect_if_real", float("nan")))
    else:
        return "%s cannot decide on signal %s (reason: %s)." % (where[0].upper() + where[1:], sig,
                                                                 res.get("reason", "not recorded"))
    if rob == "not robust":
        verb = "excluded it" if kind == "O1" else "detected a dip"
        s += ("; a processing variant with adequate sensitivity %s (variant %s); this is reported and does not change "
              "the outcome" % (verb, ", ".join(against)))
    elif rob:
        s += "; " + rob
    if kind in ("O1", "O2"):
        s += "; " + (loc_clause_A(loc) if sig == "A" else loc_clause_B(loc))
    else:
        s += "."
    return s


def main(a):
    prim = json.load(open(a.primary))
    if a.dry_run:
        if not prim.get("dry_run") or prim.get("sector") == 110:
            fail("--dry-run needs an archival dry-run primary")
        where = "the archival dry run (Sector %s)" % prim.get("sector")
    else:
        if prim.get("dry_run") or prim.get("sector") != 110 or not prim.get("registered_settings"):
            fail("the primary JSON is not a registered Sector 110 run")
        where = "Sector 110"
    check_versions(prim, "the primary JSON", "s110_final.py")
    spoc = prim["args"]["product"] == "spoc"
    applicable = ["a", "b", "e", "f"] + (["c", "d"] if spoc else [])
    variants, not_avail = {}, {}
    for v in a.variant or []:
        k, f = v.split("=", 1)
        if k not in applicable or k in variants:
            fail("variant key %r is not an applicable registered variant (or given twice)" % k)
        variants[k] = json.load(open(f))
        check_variant(k, variants[k], prim)
    for v in a.not_available or []:
        k, ef = v.split("=", 1)
        if k not in applicable or k in variants or k in not_avail:
            fail("--not-available %r is not an applicable variant without a file (or given twice)" % k)
        if not os.path.isfile(ef) or not open(ef, encoding="utf-8", errors="replace").read().strip():
            fail("--not-available %s needs the saved, non-empty error output of the attempted run (got %r)" % (k, ef))
        not_avail[k] = dict(error_file=os.path.basename(ef), sha256=sha(ef),
                            text=open(ef, encoding="utf-8", errors="replace").read())
    d_status = a.d_status if spoc else "final"
    if a.today and not a.dry_run:
        fail("--today is only for code checks with --dry-run; a registered summary uses the real UTC date")
    today = a.today or datetime.datetime.now(datetime.timezone.utc).date().isoformat()
    d_reason = None
    if spoc and d_status == "pending":
        # QLP Sector 110 not yet public: (d) is not run; a saved MAST search may be given with --not-available d=FILE
        # and is then copied into the summary as the record of that search; the qualifier stays provisional (audit 6)
        if "d" in variants:
            fail("--d-status pending with a variant (d) file: a run (d) makes the status final")
        if a.d_final_reason:
            fail("--d-final-reason is only for --d-status final")
        if today >= D_DEADLINE:
            fail("on or after %s the status of variant (d) is final (section 3)" % D_DEADLINE)
    elif spoc:
        if "d" in variants:
            d_reason = "run"
            if a.d_final_reason not in (None, "run"):
                fail("variant (d) was run: --d-final-reason must be omitted or 'run'")
        elif "d" in not_avail:
            d_reason = a.d_final_reason
            if d_reason not in ("crashed", "not-public-on-deadline"):
                fail("--d-status final without a variant (d) file needs --d-final-reason crashed (a public release "
                     "was attempted and failed; FILE is its error output) or not-public-on-deadline (FILE is the saved "
                     "MAST search of %s or later finding no QLP Sector 110 release)" % D_DEADLINE)
            if d_reason == "not-public-on-deadline" and today < D_DEADLINE:
                fail("'not public' makes the status of (d) final only on or after %s (today %s)" % (D_DEADLINE, today))
    missing = [k for k in applicable if k not in variants and k not in not_avail and not (k == "d" and d_status == "pending")]
    if missing:
        fail("variants %s are neither given nor declared not available" % missing)
    d_pending = spoc and d_status == "pending"
    psha = sha(a.primary)
    locs = {}
    ptype = prim.get("product_type") or ("spoc_2min" if spoc else "qlp")
    fallback = ptype != "spoc_2min" if not a.dry_run else not spoc
    for sig, f in (("A", a.loc_A), ("B", a.loc_B)):
        need = short(prim["results"][sig].get("outcome")) in ("O1", "O2")
        if fallback:
            if f:
                fail("no localisation file is accepted for a fallback primary product")
            locs[sig] = FALLBACK_LOC
            continue
        if need and not f:
            fail("signal %s has outcome %s: its localisation JSON is required" % (sig, prim["results"][sig]["outcome"]))
        if f:
            L = json.load(open(f))
            if L["args"].get("signal") != sig or bool(L.get("dry_run")) != bool(a.dry_run):
                fail("localisation file %s is not the %s localisation of this kind of run" % (f, sig))
            check_versions(L, "localisation file %s" % f, "s110_localise.py")
            if L["provenance"]["inputs"].get(os.path.basename(a.primary)) != psha:
                fail("localisation file %s was not computed from this primary JSON" % f)
            if L.get("file_sector") != prim.get("sector") or L.get("file_tic") not in (TARGET_TIC, str(TARGET_TIC)):
                fail("localisation file %s is not from a Sector %s pixel file of TIC %d" % (f, prim.get("sector"), TARGET_TIC))
            locs[sig] = L
    out = dict(script="s110_summary.py", script_sha256=sha(os.path.abspath(__file__)), levels_sha256=sha(LEVELS_FILE),
               inputs_sha256={os.path.basename(f): sha(f) for f in [a.primary] + [v.split("=", 1)[1] for v in a.variant or []]
                              + [v.split("=", 1)[1] for v in a.not_available or []] + [f for f in (a.loc_A, a.loc_B) if f]},
               dry_run=a.dry_run, variants=sorted(variants), not_available=not_avail, d_status=d_status,
               d_final_reason=d_reason, summary_date_utc=today, sentences={})
    for sig in ("A", "B"):
        res = prim["results"][sig]
        rob, against = robustness(sig, res.get("outcome"), variants, d_pending)
        out["sentences"][sig] = sentence(sig, res, rob, against, locs.get(sig), prim, where)
        print(out["sentences"][sig])
    out["fixed_note"] = ("No combined statement about A and B or about the system is part of this test; any such "
                         "statement is exploratory and cannot be cited for or against either outcome.")
    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--primary", required=True)
    p.add_argument("--variant", action="append")
    p.add_argument("--not-available", dest="not_available", action="append")
    p.add_argument("--loc-A", dest="loc_A")
    p.add_argument("--loc-B", dest="loc_B")
    p.add_argument("--d-status", choices=["pending", "final"], default="pending")
    p.add_argument("--d-final-reason", dest="d_final_reason", choices=["run", "crashed", "not-public-on-deadline"])
    p.add_argument("--today", help="YYYY-MM-DD, only with --dry-run (code checks of the date rule)")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--out", required=True)
    main(p.parse_args())
