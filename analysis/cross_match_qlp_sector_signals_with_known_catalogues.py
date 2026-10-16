"""Compare the vetted signals of one sector with the known TOIs, CTOIs, SPOC TCEs and TESS EBs, and count each stage.

A signal matches a known entry when the entry is on the same TIC ID and its period is 1/3, 1/2, 1, 2 or 3 times
the signal period, within 1%.

The script also measures the search against the known TOIs: it takes every TOI on a searched star with a period
inside the search range, and counts how many the search found and the vetting passed, by TFOPWG disposition.

Run from the repository root, after vet_qlp_sector_signals.py and download_known_transit_catalogues.py:
    .venv/bin/python analysis/cross_match_qlp_sector_signals_with_known_catalogues.py 98

Outputs, in results/qlp_sector_<nnnn>_known_catalogue_cross_match/:
- signals_cross_matched.csv: every vetted signal, with the matching TOIs, CTOIs, TCEs and TESS EBs;
- new_candidates.csv: the signals that pass every vetting test and match no known entry, highest MES first;
- summary.json: the number of signals at each stage, the TOI recovery, and the rank of each test star.
The script also starts the verification record of each new candidate, in results/candidate_verification_records/.
"""

import argparse
import csv
import json
from collections import defaultdict

import numpy as np

from build_qlp_sector_star_list import REPOSITORY_ROOT, sector_directory
from candidate_verification_records import append_verification_record_entry, verification_stamp
from pipeline_status import PipelineStatus
from search_qlp_sector import MAXIMUM_PERIOD_DAYS, MINIMUM_PERIOD_DAYS, output_directory as search_output_directory
from vet_qlp_sector_signals import output_directory as vetting_output_directory

CATALOGUE_DIRECTORY = REPOSITORY_ROOT / "data/known_transit_catalogues"
PERIOD_RATIOS = (1 / 3, 1 / 2, 1, 2, 3)
PERIOD_MATCH_TOLERANCE = 0.01
# Known positives to locate in the ranked list: TIC ID -> name.
TEST_STARS = {4206066: "TIC 4206066 (Rabtsevich)"}


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_known_catalogue_cross_match"


def to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return np.nan


def known_entries_by_tic_id():
    """Return {TIC ID: [(catalogue, identifier, period, toi_disposition)]}."""
    entries = defaultdict(list)
    for row in csv.DictReader(open(CATALOGUE_DIRECTORY / "toi.csv")):
        disposition = row["TFOPWG Disposition"] or "none"
        entries[int(row["TIC ID"])].append(("TOI", row["TOI"], to_float(row["Period (days)"]), disposition))
    for row in csv.DictReader(open(CATALOGUE_DIRECTORY / "ctoi.csv")):
        entries[int(row["TIC ID"])].append(("CTOI", row["CTOI"], to_float(row["Period (days)"]), ""))
    for path in sorted((CATALOGUE_DIRECTORY / "spoc_tce").glob("*.csv")):
        lines = (line for line in open(path) if not line.startswith("#"))
        for row in csv.DictReader(lines):
            entries[int(row["ticid"])].append(("TCE", row["tceid"], to_float(row["tce_period"]), ""))
    for row in csv.DictReader(open(CATALOGUE_DIRECTORY / "tess_eb.csv")):
        entries[int(row["tic_id"])].append(("TESS_EB", row["tic_id"], to_float(row["period_days"]), ""))
    return entries


def periods_match(period, known_period):
    if not np.isfinite(known_period) or known_period <= 0:
        return False
    return any(abs(period / known_period / ratio - 1) < PERIOD_MATCH_TOLERANCE for ratio in PERIOD_RATIOS)


def matching_identifiers(entries, period, catalogue):
    return [identifier for name, identifier, known_period, _ in entries
            if name == catalogue and periods_match(period, known_period)]


def search_signals_by_tic_id(sector):
    signals = defaultdict(list)
    for chunk_path in sorted(search_output_directory(sector).glob("chunk_*[0-9].csv")):
        for row in csv.DictReader(open(chunk_path)):
            signals[int(row["tic_id"])].append(row)
    return signals


def passes_vetting(row):
    return row["passes_false_alarm_tests"] == "True" and row["passes_false_positive_tests"] == "True"


def recovery_counts(toi_tests):
    return {
        "tois_testable": len(toi_tests),
        "found_by_search": sum(test["found_by_search"] for test in toi_tests),
        "passed_vetting": sum(test["passed_vetting"] for test in toi_tests),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    sector = parser.parse_args().sector
    with PipelineStatus(sector, 6, "cross_match", "cross_match_qlp_sector_signals_with_known_catalogues.py", 0,
                        item_name="signals"):
        cross_match(sector)


def cross_match(sector):
    known = known_entries_by_tic_id()
    searched = search_signals_by_tic_id(sector)
    vetted = list(csv.DictReader(open(vetting_output_directory(sector) / "vetting.csv")))

    for row in vetted:
        entries = known.get(int(row["tic_id"]), [])
        period = float(row["period_days"])
        for catalogue in ("TOI", "CTOI", "TCE", "TESS_EB"):
            row[f"matching_{catalogue.lower()}"] = " ".join(matching_identifiers(entries, period, catalogue))
        row["star_has_known_entry"] = bool(entries)
        row["is_known"] = any(row[f"matching_{name}"] for name in ("toi", "ctoi", "tce", "tess_eb"))

    passed = [row for row in vetted if passes_vetting(row)]
    passed.sort(key=lambda row: -to_float(row["MES"]))
    new_candidates = [row for row in passed if not row["is_known"]]

    # TOI recovery: TOIs on searched stars with a period inside the search range.
    # A searched star has a QLP light curve in this sector, so TESS observed it. The ExoFOP Sectors column lags
    # behind the newest sectors, so the script does not use it.
    toi_tests = []
    for tic_id, entries in known.items():
        if tic_id not in searched:
            continue
        for name, identifier, period, disposition in entries:
            if name == "TOI" and MINIMUM_PERIOD_DAYS <= period <= MAXIMUM_PERIOD_DAYS:
                found = any(row["signal_number"] != "0" and periods_match(float(row["period_days"]), period)
                            for row in searched[tic_id])
                vetted_and_passed = any(int(row["tic_id"]) == tic_id and periods_match(float(row["period_days"]), period)
                                        for row in passed)
                toi_tests.append({"toi": identifier, "disposition": disposition,
                                  "found_by_search": found, "passed_vetting": vetted_and_passed})

    test_star_ranks = {}
    for tic_id, name in TEST_STARS.items():
        ranks = [rank for rank, row in enumerate(passed, start=1) if int(row["tic_id"]) == tic_id]
        test_star_ranks[name] = {
            "ranks_by_MES_among_passed_signals": ranks,
            "periods_passed": [passed[rank - 1]["period_days"] for rank in ranks],
        }

    summary = {
        "stars_searched": len(searched),
        "stars_failed_to_process": sum(1 for rows in searched.values() if rows[0]["signal_number"] == "0"),
        "signals_vetted (SNR and transit cut)": len(vetted),
        "signals_with_vetting_error": sum(1 for row in vetted if row["error"]),
        "signals_passing_false_alarm_tests": sum(1 for row in vetted if row["passes_false_alarm_tests"] == "True"),
        "signals_passing_all_tests": len(passed),
        "signals_passing_all_tests_and_known": len(passed) - len(new_candidates),
        "signals_passing_all_tests_and_new": len(new_candidates),
        "stars_with_new_candidates": len({row["tic_id"] for row in new_candidates}),
        "toi_recovery_by_tfopwg_disposition": {
            "all": recovery_counts(toi_tests),
            **{disposition: recovery_counts([test for test in toi_tests if test["disposition"] == disposition])
               for disposition in sorted({test["disposition"] for test in toi_tests})},
        },
        "test_stars": test_star_ranks,
        "catalogue_download_date": (CATALOGUE_DIRECTORY / "download_date.txt").read_text().strip(),
    }

    directory = output_directory(sector)
    directory.mkdir(parents=True, exist_ok=True)
    columns = list(vetted[0].keys())
    for file_name, rows in (("signals_cross_matched.csv", vetted), ("new_candidates.csv", new_candidates)):
        with open(directory / file_name, "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
    stamp = verification_stamp()
    summary.update(stamp)
    with open(directory / "summary.json", "w") as file:
        json.dump(summary, file, indent=1)
    # The script appends the record rows after it writes all its output, so each record row has the output behind it.
    for row in new_candidates:
        append_verification_record_entry(
            sector, row["tic_id"], row["period_days"], "search, vetting and known-signal cross-match", "new candidate",
            f"BLS SNR {to_float(row['snr']):.1f}, MES {to_float(row['MES']):.1f}, depth {to_float(row['depth_ppm']):.0f} ppm, "
            f"{row['transit_count']} transits; passes the LEO-Vetter false alarm and false positive tests; no TOI, CTOI, "
            f"TCE or TESS EB match in the catalogues of {summary['catalogue_download_date']}",
            stamp,
        )
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
