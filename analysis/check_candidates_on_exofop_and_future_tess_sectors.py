"""Check the candidates that TRICERATOPS did not reject for new ExoFOP entries and for TESS sectors not yet used.

The candidates are the rows of each results/qlp_sector_<nnnn>_triceratops_validation/triceratops_validation.csv
with the verdict "likely planet" or "validated without imaging". For each TIC ID of these candidates, the script:
1. downloads the ExoFOP target page as JSON. ExoFOP shows the TOIs, the CTOIs and the follow-up observations that
   observers upload: imaging, spectroscopy, time series, stellar companions and files;
2. compares it with the copy from the last run, and prints each entry that is new;
3. writes the new copy in place of the old one;
4. finds the TESS sectors that observe the star with tess-point, and the sectors that the all-sector verification
   did not use. A sector in the second list has new data when QLP or SPOC release it.

The script only reads ExoFOP. It does not submit anything.

Run from the repository root, after validate_qlp_sector_candidates_with_triceratops.py:
    .venv/bin/python analysis/check_candidates_on_exofop_and_future_tess_sectors.py

Output: results/candidate_exofop_and_future_tess_sector_checks/
- exofop_target_pages/tic<id>.json, the ExoFOP target page of each TIC ID, from the last run;
- candidate_exofop_and_future_tess_sector_checks.csv, one row per TIC ID.
"""

import csv
import io
import json
import os
import re
import urllib.request
from datetime import datetime, timezone

from tess_stars2px import tess_stars2px_function_entry

from build_qlp_sector_star_list import REPOSITORY_ROOT

TRICERATOPS_VERDICTS_TO_CHECK = {"likely planet", "validated without imaging"}
EXOFOP_TARGET_JSON_URL = "https://exofop.ipac.caltech.edu/tess/target.php?id={tic_id}&json"
DOWNLOAD_TIMEOUT_SECONDS = 60
# The ExoFOP target page sections that hold catalogue entries and follow-up observations, one list of entries each.
EXOFOP_SECTIONS_TO_COMPARE = ["tois", "ctois", "imaging", "spectroscopy", "time_series", "stellar_companions", "files"]
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "results/candidate_exofop_and_future_tess_sector_checks"
EXOFOP_TARGET_PAGE_DIRECTORY = OUTPUT_DIRECTORY / "exofop_target_pages"

OUTPUT_COLUMNS = [
    "tic_id", "candidates", "confirmed_planets", "tois", "ctois",
    *[f"exofop_{section}_count" for section in EXOFOP_SECTIONS_TO_COMPARE if section not in ("tois", "ctois")],
    "new_exofop_entries", "tess_sectors", "tess_sectors_not_in_all_sector_verification", "check_date_utc",
]


def candidates_by_tic_id():
    """Return, for each TIC ID, the TRICERATOPS rows with a verdict in TRICERATOPS_VERDICTS_TO_CHECK and their sector."""
    candidates = {}
    for path in sorted(REPOSITORY_ROOT.glob("results/qlp_sector_*_triceratops_validation/triceratops_validation.csv")):
        sector = int(re.search(r"qlp_sector_(\d+)_", path.parent.name).group(1))
        for row in csv.DictReader(open(path)):
            if row["triceratops_verdict"] in TRICERATOPS_VERDICTS_TO_CHECK:
                candidates.setdefault(row["tic_id"], []).append({**row, "sector": sector})
    return candidates


def download_exofop_target_page(tic_id):
    with urllib.request.urlopen(EXOFOP_TARGET_JSON_URL.format(tic_id=tic_id),
                                timeout=DOWNLOAD_TIMEOUT_SECONDS) as response:
        return json.load(response)


def new_exofop_entries(previous_page, page):
    """Return a description of each entry of EXOFOP_SECTIONS_TO_COMPARE that is in page and not in previous_page."""
    descriptions = []
    for section in EXOFOP_SECTIONS_TO_COMPARE:
        previous_entries = {json.dumps(entry, sort_keys=True) for entry in previous_page.get(section, [])}
        for entry in page.get(section, []):
            if json.dumps(entry, sort_keys=True) not in previous_entries:
                descriptions.append(f"{section}: {json.dumps(entry, sort_keys=True)}")
    return descriptions


def write_atomically(path, text):
    """Write text to a temporary file, then replace path with it, so that a stop never leaves a partial file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_name(path.name + ".partial")
    with open(temporary_path, "w") as file:
        file.write(text)
        file.flush()
        os.fsync(file.fileno())
    os.replace(temporary_path, path)


def check_tic_id(tic_id, candidates):
    page_path = EXOFOP_TARGET_PAGE_DIRECTORY / f"tic{tic_id}.json"
    previous_page = json.loads(page_path.read_text()) if page_path.exists() else {}
    page = download_exofop_target_page(tic_id)
    # On the first run there is no earlier copy, so every entry is reported as new.
    new_entries = new_exofop_entries(previous_page, page)
    write_atomically(page_path, json.dumps(page, indent=1, sort_keys=True) + "\n")

    ra, dec = float(page["coordinates"]["ra"]), float(page["coordinates"]["dec"])
    tess_sectors = sorted(set(int(sector) for sector in tess_stars2px_function_entry(int(tic_id), ra, dec)[3]))
    sectors_used = {int(sector) for candidate in candidates for sector in candidate["sectors_used"].split()}
    return {
        "tic_id": tic_id,
        "candidates": "; ".join(f"S{candidate['sector']} P {candidate['fit_period_days']} d "
                                f"{candidate['triceratops_verdict']}" for candidate in candidates),
        "confirmed_planets": page["basic_info"]["confirmed_planets"],
        "tois": " ".join(entry["toi"] for entry in page.get("tois", [])),
        "ctois": " ".join(entry["ctoi"] for entry in page.get("ctois", [])),
        **{f"exofop_{section}_count": len(page.get(section, []))
           for section in EXOFOP_SECTIONS_TO_COMPARE if section not in ("tois", "ctois")},
        "new_exofop_entries": " | ".join(new_entries),
        "tess_sectors": " ".join(str(sector) for sector in tess_sectors),
        "tess_sectors_not_in_all_sector_verification": " ".join(
            str(sector) for sector in tess_sectors if sector not in sectors_used),
        "check_date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
    }


def main():
    candidates = candidates_by_tic_id()
    print(f"{len(candidates)} TIC IDs with a TRICERATOPS verdict of "
          f"{' or '.join(sorted(TRICERATOPS_VERDICTS_TO_CHECK))}", flush=True)
    rows = []
    for tic_id, tic_candidates in candidates.items():
        row = check_tic_id(tic_id, tic_candidates)
        rows.append(row)
        print(f"TIC {tic_id}: TOIs [{row['tois']}], CTOIs [{row['ctois']}], "
              f"TESS sectors not in the all-sector verification [{row['tess_sectors_not_in_all_sector_verification']}]",
              flush=True)
        for description in row["new_exofop_entries"].split(" | ") if row["new_exofop_entries"] else []:
            print(f"  new on ExoFOP: {description}", flush=True)

    output_path = OUTPUT_DIRECTORY / "candidate_exofop_and_future_tess_sector_checks.csv"
    text = io.StringIO()
    writer = csv.DictWriter(text, OUTPUT_COLUMNS)
    writer.writeheader()
    writer.writerows(rows)
    write_atomically(output_path, text.getvalue())
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
