"""Append the result of a verification step to the verification record of a candidate.

Each candidate has one Markdown file in results/candidate_verification_records/, named by its TIC ID and the period
of the sector search. Each step script appends one table row per candidate each time it runs: the date, the sector,
the step, the script, the version of the code, the verdict and the result. A re-run appends a new row, so the file
keeps the full history of the candidate.

The code version is the git commit of the repository when the step starts. It ends in "-modified" when analysis/ has
uncommitted changes at that time.

A step stores the date and the code version in its result row with verification_stamp, writes the row to its CSV,
and then appends the record row. A stop between the two writes leaves a CSV row with no record row. The step then
appends the record rows of its earlier rows when it continues, so each record row has a CSV row behind it.
"""

import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from build_qlp_sector_star_list import REPOSITORY_ROOT

VERIFICATION_RECORDS_DIRECTORY = REPOSITORY_ROOT / "results/candidate_verification_records"
TABLE_HEADER = (
    "| date (UTC) | sector | step | script | code version | verdict | result |\n"
    "|---|---|---|---|---|---|---|\n"
)


def verification_record_path(tic_id, period_days):
    return VERIFICATION_RECORDS_DIRECTORY / f"tic{int(tic_id)}_period_{float(period_days):.5f}_days.md"


def code_version():
    def git(*arguments):
        return subprocess.run(["git", *arguments], cwd=REPOSITORY_ROOT, capture_output=True, text=True).stdout.strip()

    commit = git("rev-parse", "--short", "HEAD") or "no commit"
    return f"{commit}-modified" if git("status", "--porcelain", "analysis") else commit


# A commit made while a step runs must not change the code version of its results. So the module reads the code
# version once, when the step script imports it at its start.
CODE_VERSION_AT_START = code_version()


def verification_stamp():
    """Return the date and the code version of a result, as columns for the result row of a step."""
    return {
        "verification_date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
        "code_version": CODE_VERSION_AT_START,
    }


def append_verification_record_entry(sector, tic_id, period_days, step, verdict, result, stamp):
    """Append one row to the verification record of the candidate. Create the record when it does not exist.

    The function does not append a row that the record already has, so a step can call it again for its earlier rows.
    """
    path = verification_record_path(tic_id, period_days)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(f"# TIC {int(tic_id)}, period {float(period_days):.5f} d\n\n{TABLE_HEADER}")
    cells = [stamp["verification_date_utc"], str(sector), step, Path(sys.argv[0]).name, stamp["code_version"],
             verdict, result]
    line = "| " + " | ".join(str(cell).replace("|", "/").replace("\n", " ") for cell in cells) + " |\n"
    if line in path.read_text():
        return
    with open(path, "a") as file:
        file.write(line)
        file.flush()
        os.fsync(file.fileno())
