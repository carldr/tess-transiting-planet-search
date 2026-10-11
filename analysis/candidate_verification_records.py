"""Append the result of a verification step to the verification record of a candidate.

Each candidate has one Markdown file in results/candidate_verification_records/, named by its TIC ID and the period
of the sector search. Each step script appends one table row per candidate each time it runs: the date, the sector,
the step, the script, the version of the code, the verdict and the result. A re-run appends a new row, so the file
keeps the full history of the candidate.

The code version is the git commit of the repository. It ends in "-modified" when analysis/ has uncommitted changes.
"""

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


def append_verification_record_entry(sector, tic_id, period_days, step, verdict, result):
    """Append one row to the verification record of the candidate. Create the record when it does not exist."""
    path = verification_record_path(tic_id, period_days)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(f"# TIC {int(tic_id)}, period {float(period_days):.5f} d\n\n{TABLE_HEADER}")
    date = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")
    cells = [date, str(sector), step, Path(sys.argv[0]).name, code_version(), verdict, result]
    with open(path, "a") as file:
        file.write("| " + " | ".join(cell.replace("|", "/").replace("\n", " ") for cell in cells) + " |\n")
