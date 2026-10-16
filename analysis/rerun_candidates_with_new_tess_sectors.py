"""Check the candidates for new ExoFOP entries and new TESS sectors, and run steps 7 to 11 again on a new sector.

A launchd job runs this script every night (README, "Check the candidates for follow-up"). The script:
1. stops when a step of the pipeline is running, so that two runs never write the same results;
2. runs check_candidates_on_exofop_and_future_tess_sectors.py;
3. for each candidate that TRICERATOPS did not reject and that has no "dropped" decision, reads the released TESS
   sectors that step 9 did not use;
4. runs steps 7 to 11 again with --restart --tic-id for each candidate with such sectors. It skips a candidate
   whose released sectors are the same as at its last rerun. Step 9 can drop a sector without a transit, so the
   same sector can stay in the list after a rerun;
5. commits the results of the check and of each rerun candidate, and pushes them.

Steps 5 and 6 write no status file, so the script cannot see them run. A rerun that reads the catalogues while step 5
downloads them can fail. The next night then tries again.

A step that fails stops the reruns. The script does not record the rerun, so the next night tries again.

Run from the repository root:
    .venv/bin/python analysis/rerun_candidates_with_new_tess_sectors.py

Output: results/candidate_exofop_and_future_tess_sector_checks/new_tess_sector_reruns.csv, one row per rerun.
"""

import csv
import subprocess
import sys
from datetime import datetime, timezone

from build_qlp_sector_star_list import REPOSITORY_ROOT
from candidate_verification_records import verification_record_path
from check_candidates_on_exofop_and_future_tess_sectors import OUTPUT_DIRECTORY, candidates_by_tic_id
from serve_candidate_board import candidate_decisions, step_statuses

PYTHON = str(REPOSITORY_ROOT / ".venv/bin/python")
TRICERATOPS_PYTHON = str(REPOSITORY_ROOT / ".venv_triceratops/bin/python")
CHECK_PATH = OUTPUT_DIRECTORY / "candidate_exofop_and_future_tess_sector_checks.csv"
RERUNS_PATH = OUTPUT_DIRECTORY / "new_tess_sector_reruns.csv"
RERUN_COLUMNS = ["tic_id", "period_days", "sector", "released_sectors", "rerun_date_utc"]
# The script and the Python environment of steps 7 to 11, and the results folder of each.
STEPS = [
    (PYTHON, "check_qlp_sector_new_candidates.py", "new_candidate_checks"),
    (PYTHON, "check_qlp_sector_candidate_pixels.py", "candidate_pixel_checks"),
    (PYTHON, "verify_qlp_sector_candidates_over_all_sectors.py", "all_sector_verification"),
    (TRICERATOPS_PYTHON, "validate_qlp_sector_candidates_with_triceratops.py", "triceratops_validation"),
    (PYTHON, "decide_qlp_sector_candidate_submissions.py", "submission_decisions"),
]


def log(message):
    print(f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S} UTC {message}", flush=True)


def run(command):
    log("run " + " ".join(command))
    return subprocess.run(command, cwd=REPOSITORY_ROOT).returncode


def last_rerun_sectors():
    if not RERUNS_PATH.exists():
        return {}
    return {(row["tic_id"], row["period_days"]): row["released_sectors"] for row in csv.DictReader(open(RERUNS_PATH))}


def record_rerun(candidate, released_sectors):
    new_file = not RERUNS_PATH.exists()
    with open(RERUNS_PATH, "a", newline="") as file:
        writer = csv.DictWriter(file, RERUN_COLUMNS)
        if new_file:
            writer.writeheader()
        writer.writerow({"tic_id": candidate["tic_id"], "period_days": candidate["period_days"],
                         "sector": candidate["sector"], "released_sectors": released_sectors,
                         "rerun_date_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M")})


def candidates_to_rerun():
    """Return the candidates with released sectors that step 9 did not use, and those sectors."""
    released_by_tic_id = {row["tic_id"]: row["released_tess_sectors_not_in_all_sector_verification"]
                          for row in csv.DictReader(open(CHECK_PATH))}
    decisions = candidate_decisions()
    last_reruns = last_rerun_sectors()
    reruns = []
    for tic_id, candidates in candidates_by_tic_id().items():
        released = released_by_tic_id.get(tic_id, "")
        for candidate in candidates:
            key = (candidate["tic_id"], candidate["period_days"])
            decision = decisions.get(key, {}).get("decision")
            if released and decision != "dropped" and last_reruns.get(key) != released:
                reruns.append((candidate, released))
    return reruns


def commit_and_push(message, candidate=None):
    """Commit the results of the check, and of the rerun of the candidate, and push them.

    The commit holds only these paths, so it never takes the results of a sector that another run is searching,
    or files that someone has staged.
    """
    paths = [OUTPUT_DIRECTORY]
    if candidate:
        paths += [REPOSITORY_ROOT / f"results/qlp_sector_{int(candidate['sector']):04d}_{folder}" for _, _, folder in STEPS]
        paths += [verification_record_path(candidate["tic_id"], candidate["period_days"])]
    paths = [str(path.relative_to(REPOSITORY_ROOT)) for path in paths if path.exists()]
    run(["git", "add", "--", *paths])
    if subprocess.run(["git", "diff", "--cached", "--quiet", "--", *paths], cwd=REPOSITORY_ROOT).returncode == 0:
        log("nothing to commit")
        return
    run(["git", "commit", "-q", "-m", message, "--", *paths])
    run(["git", "push", "-q", "origin", "HEAD"])


def main():
    running = [status for status in step_statuses() if status["state"] == "running"]
    if running:
        log("a step is running, so the check waits for the next night: "
            + ", ".join(f"sector {status['sector']} step {status['step_number']}" for status in running))
        return 0
    if run([PYTHON, "analysis/check_candidates_on_exofop_and_future_tess_sectors.py"]) != 0:
        log("the check failed")
        return 1
    reruns = candidates_to_rerun()
    log(f"{len(reruns)} candidates have new released sectors")
    for candidate, released in reruns:
        log(f"TIC {candidate['tic_id']} of sector {candidate['sector']}: new sectors {released}")
        for python, script, _ in STEPS:
            command = [python, f"analysis/{script}", str(candidate["sector"]), "--restart",
                       "--tic-id", candidate["tic_id"]]
            if run(command) != 0:
                log(f"{script} failed, so the rerun stops")
                commit_and_push(f"Rerun steps 7 to 11 of TIC {candidate['tic_id']} with new sectors {released}, "
                                f"stopped at {script}", candidate)
                return 1
        record_rerun(candidate, released)
        commit_and_push(f"Rerun steps 7 to 11 of TIC {candidate['tic_id']} with new sectors {released}", candidate)
    commit_and_push("Check the candidates on ExoFOP and for new TESS sectors")
    return 0


if __name__ == "__main__":
    sys.exit(main())
