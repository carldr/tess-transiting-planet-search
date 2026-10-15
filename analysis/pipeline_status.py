"""Write the live progress of a running step to a status file, for serve_candidate_board.py.

A step opens a PipelineStatus for its run. It calls start_candidate before each candidate and finish_candidate after
it. Code deep in a step calls report_candidate_progress with a short text, such as "MCMC step 23000", without a
reference to the PipelineStatus. The status file is results/pipeline_status/qlp_sector_<nnnn>_<step name>.json. The
script writes it atomically, so the board never reads a partial file. The file keeps the process ID, so the board
can tell a stopped step from a running one.
"""

import json
import os
from datetime import datetime, timezone

from build_qlp_sector_star_list import REPOSITORY_ROOT

PIPELINE_STATUS_DIRECTORY = REPOSITORY_ROOT / "results/pipeline_status"

# The status of the step that runs in this process, for report_candidate_progress.
current_pipeline_status = None


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


class PipelineStatus:
    def __init__(self, sector, step_number, step_name, script, candidate_count):
        self.path = PIPELINE_STATUS_DIRECTORY / f"qlp_sector_{sector:04d}_step_{step_number:02d}_{step_name}.json"
        self.status = {
            "sector": sector, "step_number": step_number, "step_name": step_name, "script": script,
            "process_id": os.getpid(), "started_utc": utc_now(), "updated_utc": utc_now(), "finished_utc": None,
            "candidate_count": candidate_count, "candidates_done": 0, "current_candidate": None,
        }

    def __enter__(self):
        global current_pipeline_status
        current_pipeline_status = self
        self.write()
        return self

    def __exit__(self, *exception):
        global current_pipeline_status
        current_pipeline_status = None
        self.status["current_candidate"] = None
        self.status["finished_utc"] = utc_now()
        self.write()

    def start_candidate(self, tic_id, period_days):
        self.status["current_candidate"] = {"tic_id": str(tic_id), "period_days": str(period_days),
                                            "started_utc": utc_now(), "progress": ""}
        self.write()

    def report_progress(self, progress):
        if self.status["current_candidate"]:
            self.status["current_candidate"]["progress"] = progress
            self.write()

    def finish_candidate(self):
        self.status["candidates_done"] += 1
        self.status["current_candidate"] = None
        self.write()

    def write(self):
        self.status["updated_utc"] = utc_now()
        PIPELINE_STATUS_DIRECTORY.mkdir(parents=True, exist_ok=True)
        temporary_path = self.path.with_name(self.path.name + ".partial")
        temporary_path.write_text(json.dumps(self.status, indent=1))
        os.replace(temporary_path, self.path)


def report_candidate_progress(progress):
    """Record the progress of the current candidate, when a step with a PipelineStatus runs in this process."""
    if current_pipeline_status is not None:
        current_pipeline_status.report_progress(progress)
