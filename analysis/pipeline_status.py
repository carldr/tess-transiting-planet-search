"""Write the live progress of a running step to a status file, for serve_candidate_board.py.

A step opens a PipelineStatus for its run, with the number of items it has to do: stars, signals or candidates.
A step of stars or signals calls set_items_done as it goes. A step of candidates calls start_candidate before each
candidate and finish_candidate after it. Code deep in a step calls report_candidate_progress with a short text, such
as "MCMC step 23000", without a reference to the PipelineStatus. The status file is
results/pipeline_status/qlp_sector_<nnnn>_step_<nn>_<step name>.json, or all_sectors_step_<nn>_<step name>.json for a
step that serves all sectors. The script writes it atomically, so the board
never reads a partial file. The file keeps the process ID, so the board can tell a stopped step from a running one.
"""

import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

# build_qlp_sector_star_list.py imports this module, so this module does not import REPOSITORY_ROOT from it.
PIPELINE_STATUS_DIRECTORY = Path(__file__).resolve().parent.parent / "results/pipeline_status"
# set_items_done writes the file at most this often, as a step of stars can do thousands of items a second.
MINIMUM_SECONDS_BETWEEN_WRITES = 2.0

# The status of the step that runs in this process, for report_candidate_progress.
current_pipeline_status = None


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


class PipelineStatus:
    def __init__(self, sector, step_number, step_name, script, item_count, item_name="candidates"):
        """sector is None for a step that serves all sectors, such as the download of the known catalogues."""
        prefix = "all_sectors" if sector is None else f"qlp_sector_{sector:04d}"
        self.path = PIPELINE_STATUS_DIRECTORY / f"{prefix}_step_{step_number:02d}_{step_name}.json"
        self.status = {
            "sector": sector, "step_number": step_number, "step_name": step_name, "script": script,
            "process_id": os.getpid(), "started_utc": utc_now(), "updated_utc": utc_now(), "finished_utc": None,
            "item_name": item_name, "item_count": item_count, "items_done": 0, "current_candidate": None,
        }
        self.last_write = 0.0

    def __enter__(self):
        global current_pipeline_status
        current_pipeline_status = self
        self.write()
        return self

    def __exit__(self, exception_type, *exception):
        global current_pipeline_status
        current_pipeline_status = None
        self.status["current_candidate"] = None
        # A step that stops with an exception has not finished, so the board shows it as stopped.
        if exception_type is None:
            self.status["finished_utc"] = utc_now()
        self.write()

    def set_items_done(self, items_done):
        self.status["items_done"] = items_done
        if time.monotonic() - self.last_write >= MINIMUM_SECONDS_BETWEEN_WRITES:
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
        self.status["items_done"] += 1
        self.status["current_candidate"] = None
        self.write()

    def write(self):
        self.last_write = time.monotonic()
        self.status["updated_utc"] = utc_now()
        PIPELINE_STATUS_DIRECTORY.mkdir(parents=True, exist_ok=True)
        temporary_path = self.path.with_name(self.path.name + ".partial")
        temporary_path.write_text(json.dumps(self.status, indent=1))
        os.replace(temporary_path, self.path)


def report_candidate_progress(progress):
    """Record the progress of the current candidate, when a step with a PipelineStatus runs in this process."""
    if current_pipeline_status is not None:
        current_pipeline_status.report_progress(progress)
