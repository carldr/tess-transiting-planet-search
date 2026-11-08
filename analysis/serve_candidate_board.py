"""Serve a live board of the candidates of all searched sectors, with the progress of the running steps.

The board has one card per new candidate of step 6. A card is in the column of the first step that has not passed it:
a step that has not checked it yet, or the step that stopped it. Step 11 proposes a candidate for submission or
for drop, or holds it, and the card goes to the column of that decision. A decision of a person in
results/candidate_decisions.csv moves the card on to "Held", "Submitted" or "Dropped". Before step 11 has decided,
a card with a decision of a person stays in the column of its step and shows the decision. The page reads the results
again every REFRESH_SECONDS, so it shows each result as soon as a step writes it. The steps write their progress to
results/pipeline_status/ (pipeline_status.py), and a card shows the progress while a step works on it.

The server only reads files. It serves the files in results/ for the figures and the verification records.

Run from the repository root:
    .venv/bin/python analysis/serve_candidate_board.py
Then open http://localhost:8000, or http://<the address of this computer>:8000 from another device on the network.
"""

import argparse
import csv
import json
import mimetypes
import os
import re
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import unquote, urlparse

from build_qlp_sector_star_list import REPOSITORY_ROOT
from candidate_verification_records import verification_record_path
from incremental_csv_output import complete_csv_rows
from pipeline_status import PIPELINE_STATUS_DIRECTORY

RESULTS_DIRECTORY = REPOSITORY_ROOT / "results"
CANDIDATE_DECISIONS_PATH = RESULTS_DIRECTORY / "candidate_decisions.csv"
CANDIDATE_DECISIONS = ("held", "submitted", "dropped")
REFRESH_SECONDS = 10
# The CSV of each step, under results/qlp_sector_<nnnn>_<directory>/, its verdict column, and its passing verdicts.
STEPS = [
    {"number": 7, "name": "Independent checks", "path": "new_candidate_checks/candidate_checks.csv",
     "verdict_column": "verdict", "passing_verdicts": {"supported", "unconfirmed"}},
    {"number": 8, "name": "Pixel check", "path": "candidate_pixel_checks/pixel_checks.csv",
     "verdict_column": "pixel_verdict", "passing_verdicts": {"on target", "unresolved"}},
    {"number": 9, "name": "All-sector fit", "path": "all_sector_verification/all_sector_verification.csv",
     "verdict_column": "verdict", "passing_verdicts": {"passed"}},
    {"number": 10, "name": "TRICERATOPS", "path": "triceratops_validation/triceratops_validation.csv",
     "verdict_column": "triceratops_verdict", "passing_verdicts": {"likely planet", "validated without imaging"}},
    {"number": 11, "name": "Submission decision", "path": "submission_decisions/submission_decisions.csv",
     "verdict_column": "decision", "passing_verdicts": {"proposed for submission", "proposed for drop", "held"}},
]
# The column of each decision of step 11.
SUBMISSION_DECISION_COLUMNS = {"proposed for submission": "proposed_submission", "proposed for drop": "proposed_drop",
                               "held": "held"}


def searched_sectors():
    sectors = []
    for path in RESULTS_DIRECTORY.glob("qlp_sector_*_known_catalogue_cross_match/new_candidates.csv"):
        match = re.fullmatch(r"qlp_sector_(\d{4})_known_catalogue_cross_match", path.parent.name)
        if match:
            sectors.append(int(match.group(1)))
    return sorted(sectors)


# The output that each of steps 1 to 6 leaves when it has finished, for a sector run before the status files existed.
SECTOR_STEPS = [
    {"number": 1, "name": "star list", "output": "data/qlp_sector_{sector:04d}/star_list.csv"},
    {"number": 2, "name": "download", "output": "data/qlp_sector_{sector:04d}/failed_downloads.csv"},
    {"number": 3, "name": "BLS search", "output": "results/qlp_sector_{sector:04d}_bls_search"},
    {"number": 4, "name": "LEO-Vetter", "output": "results/qlp_sector_{sector:04d}_leo_vetting/vetting.csv.gz"},
    {"number": 6, "name": "cross-match", "output": "results/qlp_sector_{sector:04d}_known_catalogue_cross_match/summary.json"},
]


def sector_progress(statuses):
    """Return the state of steps 1 to 6 of each sector with a status file, a star list or a cross-match result."""
    sectors = set(searched_sectors()) | {status["sector"] for status in statuses if status["sector"] is not None}
    sectors |= {int(match.group(1)) for path in (REPOSITORY_ROOT / "data").glob("qlp_sector_*/star_list.csv")
                if (match := re.fullmatch(r"qlp_sector_(\d{4})", path.parent.name))}
    progress = []
    for sector in sorted(sectors):
        steps = []
        for step in SECTOR_STEPS:
            status = next((status for status in statuses
                           if status["sector"] == sector and status["step_number"] == step["number"]), None)
            if status:
                steps.append({"number": step["number"], "name": step["name"], "state": status["state"],
                              "items_done": status.get("items_done", 0), "item_count": status.get("item_count", 0),
                              "item_name": status.get("item_name", "")})
            else:
                finished = (REPOSITORY_ROOT / step["output"].format(sector=sector)).exists()
                steps.append({"number": step["number"], "name": step["name"],
                              "state": "finished" if finished else "not started"})
        summary_path = RESULTS_DIRECTORY / f"qlp_sector_{sector:04d}_known_catalogue_cross_match/summary.json"
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        progress.append({"sector": sector, "steps": steps, "summary": {
            name: summary.get(name) for name in (
                "stars_searched", "stars_failed_to_process", "signals_vetted (SNR and transit cut)",
                "signals_passing_all_tests", "signals_passing_all_tests_and_known", "signals_passing_all_tests_and_new",
                "catalogue_download_date")}})
    return progress


def rows_by_candidate(path):
    return {(row["tic_id"], row["period_days"]): row for row in complete_csv_rows(path)}


def process_is_running(process_id):
    try:
        os.kill(process_id, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def step_statuses():
    statuses = []
    for path in sorted(PIPELINE_STATUS_DIRECTORY.glob("*.json")):
        try:
            status = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if status["finished_utc"]:
            status["state"] = "finished"
        else:
            status["state"] = "running" if process_is_running(status["process_id"]) else "stopped"
        statuses.append(status)
    return statuses


def candidate_decisions():
    return {(row["tic_id"], row["period_days"]): row for row in complete_csv_rows(CANDIDATE_DECISIONS_PATH)}


def stop_reason(step, row):
    if step["number"] in (7, 9):
        return row.get("refuted_by") or row[step["verdict_column"]]
    return row[step["verdict_column"]]


def figure_paths(sector, tic_id):
    paths = []
    for directory in sorted(RESULTS_DIRECTORY.glob(f"qlp_sector_{sector:04d}_*")):
        if directory.name.endswith("_quality_zero_mask"):
            continue
        paths += [str(path.relative_to(RESULTS_DIRECTORY)) for path in sorted(directory.glob(f"tic{tic_id}_*.png"))]
    return paths


def candidate_card(sector, new_candidate, step_rows, decisions, running):
    key = (new_candidate["tic_id"], new_candidate["period_days"])
    card = {
        "sector": sector, "tic_id": key[0], "period_days": key[1], "depth_ppm": new_candidate["depth_ppm"],
        "tmag": new_candidate["tmag"], "trail": [], "cautions": [], "column": "", "stop_reason": "",
        "error": "", "running": running.get(key), "figures": figure_paths(sector, key[0]),
    }
    record = verification_record_path(*key)
    card["record"] = str(record.relative_to(RESULTS_DIRECTORY)) if record.exists() else ""
    for step in STEPS:
        row = step_rows[step["number"]].get(key)
        if row is None:
            card["column"] = f"step_{step['number']}"
            break
        if row.get("error"):
            card["column"], card["error"] = f"step_{step['number']}", row["error"]
            break
        verdict = row[step["verdict_column"]]
        card["trail"].append({"step": step["number"], "verdict": verdict,
                              "date_utc": row.get("verification_date_utc", "")})
        card["cautions"] += [caution for caution in row.get("cautions", "").split("; ") if caution]
        if step["number"] == 7:
            card["radius_earth_radii"] = row.get("companion_radius_earth_radii", "")
        if step["number"] == 9:
            card["radius_earth_radii"] = row.get("planet_radius_earth_radii", "")
            card["fit_period_days"] = row.get("fit_period_days", "")
            card["ellipsoidal"] = (f"{row['ellipsoidal_amplitude_ppm']} +/- {row['ellipsoidal_amplitude_error_ppm']} ppm"
                                   if row.get("ellipsoidal_amplitude_ppm") else "")
        if step["number"] == 10:
            card["fpp"] = f"{row['fpp_mean']} +/- {row['fpp_standard_deviation']}"
        if step["number"] == 11:
            card["column"] = SUBMISSION_DECISION_COLUMNS[verdict]
            card["agent_summary"] = row.get("agent_summary", "")
            card["agent_reasons"] = json.loads(row["agent_reasons"] or "[]")
            if row["decisions_agree"] != "True":
                card["mechanical_decision"] = f"{row['mechanical_decision']}: {row['mechanical_reasons']}"
        if verdict not in step["passing_verdicts"]:
            card["column"], card["stop_reason"] = f"stopped_{step['number']}", stop_reason(step, row)
            break
    # A decision of a person is about a candidate that step 11 has decided. While a step checks the candidate again,
    # or after a step stops it, the card stays in the column of that step and shows the decision.
    decision = decisions.get(key)
    if decision and decision["decision"] in CANDIDATE_DECISIONS:
        card["decision"], card["decision_note"] = decision["decision"], decision.get("note", "")
        card["decision_date"] = decision.get("decision_date", "")
        if card["column"] in SUBMISSION_DECISION_COLUMNS.values():
            card["column"] = decision["decision"]
    return card


def board():
    statuses = step_statuses()
    decisions = candidate_decisions()
    cards = []
    for sector in searched_sectors():
        running = {(status["current_candidate"]["tic_id"], status["current_candidate"]["period_days"]):
                   {"step": status["step_number"], "progress": status["current_candidate"]["progress"],
                    "started_utc": status["current_candidate"]["started_utc"]}
                   for status in statuses
                   if status["sector"] == sector and status["state"] == "running" and status["current_candidate"]}
        step_rows = {step["number"]: rows_by_candidate(RESULTS_DIRECTORY / f"qlp_sector_{sector:04d}_{step['path']}")
                     for step in STEPS}
        new_candidates = complete_csv_rows(
            RESULTS_DIRECTORY / f"qlp_sector_{sector:04d}_known_catalogue_cross_match/new_candidates.csv")
        cards += [candidate_card(sector, candidate, step_rows, decisions, running) for candidate in new_candidates]
    columns = ([{"id": f"step_{step['number']}", "title": f"Waiting for or in step {step['number']}: {step['name']}"}
                for step in STEPS]
               + [{"id": "proposed_submission", "title": "Proposed for submission"},
                  {"id": "proposed_drop", "title": "Proposed for drop"}, {"id": "held", "title": "Held"},
                  {"id": "submitted", "title": "Submitted"}])
    stopped = ([{"id": f"stopped_{step['number']}", "title": f"Stopped at step {step['number']}: {step['name']}"}
                for step in STEPS] + [{"id": "dropped", "title": "Dropped"}])
    return {"columns": columns, "stopped_columns": stopped, "cards": cards, "steps": statuses,
            "sectors": sector_progress(statuses), "refresh_seconds": REFRESH_SECONDS}


class CandidateBoardRequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        path = unquote(urlparse(self.path).path)
        if path == "/":
            self.send(200, "text/html; charset=utf-8", BOARD_PAGE.encode())
        elif path == "/api/board":
            self.send(200, "application/json", json.dumps(board()).encode())
        elif path.startswith("/results/"):
            self.send_result_file(path[len("/results/"):])
        else:
            self.send(404, "text/plain", b"not found")

    def send_result_file(self, relative_path):
        # Serve only files inside results/.
        file_path = (RESULTS_DIRECTORY / relative_path).resolve()
        if RESULTS_DIRECTORY.resolve() not in file_path.parents or not file_path.is_file():
            self.send(404, "text/plain", b"not found")
            return
        content_type = ("text/plain; charset=utf-8" if file_path.suffix in (".md", ".csv")
                        else mimetypes.guess_type(file_path.name)[0] or "application/octet-stream")
        self.send(200, content_type, file_path.read_bytes())

    def send(self, code, content_type, body):
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *arguments):
        pass


BOARD_PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Candidate Board</title>
<style>
:root {
  --background: #f4f5f7; --surface: #ffffff; --text: #1d2330; --muted: #5d6677; --border: #d9dde5;
  --pass: #1f8a4c; --stop: #b3261e; --caution: #a15c00; --running: #2457c5; --running-background: #e8efff;
}
@media (prefers-color-scheme: dark) {
  :root {
    --background: #111318; --surface: #1b1e25; --text: #e6e8ee; --muted: #9aa3b5; --border: #2c313c;
    --pass: #5cc98a; --stop: #f08a82; --caution: #f0b860; --running: #8fb0ff; --running-background: #1d2a48;
  }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--background); color: var(--text);
       font: 14px/1.4 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
header { padding: 12px 16px; border-bottom: 1px solid var(--border); background: var(--surface); }
h1 { font-size: 18px; margin: 0 0 4px; }
.updated { color: var(--muted); font-size: 12px; }
.steps { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 8px; }
.step { border: 1px solid var(--border); border-radius: 6px; padding: 6px 10px; font-size: 12px; }
.step.running { border-color: var(--running); background: var(--running-background); }
.step.stopped { border-color: var(--stop); }
.bar { height: 4px; background: var(--border); border-radius: 2px; margin-top: 4px; }
.bar div { height: 4px; background: var(--running); border-radius: 2px; }
.sector { margin-top: 8px; font-size: 13px; }
.sector-step { display: inline-block; font-size: 12px; margin: 2px 4px 2px 0; padding: 1px 8px; border-radius: 9px;
               border: 1px solid var(--border); }
.sector-step.finished { border-color: var(--pass); color: var(--pass); }
.sector-step.running { border-color: var(--running); color: var(--running); background: var(--running-background); }
.sector-step.stopped { border-color: var(--stop); color: var(--stop); }
.sector-step.not-started { color: var(--muted); }
main { padding: 12px 16px; }
.board { display: flex; gap: 12px; overflow-x: auto; align-items: flex-start; padding-bottom: 8px; }
.column { flex: 0 0 260px; background: var(--surface); border: 1px solid var(--border); border-radius: 8px; }
.column h2 { font-size: 13px; margin: 0; padding: 8px 10px; border-bottom: 1px solid var(--border); }
.count { color: var(--muted); font-weight: normal; }
.cards { padding: 8px; display: flex; flex-direction: column; gap: 8px; min-height: 20px; }
.card { border: 1px solid var(--border); border-radius: 6px; padding: 8px; cursor: pointer; background: var(--surface); }
.card.running { border-color: var(--running); background: var(--running-background); }
.card .title { font-weight: 600; }
.card .line { color: var(--muted); font-size: 12px; }
.card .date { color: var(--text); font-size: 12px; }
.card .progress { color: var(--running); font-size: 12px; margin-top: 4px; }
.card .stop { color: var(--stop); font-size: 12px; }
.card .error { color: var(--stop); font-size: 12px; }
.tag { display: inline-block; font-size: 11px; padding: 1px 6px; border-radius: 9px; margin: 2px 2px 0 0;
       border: 1px solid var(--caution); color: var(--caution); }
.trail { font-size: 11px; margin-top: 4px; }
.trail span { color: var(--pass); margin-right: 6px; }
.trail span.failed { color: var(--stop); }
details.stopped { margin-top: 16px; }
details.stopped summary { cursor: pointer; font-weight: 600; }
.detail { position: fixed; inset: 0; background: rgba(0, 0, 0, 0.5); display: none; }
.detail.open { display: block; }
.detail .panel { position: absolute; top: 0; right: 0; bottom: 0; width: min(720px, 100%); background: var(--surface);
                 overflow-y: auto; padding: 16px; }
.detail pre { white-space: pre-wrap; font-size: 12px; background: var(--background); padding: 8px; border-radius: 6px; }
.dates { border-collapse: collapse; font-size: 12px; }
.dates th, .dates td { text-align: left; padding: 2px 12px 2px 0; border-bottom: 1px solid var(--border); }
.detail img { max-width: 100%; border: 1px solid var(--border); border-radius: 4px; margin-bottom: 8px; }
.detail button { float: right; font-size: 14px; }
@media (max-width: 640px) {
  .board { flex-direction: column; }
  .column { flex: 1 1 auto; width: 100%; }
}
</style>
</head>
<body>
<header>
  <h1>Candidate Board</h1>
  <div class="updated" id="updated">Loading</div>
  <div class="sectors" id="sectors"></div>
  <div class="steps" id="steps"></div>
</header>
<main>
  <div class="board" id="board"></div>
  <details class="stopped" id="stopped-section"><summary>Stopped candidates</summary>
    <div class="board" id="stopped"></div></details>
</main>
<div class="detail" id="detail"><div class="panel" id="panel"></div></div>
<script>
let latest = null;

function text(value) {
  const element = document.createElement("span");
  element.textContent = value ?? "";
  return element.innerHTML;
}

function cardHtml(card, index) {
  const lines = [`P ${text(card.fit_period_days || card.period_days)} d, depth ${text(card.depth_ppm)} ppm, Tmag ${text(card.tmag)}`];
  if (card.radius_earth_radii) lines.push(`radius ${text(card.radius_earth_radii)} Earth radii`);
  if (card.fpp) lines.push(`FPP ${text(card.fpp)}`);
  const trail = card.trail.map(item => `<span class="${item === card.trail[card.trail.length - 1] && card.stop_reason ? "failed" : ""}">${item.step}: ${text(item.verdict)}</span>`).join("");
  const last = card.trail[card.trail.length - 1];
  const date = last && last.date_utc ? `<div class="date">step ${last.step}: ${text(last.date_utc)} UTC</div>` : "";
  const running = card.running ? `<div class="progress">Step ${card.running.step} running since ${text(card.running.started_utc)} UTC${card.running.progress ? ": " + text(card.running.progress) : ""}</div>` : "";
  return `<div class="card ${card.running ? "running" : ""}" data-index="${index}">
    <div class="title">TIC ${text(card.tic_id)} <span class="line">S${card.sector}</span></div>
    ${date}
    ${lines.map(line => `<div class="line">${line}</div>`).join("")}
    <div class="trail">${trail}</div>
    ${card.stop_reason ? `<div class="stop">${text(card.stop_reason)}</div>` : ""}
    ${card.error ? `<div class="error">error: ${text(card.error)}</div>` : ""}
    ${card.agent_summary ? `<div class="stop">${text(card.agent_summary)}</div>` : ""}
    ${card.mechanical_decision ? `<div class="stop">rules disagree: ${text(card.mechanical_decision)}</div>` : ""}
    ${card.decision ? `<div class="line"><b>${text(card.decision)}</b>${card.decision_date ? " " + text(card.decision_date) : ""}${card.decision_note ? ": " + text(card.decision_note) : ""}</div>` : ""}
    <div>${card.cautions.map(caution => `<span class="tag">${text(caution)}</span>`).join("")}</div>
    ${running}</div>`;
}

function columnsHtml(columns, cards) {
  return columns.map(column => {
    const members = cards.map((card, index) => [card, index]).filter(([card]) => card.column === column.id);
    members.sort(([a], [b]) => (b.running ? 1 : 0) - (a.running ? 1 : 0));
    return `<section class="column"><h2>${text(column.title)} <span class="count">${members.length}</span></h2>
      <div class="cards">${members.map(([card, index]) => cardHtml(card, index)).join("")}</div></section>`;
  }).join("");
}

function stepsHtml(steps) {
  return steps.filter(step => step.state !== "finished").map(step => {
    const count = step.item_count ?? step.candidate_count, itemsDone = step.items_done ?? step.candidates_done;
    const done = count ? Math.round(100 * itemsDone / count) : 0;
    const current = step.current_candidate ? `, TIC ${text(step.current_candidate.tic_id)}${step.current_candidate.progress ? ": " + text(step.current_candidate.progress) : ""}` : "";
    return `<div class="step ${step.state}">${step.sector === null ? "All sectors" : "Sector " + step.sector}, step ${step.step_number} ${text(step.step_name.replaceAll("_", " "))}: ${step.state}, ${itemsDone.toLocaleString()} of ${count.toLocaleString()} ${text(step.item_name || "candidates")} done${current}
      <div class="bar"><div style="width: ${done}%"></div></div></div>`;
  }).join("") || `<div class="step">No step is running.</div>`;
}

function sectorsHtml(sectors) {
  return sectors.map(sector => {
    const steps = sector.steps.map(step => {
      const count = step.item_count ? `: ${step.items_done.toLocaleString()} of ${step.item_count.toLocaleString()} ${text(step.item_name)}` : "";
      const mark = {finished: "done", running: "running", stopped: "stopped", "not started": "not started"}[step.state];
      return `<span class="sector-step ${step.state.replace(" ", "-")}">${step.number} ${text(step.name)}: ${mark}${step.state === "finished" ? "" : count}</span>`;
    }).join("");
    const s = sector.summary;
    const summary = s.stars_searched ? `${s.stars_searched.toLocaleString()} stars searched (${s.stars_failed_to_process} failed), ${s["signals_vetted (SNR and transit cut)"].toLocaleString()} signals vetted, ${s.signals_passing_all_tests} passed vetting: ${s.signals_passing_all_tests_and_known} known, ${s.signals_passing_all_tests_and_new} new candidates. Catalogues of ${text(s.catalogue_download_date)}.` : "";
    return `<div class="sector"><b>Sector ${sector.sector}</b> ${steps}<div class="line">${summary}</div></div>`;
  }).join("");
}

async function refresh() {
  try {
    const response = await fetch("/api/board", {cache: "no-store"});
    latest = await response.json();
    document.getElementById("board").innerHTML = columnsHtml(latest.columns, latest.cards);
    document.getElementById("stopped").innerHTML = columnsHtml(latest.stopped_columns, latest.cards);
    document.getElementById("steps").innerHTML = stepsHtml(latest.steps);
    document.getElementById("sectors").innerHTML = sectorsHtml(latest.sectors);
    document.getElementById("updated").textContent = `Updated ${new Date().toLocaleTimeString()}, every ${latest.refresh_seconds} s`;
  } catch (error) {
    document.getElementById("updated").textContent = `Update failed at ${new Date().toLocaleTimeString()}: ${error}`;
  }
  setTimeout(refresh, (latest?.refresh_seconds || 10) * 1000);
}

async function openCard(index) {
  const card = latest.cards[index];
  const panel = document.getElementById("panel");
  const record = card.record ? await (await fetch("/results/" + card.record, {cache: "no-store"})).text() : "No verification record yet.";
  panel.innerHTML = `<button id="close">Close</button><h2>TIC ${text(card.tic_id)}, P ${text(card.period_days)} d, Sector ${card.sector}</h2>
    <h3>Steps</h3><table class="dates"><tr><th>step</th><th>verdict</th><th>date (UTC)</th></tr>
      ${card.trail.map(item => `<tr><td>${item.step}</td><td>${text(item.verdict)}</td><td>${text(item.date_utc)}</td></tr>`).join("")}
      ${card.decision ? `<tr><td>decision</td><td>${text(card.decision)}</td><td>${text(card.decision_date)}</td></tr>` : ""}</table>
    ${card.ellipsoidal ? `<p>Ellipsoidal variation: ${text(card.ellipsoidal)}</p>` : ""}
    ${card.agent_reasons ? `<h3>Step 11 reasons</h3><ul>${card.agent_reasons.map(reason => `<li>${text(reason)}</li>`).join("")}</ul>` : ""}
    ${card.mechanical_decision ? `<p>The rules decide otherwise: ${text(card.mechanical_decision)}</p>` : ""}
    <p><a href="https://exofop.ipac.caltech.edu/tess/target.php?id=${text(card.tic_id)}" target="_blank" rel="noopener">ExoFOP target page</a></p>
    <h3>Verification record</h3><pre>${text(record)}</pre>
    <h3>Figures</h3>${card.figures.map(path => `<div class="line">${text(path)}</div><img src="/results/${encodeURI(path)}" loading="lazy">`).join("") || "None yet."}`;
  document.getElementById("detail").classList.add("open");
  document.getElementById("close").onclick = () => document.getElementById("detail").classList.remove("open");
}

document.addEventListener("click", event => {
  const card = event.target.closest(".card");
  if (card) openCard(Number(card.dataset.index));
  if (event.target.id === "detail") event.target.classList.remove("open");
});
refresh();
</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--host", default="0.0.0.0", help="the address to listen on; 0.0.0.0 is every network")
    parser.add_argument("--port", type=int, default=8000)
    arguments = parser.parse_args()
    server = ThreadingHTTPServer((arguments.host, arguments.port), CandidateBoardRequestHandler)
    print(f"Serving the candidate board on http://{arguments.host}:{arguments.port}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
