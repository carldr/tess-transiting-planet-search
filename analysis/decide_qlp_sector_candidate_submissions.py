"""Propose each candidate of one sector that passed TRICERATOPS for submission as a CTOI or for drop, or hold it.

The candidates are the rows of triceratops_validation.csv with the verdict "likely planet" or "validated without
imaging". Steps 7 to 10 stop a candidate only on a refutation. Their cautions never stop it, so this step weighs them.
It has two parts:
1. mechanical rules give a decision with its reasons, in the order of MECHANICAL_RULES below;
2. a Claude agent, `claude -p --model fable`, gives the final decision with its reasons. Its prompt holds the results
   of steps 7 to 10, the verification record, the ExoFOP entries, the unused released sectors, the sectors that observe the
   star and are not released yet, the paths of the
   figures and the mechanical decision. The agent has the Read tool only, so it can read the files and look at the
   figures, and change nothing.

The decision is one of:
- "proposed for submission": the candidate is a plausible new planet, ready for a CTOI;
- "proposed for drop": the candidate is already known, or most likely not a planet;
- "held": the data cannot decide yet, for example until a new sector is released.
A person accepts or rejects a proposal in results/candidate_decisions.csv. This script submits nothing.

MECHANICAL_RULES:
1. proposed for drop: ExoFOP has a TOI or a CTOI for the star, or
   reference/published_transit_candidates_not_in_catalogues.csv lists the star with a period within
   PUBLISHED_PERIOD_TOLERANCE;
2. proposed for drop: the fitted radius is above COMPANION_MAXIMUM_RADIUS_EARTH_RADII and a caution in
   STELLAR_COMPANION_CAUTIONS also points to a stellar companion;
3. held: the pixel verdict is "unresolved", a caution is in HOLD_CAUTIONS, or MAST has a released light curve of a
   sector that step 9 did not use;
4. proposed for submission: any other result.

Run from the repository root, after validate_qlp_sector_candidates_with_triceratops.py:
    .venv/bin/python analysis/decide_qlp_sector_candidate_submissions.py 98

Output: results/qlp_sector_<nnnn>_submission_decisions/
- submission_decisions.csv, one row per candidate;
- tic<id>_<period>_agent_prompt.md and tic<id>_<period>_agent_reply.json, the prompt and the full reply of the agent.
The script writes each row as soon as it has it. It keeps the rows of an earlier run without an error, and decides
only the other candidates. --restart decides all candidates again.
The script also appends the decision of each candidate to its verification record.
"""

import argparse
import csv
import json
import subprocess

from tess_stars2px import tess_stars2px_function_entry

from build_qlp_sector_star_list import REPOSITORY_ROOT
from candidate_verification_records import (
    append_verification_record_entry, verification_record_path, verification_stamp,
)
from check_candidates_on_exofop_and_future_tess_sectors import download_exofop_target_page
from check_qlp_sector_new_candidates import (
    COMPANION_MAXIMUM_RADIUS_EARTH_RADII, other_sector_products, output_directory as candidate_checks_output_directory,
)
from check_qlp_sector_candidate_pixels import output_directory as pixel_checks_output_directory
from incremental_csv_output import IncrementalCsvWriter, earlier_rows_to_keep
from pipeline_status import PipelineStatus, report_candidate_progress
from validate_qlp_sector_candidates_with_triceratops import (
    output_directory as triceratops_output_directory,
)
from verify_qlp_sector_candidates_over_all_sectors import (
    MCMC_MAXIMUM_STEPS, output_directory as all_sector_verification_output_directory,
)

PUBLISHED_CANDIDATES_PATH = REPOSITORY_ROOT / "reference/published_transit_candidates_not_in_catalogues.csv"
PUBLISHED_PERIOD_TOLERANCE = 0.01
TRICERATOPS_VERDICTS_TO_DECIDE = {"likely planet", "validated without imaging"}
STELLAR_COMPANION_CAUTIONS = {"Gaia radial velocity variable", "RUWE above 1.4", "ellipsoidal variation detected",
                              "secondary eclipse detected"}
# A dip longer than a transit, or a phase 0.5 brighter than a planet can make it, can mean a sinusoid and not a
# transit. An eccentric orbit or stellar variability can also give them, so they hold the candidate.
HOLD_CAUTIONS = {"ellipsoidal variation detected", f"transit fit not converged in {MCMC_MAXIMUM_STEPS} steps",
                 "dip longer than a transit at this period", "phase 0.5 brighter than a planet can make it",
                 "star varies at half the period"}
DECISIONS = ["proposed for submission", "proposed for drop", "held"]
AGENT_MODEL = "fable"
AGENT_TIMEOUT_SECONDS = 1800
AGENT_DECISION_SCHEMA = {
    "type": "object",
    "properties": {
        "decision": {"type": "string", "enum": DECISIONS},
        "summary": {"type": "string"},
        "reasons": {"type": "array", "items": {"type": "string"}, "minItems": 1},
    },
    "required": ["decision", "summary", "reasons"],
    "additionalProperties": False,
}
AGENT_SYSTEM_PROMPT = f"""You review one TESS transiting-planet candidate and decide what to do with it.

The candidate comes from a search of one TESS sector of Quick-Look Pipeline light curves. Ten automatic steps have
checked it: the search, LEO-Vetter, a cross-match with the known catalogues, independent checks (step 7), a
difference-image check of the source of the dip (step 8), a transit fit and vetting over all TESS sectors (step 9),
and TRICERATOPS (step 10). Each step stops a candidate only when a check refutes it. A caution is a weaker warning
that never stops a candidate. Your task is to weigh the cautions and the values together.

Give one decision:
- "proposed for submission": the candidate is a plausible new planet and is ready to submit to ExoFOP as a CTOI.
  A CTOI does not need to be a confirmed planet, but it must be a candidate that the data support and that is new.
- "proposed for drop": the star already has a TOI or CTOI at this period, a publication already reports the signal,
  or the evidence shows that the signal is most likely not a planet, such as a stellar companion or a blend.
- "held": the data cannot decide yet, and more data can decide it, such as a released TESS sector that is not yet
  used, or a sector that TESS observes and has not released yet (unreleased_tess_sectors_observing_the_star).
  Give in a reason what would decide it.

A person reads your decision and accepts or rejects it. You submit nothing.

Give a summary: the main reason for the decision in one sentence of at most 15 words, with its key value.
For example: "Secondary eclipse is 49% of the transit depth; Sector 110 can test it."

Rules for the reasons:
- Give one reason per item. Cite each value that you use, with its name as the input gives it.
- Use only the values in the input and in the files it names. Do not assume a value that is not there.
- Say where you disagree with the mechanical decision, and why.
- A radius above {COMPANION_MAXIMUM_RADIUS_EARTH_RADII:.1f} Earth radii (2 Jupiter radii) is a star, but a grazing
  transit leaves the radius almost unconstrained.

You have the Read tool only. Use it to read the verification record and to look at the figures that the input names.
"""

OUTPUT_COLUMNS = [
    "tic_id", "period_days", "tmag", "pixel_verdict", "triceratops_verdict", "fpp_mean", "planet_radius_earth_radii",
    "cautions", "exofop_tois", "exofop_ctois", "published_signal", "released_tess_sectors_not_in_all_sector_verification",
    "unreleased_tess_sectors_observing_the_star", "mechanical_decision", "mechanical_reasons", "agent_decision", "agent_summary", "agent_reasons", "agent_model", "agent_cost_usd",
    "decision", "decisions_agree", "error", "verification_date_utc", "code_version",
]


def output_directory(sector):
    return REPOSITORY_ROOT / f"results/qlp_sector_{sector:04d}_submission_decisions"


def rows_by_candidate(path):
    return {(row["tic_id"], row["period_days"]): row for row in csv.DictReader(open(path))}


def published_signal(tic_id, period_days):
    """Return the source of a published signal of the star at the period, or "" when there is none."""
    for row in csv.DictReader(open(PUBLISHED_CANDIDATES_PATH)):
        if (int(row["tic_id"]) == int(tic_id)
                and abs(float(period_days) / float(row["period_days"]) - 1) < PUBLISHED_PERIOD_TOLERANCE):
            return f"P {row['period_days']} d, {row['source']}"
    return ""


def cautions_of(*rows):
    cautions = []
    for row in rows:
        cautions += [caution for caution in row.get("cautions", "").split("; ") if caution and caution not in cautions]
    return cautions


def mechanical_decision(facts):
    """Return the decision of MECHANICAL_RULES and its reasons."""
    if facts["exofop_tois"] or facts["exofop_ctois"]:
        return "proposed for drop", [f"ExoFOP has TOI [{facts['exofop_tois']}] CTOI [{facts['exofop_ctois']}]"]
    if facts["published_signal"]:
        return "proposed for drop", [f"already published: {facts['published_signal']}"]
    radius = float(facts["planet_radius_earth_radii"] or "nan")
    stellar_cautions = [caution for caution in facts["cautions"] if caution in STELLAR_COMPANION_CAUTIONS]
    if radius > COMPANION_MAXIMUM_RADIUS_EARTH_RADII and stellar_cautions:
        return "proposed for drop", [f"radius {radius:.1f} Earth radii is above "
                                     f"{COMPANION_MAXIMUM_RADIUS_EARTH_RADII:.1f}, with " + ", ".join(stellar_cautions)]
    reasons = []
    if facts["pixel_verdict"] == "unresolved":
        reasons.append("pixel verdict unresolved: the dip can be on a nearby star")
    reasons += [f"caution: {caution}" for caution in facts["cautions"] if caution in HOLD_CAUTIONS]
    if facts["released_tess_sectors_not_in_all_sector_verification"]:
        reasons.append("released sectors not used in step 9: "
                       + facts["released_tess_sectors_not_in_all_sector_verification"])
    if reasons:
        return "held", reasons
    return "proposed for submission", ["no rule holds or drops it"]


def candidate_facts(candidate, rows_of_steps_7_to_9):
    """Return the values that the mechanical rules and the agent use, from steps 7 to 10, ExoFOP and MAST."""
    candidate_checks, pixel_checks, verification = rows_of_steps_7_to_9
    page = download_exofop_target_page(candidate["tic_id"])
    sectors_used = {int(sector) for sector in candidate["sectors_used"].split()}
    ra, dec = float(page["coordinates"]["ra"]), float(page["coordinates"]["dec"])
    observing = {int(sector) for sector in tess_stars2px_function_entry(int(candidate["tic_id"]), ra, dec)[3]}
    # Sector 0 does not exist, so other_sector_products returns every sector with a light curve.
    released = {sector for sector, _, _ in other_sector_products(int(candidate["tic_id"]), 0)}
    return {
        "tic_id": candidate["tic_id"], "period_days": candidate["period_days"], "tmag": candidate["tmag"],
        "pixel_verdict": candidate["pixel_verdict"], "triceratops_verdict": candidate["triceratops_verdict"],
        "fpp_mean": candidate["fpp_mean"], "planet_radius_earth_radii": verification["planet_radius_earth_radii"],
        "cautions": cautions_of(candidate_checks, verification),
        "exofop_tois": " ".join(entry["toi"] for entry in page.get("tois", [])),
        "exofop_ctois": " ".join(entry["ctoi"] for entry in page.get("ctois", [])),
        "published_signal": published_signal(candidate["tic_id"], candidate["fit_period_days"]),
        "released_tess_sectors_not_in_all_sector_verification": " ".join(
            str(sector) for sector in sorted(released - sectors_used)),
        "unreleased_tess_sectors_observing_the_star": " ".join(
            str(sector) for sector in sorted(observing - released - sectors_used)),
    }


def figure_paths(sector, tic_id):
    paths = []
    for directory in sorted(REPOSITORY_ROOT.glob(f"results/qlp_sector_{sector:04d}_*")):
        if not directory.name.endswith("_quality_zero_mask"):
            paths += sorted(directory.glob(f"tic{tic_id}_*.png"))
    return paths


def agent_prompt(sector, candidate, rows_of_steps_7_to_9, facts, mechanical):
    record_path = verification_record_path(candidate["tic_id"], candidate["period_days"])
    step_rows = dict(zip(["step 7, independent checks", "step 8, pixel check", "step 9, all-sector verification"],
                         rows_of_steps_7_to_9))
    step_rows["step 10, TRICERATOPS"] = candidate
    sections = [
        f"# TIC {candidate['tic_id']}, Sector {sector}, period {candidate['period_days']} d",
        "## Summary\n\n```json\n" + json.dumps({key: value for key, value in facts.items()}, indent=1) + "\n```",
        f"## Mechanical decision\n\n{mechanical[0]}\n\n" + "\n".join(f"- {reason}" for reason in mechanical[1]),
        *[f"## Result row of {name}\n\n```json\n" + json.dumps(row, indent=1) + "\n```"
          for name, row in step_rows.items()],
        f"## Files to read\n\n- Verification record: {record_path}\n"
        + "\n".join(f"- Figure: {path}" for path in figure_paths(sector, candidate["tic_id"])),
    ]
    return "\n\n".join(sections) + "\n"


def run_agent(prompt):
    """Run the agent on the prompt, and return its reply as the JSON of `claude -p --output-format json`.

    The agent loads no settings, memory, MCP servers or skills of this computer, so only the system prompt and the
    prompt instruct it.
    """
    command = [
        "claude", "-p", "--model", AGENT_MODEL, "--system-prompt", AGENT_SYSTEM_PROMPT,
        "--tools", "Read", "--permission-mode", "dontAsk",
        "--setting-sources", "", "--settings", json.dumps({"autoMemoryEnabled": False}),
        "--strict-mcp-config", "--disable-slash-commands", "--no-session-persistence",
        "--output-format", "json", "--json-schema", json.dumps(AGENT_DECISION_SCHEMA),
    ]
    completed = subprocess.run(command, input=prompt, capture_output=True, text=True, cwd=REPOSITORY_ROOT,
                               timeout=AGENT_TIMEOUT_SECONDS)
    if completed.returncode != 0:
        raise RuntimeError(f"claude exited with code {completed.returncode}: {completed.stderr.strip()[:500]}")
    reply = json.loads(completed.stdout)
    if reply.get("is_error") or not reply.get("structured_output"):
        raise RuntimeError(f"claude gave no decision: {reply.get('subtype')}: {str(reply.get('result'))[:500]}")
    return reply


def decide(sector, candidate, rows_of_steps_7_to_9, directory):
    facts = candidate_facts(candidate, rows_of_steps_7_to_9)
    mechanical = mechanical_decision(facts)
    prompt = agent_prompt(sector, candidate, rows_of_steps_7_to_9, facts, mechanical)
    stem = f"tic{candidate['tic_id']}_{float(candidate['period_days']):.5f}"
    (directory / f"{stem}_agent_prompt.md").write_text(prompt)
    report_candidate_progress(f"agent {AGENT_MODEL} deciding; mechanical decision {mechanical[0]}")
    reply = run_agent(prompt)
    (directory / f"{stem}_agent_reply.json").write_text(json.dumps(reply, indent=1) + "\n")
    decision = reply["structured_output"]
    return {
        **facts, "cautions": "; ".join(facts["cautions"]),
        "mechanical_decision": mechanical[0], "mechanical_reasons": "; ".join(mechanical[1]),
        "agent_decision": decision["decision"], "agent_summary": decision["summary"], "agent_reasons": json.dumps(decision["reasons"]),
        "agent_model": " ".join(reply.get("modelUsage", {})) or AGENT_MODEL,
        "agent_cost_usd": f"{reply.get('total_cost_usd', 0):.3f}",
        "decision": decision["decision"], "decisions_agree": decision["decision"] == mechanical[0],
    }


def verification_record_result(row):
    reasons = "; ".join(json.loads(row["agent_reasons"]))
    return (f"agent {row['agent_model']}: {row['agent_decision']}: {row['agent_summary']} Reasons: {reasons}. "
            f"Mechanical: {row['mechanical_decision']}: {row['mechanical_reasons']}")


def append_verification_record(sector, row):
    """Append the row of a candidate, as written to submission_decisions.csv, to its verification record.

    A row with an error gets no record row: the next run decides that candidate again.
    """
    if row["error"]:
        return
    append_verification_record_entry(sector, row["tic_id"], row["period_days"], "submission decision",
                                     row["decision"], verification_record_result(row), row)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sector", type=int)
    parser.add_argument("--restart", action="store_true",
                        help="decide all candidates again, not only those without a row in submission_decisions.csv")
    parser.add_argument("--tic-id", dest="tic_ids", action="append", default=[],
                        help="decide only the candidates of this TIC ID; give it more than once for more TIC IDs")
    arguments = parser.parse_args()
    sector = arguments.sector

    directory = output_directory(sector)
    directory.mkdir(parents=True, exist_ok=True)
    output_path = directory / "submission_decisions.csv"
    validations = csv.DictReader(open(triceratops_output_directory(sector) / "triceratops_validation.csv"))
    passed = [row for row in validations if row["triceratops_verdict"] in TRICERATOPS_VERDICTS_TO_DECIDE]
    step_rows = [rows_by_candidate(candidate_checks_output_directory(sector) / "candidate_checks.csv"),
                 rows_by_candidate(pixel_checks_output_directory(sector) / "pixel_checks.csv"),
                 rows_by_candidate(all_sector_verification_output_directory(sector) / "all_sector_verification.csv")]
    earlier_rows = earlier_rows_to_keep(
        output_path, arguments.restart, {(row["tic_id"], row["period_days"]) for row in passed},
        tic_ids=arguments.tic_ids)
    decided = {(row["tic_id"], row["period_days"]) for row in earlier_rows}
    candidates = [row for row in passed if (row["tic_id"], row["period_days"]) not in decided
                  and (not arguments.tic_ids or row["tic_id"] in arguments.tic_ids)]
    print(f"Sector {sector}: {len(passed)} candidates passed TRICERATOPS, {len(earlier_rows)} decided in an "
          f"earlier run", flush=True)

    # A stop between a CSV row and its record row leaves the record row missing. This appends it.
    for row in earlier_rows:
        append_verification_record(sector, row)

    status = PipelineStatus(sector, 11, "submission_decision", "decide_qlp_sector_candidate_submissions.py",
                            len(candidates))
    with IncrementalCsvWriter(output_path, OUTPUT_COLUMNS, earlier_rows) as writer, status:
        for count, candidate in enumerate(candidates, start=1):
            status.start_candidate(candidate["tic_id"], candidate["period_days"])
            key = (candidate["tic_id"], candidate["period_days"])
            row = {column: candidate.get(column, "") for column in OUTPUT_COLUMNS}
            # The row of step 10 also has the column error. This step gives its own.
            row["error"] = ""
            try:
                row.update(decide(sector, candidate, [rows[key] for rows in step_rows], directory))
            except Exception as error:
                row["error"] = f"{type(error).__name__}: {error}"
            row.update(verification_stamp())
            row = writer.write_row(row)
            append_verification_record(sector, row)
            status.finish_candidate()
            print(f"  {count} of {len(candidates)}: TIC {row['tic_id']} {row['decision'] or 'error'}: "
                  f"{row['error'] or verification_record_result(row)}", flush=True)
    print(f"Wrote {output_path.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()
