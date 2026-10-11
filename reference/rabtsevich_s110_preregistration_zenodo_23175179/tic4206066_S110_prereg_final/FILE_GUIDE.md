# File guide

Every file of this deposit and what it is for. Upper-case words mark the files that carry the predictions,
the rules (code) and the reproduction. `SHA256SUMS` lists the sha256 of every file except itself.

## Read first

| file | purpose |
|---|---|
| `CHANGELOG.md` | history of the rules from the first internal draft to this version |
| `CITATION.cff` | how to cite (machine-readable) |
| `FILE_GUIDE.md` | this list |
| `LICENSE` | licences: MIT for *.py and *.sh, CC BY 4.0 for everything else |
| `PREREG_S110_FINAL.pdf` | THE PRE-REGISTRATION: question, how the signals were found, predictions, decision rules (Sections 1 to 7), and the full rules and commitments (Appendices A to N) |
| `PREREG_S110_SHORT.pdf` | short reader-oriented overview (three pages, with one figure); start here, then use the full protocol for execution |
| `protocol_source/` | LaTeX source of both PDFs (main files `PREREG_S110_FINAL.tex` and `PREREG_S110_SHORT.tex`, both use class `aastex702.cls`; the full protocol also uses the bibliography `references.bib`); same text as the PDF |
| `protocol_source/S110_transits.pdf`, `protocol_source/make_figure1.py` | Figure 1 of the summary (archival transits of A and B per sector, folded on the frozen ephemeris, and the empty Sector 110 panel) and the script that draws it from `calib/expected_depths.json`, `ephemeris_final.json` and the three archival light curves with the frozen detrending (inputs checked against their sha256) |
| `PRF_SHA256SUMS` | sha256 of the 25 public TESS PRF model files that the localisation must use (files not included) |
| `README.md` | start here: what this is, what to read, how to reproduce |
| `SHA256SUMS` | sha256 of every other file of the deposit |
| `ephemeris_final.json` | PREDICTIONS: frozen periods, reference times, durations, covariance and the 12 predicted Sector 110 mid-times with their uncertainties |
| `requirements.txt` | binding software pin of every registered run (Python 3.12.13, package versions, arm64, macOS) |
| `requirements_lock.txt` | full environment with distribution hashes, in which all outputs were reproduced on 2026-10-03 |
| `s110_final.py` | CODE: the photometric test (counting rule, depth statistic, noise model, thresholds, outcomes O1 to O4, timing) |
| `s110_localise.py` | CODE: pixel localisation of the dip (labels L1 to L3 for signal A) |
| `s110_summary.py` | CODE: turns the output files into the only sentences allowed to describe the outcomes; refuses incomplete inputs |

## calib/: predictions and their inputs, calibrations

| file | purpose |
|---|---|
| `calib/ephemeris_inputs.json` | PREDICTION INPUT: how ephemeris_final.json was built (A fit, T14, B ephemeris and sigmas) |
| `calib/ephemeris_sources/README_sources.md` | explains the files of this folder and how each ephemeris number follows from them |
| `calib/ephemeris_sources/chain_B_pipe_mann_frozen_w0.75.npz` | posterior samples of the joint fit of B (source of the B ephemeris and its sigmas) |
| `calib/ephemeris_sources/mcmc_samples.npz` | posterior samples of the transit model of A (source of T14) |
| `calib/ephemeris_sources/timing_injections.csv` | timing injections per transit of A (source of the transit-time errors) |
| `calib/ephemeris_sources/timing_injections_deltas.npz` | every timing injection of A (recovered minus injected time) |
| `calib/expected_depths.json` | PREDICTIONS: expected depths D (Sectors 6 + 32) and D98 (Sector 98), made by s110_final.py --calibrate |
| `calib/false_alarm_A.json` | false-alarm trials, A-like fake ephemerides, Sectors 6, 32, 98 |
| `calib/false_alarm_A_S32.json` | false-alarm trials, A-like, SPOC 2-min Sector 32 (the primary product type) |
| `calib/false_alarm_A_S32_a3.json` | false-alarm trials, A-like, the Sector 32 set first run by audit 3 (kept so that every set run is pooled) |
| `calib/false_alarm_B.json` | false-alarm trials, B-like fake ephemerides, Sectors 32, 98 |
| `calib/false_alarm_B_S32.json` | false-alarm trials, B-like, SPOC 2-min Sector 32 |
| `calib/false_alarm_B_S32_a3.json` | false-alarm trials, B-like, the Sector 32 set first run by audit 3 |
| `calib/false_alarm_levels.json` | measured false-detection levels per product and number of transits (quoted by s110_summary.py) |
| `calib/false_alarm_summary.txt` | console summary of the six false-alarm trial sets |
| `calib/localisation_archival.json` | archival pixel positions of A (Sectors 32 and 98) for the descriptive combined position |
| `calib/localisation_coverage.json` | coverage check of the empirical localisation p (made by build/localisation_coverage.py) |
| `calib/transit_times_A.csv` | PREDICTION INPUT: the 22 archival transit times of A with their errors |

## deps/: frozen helper code

| file | purpose |
|---|---|
| `deps/lcio.py` | CODE: frozen light-curve loader and detrender used by s110_final.py |
| `deps/pixloc_lib.py` | CODE: frozen pixel-localisation library used by s110_localise.py (its default file paths point next to the deposit and are replaced at run time) |

## build/: how the outputs were made, reproduction

| file | purpose |
|---|---|
| `build/calibrate.log` | console output of s110_final.py --calibrate |
| `build/compare_outputs.py` | REPRODUCTION: compares a rebuilt copy with the deposited outputs |
| `build/false_alarm_check.py` | makes the false-alarm trial files |
| `build/false_alarm_levels.py` | makes calib/false_alarm_levels.json from the trial files |
| `build/fetch_archival_inputs.sh` | REPRODUCTION: downloads the public TESS files and PRF models and checks their sha256 |
| `build/localisation_coverage.log` | console output of build/localisation_coverage.py |
| `build/localisation_coverage.py` | makes calib/localisation_coverage.json |
| `build/make_ephemeris.py` | how ephemeris_final.json was written (reads the author's archival folders; the same inputs are in calib/) |
| `build/reproduce.sh` | REPRODUCTION: one command that copies the deposit, runs build/run_all.sh and compares every output |
| `build/run_all.log` | console output of build/run_all.sh |
| `build/run_all.sh` | rebuilds calib/expected_depths.json, calib/localisation_coverage.json and dryrun/ (expects the archival inputs next to the deposit, as made by build/fetch_archival_inputs.sh; build/reproduce.sh runs it in a copy) |

## dryrun/: the frozen code run on archival data (validation only)

| file | purpose |
|---|---|
| `dryrun/loc_s32_B.json` | localisation of B on Sector 32: not possible (B not detected there) |
| `dryrun/loc_s32_B.log` | console output of dryrun/loc_s32_B.json |
| `dryrun/loc_s32_coverage_seed5005.json` | localisation dry run with another calibration seed (input of the coverage check) |
| `dryrun/loc_s32_final_pipeline.json` | localisation dry run of A on Sector 32 pixels, final pipeline (would be L1) |
| `dryrun/loc_s32_final_pipeline.log` | console output of dryrun/loc_s32_final_pipeline.json |
| `dryrun/loc_s32_items10_fidelity.json` | localisation run that reproduces the archival item 10/11 result exactly |
| `dryrun/loc_s32_items10_fidelity.log` | console output of dryrun/loc_s32_items10_fidelity.json |
| `dryrun/loc_s32_null_gate.json` | localisation after the null test: not possible (A not detected) |
| `dryrun/loc_s32_null_gate.log` | console output of dryrun/loc_s32_null_gate.json |
| `dryrun/loc_s32_refused.log` | REFUSAL CHECK: a labelled localisation on archival files must be refused (it is) |
| `dryrun/loc_s6_ccd1_refused.log` | REFUSAL CHECK: a file from another sector and CCD must be refused before anything is written (it is) |
| `dryrun/na_not_run.txt` | TEST FIXTURE: text used to declare variants not available in the dry-run code checks |
| `dryrun/phot_s32_refused.log` | REFUSAL CHECK: a registered photometric run on an archival file must be refused (it is) |
| `dryrun/s32_desat1.json` | dry run, Sector 32, variant (c) desaturation cadences removed |
| `dryrun/s32_desat1.log` | console output of dryrun/s32_desat1.json |
| `dryrun/s32_null.json` | dry run, Sector 32, ephemeris shifted by half a period, real transits cut (null test: O4) |
| `dryrun/s32_null.log` | console output of dryrun/s32_null.json |
| `dryrun/s32_real.json` | dry run, Sector 32, primary settings (A: O1) |
| `dryrun/s32_real.log` | console output of dryrun/s32_real.json |
| `dryrun/s32_sap.json` | dry run, Sector 32, variant (e) SAP flux |
| `dryrun/s32_sap.log` | console output of dryrun/s32_sap.json |
| `dryrun/s32_sap_other_version.json` | TEST FIXTURE: copy of s32_sap.json with its script hash set to zeros on purpose (the summary must refuse it) |
| `dryrun/s32_v2rule.json` | dry run, Sector 32, variant (f) old counting rule |
| `dryrun/s32_v2rule.log` | console output of dryrun/s32_v2rule.json |
| `dryrun/s32_w05.json` | dry run, Sector 32, variant (a) window 0.5 d |
| `dryrun/s32_w05.log` | console output of dryrun/s32_w05.json |
| `dryrun/s32_w10.json` | dry run, Sector 32, variant (b) window 1.0 d |
| `dryrun/s32_w10.log` | console output of dryrun/s32_w10.json |
| `dryrun/s98_real.json` | dry run, Sector 98, QLP fallback settings (A and B: O1) |
| `dryrun/s98_real.log` | console output of dryrun/s98_real.json |
| `dryrun/s98_v2rule.json` | dry run, Sector 98, variant (f) old counting rule |
| `dryrun/s98_v2rule.log` | console output of dryrun/s98_v2rule.json |
| `dryrun/summary_d_final_refused.log` | REFUSAL CHECK: the summary must refuse 'QLP not public' as final before 2027-02-28 (it does) |
| `dryrun/summary_na_refused.log` | REFUSAL CHECK: the summary must refuse 'not available' without a saved error output (it does) |
| `dryrun/summary_refused.log` | REFUSAL CHECK: the summary must refuse when a variant is silently left out (it does) |
| `dryrun/summary_s32.json` | fixed sentences produced from the Sector 32 dry runs |
| `dryrun/summary_s32.log` | console output of dryrun/summary_s32.json |
| `dryrun/summary_s32_null.json` | fixed sentences produced from the null test |
| `dryrun/summary_s32_null.log` | console output of dryrun/summary_s32_null.json |
| `dryrun/summary_s98.json` | fixed sentences produced from the Sector 98 dry run |
| `dryrun/summary_s98.log` | console output of dryrun/summary_s98.json |
| `dryrun/summary_version_refused.log` | REFUSAL CHECK: the summary must refuse a variant made by another script version (it does) |

## Not in this deposit

The unpublished internal drafts v1 (2026-09-27) and v2 (2026-09-30) are not part of this deposit. They are not needed to reproduce anything in this deposit, and they never define the
rules. `calib/ephemeris_inputs.json` still notes that the B ephemeris is also in v1; the values are in that file itself.
