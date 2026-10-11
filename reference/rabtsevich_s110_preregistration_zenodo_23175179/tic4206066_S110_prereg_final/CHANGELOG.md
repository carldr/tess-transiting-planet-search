# CHANGELOG: TIC 4206066 Sector 110 pre-registration

The internal/independent package checks below were performed by separate AI review sessions.
This does not concern outside astronomers who supplied data or comments on the target.

- **v1, 2026-09-27** (internal draft, not published): predicted times, box statistic, expected depths, fixed
  thresholds from Sector 98 noise, two labels (DETECT / REFUTE), QLP as the decisive product. Not published.
- **v2, 2026-09-30** (internal draft, not published): four mutually exclusive outcomes; SPOC 2-min primary; one quality
  mask per product; thresholds from the analysed data; sensitivity requirement for "excluded"; timing and localisation
  rules; processing variants and a robustness rule; frozen script `s110_test.py`. Revised the same day after a ten-pass
  internal audit. Not published.
- **2026-10-01**: final re-derivation of all archival numbers from raw data. New ephemeris of A (about 3.0 min earlier
  in Sector 110, calibrated prediction sigma 4.1 min); B unchanged.
- **2026-10-02**: independent read-only check of v2 against the re-derivation. Found: pooled vs per-transit depth
  estimators, no protection against edge transits, sigma_meas = 0 for one transit, individual times and input hashes
  promised but not produced, localisation code not frozen by the package, a manifest path mismatch.
- **2026-10-02, adversarial audit of a first draft of the final version**: the noise windows were not processed like
  the real windows (stated 0.3 % level really about 1 to 6 % on archival data; v1 and v2 had the same flaw); product
  dates, crash handling, localisation settings, the p convention and the environment left choices after the data.
  Fixed: noise pool from fake ephemerides through the full pipeline; edge zones of W/2 at the ends of every continuous
  segment (transits there are "not tested"); fixed product dates and procedure; the localisation script refuses
  non-registered settings; pinned environment.
- **2026-10-02, adversarial audit 2**: realised false-alarm level still above the stated 0.3 %; edge-zone segments not
  the detrender's own; dry-run bypass of the localisation gates; crashes on O4 paths; W defined inconsistently.
  Fixed: det / lo / hi at 99.9 / 0.1 / 99.9 percent with the realised levels measured and reported; segments taken from
  the series the detrender is given; gates and crash paths closed; "code decides" clause.
- **2026-10-02, adversarial audit 3**: localisation fake stacks smaller than the real stack (biased toward
  "compatible"); fixed wording incomplete; false-detection level of A averaged over products. Fixed: equal stacks;
  `s110_summary.py` produces the only allowed sentences, with fixed O4 reasons and a computed robustness qualifier;
  realised levels per product and per number of transits; environment and input checks for every registered run;
  QLP fallback treated as a registered run; exact PRF file set; reporting deadlines.
- **2026-10-02, adversarial audit 4**: the equal-stack rule kept only phases where enough fakes fitted, forcing noisy
  fakes next to gaps into the stacks; the A false-detection level came from one trial set and a grouping of n chosen
  after seeing the counts. Fixed: fakes obey the photometric edge-zone rule and are completed from further phases, with
  a floor of 15 phase slots; all trial sets pooled, levels overall and for the exact n with upper bounds;
  `s110_summary.py` checks its inputs (registered variants, localisation computed from the primary) and refuses
  incomplete sets; file-name checks; repeat-once rule for interrupted runs; deadlines for variants and localisation.
- **2026-10-02, adversarial audit 5**: the Gaussian tail p = exp(-d^2/2) of the localisation label was not calibrated
  (12 of 300 on-target archival injections below 0.0027); the camera/CCD rule of the localisation ran before the origin
  checks, so an archival file from another sector and CCD wrote an official "not possible"; `s110_summary.py` did not
  check the sector of a localisation file and accepted "not available" with any text. Fixed: empirical p from 3000
  further on-target injections with unchanged thresholds, coverage checked; every origin check before the CCD rule;
  sector and star checked in the summary; "not available" needs the saved error output.
- **2026-10-02, audit 6**: a saved "QLP not yet public" record made the robustness qualifier final, although the
  rule keeps it provisional until 2027-02-28; the summary accepted a variant carrying another script, ephemeris,
  expected-depth or helper-module hash; four localisation dry-run files carried the hash of an earlier text of
  `s110_localise.py`; the false-alarm trials used one seed where the Sector 110 run uses the median of 20 (one archival
  phase changed class); the fixed sentence called a binomial bound of dependent trials "at most". Fixed: "not public"
  keeps the qualifier provisional, and a final status needs a run (d), a failed attempt on a public release, or a MAST
  search on or after 2027-02-28; every result must carry this package's hashes; the four localisation files re-run with
  the final script (all numbers identical); every trial set re-run with the Sector 110 thresholds (14 of 17 835
  classifications changed; SPOC 2-min A 0.47 %, B 0.49 %); the bound is called approximate in the sentences and in
  section 4.
- **2026-10-03, final independent check**: wording and documentation only; no code, number, threshold or rule of
  computation changed, and every output was reproduced exactly from a clean copy with freshly downloaded MAST files.
  Fixed: the predicted-time column headed "UTC" gave the barycentric time as a calendar date, 5.6 to 7.2 min later
  than the arrival time; the table now gives both. A fallback primary needs a saved MAST search of 2027-02-28 or
  later, and "not available" needs the unedited error output of the registered command (a refusal does not count);
  readers compare posted hashes with `SHA256SUMS`. The O4 row of section 4 now says "d <= det" as the code does.
  Added `calib/transit_times_A.csv` and `calib/ephemeris_inputs.json` (inputs of the ephemeris), `LICENSE`,
  acknowledgements and a tested reproduction recipe in `README.md`.
- **2026-10-03, packaging for Zenodo**: no code, number, threshold or rule changed. Added a plain-language opening
  section, a guide to the sections and a list of terms to `PREREG_S110_FINAL.md`, the public reference of the DDT
  approval (DDT 100), and a precision to point 7 of section 9 (an adaptation under point 4 has its own posted hash);
  rewrote `README.md`; added `PREREG_S110_FINAL.pdf` (rendering of the .md), `FILE_GUIDE.md`, `CITATION.cff`,
  `requirements_lock.txt`, `calib/ephemeris_sources/` (sources of the ephemeris: timing injections and the posterior
  samples of A and B) and three reproduction helpers in `build/` (`fetch_archival_inputs.sh`, `reproduce.sh`,
  `compare_outputs.py`), tested from a fresh copy with freshly downloaded MAST files.
- **2026-10-05, rules the code does not enforce**: no code, number, threshold, prediction or dry-run output changed.
  An adversarial check of the deposit (2026-10-03) showed that the summary accepts any non-empty file as the error
  output of a variant declared not available, that failures leaving no JSON had no path to a fixed sentence, that a
  reading adaptation did not say under which file name it runs, that "robust" also covers a sensitive variant
  returning O2, and that the code checks neither the choice of a fallback product nor the time system of an input.
  Added section 9c (invalid results; product search record; registered time-system header values; attempt records
  and "not available"; fixed sentences for failures without a JSON; adapted reader keeps its file name; exact meaning
  of "robust" and a fixed line of variant outcomes), cross-references in sections 0, 2, 3, 9 and 9b, row 9h of
  section 1, and matching sentences in `README.md`. The code does not enforce these rules; the public record does.
- **2026-10-05, after two independent checks of the packaged deposit**: no code, number, threshold, prediction or
  dry-run output changed. Section 9c.1 now states that the time limits of section 9 points 1 and 3 (30, 7 and 14
  days) are administrative (no outcome depends on them; the only date any script compares is 2027-02-28), while the
  rules of order, 2027-02-28 and the 30-day limit for an adapted reader stay binding. Wording made consistent: the
  file given for a variant not available is its attempt record (sections 3, 9, 9c.4); a reading failure needs one
  record; a (d) file failing the time check; the helper modules are never adapted; listing and hashing files is not
  an opening; attempt records written with `script` in UTC, variants on one UTC day; the variant line prints
  P(detect | real) where computed; (d) is attempted with the other variants only if a QLP release is public that
  day; a broken rule of order cannot be repaired by repeating a run; the localisation part of a hand-written failure
  sentence is written like section 9b. `CITATION.cff`: reference type `data` (CFF 1.2.0). The deposit no longer includes
  the unpublished drafts v1 and v2; they are not needed to reproduce anything here.
- **2026-10-05, neutral paths and comments; last wording fixes**: no number, threshold, prediction or rule of
  computation changed. Comments and default file paths that named the author's computer or private working folders
  were made neutral in `deps/pixloc_lib.py`, `deps/lcio.py`, `s110_final.py` (folder of the archival inputs for
  `--calibrate`, now `../archival_inputs/mast`), `s110_localise.py`, `build/false_alarm_check.py` (finds its own
  folder), `build/run_all.sh`, `build/fetch_archival_inputs.sh`, `build/reproduce.sh`, `build/make_ephemeris.py`, and
  in the notes of `calib/ephemeris_inputs.json` and `calib/localisation_archival.json`. Because the code hashes
  changed, `calib/expected_depths.json`, `calib/localisation_coverage.json`, every `dryrun/` output and the build logs
  were rebuilt with `build/run_all.sh`; every number, trial list and fixed sentence is identical to the previous
  outputs (`build/compare_outputs.py`). The false-alarm trial files were not re-run (hours of CPU; the script changed
  only in its path line and comments, and the files carry no code hash). Section 9c: one continuous session for the
  variants, with the date of a late (d); the two kinds of 30-day limit (administrative run limit, binding reading
  limits) named and separated; the remedy when a date has passed; the form of the variant line; where (d)'s search
  record is posted; the qualifier in the hand-written failure cases. `README.md`: `SHA256SUMS.txt` lists the four
  other uploaded files.
- **Final, 2026-10-02, wording clarified 2026-10-03 and 2026-10-05** (this deposit, `PREREG_S110_FINAL.md`): the only public
  pre-registration. All changes from v2
  and their reasons are in its section 1. No Sector 110 data existed at any point.
- **2026-10-06, version 1.1: text rewritten for outside readers** (this deposit: `PREREG_S110_FINAL.pdf`, the
  authoritative text, and `PREREG_S110_SHORT.pdf`, its Sections 1 to 7 alone; LaTeX source in `protocol_source/`).
  New title. No threshold, window, ephemeris, seed, order of runs, fallback, outcome rule or code changed; every file
  other than the text files is byte-identical to the deposit of 2026-10-05. That deposit (the earlier Markdown text,
  `PREREG_S110_FINAL.md`, described by the previous entry as "the only public pre-registration") was published on
  Zenodo on 2026-10-05 as record 23171456 and deleted by the author the same day; version 1.1 replaces it and does not
  include it. Precedence: version 1.1 is authoritative; within it, the formal rules of the appendices govern over the
  summaries of Sections 1 to 7; the frozen code still defines how an outcome is computed from valid inputs.
  Corrections without change of rule: the table of measured false-detection levels is completed from
  `calib/false_alarm_levels.json` (SPOC 2-min A n = 9: 0 of 24; QLP A per n; the TESS-SPOC 30-min set of 618 trials;
  every sum checked, no value changed); the coded POS_CORR condition of the localisation is stated next to the earlier
  wording; two cases where the code words a sentence differently from the tables are documented (Appendix H); factual
  wording corrected in the limits ("never more lenient" replaced by "makes detection stricter and widens the
  consistency band slightly"; the meaning column of O1 and O3; "tighter" archival tolerances); the section numbers cited
  by the code are mapped in Appendix N. Added: how the signals were found and why their selection was post hoc
  (Section 2), the archival localisation of A, ephemeris uncertainties of A from the frozen covariance, abstract,
  references, acknowledgments, the placeholders of the fixed sentences, the statement that the internal reviews were
  made by separate AI sessions, and commands split with shell line continuation.

- **2026-10-06, completion of independent checks and the short version**: Sections 4 and 7 now explicitly state
  the existing counting conditions (80%, two baseline cadences per side, segment gaps >0.3 d), QLP flux-column
  order, fallback priority, 20 bootstrap seeds and medians, and the main order, failure and deadline contingencies
  already registered in the appendices. The O4 QLP cross-reference now points to Appendix I.
  Sections 5 and E.4 clarify an inherited text/code discrepancy: the frozen empirical-p denominator is
  N_cal + 1 for the successful calibration trials (2700 to 3000), not always 3001. This follows the existing
  code-precedence rule and changes no code, threshold or computed outcome. Appendix H clarifies that the
  NOT_POSSIBLE label in the non-finite-distance case is assigned to A, whereas B receives no label.
  README acknowledges registered hand-written failure sentences. Codex (OpenAI) provided editorial and
  verification assistance. No frozen analysis file, calibration, dry-run output or scientific rule was changed.
  The empirical-p resolution in Appendix K now uses the same N_cal + 1 denominator; the historical
  description in Appendix L distinguishes stricter detection from the wider consistency band.
  The author confirmed that internal/independent package checks were separate AI review sessions, not human
  peer review; this says nothing about outside astronomers' data or comments on the target.
  Cross-document LaTeX labels are prefixed to avoid duplicates; Table 3 spacing and DOI labels were adjusted
  to prevent overflow. These presentation edits do not change numerical values.
  Sections 4.5 and D.2 clarify existing summary-code behavior: a qualifying contradiction produces
  not robust immediately, even while variant (d) is pending; robust and robustness not testable remain
  provisional. Table H's robust condition explicitly requires at least one qualifying variant.
  In the short PDF, references to appendices/tables in the full protocol print their correct numbers
  without a dangling local hyperlink; internal references in each PDF remain clickable.

- **2026-10-06, four-page reader's summary**: replaced the initial eight-page extract of Sections 1–7 with
  a standalone four-page summary. The full v1.1 protocol and its appendices remain authoritative.
  Predictions, definitions, decision boundaries, timing, localisation, robustness and the main procedural
  commitments are retained in the summary; execution details remain in the full protocol. Removed repeated
  editorial chronology from the main text and summary; retained scientifically relevant post-hoc changes
  in Appendix L and detailed chronology here. Condensed AI disclosure, placed at the end of the summary;
  separate AI review sessions are not human peer review. No frozen analysis code, calibration, prediction
  or archival output changed.

- **2026-10-06, reader-oriented overview**: rewrote the Summary as a two-page introduction to the
  experiment. Replaced calculation formulas and execution details with plain-English descriptions.
  Preserved the predictions, primary/fallback principle, fixed-time detection, empirical noise and
  calibration caveat, sensitivity gate, O1–O4, timing, localisation, narrow robustness rule and limits.
  Prediction values displayed in the Summary are explicitly rounded; operational predictions remain
  unchanged in the full protocol and JSON files. The full protocol's two file descriptions now call
  the Summary a reader-oriented overview without specifying its length. No scientific rule, threshold,
  prediction, frozen code, calibration or archival output changed.

- **2026-10-06, entry-point simplification**: applied the author's additional reader-facing brief.
  The prediction table now gives period, primary expected box depth and Sector110 event count only.
  Data and test is three paragraphs; timing, localisation and robustness are two or three sentences each.
  Removed localisation numerical boundaries and the variant list/qualifier mechanics from the Summary;
  all their exact definitions remain in the binding full protocol and frozen code. Preserved every
  O1–O4 gate, sensitivity distinction, archival selection caveat, calibration caveat and interpretation
  limits. The detailed full PDF and all its scientific sources are unchanged in this edit.
