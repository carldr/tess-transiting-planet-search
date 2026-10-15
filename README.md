# TESS: search for new transiting planet candidates

This project searches TESS light curves for transiting planets that nobody has reported yet.
It follows the method of Pavel Rabtsevich, who found two candidate signals around the star TIC 4206066.
His star is the known-positive test: the search must find his 3.18 d and 11.13 d signals in Sector 98 before we trust it on other stars.

## Layout

| path | content |
|---|---|
| `analysis/` | The search code. |
| `data/tic4206066_mast/` | TESS light curves and target pixel files of TIC 4206066 for Sectors 6, 32 and 98. |
| `data/tic4206066_tesscut/` | 21x21-pixel TESScut cutouts of the full-frame images of TIC 4206066 for Sectors 6, 32 and 98. |
| `data/qlp_bulk_download_scripts/` | The MAST bulk download scripts of QLP light curves, one per sector. The search code reads the TIC IDs and URLs from them and never runs them. |
| `data/qlp_sector_<nnnn>/` | The star list and the downloaded light curves of one sector. |
| `data/known_transit_catalogues/` | The TOI, CTOI, SPOC TCE and TESS EB lists, to tell new signals from known ones. |
| `data/candidate_other_sector_light_curves/` | The light curves of the new candidates from other TESS sectors. |
| `data/candidate_tesscut_cutouts/` | TESScut cutouts of the full-frame images around the new candidates, one folder per TIC ID. |
| `results/` | Outputs of the search code, one folder per script. The folders that end in `_quality_zero_mask` hold the first Sector 98 run, which also removed the cadences with QUALITY bit 30. |
| `results/candidate_verification_records/` | One file per candidate, with a row for each verification step that ran on it. |
| `results/candidate_decisions.csv` | The decisions that no step makes, such as a candidate held for a later sector. |
| `results/pipeline_status/` | The live progress of the running steps, for the candidate board. Not in the repository. |
| `results/ctoi_drafts/` | One file per CTOI draft, with the values to submit and the checks made before submission. |
| `reference/rabtsevich_preprint_zenodo_22967456/` | Rabtsevich's research note, `manuscript.pdf` (doi:10.5281/zenodo.22967456). |
| `reference/rabtsevich_data_zenodo_22999133/` | His supporting data and code, version 1.0.1 (doi:10.5281/zenodo.22999133). `unpacked/` holds the contents of the three ZIP files. |
| `reference/rabtsevich_s110_preregistration_zenodo_23175179/` | His pre-registered test of TESS Sector 110 (doi:10.5281/zenodo.23175179). `tic4206066_S110_prereg_final/` holds the contents of the ZIP file. |
| `.venv/` | The Python 3.12.13 environment. `requirements.txt` lists its packages. |
| `.venv_triceratops/` | A separate Python 3.12.13 environment for TRICERATOPS, which needs numpy 1 and an older scipy. `requirements_triceratops.txt` lists its packages. |

The files under `reference/` and `data/` are downloads. Do not edit them.
The repository does not hold `data/` or the Python environments.
Step 10 moves the Gaia background star list that TRICERATOPS writes to `data/triceratops_gaia_background/`.

## Set up

1. Create the Python environment:

   ```
   uv venv --python 3.12.13 .venv
   uv pip install --python .venv/bin/python -r requirements.txt
   ```

   Create the TRICERATOPS environment:

   ```
   uv venv --python 3.12.13 .venv_triceratops
   uv pip install --python .venv_triceratops/bin/python -r requirements_triceratops.txt
   ```

2. Download the archival TESS files of TIC 4206066:

   ```
   .venv/bin/python analysis/download_tic4206066_archival_tess_files.py
   ```

   The script reads the file list from Rabtsevich's `01_inputs/INPUT_PRODUCTS.csv`, so the test star uses the same files he did.
   It downloads about 320 MB and checks each file against the sha256 in that list.
   The three TESScut cutouts always fail this check, because MAST makes each cutout on request and the file headers change.
   Every other file must match.

## Run

Run each script from the repository root.

### Search one sector

Run these steps in order. The examples use Sector 98, the sector that Rabtsevich searched.

Steps 2–4 and 7–10 continue a stopped run by default. A stop, such as a reboot or a power cut, loses at most the item in progress:

- Steps 2 and 3 skip the stars and chunks that are already done.
- Steps 4 and 7–10 write each result to disk as soon as they have it. They keep the earlier results without an error, and do the other items.
  They drop the earlier results of items that are no longer in their input, such as a candidate that now matches a TOI.
  They do the items with an error again, so they retry a temporary failure, such as an archive outage.

Run the same command again to continue.
Steps 1, 5 and 6 write their output only at the end, and always start from the beginning.

Add `--restart` to each later step after you change the code or the output of an earlier step.
Otherwise the later steps keep their results from the old input.
Add `--tic-id <TIC ID>` to steps 7 to 10 to do only the candidates of that star. Give it more than once for more stars.
With `--restart`, the step then does only those candidates again, and keeps the results of the other candidates.

Commit the output of each step when the step is done. The repository does not hold `data/`, so a disk failure loses all uncommitted results.

1. List the stars with a QLP light curve in the sector, with their TIC magnitude, radius, mass and temperature:

   ```
   .venv/bin/python analysis/build_qlp_sector_star_list.py 98
   ```

   The script writes `data/qlp_sector_0098/star_list.csv`. Sector 98 has 628,730 stars; the TIC query takes about 15 min.

2. Download the light curves of the stars at or brighter than TESS magnitude 11.5:

   ```
   .venv/bin/python analysis/download_qlp_sector_light_curves.py 98 --max-tmag 11.5
   ```

   The script keeps only time, SAP flux and quality, in `data/qlp_sector_0098/light_curves/`.
   For Sector 98 it downloads about 185 GB and keeps about 30 GB.
   It skips stars that are already saved, so run it again to continue or to retry failed downloads.

3. Search every downloaded light curve with box least squares (BLS):

   ```
   .venv/bin/python analysis/search_qlp_sector.py 98
   ```

   The script removes flagged cadences, but keeps the cadences flagged only with QLP QUALITY bit 30 (2^29).
   QLP does not document bit 30, and those cadences are nearly as clean as unflagged ones.
   The script detrends each light curve and searches periods of 0.5 to 14 d.
   It masks each signal it finds and searches again, up to three signals per star.
   It writes one row per signal to `results/qlp_sector_0098_bls_search/`, in chunks of 1000 stars.
   It skips chunks that are already done, so run it again to continue.

4. Vet the signals with LEO-Vetter:

   ```
   .venv/bin/python analysis/vet_qlp_sector_signals.py 98
   ```

   The script vets each signal with a BLS SNR of at least 7 and at least 3 transits.
   It runs the LEO-Vetter flux-level tests and checks them against the LEO-Vetter default thresholds.
   The false alarm tests reject noise and systematics. The false positive tests reject eclipsing binaries.
   A signal that passes both groups is a planet candidate.
   The script writes `results/qlp_sector_0098_leo_vetting/vetting.csv`, with the names of the tests that each signal failed.
   It does not run the LEO-Vetter pixel-level tests.
   It uses all but two CPU cores. Add `--workers 9` to use 9 processes.

5. Download the catalogues of known transit signals:

   ```
   .venv/bin/python analysis/download_known_transit_catalogues.py
   ```

   The script downloads the TOI and CTOI lists from ExoFOP, every SPOC TCE table from MAST, and the TESS eclipsing binary (EB) catalogue from VizieR, about 300 MB.
   It writes them to `data/known_transit_catalogues/`, with the download date in `download_date.txt`.
   Run it again before step 6 when the lists are more than a few days old.

6. Compare the vetted signals with the known signals, and count each stage:

   ```
   .venv/bin/python analysis/cross_match_qlp_sector_signals_with_known_catalogues.py 98
   ```

   A signal is known when a TOI, CTOI, TCE or TESS EB is on the same star with a period of 1/3, 1/2, 1, 2 or 3 times the signal period, within 1%.
   The script writes three files to `results/qlp_sector_0098_known_catalogue_cross_match/`:

   | file | content |
   |---|---|
   | `new_candidates.csv` | The signals that pass every vetting test and match no known signal, highest MES first. |
   | `signals_cross_matched.csv` | Every vetted signal, with the matching TOIs, CTOIs, TCEs and TESS EBs. |
   | `summary.json` | The number of signals at each stage, the TOI recovery, and the rank of TIC 4206066. |

   The TOI recovery measures the search: of the TOIs on searched stars, how many the search found and the vetting passed.
   `summary.json` gives the counts for all TOIs and for each TFOPWG disposition, such as confirmed planet (CP) or false positive (FP).

7. Check each new candidate against independent data:

   ```
   .venv/bin/python analysis/check_qlp_sector_new_candidates.py 98
   ```

   The script runs these checks on each candidate in `new_candidates.csv`:

   | check | data | refutes the candidate when |
   |---|---|---|
   | Other TESS sectors | The star's SPOC, TESS-SPOC or QLP light curve of each other sector, searched near the candidate period | The transit should be clear in a sector, but is not there. |
   | Period pile-up | The periods of all the signals of the sector search | Many more stars have a signal at this period than at nearby periods. The period is a systematic of the sector. |
   | Ephemeris match | The signals of the sector search, the TOIs and the TESS EBs within 10 arcmin | A nearby star has a dip at the same period that removes more light. That star is the likely source. |
   | Gaia DR3 target | RUWE, non-single-star solution, eclipsing binary classification | Gaia classifies the star as an eclipsing binary. |
   | Gaia DR3 neighbours | The sources within 60 arcsec | A neighbour that could cause the dip is a Gaia eclipsing binary. |
   | VSX | The known variable stars within 60 arcsec | One of them is an eclipsing or ellipsoidal binary. |
   | Companion | The LEO-Vetter transit fit and the TIC stellar radius | The companion is larger than 2 Jupiter radii, so it is a star. |

   The script gives each candidate a verdict.
   "refuted" means that a check refutes it.
   "supported" means that the joint search detects the transit and no check refutes it.
   The joint search searches all sectors together, near the candidate period.
   A real transit stays in phase over the years between sectors, so its SNR grows.
   The script also searches the same light curves at ten wrong periods.
   The joint search detects the transit when its SNR is 1.5 times the highest SNR at a wrong period.
   "unconfirmed" means neither.
   The column `cautions` lists weaker warnings, such as a Gaia RUWE above 1.4, a grazing transit, or a neighbour bright enough to cause the dip.

   The script saves the other-sector light curves in `data/candidate_other_sector_light_curves/`.
   It writes `results/qlp_sector_0098_new_candidate_checks/candidate_checks.csv`, one row per candidate.

8. Find which star causes the dip of each candidate that step 7 did not refute:

   ```
   .venv/bin/python analysis/check_qlp_sector_candidate_pixels.py 98
   ```

   TESS pixels are 21 arcsec wide, so a dip on a nearby star can look like a transit on the target.
   The script checks every sector in which the star falls on the science pixels and TESScut has full-frame images.
   In each sector it makes a difference image from a TESScut cutout: the image out of transit minus the image in transit.
   The pixels that lose light in the transit are bright in it.
   The script fits the difference image with the TESS pixel response function at the target and at each Gaia star bright enough to cause the dip.
   The star that fits best is the likely source.
   The fit cannot separate stars less than half a pixel apart, so "target" below includes the stars within half a pixel of it.

   | `pixel_verdict` | meaning |
   |---|---|
   | on target | The target fits clearly better than every other star. |
   | on neighbour | A nearby star fits clearly better than the target. The candidate is a false positive. |
   | unresolved | The target and a nearby star fit about equally well. |
   | too faint | The difference image is too noisy to locate the dip. |

   The verdict of the candidate is "on neighbour" when any sector gives "on neighbour".
   Otherwise it is the first of "on target", "unresolved" and "too faint" that a sector gives.

   The script saves the cutouts in `data/candidate_tesscut_cutouts/`, about 200 MB per star and sector.
   It writes two files and one image per candidate and sector to `results/qlp_sector_0098_candidate_pixel_checks/`:

   | file | content |
   |---|---|
   | `pixel_checks.csv` | One row per candidate, with the combined verdict. |
   | `pixel_checks_by_sector.csv` | One row per candidate and sector. |

9. Verify each candidate that step 8 puts on the target or leaves unresolved, with the light curves of all TESS sectors:

   ```
   .venv/bin/python analysis/verify_qlp_sector_candidates_over_all_sectors.py 98
   ```

   The script fits one transit model to all sectors, and gives the values that a CTOI needs: period, epoch, depth, duration and planet radius, each with its error.
   It then runs these tests:

   | test | refutes the candidate when |
   |---|---|
   | Odd and even transits | Their depths differ by 3 sigma or more. The dip is an eclipsing binary at twice the period. |
   | Secondary eclipse at phase 0.5 | It is deeper than a planet can give by 3 sigma or more. |
   | Planet radius | It is larger than 2 Jupiter radii. |
   | Transit fit | The MCMC chain does not converge in 100,000 steps. The data do not constrain the transit shape. |
   | Stellar density from the fit | Never. A density that does not fit the TIC star is a caution. |
   | Depth in each sector, with the same apertures in the TESScut cutout of each sector | Never. Different depths are a caution. |
   | Gaia DR3 radial velocity | Never. A variable radial velocity means a stellar companion, and is a caution. For a bright star, the script gives the largest companion mass that the Gaia epochs allow. |
   | LAMOST DR7 spectra within 15 arcsec | Never. The script lists the class and the radial velocity of each spectrum, and the star in its fibre. A fibre with no Gaia star in it recorded the light of a near star, often the target. |
   | Ellipsoidal variation at half the period | Never. A variation of 3 sigma or more is a caution, and the script gives the companion mass that it needs. A planet gives a few ppm at most, but a low-mass star gives hundreds of ppm. |

   The verdict is "refuted" or "passed".
   Some known planets also show an ellipsoidal variation of 3 sigma, so this test does not refute a candidate.
   To measure how often the test gives a false detection, run the test on known planets of hot stars:

   ```
   .venv/bin/python analysis/calibrate_ellipsoidal_variation_test_on_known_planets.py
   ```

   The script prints the spread of the significance for the planets that cannot cause a detectable variation, and writes `results/ellipsoidal_variation_test_calibration/`.
   The script writes `results/qlp_sector_0098_all_sector_verification/all_sector_verification.csv`, and one image per candidate with the folded transit of each sector.

10. Compute the false positive probability of each candidate that passed step 9:

    ```
    .venv_triceratops/bin/python analysis/validate_qlp_sector_candidates_with_triceratops.py 98
    ```

    The script runs TRICERATOPS in its own environment, on the light curves of all sectors, with the period, epoch and depth from step 9.
    TRICERATOPS compares the folded transit with each false positive scenario, such as an eclipsing binary on the target, on a companion or on a nearby star.
    It gives the false positive probability (FPP) and the nearby false positive probability (NFPP).
    TRICERATOPS runs 5 times for each candidate, and one candidate can take many hours.
    The script prints a line with the time at the start of the first run and at the end of each run.
    For an "on target" candidate, the script removes the scenarios on the stars that step 8 rules out.

    | `triceratops_verdict` | condition |
    |---|---|
    | validated without imaging | FPP < 0.015 and NFPP < 0.001 |
    | likely planet | FPP < 0.5 and NFPP < 0.001 |
    | not validated | any other result |

    A full validation also needs high-resolution imaging of the star. This project does not have it.
    The script writes `results/qlp_sector_0098_triceratops_validation/triceratops_validation.csv`.

### Verification sequence and record

A candidate goes to the next step only when it meets the condition of that step:

| step | script | condition to go to the next step |
|---|---|---|
| 6 | `cross_match_qlp_sector_signals_with_known_catalogues.py` | It passes every LEO-Vetter test and matches no TOI, CTOI, TCE or TESS EB. |
| 7 | `check_qlp_sector_new_candidates.py` | The verdict is "supported" or "unconfirmed". |
| 8 | `check_qlp_sector_candidate_pixels.py` | The pixel verdict is "on target" or "unresolved". |
| 9 | `verify_qlp_sector_candidates_over_all_sectors.py` | The verdict is "passed". |
| 10 | `validate_qlp_sector_candidates_with_triceratops.py` | The verdict is "likely planet" or "validated without imaging". |

Before a CTOI submission:

1. Download the catalogues again (step 5), and run step 6 again.
2. Make sure that the candidate is still in `new_candidates.csv`. If it is not, a TOI, CTOI, TCE or TESS EB now matches it.
3. Find the code version of the candidate's rows for steps 7 to 10 in its verification record.
4. Run `git diff --stat <code version> HEAD -- analysis/` to list the changes to the code since then.
5. If a change can change a result of a step, run that step again with `--restart --tic-id <TIC ID>`.
6. Run a later step again only if it uses a value that the rerun changed. Step 10 uses the sectors, period, epoch, depth and duration of step 9.

A CTOI is published under the name of the person who submits it.

Each of steps 6 to 10 appends one row per candidate to the verification record of the candidate, `results/candidate_verification_records/tic<TIC ID>_period_<period>_days.md`.
The row gives the date, the sector, the step, the script, the git commit of the code, the verdict and the result.
The commit is the one checked out when the step starts, so a commit made during a long run does not change it.
The commit ends in `-modified` when `analysis/` has uncommitted changes at the start.
Each step writes its result to its CSV first, with the same date and code version, and then appends the record row.
So each record row has a result row in `results/` behind it.
A result with an error gets no record row, because the next run does that candidate again.
A run with `--restart` adds a new row, so the record keeps the full history of the candidate.
Commit the records with the results.

### Check the candidates for follow-up

Run this command after step 10, and again at intervals, such as once a month:

```
.venv/bin/python analysis/check_candidates_on_exofop_and_future_tess_sectors.py
```

The script checks each TIC ID with a step 10 verdict of "likely planet" or "validated without imaging", in all sectors.
It downloads the ExoFOP target page of the star and prints each TOI, CTOI, follow-up observation or file that the last run did not have.
Observers do not publish their follow-up schedules, so ExoFOP shows an observation only after the observer uploads it.
The script also lists the TESS sectors that observe the star and that step 9 did not use.
The column `released_tess_sectors_not_in_all_sector_verification` gives those of them with a light curve on MAST.
Steps 7 to 10 find the other sectors of a candidate on MAST each time they check it.
So when that column has a sector, run steps 7 to 10 again with `--restart --tic-id <TIC ID>` on the sector of the candidate.
The script only reads ExoFOP.

The script writes `results/candidate_exofop_and_future_tess_sector_checks/candidate_exofop_and_future_tess_sector_checks.csv`, one row per TIC ID.
It keeps the last ExoFOP page of each star in `exofop_target_pages/` in the same folder, to compare with the next run.
Commit the folder after each run.

### Watch the candidates on the board

```
.venv/bin/python analysis/serve_candidate_board.py
```

Open http://localhost:8000, or http://<the address of this computer>:8000 from another device on the same network.
The server listens on every network of the computer. Give `--host 127.0.0.1` to limit it to this computer.

The board has one card per new candidate of step 6, in all searched sectors.
A card is in the column of the step that checks it next, or of the step that stopped it.
The page reads the results again every 10 seconds, so a card moves as soon as a step writes its result.
Steps 7 to 10 write their progress to `results/pipeline_status/`, and the card of the candidate that a step works on shows that progress.
Click a card to see its verification record and its figures.

`results/candidate_decisions.csv` holds the decisions that no step makes, one row per candidate.
The `decision` column is `held`, `submitted` or `dropped`, and the card goes to the column of that decision.
Write the reason in the `note` column, and commit the file after each change.

### Test star

| command | what it does |
|---|---|
| `.venv/bin/python analysis/search_tic4206066_tls_per_sector.py` | Searches each archival sector of TIC 4206066 with Transit Least Squares (TLS), without Rabtsevich's periods. Writes `results/tls_per_sector/`. |

`analysis/tic4206066_light_curves.py` loads and detrends the light curves for the scripts.
It removes every cadence with a non-zero QUALITY flag, and removes trends with a wotan biweight filter with a 0.5 d window.

## Check the Rabtsevich downloads

Run each command from the folder that holds the checksum file. Each command prints `OK` for every file.

```
cd reference/rabtsevich_data_zenodo_22999133/unpacked && shasum -a 256 -c SHA256SUMS
cd reference/rabtsevich_s110_preregistration_zenodo_23175179 && shasum -a 256 -c SHA256SUMS.txt
```
