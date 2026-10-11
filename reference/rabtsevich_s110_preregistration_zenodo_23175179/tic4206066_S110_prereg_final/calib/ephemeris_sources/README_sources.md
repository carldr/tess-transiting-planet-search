# Sources of the ephemeris (added 2026-10-03)

These files are the deeper inputs behind `ephemeris_final.json`, copied unchanged from the author's archival analyses so
that every number of the ephemeris can be traced inside this deposit. They do not change anything; the ephemeris was
built from them on 2026-10-01 and is frozen.

| file | what it is | origin | used for |
|---|---|---|---|
| `../transit_times_A.csv` | the 22 archival mid-times of signal A (BTJD) with formal and injection-calibrated errors | final re-derivation, item 03 (`transit_times_calibrated.csv`, 2026-10-01) | weighted least-squares fit of T0, P and their covariance |
| `timing_injections.csv` | per-transit summary of the timing injections (200 per transit) | final re-derivation, item 03 (2026-10-01) | `sig68_min` = the injection-calibrated error `tc_err_inj_min` |
| `timing_injections_deltas.npz` | recovered minus injected mid-time (min) of every injection, keys `S<sector>_<epoch>_delta` | final re-derivation, item 03 (2026-10-01) | half the 16-84 percentile width of each array is `tc_err_inj_min` |
| `mcmc_samples.npz` | 120 000 posterior samples of the transit model of A (names in the array `names`: P, T0, k, b, aRs, ...) | final re-derivation, item 03 (2026-10-01) | T14 of A = median of P/pi * arcsin(sqrt((1+k)^2 - b^2) / (aRs sin i)), cos i = b/aRs |
| `chain_B_pipe_mann_frozen_w0.75.npz` | 13 152 posterior samples of the joint transit fit of B (names in `names`) | archival joint fit of 2026-09-25, used unchanged by the drafts v1 and v2 | P and T0 of B (medians); sigma_pred of B = 1.5 x the standard deviation of T0 + n P over the samples |

Checked on 2026-10-03: the errors in `../transit_times_A.csv` equal the injection widths to 4e-5 min; the weighted fit
gives T0 and P of `ephemeris_final.json` to 1e-9 d; the T14 formula above gives 1.870548542043 h (frozen
1.870548542043308 h); the B samples give P and T0 exactly and sigma_pred within 0.0005 min (the frozen values were
rounded to 0.001 min before the factor 1.5).

Not part of this deposit: the code that measured the transit times and ran the two posterior fits on the TESS light
curves (the author's archival analysis). Those steps are inputs here, not reproduced by `build/reproduce.sh`.
