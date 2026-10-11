#!/bin/zsh
# Rebuilds calib/expected_depths.json, calib/localisation_coverage.json and every dryrun/ output with the frozen scripts.
set -e
cd "${0:A:h}/.."
PY=python   # a Python from requirements.txt; build/reproduce.sh sets this line
M=../archival_inputs/mast
S32=$M/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits
S98=$M/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits
TPF=$M/tess2020324010417-s0032-0000000004206066-0200-s_tp.fits
TPF6=$M/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits
PRF=../tess_prf/start_s0004
rm -rf dryrun; mkdir dryrun
$PY -u s110_final.py --calibrate --out calib/expected_depths.json > build/calibrate.log 2>&1
run(){ n=$1; shift; $PY -u s110_final.py --dry-run "$@" --out dryrun/$n.json > dryrun/$n.log 2>&1; }
run s32_real --lc $S32 --product spoc & run s32_null --lc $S32 --product spoc --t0-shift-frac 0.5 --cut-true &
run s98_real --lc $S98 --product qlp & run s98_v2rule --lc $S98 --product qlp --count-rule v2 & wait
run s32_sap --lc $S32 --product spoc --flux SAP_FLUX & run s32_w05 --lc $S32 --product spoc --window 0.5 &
run s32_w10 --lc $S32 --product spoc --window 1.0 & run s32_desat1 --lc $S32 --product spoc --drop-desat-hours 1 &
run s32_v2rule --lc $S32 --product spoc --count-rule v2 &
OMP_NUM_THREADS=2 $PY -u s110_localise.py --dry-run --tpf $TPF --prf-dir $PRF --ephem items10 --events own --out dryrun/loc_s32_items10_fidelity.json > dryrun/loc_s32_items10_fidelity.log 2>&1 &
wait
OMP_NUM_THREADS=2 $PY -u s110_localise.py --dry-run --tpf $TPF --prf-dir $PRF --photometry dryrun/s32_real.json --out dryrun/loc_s32_final_pipeline.json > dryrun/loc_s32_final_pipeline.log 2>&1 &
OMP_NUM_THREADS=2 $PY -u s110_localise.py --dry-run --tpf $TPF --prf-dir $PRF --photometry dryrun/s32_real.json --signal B --out dryrun/loc_s32_B.json > dryrun/loc_s32_B.log 2>&1 &
$PY -u s110_localise.py --dry-run --tpf $TPF --prf-dir $PRF --photometry dryrun/s32_null.json --out dryrun/loc_s32_null_gate.json > dryrun/loc_s32_null_gate.log 2>&1 &
wait
# coverage of the empirical localisation p (audit 5): the same dry run with another calibration seed
OMP_NUM_THREADS=2 $PY -u build/localisation_coverage.py --tpf $TPF --prf-dir $PRF --photometry dryrun/s32_real.json \
    --reference dryrun/loc_s32_final_pipeline.json --out-loc dryrun/loc_s32_coverage_seed5005.json \
    --out calib/localisation_coverage.json > build/localisation_coverage.log 2>&1
# registered (non-dry) runs on archival files must be refused
$PY -u s110_localise.py --tpf $TPF --prf-dir $PRF --photometry dryrun/s32_real.json --out dryrun/loc_s32_refused.json > dryrun/loc_s32_refused.log 2>&1 || true
# audit 5: a file from another sector on another CCD (S6, camera 1 / CCD 1) is refused before the CCD rule can write
$PY -u s110_localise.py --tpf $TPF6 --prf-dir $PRF --n-trials 1 --no-clear --out dryrun/loc_s6_ccd1_refused.json > dryrun/loc_s6_ccd1_refused.log 2>&1 || true
$PY -u s110_final.py --lc $S32 --product spoc --out dryrun/phot_s32_refused.json > dryrun/phot_s32_refused.log 2>&1 || true
# fixed sentences from the dry-run outputs (code check of s110_summary.py; needs calib/false_alarm_levels.json)
# --not-available needs a saved error output (audit 5); in these code checks the variants were simply not run
NAF=dryrun/na_not_run.txt
echo "dry run: variant not run (code check of s110_summary.py only; in Sector 110 this file is the saved error output of the attempted run)" > $NAF
NA="--not-available a=$NAF --not-available b=$NAF --not-available c=$NAF --not-available e=$NAF --not-available f=$NAF"
$PY -u s110_summary.py --dry-run --primary dryrun/s32_real.json --variant a=dryrun/s32_w05.json --variant b=dryrun/s32_w10.json \
    --variant c=dryrun/s32_desat1.json --variant e=dryrun/s32_sap.json --variant f=dryrun/s32_v2rule.json \
    --loc-A dryrun/loc_s32_final_pipeline.json --loc-B dryrun/loc_s32_B.json --out dryrun/summary_s32.json > dryrun/summary_s32.log 2>&1
$PY -u s110_summary.py --dry-run --primary dryrun/s32_null.json ${=NA} --out dryrun/summary_s32_null.json > dryrun/summary_s32_null.log 2>&1
$PY -u s110_summary.py --dry-run --primary dryrun/s98_real.json --not-available a=$NAF --not-available b=$NAF \
    --not-available e=$NAF --not-available f=$NAF --out dryrun/summary_s98.json > dryrun/summary_s98.log 2>&1
# the summary must refuse when a registered variant is silently left out
$PY -u s110_summary.py --dry-run --primary dryrun/s32_real.json --variant a=dryrun/s32_w05.json \
    --loc-A dryrun/loc_s32_final_pipeline.json --out dryrun/summary_refused.json > dryrun/summary_refused.log 2>&1 || true
# audit 5: a variant declared not available without a saved error output is refused
$PY -u s110_summary.py --dry-run --primary dryrun/s32_null.json --not-available a=dryrun/no_such_error_file.txt \
    --out dryrun/summary_na_refused.json > dryrun/summary_na_refused.log 2>&1 || true
# audit 6: a variant that carries another script hash is refused; "(d) not public" cannot make the status final before 2027-02-28
$PY -c "import json; j = json.load(open('dryrun/s32_sap.json')); j['provenance']['script_sha256'] = '0' * 64; json.dump(j, open('dryrun/s32_sap_other_version.json', 'w'), indent=1)"
$PY -u s110_summary.py --dry-run --primary dryrun/s32_real.json --variant a=dryrun/s32_w05.json --variant b=dryrun/s32_w10.json \
    --variant c=dryrun/s32_desat1.json --variant e=dryrun/s32_sap_other_version.json --variant f=dryrun/s32_v2rule.json \
    --loc-A dryrun/loc_s32_final_pipeline.json --loc-B dryrun/loc_s32_B.json \
    --out dryrun/summary_version_refused.json > dryrun/summary_version_refused.log 2>&1 || true
$PY -u s110_summary.py --dry-run --primary dryrun/s32_real.json --variant a=dryrun/s32_w05.json --variant b=dryrun/s32_w10.json \
    --variant c=dryrun/s32_desat1.json --variant e=dryrun/s32_sap.json --variant f=dryrun/s32_v2rule.json \
    --not-available d=$NAF --d-status final --d-final-reason not-public-on-deadline \
    --loc-A dryrun/loc_s32_final_pipeline.json --loc-B dryrun/loc_s32_B.json \
    --out dryrun/summary_d_final_refused.json > dryrun/summary_d_final_refused.log 2>&1 || true
echo ALL DONE
