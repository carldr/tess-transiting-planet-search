#!/bin/zsh
# One-command rebuild check (added 2026-10-03). It never writes into this deposit.
#   1. zsh build/fetch_archival_inputs.sh          (public TESS files and PRF models, checked by sha256)
#   2. zsh build/reproduce.sh /path/to/python [WORKDIR]
# The python must come from an environment built from requirements.txt (or requirements_lock.txt).
# A copy of the deposit is made in WORKDIR (default: ../s110_rebuild next to the deposit); build/run_all.sh is run
# there with only its PY= and PRF= lines set; then build/compare_outputs.py compares every rebuilt output with this
# deposit (numbers, trial lists and fixed sentences; run timestamps and output hashes are expected to differ).
set -e
DEP="${0:A:h:h}"
PY="${1:?usage: zsh build/reproduce.sh /path/to/python [WORKDIR]}"
PY="${PY:a}"   # absolute path; symlinks of a venv must not be resolved
WORK="${2:-$DEP/../s110_rebuild}"
IN="${DEP:h}/archival_inputs"
PRF="${DEP:h}/tess_prf/start_s0004"
if [ ! -d "$IN/mast" ] || [ ! -d "$PRF/cam1_ccd2" ]; then
  echo "archival inputs missing: run  zsh build/fetch_archival_inputs.sh  first"; exit 1
fi
mkdir -p "$WORK"; WORK="${WORK:A}"
C="$WORK/tic4206066_S110_prereg_final"
rm -rf "$C" "$WORK/archival_inputs"
cp -Rp "$DEP" "$C"
ln -s "$IN" "$WORK/archival_inputs"
sed -e "s|^PY=.*|PY=\"$PY\"|" -e "s|^PRF=.*|PRF=\"$PRF\"|" "$C/build/run_all.sh" > "$C/build/run_all_local.sh"
echo "rebuilding in $C (25 to 40 min on an Apple Silicon Mac)"
zsh "$C/build/run_all_local.sh"
"$PY" "$DEP/build/compare_outputs.py" "$DEP" "$C"
