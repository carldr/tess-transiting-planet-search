#!/bin/zsh
# Downloads the public archival TESS files and PRF models used by the rebuild, puts them where the frozen scripts expect
# them, and checks their sha256. Run from anywhere: zsh build/fetch_archival_inputs.sh
# Layout created next to the deposit folder:
#   ../archival_inputs/mast/                     (path fixed inside s110_final.py --calibrate and build/run_all.sh)
#   ../tess_prf/start_s0004/cam1_ccd2/           (passed as --prf-dir ../tess_prf/start_s0004)
set -e
DEP="${0:A:h:h}"
M="$DEP/../archival_inputs/mast"
PRF="$DEP/../tess_prf/start_s0004/cam1_ccd2"
mkdir -p "$M" "$PRF"
TID=https://archive.stsci.edu/missions/tess/tid/s0032/0000/0000/0420/6066
SPOC6=https://archive.stsci.edu/hlsps/tess-spoc/s0006/target/0000/0000/0420/6066
QLP98=https://archive.stsci.edu/hlsps/qlp/s0098/0000/0000/0420/6066
cd "$M"
for u in $SPOC6/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits \
         $SPOC6/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits \
         $TID/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits \
         $TID/tess2020324010417-s0032-0000000004206066-0200-s_tp.fits \
         $QLP98/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits; do
  [ -f "${u:t}" ] || curl -sf -O "$u"
done
shasum -a 256 -c - <<'EOF'
ac7736ea403e33360c53f88b5f22bb4e41c8235fe17ff2cdb934ba474037882e  hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits
7a03d3baf07e30ac2a57160368077dc6bcc9ea6862085b336d7d7f1f8ef7003d  hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits
764dad34b3bbd467eef973a1d09f640a97321824e6739ac6a0c69cf40c31b0e7  tess2020324010417-s0032-0000000004206066-0200-s_lc.fits
f1f76576cc088d134cd94b8b004cbd7af27be161a9a80578bbe96b7335f1b2d1  tess2020324010417-s0032-0000000004206066-0200-s_tp.fits
a034136697b4f75765a91503107f55a455b365d19cd270ed633845efcdcb188e  hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits
EOF
cd "$PRF"
grep -v '^#' "$DEP/PRF_SHA256SUMS" | while read h f; do
  [ -f "${f:t}" ] || curl -sf -O "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0004/$f"
done
grep -v '^#' "$DEP/PRF_SHA256SUMS" | sed 's|cam1_ccd2/||' | shasum -a 256 -c -
echo "archival inputs ready: ${M:A} and ${PRF:A:h}"
