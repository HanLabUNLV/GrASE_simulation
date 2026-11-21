#!/bin/bash
# MAJIQ on the strand-split BAMs without de novo junctions, then its unit ground
# truth: majiq/run_majiq_stranded.sh (build --disable-denovo, deltapsi, voila tsv),
# make_lsv_junction_keys.R and infer_lsv_junction_gt.R on the stranded build.
# The evaluation (30_merge_downstream.sh steps 3-5) reads these afterwards.
#
# Usage: bash scripts/42_majiq_nodenovo.sh
set -euo pipefail
cd /mnt/data1/home/mirahan/GrASE_simulation
export STRANDED=1
echo "=== MAJIQ  $(date '+%F %T') ==="
bash majiq/run_majiq_stranded.sh
echo "=== LSV junction keys  $(date '+%F %T') ==="
Rscript scripts/make_lsv_junction_keys.R
echo "=== LSV ground truth  $(date '+%F %T') ==="
Rscript scripts/infer_lsv_junction_gt.R
echo "=== done  $(date '+%F %T') ==="
