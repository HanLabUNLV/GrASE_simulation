#!/bin/bash
# Post-hoc lfc_diff_net filter evaluation on the STRANDED MERGED EBapprox
# results -- the lane every other figure uses. run_evaluate_posthoc_lfc_filter.sh
# reads bipartition.test/, the unstranded pre-merge exonic arm.
#
# Note the merged arm names its outputs test_bipartition.merged_*, not
# test_bipartition.{internal,TSSTTS}_*.
#
# Usage: bash scripts/run_posthoc_lfc_stranded.sh [delta]
set -eu
SIM_DIR=/mnt/data1/home/mirahan/GrASE_simulation
SCRIPT_DIR=/mnt/data1/home/mirahan/GrASE_simulation/scripts
MODEL=betabinom_EBapprox
DELTA=${1:-0}

TEST_FILES="${SIM_DIR}/bipartition.merged.stranded.test.EBapprox/test_bipartition.merged_${MODEL}.annotated.txt,${SIM_DIR}/bipartition.merged.TSSTTS.stranded.test.EBapprox/test_bipartition.merged_${MODEL}.annotated.txt"
GT_DIR="${SIM_DIR}/results/sim_exon_info"
SIM_RDA="${SIM_DIR}/swimdown/simulate/data/simulate.rda"
OUT_DIR="${SIM_DIR}/posthoc_lfc_filter.stranded_merged.${MODEL}.delta${DELTA}"

Rscript "${SCRIPT_DIR}/evaluate_posthoc_lfc_filter.R" \
  "${TEST_FILES}" "${GT_DIR}" "${OUT_DIR}" "${SIM_RDA}" "${DELTA}"
echo "=== eval done -> ${OUT_DIR} ==="

Rscript "${SCRIPT_DIR}/plot_delta_roc.R" \
  "${OUT_DIR}/roc_data.full.txt" \
  "${SIM_DIR}/scripts/plots/delta_roc.stranded_merged.pdf"
echo "=== plot done ==="
