#!/bin/bash
# Precision-recall and partial ROC laid out as DTE and DTU ROWS, one file per
# universe. The sweep's own panels split the categories across separate files
# with the universes in rows, which makes the DTE-against-DTU contrast -- the
# comparison the results text actually makes -- impossible to see without
# flipping between images. These transpose it.
#
# Both scripts READ results/eval_metricfair/pr_three_levels.GT_rule.txt and
# recompute nothing, so run 30_merge_downstream.sh stranded first if that table
# is stale. Curves start at padj 0.001, matching the other plotting scripts; see
# the note in README.md.
#
# Outputs: plots/pr_dte_dtu.{full,restricted}.{pdf,png}
#          plots/pr_all.{full,restricted}.{pdf,png}      (written by the same script)
#          plots/roc_dte_dtu.{full,restricted}.{pdf,png}
#
# Usage: bash scripts/86_plots_dte_dtu.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation

TAB=results/eval_metricfair/pr_three_levels.GT_rule.txt
[ -f "$TAB" ] || { echo "missing $TAB -- run 30_merge_downstream.sh stranded" >&2; exit 1; }
echo "=== input: $TAB  ($(date -r "$TAB" '+%F %H:%M')) ==="

echo "=== precision-recall ==="
Rscript scripts/pr_dte_dtu_by_universe.R

echo "=== partial ROC ==="
Rscript scripts/roc_dte_dtu_by_universe.R

ls -l --time-style=long-iso plots/pr_dte_dtu.*.pdf plots/pr_all.*.pdf plots/roc_dte_dtu.*.pdf
