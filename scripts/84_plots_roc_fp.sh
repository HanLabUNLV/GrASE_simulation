#!/bin/bash
# The two manuscript figure producers that had no driver: the partial ROC and
# the FP-location panel.
#
# scripts/visualize_eval.py names four current producers for the benchmark
# figures and tables. Two were already driven -- pr_curves_three_levels_gtrule.R
# by 30_merge_downstream.sh, benchmark_tables.R by 85_tables_manuscript.sh -- and
# these two were not, so their outputs silently fell behind the lane by a week.
#
# Both read results/eval_metricfair/pr_calls_cache.GT_rule.rds, so run
# 30_merge_downstream.sh stranded first.
#
# Usage: bash scripts/84_plots_roc_fp.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation

CACHE=results/eval_metricfair/pr_calls_cache.GT_rule.rds
[ -f "$CACHE" ] || { echo "missing $CACHE -- run 30_merge_downstream.sh stranded" >&2; exit 1; }
echo "=== input: $CACHE  ($(date -r "$CACHE" '+%F %H:%M')) ==="

echo "=== partial ROC ==="
Rscript scripts/plot_roc_partial.R

echo "=== FP location ==="
Rscript scripts/plot_fp_universe_shift.R

echo "=== done ==="
