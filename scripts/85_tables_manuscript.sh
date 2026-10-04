#!/usr/bin/env bash
## Manuscript tables 2, 3, 4 and 5, regenerated from the evaluation outputs.
##
## This was previously reachable from no driver at all -- benchmark_tables.R had
## to be run by hand, so "regenerate the manuscript tables" had no documented
## path and it was possible for the tables to be built from a different pipeline
## state than the figures.
##
## Inputs, all written by 30_merge_downstream.sh:
##   results/eval_metricfair/pr_three_levels.GT_rule.txt
##   results/eval_metricfair/transcript_level_metrics.txt
##   results/eval_metricfair/pr_calls_cache.GT_rule.rds
##
## Outputs: markdown to stdout, plus four .txt tables in results/eval_metricfair.
##
## ALL REPORTED RESULTS ARE STRANDED. Run 30_merge_downstream.sh stranded first.
##
## Usage: bash scripts/85_tables_manuscript.sh [outfile]
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
OUT=${1:-results/eval_metricfair/manuscript_tables.md}
EV=results/eval_metricfair
for f in pr_three_levels.GT_rule.txt transcript_level_metrics.txt pr_calls_cache.GT_rule.rds; do
  [ -s "$EV/$f" ] || { echo "missing $EV/$f -- run 30_merge_downstream.sh stranded first" >&2; exit 1; }
done
echo "=== inputs ==="
for f in pr_three_levels.GT_rule.txt transcript_level_metrics.txt pr_calls_cache.GT_rule.rds; do
  echo "  $f  $(date -r "$EV/$f" +'%Y-%m-%d %H:%M')"
done
Rscript scripts/benchmark_tables.R | tee "$OUT"
echo "=== wrote $OUT ==="
