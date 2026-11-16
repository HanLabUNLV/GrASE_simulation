#!/bin/bash
# Rerun the benchmark with the GrASE_merged_default curve (exontest.R's default
# call rule on plain BH), copying current outputs to *.$STAMP first.
#
# Usage: bash scripts/91_default_gate_rerun.sh [backup suffix, default bak_nodefault]
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
EV=results/eval_metricfair
STAMP=${1:-bak_nodefault}
export STRANDED=1 GT_RULE_STRANDED=1 SJAWARE=1
for f in $EV/pr_three_levels.GT_rule.txt $EV/transcript_level_metrics.txt $EV/tsstts_attribution.txt \
         $EV/delta_sweep.GT_rule.txt $EV/delta_sweep.GT_rule.table.txt \
         $EV/native_unit_table.txt $EV/structural_reach_table.txt $EV/transcript_level_table.txt \
         $EV/gene_level_table_fmt.txt \
         plots/pr_three_levels.GT_rule*.pdf plots/roc_partial.GT_rule*.pdf plots/delta_sweep.GT_rule.*.pdf \
         plots/confusion_three_levels.GT_rule*.pdf; do
  [ -e "$f" ] && [ ! -e "$f.$STAMP" ] && cp -p "$f" "$f.$STAMP"
done
[ -f $EV/pr_calls_cache.GT_rule.rds ] && mv $EV/pr_calls_cache.GT_rule.rds $EV/pr_calls_cache.GT_rule.rds.$STAMP
echo "=== PR sweep ===";        Rscript scripts/pr_curves_three_levels_gtrule.R
echo "=== ROC ===";             Rscript scripts/plot_roc_partial.R
echo "=== transcript ===";      Rscript scripts/transcript_level_metrics.R
echo "=== lfc delta sweep ==="; Rscript scripts/pr_delta_sweep.R
echo "=== tables ===";          Rscript scripts/benchmark_tables.R
echo "=== Done ==="
