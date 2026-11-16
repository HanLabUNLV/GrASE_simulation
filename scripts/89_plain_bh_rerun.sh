#!/bin/bash
# Rerun the PR/ROC sweep and the delta_pi (effect-size) sweep with GrASE padj as
# plain BH (PADJ_METHOD default), copying the current outputs to *.$STAMP first.
#
# Usage: bash scripts/89_plain_bh_rerun.sh [backup suffix, default bak_nestedBH]
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
EV=results/eval_metricfair
STAMP=${1:-bak_nestedBH}
export STRANDED=1 GT_RULE_STRANDED=1 SJAWARE=1
for f in $EV/pr_three_levels.GT_rule.txt $EV/effect_sweep.GT_rule.txt $EV/effect_sweep.GT_rule.table.txt \
         plots/pr_three_levels.GT_rule*.pdf plots/roc_partial.GT_rule*.pdf plots/effect_sweep.GT_rule.*.pdf; do
  [ -e "$f" ] && [ ! -e "$f.$STAMP" ] && cp "$f" "$f.$STAMP"
done
[ -f $EV/pr_calls_cache.GT_rule.rds ] && mv $EV/pr_calls_cache.GT_rule.rds $EV/pr_calls_cache.GT_rule.rds.$STAMP
echo "=== PR sweep ==="
Rscript scripts/pr_curves_three_levels_gtrule.R
echo "=== ROC ==="
Rscript scripts/plot_roc_partial.R
echo "=== delta_pi sweep ==="
Rscript scripts/pr_effect_sweep.R
echo "=== Done ==="
