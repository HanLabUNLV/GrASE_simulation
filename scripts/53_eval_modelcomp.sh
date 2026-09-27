#!/bin/bash
# Per-model evaluations on the STRANDED MERGED results, in the sim_type-stratified
# framework (Background / DGE / DTE / DTU). The GT_rule table labels only DTE and
# DTU, so the Background and DGE false-positive counts quoted in the results text
# can only come from here.
#
# Usage: bash scripts/53_eval_modelcomp.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
EV=scripts/evaluate_bipartition_test.R
GT=results/sim_exon_info
RDA=swimdown/simulate/data/simulate.rda
INT=bipartition.merged.stranded.test.modelcomp
TSS=bipartition.merged.TSSTTS.stranded.test.modelcomp

for m in betabinom_EBapprox betabinom_EBmap betabinom_MLE wilcoxon; do
  out=results/eval_modelcomp_stranded_${m}
  echo "=== $m -> $out  $(date +%H:%M:%S) ==="
  Rscript $EV \
    "$INT/test_bipartition.merged_${m}.annotated.txt,$TSS/test_bipartition.merged_${m}.annotated.txt" \
    "$GT" "$out" "$RDA" &> eval_modelcomp.${m}.log
  echo "  -> $(ls $out 2>/dev/null | wc -l) files"
done
echo "=== done ==="
