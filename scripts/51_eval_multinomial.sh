#!/bin/bash
# Evaluate the STRANDED multinomial EBplugin results (run 09-19 04:29),
# replacing the 2026-07-13 unstranded eval that the cross panel still reads.
#
# Ground truth stays the part-level framework (results/sim_exon_info, gtI/II/III),
# NOT GT_rule: no GT_rule table exists for multinomial, and its long-format
# counts have no distinct/reference split for the rule to apply to. The cross
# panel therefore scores all three comparison structures the same way.
#
# NOTE --min_reads was not applied to the multinomial tests (exontest.R skips
# the read floor when diff1/diff2 are absent), so this arm is unfiltered while
# bipartition and n_choose_2 are floored at 10. State that in the caption.
#
# Usage: bash scripts/51_eval_multinomial.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
EVAL=scripts/evaluate_multinomial_test.R
GT=results/sim_exon_info
SIM=swimdown/simulate/data/simulate.rda
T=multinomial.stranded.test.EBplugin
OUT=results/eval_multinomial_EBplugin

f1="$T/test_multinomial.internal_dirmult_EBplugin.annotated.txt"
f2="$T/test_multinomial.TSSTTS_dirmult_EBplugin.annotated.txt"
for f in "$f1" "$f2"; do
  [ -f "$f" ] || { echo "MISSING $f"; exit 1; }
done

if [ -d "$OUT" ] && [ ! -d "${OUT}.pre_stranded" ]; then
  mv -f "$OUT" "${OUT}.pre_stranded"
  echo "backed up $OUT -> ${OUT}.pre_stranded"
fi

echo "=== multinomial eval $(date +%H:%M:%S) ==="
Rscript "$EVAL" "$f1,$f2" "$GT" "$OUT" "$SIM"
echo "=== done $(date +%H:%M:%S) -> $(ls "$OUT" 2>/dev/null | wc -l) files ==="
