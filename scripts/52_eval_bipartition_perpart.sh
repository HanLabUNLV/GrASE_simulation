#!/bin/bash
# Evaluate the per-part, exon-only bipartition results (02_stranded_pipeline.sh),
# the bipartition arm of the comparison of structures, with the evaluator, ground
# truth and simulate.rda that 50_eval_nc2.sh and 51_eval_multinomial.sh use, so
# the three structures are scored the same way (part-level, results/sim_exon_info).
# visualize_eval.py (73_cross_confusion.sh) reads the output.
#
# Usage: bash scripts/52_eval_bipartition_perpart.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
EVAL=scripts/evaluate_bipartition_test.R
GT=results/sim_exon_info
SIM=swimdown/simulate/data/simulate.rda
f1=bipartition.internal.stranded.test.EBapprox.perpart/test_bipartition.ex.internal_betabinom_EBapprox.annotated.txt
f2=bipartition.TSSTTS.stranded.test.EBapprox.perpart/test_bipartition.ex.terminal_betabinom_EBapprox.annotated.txt
OUT=results/eval_bipartition_EBapprox
for f in "$f1" "$f2"; do
  [ -f "$f" ] || { echo "MISSING $f -- run 02_stranded_pipeline.sh" >&2; exit 1; }
done
[ ! -e "$OUT" ] || { echo "$OUT exists; move it aside first" >&2; exit 1; }
echo "=== bipartition eval $(date +%H:%M:%S) ==="
Rscript "$EVAL" "$f1,$f2" "$GT" "$OUT" "$SIM"
echo "=== done $(date +%H:%M:%S) -> $(ls "$OUT" 2>/dev/null | wc -l) files ==="
