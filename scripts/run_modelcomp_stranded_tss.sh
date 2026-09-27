#!/bin/bash
# The TSS/TTS half of the stranded-merged model comparison. The internal half is
# scripts/run_modelcomp_stranded.sh; evaluate_bipartition_test.R takes the two
# arms comma-separated, so both are needed before the PR-across-models figure
# can be regenerated.
#
# Same principle as the internal half: ONE phi, estimated from these counts
# during the stranded run, shared by EBapprox and EBmap. MLE and wilcoxon do
# not use phi. A per-model phi would confound the model comparison with the
# dispersion estimate.
#
# EBapprox already exists on this arm and is copied in rather than refit.
# Order is wilcoxon -> MLE -> EBmap, cheapest first, so a failure surfaces early
# rather than after the multi-hour run.
#
# Usage: bash scripts/run_modelcomp_stranded_tss.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
CNT=bipartition.merged.TSSTTS.stranded.counts
SRC=bipartition.merged.TSSTTS.stranded.test.EBapprox
OUT=bipartition.merged.TSSTTS.stranded.test.modelcomp
mkdir -p $OUT

for f in phi.merged.TSSTTS.stranded.txt phi.merged.TSSTTS.stranded.approx.moderated.txt \
         test_bipartition.merged_betabinom_EBapprox.txt; do
  [ -e "$OUT/$f" ] || cp "$SRC/$f" "$OUT/$f"
done

run () {
  local m=$1; shift
  echo "=== $m  $(date +%H:%M:%S) ==="
  Rscript $S/exontest.R \
    --file bipartition.merged.exoncnt.combined.txt \
    --countdir $CNT --outdir $OUT \
    --splittype bipartition --model "$m" \
    --cond1 group1 --cond2 group2 "$@" \
    &> modelcomp.tss.$m.log
  echo "  -> $(stat -c %s $OUT/test_bipartition.merged_${m}.txt 2>/dev/null || echo FAILED) bytes"
}

run wilcoxon
run betabinom_MLE
run betabinom_EBmap --phi phi.merged.TSSTTS.stranded.txt --use_phi_loess --independent_filtering
echo "=== done ==="
