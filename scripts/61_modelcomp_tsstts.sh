#!/bin/bash
# The TSS/TTS half of the stranded-merged model comparison. The internal half is
# scripts/60_modelcomp_internal.sh; evaluate_bipartition_test.R takes the two
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
# Usage: [MC_CORES=24] bash scripts/61_modelcomp_tsstts.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
CNT=bipartition.merged.TSSTTS.stranded.counts
SRC=bipartition.merged.TSSTTS.stranded.test.EBapprox
OUT=bipartition.merged.TSSTTS.stranded.test.modelcomp
mkdir -p $OUT

# Same as the internal half: copy the EBapprox run's phi under the name this
# script passes to --phi, and bring the annotated EBapprox result along.
cp -n "$SRC/phi.merged.TSSTTS.stranded.txt" \
      "$OUT/phi.merged.TSSTTS.stranded.txt"
cp -n "$SRC/phi.merged.TSSTTS.stranded.approx.moderated.txt" \
      "$OUT/phi.merged.TSSTTS.stranded.approx.moderated.txt"
for f in test_bipartition.merged_betabinom_EBapprox.txt \
         test_bipartition.merged_betabinom_EBapprox.annotated.txt; do
  cp -n "$SRC/$f" "$OUT/$f"
done

# per-gene DEXSeq GFFs and mate length: exontest.R writes per-base pi, which the
# default call rule needs (--dpi_gate both)
GFF=/mnt/data1/home/mirahan/GrASE_simulation/dexseq.gff
run () {
  local m=$1; shift
  echo "=== $m  $(date +%H:%M:%S) ==="
  Rscript $S/exontest.R \
    --prefix bipartition.exsj.terminal \
    --countdir $CNT --outdir $OUT \
    --splittype bipartition --model "$m" \
    --cond1 group1 --cond2 group2 --gff_dir $GFF --read_length 100 --sj_overhang 1 \
    --mc_cores "${MC_CORES:-24}" "$@" \
    &> modelcomp.tss.$m.log
  for f in $OUT/test_bipartition.exsj.terminal_$m*; do
    ln -sf "$(basename "$f")" "$OUT/$(basename "$f" | sed "s/^test_bipartition\.exsj\.terminal_/test_bipartition.merged_/")"
  done
  echo "  -> $(stat -L -c %s $OUT/test_bipartition.merged_${m}.txt 2>/dev/null || echo FAILED) bytes"
}

run wilcoxon
run betabinom_MLE
run betabinom_EBmap --phi phi.merged.TSSTTS.stranded.txt --use_phi_loess --padj_threshold=0.01
echo "=== done ==="
