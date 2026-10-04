#!/bin/bash
# Model comparison (EBapprox / EBmap / MLE / wilcoxon) on the STRANDED MERGED
# internal arm, so plot_model_comparison.R matches the lane every other figure
# and table uses. It previously read bipartition.test/, which is the UNSTRANDED
# PRE-MERGE arm -- the p-values were valid but not comparable with the rest.
#
# All four run into one outdir so the plot reads a single directory. The
# existing EBapprox result and its phi tables are copied in rather than refit;
# phi is estimated once from these counts and reused, which is what makes the
# four models comparable (a different phi per model would confound the
# comparison with the dispersion estimate).
#
# --min_reads is left at its default 10, matching the rest of the manuscript.
# That is apples-to-apples across models because d_support does not depend on
# the model, so the same sides are dropped from all four. Pass --min_reads 0 if
# you want the low-support regime, where the models differ most.
#
# Usage: bash scripts/60_modelcomp_internal.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
CNT=bipartition.merged.stranded.counts
SRC=bipartition.merged.stranded.test.EBapprox
OUT=bipartition.merged.stranded.test.modelcomp
mkdir -p $OUT

# phi estimated from these counts, plus the already-computed EBapprox result.
# The EBapprox run tags its phi output with the model name, so the files are
# copied under the untagged names this script passes to --phi. Both EBapprox and
# EBmap then read the SAME phi -- a per-model phi would confound the model
# comparison with the dispersion estimate.
cp -n "$SRC/phi.merged.stranded.EBapprox.txt" \
      "$OUT/phi.merged.stranded.txt"
cp -n "$SRC/phi.merged.stranded.EBapprox.approx.moderated.txt" \
      "$OUT/phi.merged.stranded.approx.moderated.txt"
# EBapprox itself is not refit; eval_modelcomp.sh reads the annotated form.
for f in test_bipartition.merged_betabinom_EBapprox.txt \
         test_bipartition.merged_betabinom_EBapprox.annotated.txt; do
  cp -n "$SRC/$f" "$OUT/$f"
done

run () {   # run <model> [extra args]
  local m=$1; shift
  echo "=== $m ==="
  Rscript $S/exontest.R \
    --file bipartition.merged.stranded.internal.exoncnt.combined.txt \
    --countdir $CNT --outdir $OUT \
    --splittype bipartition --model "$m" \
    --cond1 group1 --cond2 group2 "$@" \
    &> modelcomp.$m.log
  echo "  -> $(ls -la $OUT/test_bipartition.merged_${m}.txt 2>/dev/null | awk '{print $5}') bytes"
}

run betabinom_EBmap --phi phi.merged.stranded.txt --use_phi_loess --independent_filtering
run betabinom_MLE
run wilcoxon
echo "=== done ==="
