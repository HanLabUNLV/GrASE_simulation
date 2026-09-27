#!/bin/bash
# Re-run the four STRANDED simulation exontest arms under the new pre-testing
# read-support filter (--min_reads 10).
#
# The filter sets p.value <- NA for any side whose tested distinct set has
# fewer than min_reads in BOTH contrasted groups, BEFORE p-value adjustment,
# so those sides never enter the FDR set. The row is kept with padj = NA,
# matching independent filtering's convention.
#
# The LRT checkpoints (ckpt_lrt_*.rds) and the phi tables already exist in each
# outdir, so nothing is refit -- this re-runs adjustment + annotation only.
# Do NOT delete the ckpt files or this becomes a multi-day job.
#
# Usage: bash scripts/20_tests_bipartition_minreads.sh
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
R=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
cd $BASE

# grase resolves via ~/.Renviron to the shared library at
# /mnt/data1/home/jaquino/R/... -- re-install there after any edit to
# Rpkg/R/*.R, or the previous build's signatures are used silently (that is how
# an add_significant() without min_dpi_sj/min_reads survived several runs):
#   R CMD INSTALL --library=/mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2 \
#                 /mnt/data1/home/mirahan/GrASE/Rpkg

STAMP=bak_premin_reads

run_arm () {
  local file=$1 cntdir=$2 outdir=$3 phi=$4 log=$5
  echo "=== $outdir ==="
  # snapshot the pre-filter outputs once, so before/after stays comparable
  for f in $outdir/test_*.txt; do
    [ -e "$f" ] || continue
    [ -e "$f.$STAMP" ] || cp "$f" "$f.$STAMP"
  done
  Rscript $R/exontest.R \
    --file $file \
    --countdir $cntdir \
    --outdir $outdir \
    --splittype bipartition \
    --model betabinom_EBapprox \
    --phi $phi \
    --cond1 group1 --cond2 group2 \
    --use_phi_loess --independent_filtering \
    --min_reads 10 &> $log
  echo "--- support filter lines:"
  grep -c "Support filter" $log || true
  grep "Support filter" $log | head -4 || true
}

run_arm bipartition.internal.exoncnt.combined.txt \
        bipartition.internal.stranded.counts \
        bipartition.internal.stranded.test.EBapprox \
        phi.internal.stranded.txt \
        minreads.internal.stranded.log

run_arm bipartition.TSSTTS.exoncnt.combined.txt \
        bipartition.TSSTTS.stranded.counts \
        bipartition.TSSTTS.stranded.test.EBapprox \
        phi.TSSTTS.stranded.txt \
        minreads.TSSTTS.stranded.log

run_arm bipartition.merged.exoncnt.combined.txt \
        bipartition.merged.stranded.counts \
        bipartition.merged.stranded.test.EBapprox \
        phi.merged.stranded.txt \
        minreads.merged.internal.stranded.log

run_arm bipartition.merged.exoncnt.combined.txt \
        bipartition.merged.TSSTTS.stranded.counts \
        bipartition.merged.TSSTTS.stranded.test.EBapprox \
        phi.merged.TSSTTS.stranded.txt \
        minreads.merged.TSSTTS.stranded.log

echo "=== Done ==="
