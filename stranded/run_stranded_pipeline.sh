#!/bin/bash
# GrASE pipeline on the strand-reconstructed DEXSeq counts.
#
# The existing counts were made with dexseq_count --stranded reverse on an
# UNSTRANDED library, which discarded ~half of every bin at random and let
# antisense neighbours' reads into a gene's bins. DEXSeq/count_files_stranded/
# fixes both (see scripts/run_stranded_recount.sh). Everything here writes to
# NEW *.stranded.* directories so the current benchmark stays reproducible.
#
# phi MUST be re-estimated: exontest.R skips estimation when a phi file already
# exists in --outdir, and at 2x depth the dispersions move materially. Fresh
# outdirs, fresh --phi basenames, nothing pre-seeded.
#
# Only the MERGED tests are run. merge_exon_sj_counts.R substitutes a junction
# count only on a side whose exonic distinct set is empty and leaves every other
# side byte-identical, so the merged test subsumes the exon-only test. The
# exonic counts are still built -- they are the merge input and the reference
# side is never substituted.
#
# Usage: bash scripts/run_stranded_pipeline.sh
set -eux
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
R=$HOME/GrASE/Rpkg/scripts
CNT=$BASE/DEXSeq/count_files_stranded

echo "=== 0. split-read counts on bipartition intronic edges ==="
# Idempotent: bipartition_sjcnt.R skips any gene whose output already exists, so
# a rerun costs only the SJ matrix build. These counts come from STAR's
# SJ.out.tab and are NOT affected by the strand reconstruction -- STAR applied
# no strand filter (its only strand option here is --outSAMstrandField
# intronMotif, which adds an XS tag and drops nothing), so junctions were always
# at full depth. Only the exonic side was halved by dexseq_count --stranded
# reverse, which is what step 1 fixes.
for T in internal TSSTTS; do
  bash $R/run_bipartition_sjcnt.sh $T
done

echo "=== 1. exonic part counts ==="
for A in internal TSS; do
  TAG=$([ "$A" = internal ] && echo internal || echo TSSTTS)
  OUT=$BASE/bipartition.$TAG.stranded.counts
  [ -s "$OUT/bipartition.$TAG.exoncnt.combined.txt" ] || \
  Rscript $R/exoncnt.R -c "$CNT" -t bipartition --cond1=group1 --cond2=group2 \
      -a $A -i $BASE/bipartition.filtered -o "$OUT"
done

echo "=== 2. merge in the exclusive junction where the exonic side is empty ==="
Rscript $R/merge_exon_sj_counts.R --exon_counts bipartition.internal.stranded.counts \
        --sj_counts sjcnt        --output bipartition.merged.stranded.counts
Rscript $R/merge_exon_sj_counts.R --exon_counts bipartition.TSSTTS.stranded.counts \
        --sj_counts sjcnt.TSSTTS --output bipartition.merged.TSSTTS.stranded.counts

echo "=== 3. merged tests (phi estimated fresh from these counts) ==="
run_test () {   # countdir combined outdir phibase annotdir
  local CD=$1 CF=$2 OD=$3 PB=$4 AD=$5
  mkdir -p "$OD"
  if [ ! -s "$CD/$CF" ]; then
    head -1 "$(ls $CD/*.exoncnt.txt | head -1)" > "$CD/$CF"
    for f in $CD/*.exoncnt.txt; do tail -n +2 "$f" >> "$CD/$CF"; done
  fi
  for f in $AD/*.bipartition.txt; do ln -sf "$(realpath "$f")" "$CD/$(basename "$f")"; done
  Rscript $R/exontest.R --file "$CF" --countdir "$CD" --outdir "$OD" \
      --splittype bipartition --model betabinom_EBapprox --phi "$PB" \
      --cond1 group1 --cond2 group2 --use_phi_loess --independent_filtering
}
## MERGED ONLY. The merged counts are a strict superset of the exonic ones --
## merge_exon_sj_counts.R substitutes a junction count ONLY on a side whose
## exonic distinct set is empty (is.na(setdiff) & !is.na(intron_distinct)) and
## leaves every other side byte-identical. So the merged test subsumes the
## exon-only test; running both doubles the compute for no extra result.
## The exonic counts from step 1 are still required -- they are the merge input,
## and the reference side is never substituted.
run_test bipartition.merged.stranded.counts   bipartition.merged.exoncnt.combined.txt \
         bipartition.merged.stranded.test.EBapprox   phi.merged.stranded.txt \
         bipartition.internal.counts
run_test bipartition.merged.TSSTTS.stranded.counts bipartition.merged.exoncnt.combined.txt \
         bipartition.merged.TSSTTS.stranded.test.EBapprox phi.merged.TSSTTS.stranded.txt \
         bipartition.TSSTTS.counts

echo "=== Done ==="
