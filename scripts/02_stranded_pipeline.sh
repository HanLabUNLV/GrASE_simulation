#!/bin/bash
# GrASE pipeline on the strand-reconstructed DEXSeq counts: filter, split-read
# counts, exonic counts, merge, merged tests.
#
# The existing counts were made with dexseq_count --stranded reverse on an
# UNSTRANDED library, which discarded ~half of every bin at random and let
# antisense neighbours' reads into a gene's bins. DEXSeq/count_files_stranded/
# fixes both (see scripts/01_stranded_recount.sh).
#
# Each output directory is moved to <dir>.olddedup the first time this runs
# (the bipartitions were regenerated with the corrected dedup), and a rerun
# never touches that backup. Every stage therefore starts from an empty
# directory: bipartition_sjcnt.R skips genes whose output exists and
# merge_exon_sj_counts.R does not clear its output, so stale per-gene files
# would otherwise survive.
#
# phi MUST be re-estimated: exontest.R skips estimation when a phi file already
# exists in --outdir. Fresh outdirs, nothing pre-seeded.
#
# Only the MERGED tests are run. merge_exon_sj_counts.R substitutes a junction
# count only on a side whose exonic distinct set is empty and leaves every other
# side byte-identical, so the merged test subsumes the exon-only test. The
# exonic counts are still built -- they are the merge input and the reference
# side is never substituted.
#
# Directory names are the ones the evaluation scripts read; the files inside
# carry the package's current names (<gene>.bipartition.{ex,exsj}.<class>.cnt.txt,
# test_bipartition.exsj.<class>_*). The evaluation scripts read the merged
# results as test_bipartition.merged_*, so each result is also linked under
# that name.
#
# Usage: bash scripts/02_stranded_pipeline.sh
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
R=$HOME/GrASE/Rpkg/scripts
CNT=$BASE/DEXSeq/count_files_stranded
BIPART=bipartition.v28.dedupfix

backup () { for d in "$@"; do [ -e $d.olddedup ] || { [ ! -e $d ] || mv $d $d.olddedup; }; done; }
cdir () { [ $1 = internal ] && echo $2 || echo $3; }   # class internal-name terminal-name

echo "=== 0. split bipartitions into internal and terminal ==="
backup bipartition.filtered
Rscript $R/filterTerminal.R --split_dir=$BIPART --outdir=bipartition.filtered --split_type=bipartition

echo "=== 1. split-read counts on bipartition intronic edges ==="
# These counts come from STAR's SJ.out.tab and are NOT affected by the strand
# reconstruction -- STAR applied no strand filter (its only strand option here
# is --outSAMstrandField intronMotif, which adds an XS tag and drops nothing),
# so junctions were always at full depth.
for T in internal terminal; do
  OUT=$(cdir $T sjcnt sjcnt.TSSTTS)
  backup $OUT
  Rscript $R/bipartition_sjcnt.R --inputdir $BASE/bipartition.filtered \
      --graphmldir $BASE/graphml.v28 --gff $BASE/ref/gencode.v28.dexseq.bygene.gff \
      --sjdir $BASE/STAR --cond1 group1 --cond2 group2 \
      --output $BASE/$OUT --type $T --cores 32
done

echo "=== 2. exonic part counts ==="
for T in internal terminal; do
  OUT=$(cdir $T bipartition.internal.stranded.counts bipartition.TSSTTS.stranded.counts)
  backup $OUT
  Rscript $R/exoncnt.R -c "$CNT" -t bipartition --cond1=group1 --cond2=group2 \
      -a $T -i $BASE/bipartition.filtered -o $BASE/$OUT
done

echo "=== 3. merge in the exclusive junction where the exonic side is empty ==="
backup bipartition.merged.stranded.counts bipartition.merged.TSSTTS.stranded.counts
Rscript $R/merge_exon_sj_counts.R --exon_counts bipartition.internal.stranded.counts \
        --sj_counts sjcnt        --output bipartition.merged.stranded.counts
Rscript $R/merge_exon_sj_counts.R --exon_counts bipartition.TSSTTS.stranded.counts \
        --sj_counts sjcnt.TSSTTS --output bipartition.merged.TSSTTS.stranded.counts

echo "=== 4. merged tests (phi estimated fresh from these counts) ==="
run_test () {   # class countdir outdir phibase annotdir
  local T=$1 CD=$2 OD=$3 PB=$4 AD=$5
  local CF=bipartition.exsj.$T.cnt.combined.txt
  backup "$OD"; mkdir -p "$OD"
  # combine_exoncnt.py aligns columns by NAME: per-gene merged files do not all
  # share one column order, and a positional concatenation corrupts the master.
  [ -s "$CD/$CF" ] || python3 $R/combine_exoncnt.py "$CD" "$CD/$CF"
  for f in $AD/*.bipartition.txt; do ln -sf "$(realpath "$f")" "$CD/$(basename "$f")"; done
  Rscript $R/exontest.R --file "$CF" --countdir "$CD" --outdir "$OD" \
      --splittype bipartition --model betabinom_EBapprox --phi "$PB" \
      --cond1 group1 --cond2 group2 --use_phi_loess --independent_filtering \
      --min_reads 10
  for f in $OD/test_bipartition.exsj.${T}_*; do
    ln -sf "$(basename "$f")" "$OD/$(basename "$f" | sed "s/^test_bipartition\.exsj\.${T}_/test_bipartition.merged_/")"
  done
}
run_test internal bipartition.merged.stranded.counts \
         bipartition.merged.stranded.test.EBapprox phi.merged.stranded.txt \
         bipartition.internal.stranded.counts & P1=$!
run_test terminal bipartition.merged.TSSTTS.stranded.counts \
         bipartition.merged.TSSTTS.stranded.test.EBapprox phi.merged.TSSTTS.stranded.txt \
         bipartition.TSSTTS.stranded.counts & P2=$!
wait $P1 || { echo "internal test FAILED" >&2; exit 1; }
wait $P2 || { echo "terminal test FAILED" >&2; exit 1; }
echo "=== Done ==="
