#!/bin/bash
# GrASE pipeline on read-signature counts: the counterpart of 02_stranded_pipeline.sh
# with each read counted once per set (scripts/03_signature_counts.sh) in place of
# per-part DEXSeq counts plus STAR junction counts, and the bipartitions of
# bipartition.v28.features (intron edges in the comparison of paths, partitions
# read off the features; 00_bipartition_v28_features.sh).
#
#   0. split the bipartitions into internal and terminal events
#   1. junction sets of every bipartition (distinct and shared introns)
#   2. per-side counts from the read signatures (sigcnt.py)
#   3. tests, default count mode (no --count_per_part), --sj_overhang 6 as STAR pass 2
#
# Writes the directories 30_merge_downstream.sh, 60_, 61_ and the evaluation
# scripts read. Their current contents (per-part counts) are moved to
# <dir>.before_signature the first time this runs; a rerun never touches that
# backup. 02_stranded_pipeline.sh stays the per-part route.
#
# Needs: counts_sig/raw/<group>_out_<i>_counts.txt.sig for the 24 samples
# (03_signature_counts.sh) and bipartition.v28.features (00_bipartition_v28_features.sh).
#
# Usage: bash scripts/04_signature_pipeline.sh [cores, default 32]
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
R=$HOME/GrASE/Rpkg/scripts
CORES=${1:-32}
BIPART=bipartition.v28.features
SIG=$BASE/counts_sig
GFF=$BASE/dexseq.gff

backup () { for d in "$@"; do [ -e $d.before_signature ] || { [ ! -e $d ] || mv $d $d.before_signature; }; done; }
cdir () { [ $1 = internal ] && echo $2 || echo $3; }   # class internal-name terminal-name

for G in group1 group2; do for i in $(seq 1 12); do
  [ -s $SIG/raw/${G}_out_${i}_counts.txt.sig ] || { echo "missing signatures: ${G}_out_$i" >&2; exit 1; }
done; done
[ -d $BIPART ] && [ ! -e $BIPART/errors.log ] || { echo "$BIPART missing or has errors" >&2; exit 1; }

echo "=== 0. split bipartitions into internal and terminal  $(date '+%F %T') ==="
backup bipartition.filtered
Rscript $R/filterTerminal.R --split_dir=$BIPART --outdir=bipartition.filtered --split_type=bipartition

echo "=== 1. junction sets  $(date '+%F %T') ==="
for T in internal terminal; do
  rm -rf $SIG/junctions.$T
  Rscript $R/bipartition_junction_sets.R --inputdir $BASE/bipartition.filtered --graphmldir $BASE/graphml.v28 \
    --gff $BASE/ref/gencode.v28.dexseq.bygene.gff --type $T --output $SIG/junctions.$T \
    --cores $CORES &> $SIG/junctions.$T.log
done

echo "=== 2. per-side counts from the read signatures  $(date '+%F %T') ==="
for G in group1 group2; do for i in $(seq 1 12); do
  printf "%s\t%s\t%s\n" $SIG/raw/${G}_out_${i}_counts.txt.sig ${G}_out_${i} $G
done; done > $SIG/samples.tsv
rm -rf $SIG/buckets
python3 $R/sigcnt.py split --samples $SIG/samples.tsv --outdir $SIG/buckets &> $SIG/split.log
for T in internal terminal; do
  OUT=$(cdir $T bipartition.merged.stranded.counts bipartition.merged.TSSTTS.stranded.counts)
  backup $OUT
  python3 $R/sigcnt.py count --buckets $SIG/buckets --samples $SIG/samples.tsv \
    --bipartdir $BASE/bipartition.filtered --junctiondir $SIG/junctions.$T --type $T \
    --outdir $BASE/$OUT --cores $CORES &> $SIG/count.$T.log
  printf "  %-8s genes %6s  junction-measured sides %s\n" $T \
    "$(ls $OUT | grep -c '\.cnt\.txt$')" "$(($(wc -l < $OUT/measured_junctions.txt) - 1))"
done

echo "=== 3. tests (phi estimated fresh from these counts)  $(date '+%F %T') ==="
run_test () {   # class countdir outdir phibase
  local T=$1 CD=$2 OD=$3 PB=$4
  backup "$OD"; mkdir -p "$OD"
  Rscript $R/exontest.R --prefix "bipartition.exsj.$T" --countdir "$CD" --outdir "$OD" \
      --splittype bipartition --model betabinom_EBapprox --phi "$PB" \
      --cond1 group1 --cond2 group2 --use_phi_loess --padj_threshold=0.01 \
      --min_reads 10 --gff_dir $GFF --read_length 100 --sj_overhang 6 &> "$OD/test.log"
  for f in $OD/test_bipartition.exsj.${T}_*; do
    ln -sf "$(basename "$f")" "$OD/$(basename "$f" | sed "s/^test_bipartition\.exsj\.${T}_/test_bipartition.merged_/")"
  done
}
run_test internal bipartition.merged.stranded.counts \
         bipartition.merged.stranded.test.EBapprox phi.merged.stranded.txt & P1=$!
run_test terminal bipartition.merged.TSSTTS.stranded.counts \
         bipartition.merged.TSSTTS.stranded.test.EBapprox phi.merged.TSSTTS.stranded.txt & P2=$!
wait $P1 || { echo "internal test FAILED" >&2; exit 1; }
wait $P2 || { echo "terminal test FAILED" >&2; exit 1; }
echo "=== done  $(date '+%F %T') ==="
