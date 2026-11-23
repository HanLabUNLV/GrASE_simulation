#!/bin/bash
# Per-part, exon-only route of the GrASE pipeline: the bipartition arm of the
# comparison of structures (bipartition, n_choose_2, multinomial). n_choose_2 and
# multinomial see exonic parts only, so this arm does too:
#   - partitions from bipartition.v28.dedupfix, enumerated with exonic parts only
#     (no intron edges in the comparison of paths, 00_bipartition_v28_dedupfix.sh);
#   - counts from the stranded DEXSeq bins (DEXSeq/count_files_stranded,
#     01_stranded_recount.sh), one count per part;
#   - no junction counts and no junction substitution: a side with no exonic
#     distinct part is not tested.
# The main results use the read-signature route (04_signature_pipeline.sh) instead.
#
# Every output directory ends in .perpart, so nothing of the signature route is
# touched. An existing output directory stops the run, so stale per-gene files
# cannot survive and phi is always estimated fresh (exontest.R skips estimation
# when a phi file is already in --outdir).
#
# Usage: bash scripts/02_stranded_pipeline.sh
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
R=$HOME/GrASE/Rpkg/scripts
CNT=$BASE/DEXSeq/count_files_stranded
BIPART=bipartition.v28.dedupfix
SPLIT=bipartition.filtered.perpart
GFF=$BASE/dexseq.gff

fresh () { for d in "$@"; do [ ! -e $d ] || { echo "$d exists; move it aside first" >&2; exit 1; }; done; }
cdir () { [ $1 = internal ] && echo $2 || echo $3; }   # class internal-name terminal-name

[ -d $BIPART ] || { echo "missing $BIPART (00_bipartition_v28_dedupfix.sh)" >&2; exit 1; }
fresh $SPLIT bipartition.internal.stranded.counts.perpart bipartition.TSSTTS.stranded.counts.perpart \
      bipartition.internal.stranded.test.EBapprox.perpart bipartition.TSSTTS.stranded.test.EBapprox.perpart

echo "=== 0. split bipartitions into internal and terminal  $(date '+%F %T') ==="
Rscript $R/filterTerminal.R --split_dir=$BIPART --outdir=$SPLIT --split_type=bipartition

echo "=== 1. exonic part counts  $(date '+%F %T') ==="
for T in internal terminal; do
  OUT=$(cdir $T bipartition.internal.stranded.counts.perpart bipartition.TSSTTS.stranded.counts.perpart)
  Rscript $R/exoncnt.R -c "$CNT" -t bipartition --cond1=group1 --cond2=group2 \
      -a $T -i $BASE/$SPLIT -o $BASE/$OUT
done

echo "=== 2. tests (phi estimated fresh from these counts)  $(date '+%F %T') ==="
run_test () {   # class countdir outdir phibase
  local T=$1 CD=$2 OD=$3 PB=$4
  mkdir -p "$OD"
  Rscript $R/exontest.R --prefix "bipartition.ex.$T" --countdir "$CD" --outdir "$OD" \
      --splittype bipartition --model betabinom_EBapprox --phi "$PB" \
      --cond1 group1 --cond2 group2 --use_phi_loess --padj_threshold=0.01 \
      --min_reads 10 --gff_dir $GFF --read_length 100 --count_per_part &> "$OD/test.log"
}
run_test internal bipartition.internal.stranded.counts.perpart \
         bipartition.internal.stranded.test.EBapprox.perpart phi.internal.stranded.txt & P1=$!
run_test terminal bipartition.TSSTTS.stranded.counts.perpart \
         bipartition.TSSTTS.stranded.test.EBapprox.perpart phi.TSSTTS.stranded.txt & P2=$!
wait $P1 || { echo "internal test FAILED" >&2; exit 1; }
wait $P2 || { echo "terminal test FAILED" >&2; exit 1; }
echo "=== done  $(date '+%F %T') ==="
