#!/bin/bash
# Rerun the simulation from the split-read counts down, after the junction
# distinct set became the INTERSECTION over a side's routes (introns every route
# crosses; NA when none). The per-route sum counted a junction once per route
# that crosses it.
#
# Every output this touches is moved to <name>.$STAMP first, so each stage
# starts empty (bipartition_sjcnt.R skips genes whose output exists, and
# 02_stranded_pipeline.sh's backup only moves a directory the first time).
# Exonic counts are unchanged and reused. Tests add --sj_overhang 6 (STAR
# --alignSJDBoverhangMin 6, STAR/08_pass2.sh) so junction-sourced sides get a
# per-base pi; the call rule is the package default (no delta_pi gate).
#
# Usage: bash scripts/92_sj_intersection_rerun.sh [backup suffix, default bak_routesum]
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
R=$HOME/GrASE/Rpkg/scripts
STAMP=${1:-bak_routesum}
GFF=$BASE/dexseq.gff

for f in sjcnt sjcnt.TSSTTS bipartition.merged.stranded.counts bipartition.merged.TSSTTS.stranded.counts \
         bipartition.merged.stranded.test.EBapprox bipartition.merged.TSSTTS.stranded.test.EBapprox \
         results/gt/count_source.internal.stranded.txt results/gt/count_source.TSSTTS.stranded.txt \
         results/gt/bipartition_merged_gt.stranded.txt; do
  if [ -e $f ]; then
    [ ! -e $f.$STAMP ] || { echo "$f.$STAMP already exists; pick another suffix" >&2; exit 1; }
    mv $f $f.$STAMP
  fi
done

echo "=== 1. split-read counts  $(date '+%F %T') ==="
for T in internal terminal; do
  [ $T = internal ] && OUT=sjcnt || OUT=sjcnt.TSSTTS
  Rscript $R/bipartition_sjcnt.R --inputdir $BASE/bipartition.filtered \
      --graphmldir $BASE/graphml.v28 --gff $BASE/ref/gencode.v28.dexseq.bygene.gff \
      --sjdir $BASE/STAR --cond1 group1 --cond2 group2 \
      --output $BASE/$OUT --type $T --cores 32
done

echo "=== 2. merge  $(date '+%F %T') ==="
Rscript $R/merge_exon_sj_counts.R --exon_counts bipartition.internal.stranded.counts \
        --sj_counts sjcnt        --output bipartition.merged.stranded.counts
Rscript $R/merge_exon_sj_counts.R --exon_counts bipartition.TSSTTS.stranded.counts \
        --sj_counts sjcnt.TSSTTS --output bipartition.merged.TSSTTS.stranded.counts

echo "=== 3. merged tests  $(date '+%F %T') ==="
run_test () {   # class countdir outdir phibase annotdir
  local T=$1 CD=$2 OD=$3 PB=$4 AD=$5
  local CF=bipartition.exsj.$T.cnt.combined.txt
  mkdir -p "$OD"
  [ -s "$CD/$CF" ] || python3 $R/combine_exoncnt.py "$CD" "$CD/$CF"
  for f in $AD/*.bipartition.txt; do ln -sfr "$(realpath "$f")" "$CD/$(basename "$f")"; done
  Rscript $R/exontest.R --file "$CF" --countdir "$CD" --outdir "$OD" \
      --splittype bipartition --model betabinom_EBapprox --phi "$PB" \
      --cond1 group1 --cond2 group2 --use_phi_loess --padj_threshold=0.01 \
      --min_reads 10 --gff_dir $GFF --read_length 100 --count_per_part --sj_overhang 6 &> "$OD/test.log"
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

echo "=== 4. count source, merged ground truth, sweep  $(date '+%F %T') ==="
bash scripts/30_merge_downstream.sh stranded
echo "=== 5. ROC, transcript, lfc sweep, tables  $(date '+%F %T') ==="
bash scripts/91_default_gate_rerun.sh $STAMP
echo "=== Done  $(date '+%F %T') ==="
