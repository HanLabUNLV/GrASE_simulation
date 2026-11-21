#!/bin/bash
# Read-signature counts of the simulation, from the strand-split BAMs that
# 01_stranded_recount.sh counts: dexseq_count_sig.py with the same flags and
# strand-split GFFs, writing per-part counts plus read signatures (parts a read
# pair touches, introns it crosses). Plus and minus are combined per sample into
# counts_sig/raw/<group>_<sample>_counts.txt and .txt.sig, as the DICE counts.
# The paired-end branch reproduces the stock per-part counts exactly (chr21,
# group1 out_1 plus). Steps 2-3 (junction sets, per-side counts) need the
# filtered bipartitions and are run after the split.
#
# Usage: bash scripts/03_signature_counts.sh [n_parallel, default 48]
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
NPAR=${1:-48}
mkdir -p counts_sig/raw counts_sig/log

one() {
  set -eu
  BASE=/mnt/data1/home/mirahan/GrASE_simulation
  PY=$HOME/miniconda3/envs/smartSim/bin/python
  DC=$HOME/GrASE/Rpkg/inst/python/dexseq_count_sig.py
  G=$1; S=$2; STR=$3
  PFX=$BASE/strand_recon/bam/${G}_${S}
  OUT=$BASE/counts_sig/raw/${G}_${S}.${STR}
  [ -s "$OUT.counts.sig" ] && return 0
  $PY $DC --format bam --order pos --paired yes --stranded no \
      $BASE/ref/gencode.v28.dexseq.bygene.$STR.gff ${PFX}.${STR}.bam $OUT.part \
      > $BASE/counts_sig/log/${G}_${S}.${STR}.log 2>&1 || { echo "FAILED ${G}_${S} $STR"; return 1; }
  mv $OUT.part.sig $OUT.counts.sig; mv $OUT.part $OUT.counts
}
export -f one
for G in group1 group2; do for i in $(seq 1 12); do for STR in plus minus; do
  echo "$G out_$i $STR"; done; done; done | xargs -P "$NPAR" -L1 bash -c 'one "$@"' _

## combine strands per sample: per-part counts as 01_stranded_recount.sh does, signatures under one header
n=0
for G in group1 group2; do for i in $(seq 1 12); do
  P=counts_sig/raw/${G}_out_${i}
  [ -s $P.plus.counts.sig ] && [ -s $P.minus.counts.sig ] || { echo "missing strand for ${G}_out_$i"; continue; }
  cat $P.plus.counts $P.minus.counts | grep -v "^_" | sed 's/"//g' | sort -k1,1 > ${P}_counts.txt
  { head -1 $P.plus.counts.sig; tail -n +2 $P.plus.counts.sig; tail -n +2 $P.minus.counts.sig; } > ${P}_counts.txt.sig
  n=$((n+1))
done; done
echo "samples combined: $n of 24"
