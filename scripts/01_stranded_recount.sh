#!/bin/bash
# Re-count DEXSeq bins with the strand recovered from each read's true origin.
#
# WHY: the library is unstranded (polyester default), yet the original count script
# ran with --stranded reverse. On unstranded data that rule discards about half
# of every bin's reads at random and does NOT stop an antisense neighbour's reads
# landing in a gene's bins. Validation on chr2: the contaminated bins of
# ENSG00000138382.14 went 1491 -> 0 and 1278 -> 130 while its reference bins
# doubled, and the global new/old ratio was 1.989 across 6,953 bins.
#
# HOW: polyester writes the source transcript into the read name, so the true
# strand is recoverable. That label is carried as BAM membership (not by
# rewriting FLAGs) and counted against a strand-split GFF with --stranded no,
# which is exactly the discrimination a stranded protocol gives. Strand only,
# never gene -- same-strand overlaps stay ambiguous as they would in reality.
#
# Usage: bash scripts/01_stranded_recount.sh [n_parallel]   (default 12)
# STAR_BAM_FMT: printf pattern of the pass 2 BAM directory per group [STAR/%s/bam];
#   STAR_align/bam/%s for 05_star_align.sh output
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
NPAR=${1:-12}
PY=$HOME/miniconda3/envs/smartSim/bin/python      # the env that has HTSeq
DC=$HOME/GrASE/Rpkg/inst/python/dexseq_count.py

mkdir -p strand_recon/bam strand_recon/log
[ -s ref/gencode.v28.dexseq.bygene.plus.gff ] || bash scripts/split_gff_by_strand.sh

work() {
  set -eu
  BASE=/mnt/data1/home/mirahan/GrASE_simulation
  PY=$HOME/miniconda3/envs/smartSim/bin/python
  DC=$HOME/GrASE/Rpkg/inst/python/dexseq_count.py
  G=$1; S=$2; SN=$3                      # group, sample (out_N), sample_0N tag
  L=$BASE/strand_recon/log/${G}_${S}.log
  OUT=$BASE/DEXSeq/count_files_stranded/$G/${S}_counts.txt
  mkdir -p "$(dirname "$OUT")"
  [ -s "$OUT" ] && { echo "skip $G/$S (done)" >> "$L"; return 0; }

  NPY=$BASE/strand_recon/${G}_${S}.npy
  [ -s "$NPY" ] || $PY $BASE/scripts/build_read_strand.py \
      $BASE/STAR/$G/${S}_${SN}_1_shuffled.fa.gz "$NPY" >> "$L" 2>&1

  PFX=$BASE/strand_recon/bam/${G}_${S}
  [ -s "${PFX}.minus.bam" ] || $PY $BASE/scripts/split_bam_by_origin.py \
      $BASE/$(printf "${STAR_BAM_FMT:-STAR/%s/bam}" "$G")/pass2_${S}.Aligned.sortedByCoord.out.bam "$NPY" "$PFX" >> "$L" 2>&1

  for STR in plus minus; do
    $PY $DC --format bam --order pos --paired yes --stranded no \
        $BASE/ref/gencode.v28.dexseq.bygene.$STR.gff \
        ${PFX}.${STR}.bam ${PFX}.${STR}.counts >> "$L" 2>&1
  done
  cat ${PFX}.plus.counts ${PFX}.minus.counts | grep -v "^_" | sed 's/"//g' \
      | sort -k1,1 > "$OUT"
  echo "done $G/$S -> $OUT ($(wc -l < "$OUT") bins)" >> "$L"
}
export -f work; export STAR_BAM_FMT

for G in group1 group2; do
  SN=$([ "$G" = group1 ] && echo sample_01 || echo sample_02)
  for i in $(seq 1 12); do echo "$G out_$i $SN"; done
done | xargs -P "$NPAR" -n 3 bash -c 'work "$@"' _

echo "=== all samples done ==="
ls DEXSeq/count_files_stranded/group1 DEXSeq/count_files_stranded/group2 | wc -l
