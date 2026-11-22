#!/bin/bash
# Two-pass STAR alignment of the 24 simulated samples against the primary assembly
# (~/GrASE/scripts/star_two_pass.sh), with the settings of the original run except the overhang:
# paired 100 bp reads, GENCODE v28 for both passes, --alignSJDBoverhangMin 1 (as
# DICE; the original run used 6, and the tests' --sj_overhang follows this value),
# junctions filtered per group. Get the genome first with
#   bash ~/GrASE/scripts/primary_assembly.sh ref 28
# Then run 01_stranded_recount.sh with STAR_BAM_FMT=STAR_align/bam/%s.
# Usage: bash scripts/05_star_align.sh [jobs] [threads per job]
# Output: STAR_align/ (see star_two_pass.sh)
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
cd $BASE
JOBS=${1:-12}; THREADS=${2:-4}
STAR_BIN=${STAR_BIN:-$HOME/miniconda3/envs/star2711b/bin/STAR}
SAMTOOLS=${SAMTOOLS:-$HOME/miniconda3/envs/cnr_bam/bin/samtools}
GENOME_FASTA=${GENOME_FASTA:-$BASE/ref/GRCh38.primary_assembly.genome.fa}
OH_MIN=${OH_MIN:-1}
mkdir -p STAR_align
for G in group1 group2; do
  SN=$([ "$G" = group1 ] && echo sample_01 || echo sample_02)
  for i in $(seq 1 12); do
    printf "out_%s\t%s\t%s\t%s\n" $i $G $BASE/STAR/$G/out_${i}_${SN}_1_shuffled.fa.gz $BASE/STAR/$G/out_${i}_${SN}_2_shuffled.fa.gz
  done
done > STAR_align/samples.tsv
bash $HOME/GrASE/scripts/star_two_pass.sh --samples STAR_align/samples.tsv --fasta "$GENOME_FASTA" \
  --gtf ref/gencode.v28.annotation.gtf --outdir STAR_align --sjdb_overhang 99 --overhang_min "$OH_MIN" \
  --star "$STAR_BIN" --samtools "$SAMTOOLS" --jobs "$JOBS" --threads "$THREADS"
