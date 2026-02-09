#!/bin/bash
set -eux

# Set input and output directories
PASS2_GENOME_DIR=$HOME"/GrASE_simulation/STAR/STAR_index_2.7.10b_pass2"
FASTA_FILE=$HOME"/GrASE_simulation/ref/GRCh38.p13.genome.fa"
GTF_FILE=$HOME"/GrASE_simulation/ref/gencode.v28.annotation.gtf"
SJDB_FILE=$HOME"/GrASE_simulation/STAR/all_celltypes.final.merged.SJ.tab"

# Create output directory if it doesn't exist

STAR \
    --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $PASS2_GENOME_DIR \
    --genomeFastaFiles $FASTA_FILE \
    --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang 99 \
    --sjdbFileChrStartEnd $SJDB_FILE

echo "Done generating pass2 genome index"


