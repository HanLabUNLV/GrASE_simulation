#!/bin/bash
set -eux

# Set input and output directories
GENOME_DIR=$HOME"/GrASE_simulation/STAR/STAR_index_2.7.10b"
FASTA_FILE=$HOME"/GrASE_simulation/ref/GRCh38.p13.genome.fa"
GTF_FILE=$HOME"/GrASE_simulation/ref/gencode.v34.annotation.gtf"

# Create output directory if it doesn't exist


STAR \
    --runThreadN 24 \
    --runMode genomeGenerate \
    --genomeDir $GENOME_DIR \
    --genomeFastaFiles $FASTA_FILE \
    --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang 99 \

echo "Done generating pass2 genome index"


