#!/bin/bash
set -eux

# Set input and output directories
#GENOME_DIR=$HOME"/GrASE_simulation/STAR/STAR_index_2.7.10b"
#GTF_FILE=$HOME"/GrASE_simulation/ref/gencode.v34.annotation.gtf"
GENOME_DIR=$HOME"/GrASE_simulation/STAR/STAR_index_2.7.10b_pass2"
GTF_FILE=$HOME"/GrASE_simulation/ref/gencode.v28.annotation.gtf"


STAR \
    --runThreadN 2 \
    --genomeDir "$GENOME_DIR" \
    --genomeLoad LoadAndExit 



