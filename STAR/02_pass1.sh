#!/bin/bash
set -eux

CELLTYPE=$1
SRRID=$2
# Set input and output directories
FASTQ_DIR=$HOME"/GrASE_simulation/"${CELLTYPE}
SJ_DIR=${FASTQ_DIR}"/sj/"
GENOME_DIR=$HOME"/GrASE_simulation/STAR/STAR_index_2.7.10b"
GTF_FILE=$HOME"/GrASE_simulation/ref/gencode.v34.annotation.gtf"

# Create output directory if it doesn't exist
mkdir -p "$SJ_DIR"

echo "Processing $SRRID..."

STAR \
    --runThreadN 4 \
    --genomeDir "$GENOME_DIR" \
    --genomeLoad LoadAndKeep \
    --alignIntronMax 299999 \
    --alignSJDBoverhangMin 6 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNmax 3 \
    --chimSegmentMin 2 \
    --readFilesIn "${FASTQ_DIR}/${SRRID}_sample_02_1_shuffled.fa.gz" "${FASTQ_DIR}/${SRRID}_sample_02_2_shuffled.fa.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${SJ_DIR}/pass1_${SRRID}." \
    --outTmpDir "pass1_${SRRID}_${$}_STARtmp" \
    --outSAMtype None \

echo "Done processing $SRRID"


