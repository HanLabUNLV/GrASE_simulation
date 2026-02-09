#!/bin/bash
set -eux

CELLTYPE=$1
SRRID=$2

if [ "$CELLTYPE" == "group1" ]; then
    SAMPLE_NUM="01"
elif [ "$CELLTYPE" == "group2" ]; then
    SAMPLE_NUM="02"
else
    echo "Error: Unknown CELLTYPE $CELLTYPE"
    exit 1
fi

# Set input and output directories
FASTQ_DIR=$HOME"/GrASE_simulation/"${CELLTYPE}
BAM_DIR=${FASTQ_DIR}"/bam/"
GENOME_DIR=$HOME"/GrASE_simulation/STAR/STAR_index_2.7.10b_pass2"
GTF_FILE=$HOME"/GrASE_simulation/ref/gencode.v28.annotation.gtf"

# Create output directory if it doesn't exist
mkdir -p "$BAM_DIR"
LOG_FILE="${BAM_DIR}/pass2_${SRRID}.Log.final.out"
#BAM_FILE="${BAM_DIR}/pass2_${SRRID}.Aligned.out.bam"
BAM_FILE="${BAM_DIR}/pass2_${SRRID}.Aligned.sortedByCoord.out.bam"


# Check if STAR already ran successfully
if [[ -f "$LOG_FILE" && -s "$BAM_FILE" ]]; then
    echo "STAR already run for $SRRID, skipping..."
else
  echo "Processing $SRRID..."
  STAR \
      --runThreadN 4 \
      --genomeDir "$GENOME_DIR" \
      --genomeLoad LoadAndKeep \
      --limitBAMsortRAM 32000000000 \
      --alignIntronMax 299999 \
      --alignSJDBoverhangMin 6 \
      --alignEndsType EndToEnd \
      --outFilterMismatchNmax 3 \
      --chimSegmentMin 2 \
      --readFilesIn "${FASTQ_DIR}/${SRRID}_sample_${SAMPLE_NUM}_1_shuffled.fa.gz" "${FASTQ_DIR}/${SRRID}_sample_${SAMPLE_NUM}_2_shuffled.fa.gz" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${BAM_DIR}/pass2_${SRRID}." \
      --outTmpDir "pass2_${SRRID}_${$}_STARtmp" \
      --outSAMstrandField intronMotif \
      --outSAMtype BAM SortedByCoordinate 

  echo "Done processing $SRRID"
fi
