#!/bin/bash
set -eux


CELLTYPE=$1
# Set input and output directories
FASTQ_DIR=$HOME"/GrASE_simulation/"${CELLTYPE}
SJ_DIR=${FASTQ_DIR}"/sj/"

cat ${SJ_DIR}/*SJ.out.tab \
  | awk '$5 > 0' \
  | sort -k1,1 -k2,2n -k3,3n \
  | uniq > ${CELLTYPE}.raw.SJ.tab
