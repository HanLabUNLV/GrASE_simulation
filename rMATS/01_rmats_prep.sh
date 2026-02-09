#!/bin/bash
set -eux

CELLTYPE=$1
RMATS_DIR=$HOME'/miniconda3/envs/rmats.4.2.0/rmats_turbo_v4_2_0'
GTF_FILE=$HOME'/GrASE_simulation/ref/gencode.v28.annotation.gtf'

mkdir -p $HOME/GrASE_simulation/rMATS/rmats_prep/
ls $HOME/GrASE_simulation/STAR/${CELLTYPE}/bam/*.Aligned.sortedByCoord.out.bam | sed -n -e ':a;N;$!ba;s/\n/,/g;p' > $HOME/GrASE_simulation/rMATS/rmats_prep/${CELLTYPE}.prep.txt

mkdir -p $HOME/GrASE_simulation/rMATS/rmats_prep/${CELLTYPE}
mkdir -p $HOME/GrASE_simulation/rMATS/rmats_prep/tmp_${CELLTYPE}
python $RMATS_DIR/rmats.py --b1 $HOME/GrASE_simulation/rMATS/rmats_prep/${CELLTYPE}.prep.txt --gtf $GTF_FILE -t single --readLength 50 --nthread 16 --od $HOME/GrASE_simulation/rMATS/rmats_prep/${CELLTYPE} --tmp $HOME/GrASE_simulation/rMATS/rmats_prep/tmp_${CELLTYPE} --task prep

mkdir -p $HOME/GrASE_simulation/rMATS/rmats_tmp_post
python $RMATS_DIR/cp_with_prefix.py prep_${CELLTYPE}_ $HOME/GrASE_simulation/rMATS/rmats_tmp_post/ $HOME/GrASE_simulation/rMATS/rmats_prep/tmp_${CELLTYPE}/*.rmats



