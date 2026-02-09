#!/bin/bash
set -eux

CELLTYPE1=$1
CELLTYPE2=$2
RMATS_DIR=$HOME'/miniconda3/envs/rmats.4.2.0/rmats_turbo_v4_2_0'
GTF_FILE=$HOME'/GrASE_simulation/ref/gencode.v28.annotation.gtf'

mkdir -p $HOME/GrASE_simulation/rMATS/rmats_post_${CELLTYPE1}_${CELLTYPE2}
python $RMATS_DIR/rmats.py --b1 $HOME/GrASE_simulation/rMATS/rmats_prep/${CELLTYPE1}.prep.txt --b2  $HOME/GrASE_simulation/rMATS/rmats_prep/${CELLTYPE2}.prep.txt --gtf $GTF_FILE -t single --readLength 50 --nthread 16 --od $HOME/GrASE_simulation/rMATS/rmats_post_${CELLTYPE1}_${CELLTYPE2}/ --tmp $HOME/GrASE_simulation/rMATS/rmats_tmp_post/ --task post --novelSS 

