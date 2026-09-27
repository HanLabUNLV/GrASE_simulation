#!/bin/bash
# Multinomial arms only. n_choose_2 already completed on 09-18; the previous
# driver aborted here because exontest.R built its read-floor support table
# from diff1/diff2, which the multinomial long-format counts do not have.
#
# exontest.R now guards that block, so multinomial runs with --min_reads NOT
# applied (the floor has no definition without a distinct/reference split).
#
# Usage: bash scripts/run_stranded_tests_multi.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts

for a in internal TSSTTS; do
  echo "=== multinomial $a  $(date +%H:%M:%S) ==="
  Rscript $S/exontest.R \
    --file=multinomial.${a}.exoncnt.combined.txt \
    --countdir=multinomial.${a}.stranded.counts/ \
    --outdir=multinomial.stranded.test.EBplugin \
    --splittype=multinomial --model=dirmult_EBplugin \
    --prec=prec_multinomial.${a}.stranded.txt \
    --min_reads=10 --cond1=group1 --cond2=group2 \
    &> multi.${a}.stranded.log
  echo "  -> $(stat -c %s multinomial.stranded.test.EBplugin/test_multinomial.${a}_dirmult_EBplugin.txt 2>/dev/null || echo FAILED)"
done
echo "=== done $(date +%H:%M:%S) ==="
