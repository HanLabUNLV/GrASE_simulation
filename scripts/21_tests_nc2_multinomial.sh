#!/bin/bash
# Tests for n_choose_2 and multinomial on the STRANDED counts, so the three
# comparison structures share a lane and a read floor. Without this the cross
# panel compares bipartition (stranded, merged, read floor) against two arms
# that are unstranded, unmerged and unfiltered -- all three differences push
# bipartition's FP count down, and Background FPs read 34 vs 2,049 and 870.
#
# --min_reads 10 matches the bipartition arm. phi (n_choose_2) and prec
# (multinomial) are NOT copied from the unstranded runs: they must be estimated
# from these counts, since the counts are what changed.
#
# One axis stays unmatched: bipartition uses the merged exon+SJ unit and no
# merged counts exist for the other two structures. State that in the caption.
#
# Usage: bash scripts/21_tests_nc2_multinomial.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts

for a in internal TSSTTS; do
  echo "=== n_choose_2 $a  $(date +%H:%M:%S) ==="
  Rscript $S/exontest.R \
    --file=n_choose_2.${a}.exoncnt.combined.txt \
    --countdir=n_choose_2.${a}.stranded.counts/ \
    --outdir=n_choose_2.stranded.test.EBapprox \
    --splittype=n_choose_2 --model=betabinom_EBapprox \
    --phi=phi.nc2.${a}.stranded.txt --use_phi_loess --independent_filtering \
    --min_reads=10 --cond1=group1 --cond2=group2 \
    &> nc2.${a}.stranded.log
  echo "  -> $(stat -c %s n_choose_2.stranded.test.EBapprox/test_n_choose_2.${a}_betabinom_EBapprox.txt 2>/dev/null || echo FAILED)"
done

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
echo "=== done ==="
