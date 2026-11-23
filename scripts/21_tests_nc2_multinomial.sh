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
# A structure and class whose annotated result exists is skipped; an interrupted
# test resumes from exontest.R's phi file and test checkpoints in the outdir.
#
# Usage: [MC_CORES=48] bash scripts/21_tests_nc2_multinomial.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
MC=${MC_CORES:-48}
done_ () { [ -s "$1/test_$2.ex.$3_$4.annotated.txt" ] && echo "  $2 $3: result exists, skipping"; }   # outdir split class model

# exoncnt.R names the combined counts <split>.ex.<internal|terminal>.cnt.combined.txt and
# exontest.R names its results after them; each result is linked under the
# test_<split>.<internal|TSSTTS>_* name the evaluation scripts (50_, 51_) read.
link_results () {   # link_results <outdir> <split> <class> <tag>
  for f in $1/test_$2.ex.$3_*; do
    ln -sf "$(basename "$f")" "$1/$(basename "$f" | sed "s/^test_$2\.ex\.$3_/test_$2.$4_/")"
  done
}
for a in internal terminal; do
  tag=$([ "$a" = internal ] && echo internal || echo TSSTTS)
  echo "=== n_choose_2 $a  $(date +%H:%M:%S) ==="
  done_ n_choose_2.stranded.test.EBapprox n_choose_2 $a betabinom_EBapprox && continue
  Rscript $S/exontest.R \
    --file=n_choose_2.ex.${a}.cnt.combined.txt \
    --countdir=n_choose_2.${tag}.stranded.counts/ \
    --outdir=n_choose_2.stranded.test.EBapprox \
    --splittype=n_choose_2 --model=betabinom_EBapprox \
    --phi=phi.nc2.${tag}.stranded.txt --use_phi_loess --padj_threshold=0.01 \
    --min_reads=10 --cond1=group1 --cond2=group2 --mc_cores=$MC \
    &> nc2.${tag}.stranded.log
  link_results n_choose_2.stranded.test.EBapprox n_choose_2 $a $tag
  echo "  -> $(stat -L -c %s n_choose_2.stranded.test.EBapprox/test_n_choose_2.${tag}_betabinom_EBapprox.txt 2>/dev/null || echo FAILED)"
done

for a in internal terminal; do
  tag=$([ "$a" = internal ] && echo internal || echo TSSTTS)
  echo "=== multinomial $a  $(date +%H:%M:%S) ==="
  done_ multinomial.stranded.test.EBplugin multinomial $a dirmult_EBplugin && continue
  Rscript $S/exontest.R --padj_threshold=0.01 \
    --file=multinomial.ex.${a}.cnt.combined.txt \
    --countdir=multinomial.${tag}.stranded.counts/ \
    --outdir=multinomial.stranded.test.EBplugin \
    --splittype=multinomial --model=dirmult_EBplugin \
    --prec=prec_multinomial.${tag}.stranded.txt \
    --min_reads=10 --cond1=group1 --cond2=group2 --mc_cores=$MC \
    &> multi.${tag}.stranded.log
  link_results multinomial.stranded.test.EBplugin multinomial $a $tag
  echo "  -> $(stat -L -c %s multinomial.stranded.test.EBplugin/test_multinomial.${tag}_dirmult_EBplugin.txt 2>/dev/null || echo FAILED)"
done
echo "=== done ==="
