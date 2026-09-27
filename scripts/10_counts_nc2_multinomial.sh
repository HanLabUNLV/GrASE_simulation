#!/bin/bash
# Stranded exonic-part counts for n_choose_2 and multinomial, so the three
# comparison structures can be compared on MATCHED inputs.
#
# Why: the cross panel currently has the bipartition arm on stranded merged
# counts with the read floor, against n_choose_2 and multinomial on unstranded
# pre-merge counts with no floor. All three differences push bipartition's FP
# count DOWN -- Background FPs read 34 for bipartition vs 2,049 and 870 -- so
# the panel measures lane and filter, not comparison structure. Antisense
# contamination in the unstranded counts is a direct FP source: reads from an
# overlapping opposite-strand gene inflate exonic counts and produce a ratio
# shift in a gene with no real splicing change.
#
# Only -c (count files) and -o (output) differ from the original unstranded
# invocations in scripts/command.txt.
#
# Usage: bash scripts/10_counts_nc2_multinomial.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
S=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
CNT=/mnt/data1/home/mirahan/GrASE_simulation/DEXSeq/count_files_stranded

for t in n_choose_2 multinomial; do
  for a in internal TSS; do
    tag=$([ "$a" = internal ] && echo internal || echo TSSTTS)
    out=${t}.${tag}.stranded.counts
    if [ -s "$out/${t}.${tag}.exoncnt.combined.txt" ]; then
      echo "=== $out already built, skipping ==="; continue
    fi
    echo "=== $t $a -> $out  $(date +%H:%M:%S) ==="
    Rscript $S/exoncnt.R -c $CNT -t $t --cond1=group1 --cond2=group2 \
      -a $a -i /mnt/data1/home/mirahan/GrASE_simulation/${t}.filtered \
      -o /mnt/data1/home/mirahan/GrASE_simulation/$out \
      &> counts.${t}.${tag}.stranded.log
    echo "  -> $(ls $out | wc -l) files"
  done
done
echo "=== done ==="
