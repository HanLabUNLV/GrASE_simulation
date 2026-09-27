#!/bin/bash
# Everything downstream of the merged (exon + split-read) exontest runs.
#
# Usage: bash scripts/30_merge_downstream.sh [stranded|unstranded]
#        (default: stranded -- the lane the manuscript reports)
#
# LANE IS NOT COSMETIC. Three different env vars select it, with three
# different names, and a stage that misses one silently reads the other lane's
# files and produces a table that looks valid:
#   pr_curves_three_levels_gtrule.R   STRANDED=1
#   infer_bipartition_gt.R            GT_RULE_STRANDED=1
#   transcript_level_metrics.R        neither -- it reads BOTH lanes and emits
#                                     GrASE_merged_* and GrASE_merged_stranded_*
#                                     as separate tools
# This script existed from 2026-09-10 with none of them set, so it ran entirely
# unstranded while every reported stranded number came from manual
# STRANDED=1 invocations (see sweep_*.log). Do not drop the exports.
#
#   1. per-event provenance (which side used a junction, and which junction)
#   2. merged unit ground truth -- same GT_rule, junction-sourced sides included
#   3. PR sweep (three levels)   -> results/eval_metricfair/pr_three_levels.GT_rule.txt
#   4. transcript-level metrics  -> transcript_level_metrics.txt
#   5. gene-level tables         -> gene_level_table.{restricted,full}.txt
#
# The sweep caches its call lists in pr_calls_cache.GT_rule.rds and SKIPS every
# call-building block when that file exists, so a new tool or a re-run test
# never appears until the cache is moved aside. Step 3 does that.
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
EV=results/eval_metricfair
LANE=${1:-stranded}

# grase resolves via ~/.Renviron to the shared library at
# /mnt/data1/home/jaquino/R/... -- re-install there after any edit to
# Rpkg/R/*.R, or the previous build's signatures are used silently:
#   R CMD INSTALL --library=/mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2 \
#                 /mnt/data1/home/mirahan/GrASE/Rpkg

if [ "$LANE" = "stranded" ]; then
  export STRANDED=1 GT_RULE_STRANDED=1
  CNT_INT=bipartition.merged.stranded.counts
  CNT_TSS=bipartition.merged.TSSTTS.stranded.counts
  GT=results/gt/bipartition_merged_gt.stranded.txt
else
  unset STRANDED GT_RULE_STRANDED || true
  CNT_INT=bipartition.merged.counts
  CNT_TSS=bipartition.merged.TSSTTS.counts
  GT=results/gt/bipartition_merged_gt.txt
fi
# SJAWARE gates the sweep's dpi floor for junction-sourced sides (0.05 vs the
# exonic 0.1). Set to match the shipped add_significant() default, which is
# source-aware as of 2026-09-14. Unset it here to sweep at a uniform floor.
export SJAWARE=1
echo "=== lane: $LANE  (STRANDED=${STRANDED:-unset} SJAWARE=${SJAWARE:-unset}) ==="

while pgrep -f "exec/R.*[e]xontest" > /dev/null; do sleep 60; done

# 1. Provenance depends only on the merged COUNTS, not on p-values, so a re-run
#    that changed only the tests does not invalidate it. Skipped when present.
#    NOTE: infer_bipartition_gt.R reads merge_provenance.<type>.txt with no
#    stranded variant, so the stranded GT is built against the unstranded
#    sidecar. Pre-existing; changing it would move every stranded GT number.
echo "=== 1. provenance ==="
for t in internal TSSTTS; do
  [ "$t" = internal ] && C=$CNT_INT || C=$CNT_TSS
  if [ -s results/gt/merge_provenance.$t.txt ]; then
    echo "  merge_provenance.$t.txt present, skipping"
  else
    python3 scripts/extract_merge_provenance.py $C $t results/gt/merge_provenance.$t.txt
  fi
done

echo "=== 2. merged unit ground truth ==="
Rscript scripts/infer_bipartition_gt.R --merged --cores 16
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} {n++; if($h["gt_positive"]=="TRUE")p++}
            END{printf "  GT: %d rows, %d positives\n", n, p}' $GT

echo "=== 3. PR sweep (cache invalidated so the re-run tests are picked up) ==="
[ -f $EV/pr_calls_cache.GT_rule.rds ] && \
  mv $EV/pr_calls_cache.GT_rule.rds $EV/pr_calls_cache.GT_rule.rds.prev_$LANE
[ -f $EV/pr_three_levels.GT_rule.txt ] && \
  cp $EV/pr_three_levels.GT_rule.txt $EV/pr_three_levels.GT_rule.txt.prev_$LANE
Rscript scripts/pr_curves_three_levels_gtrule.R

# transcript metrics FIRST: gene_level_table.py reads n_sig_calls out of
# transcript_level_metrics.txt, so running it after leaves any newly added tool
# with n_sig_calls = NA in the gene tables.
echo "=== 4. transcript-level metrics ==="
Rscript scripts/transcript_level_metrics.R

echo "=== 5. gene-level tables ==="
python3 scripts/gene_level_table.py

echo "=== Done ($LANE) ==="
