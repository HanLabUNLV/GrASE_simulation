#!/bin/bash
# Everything downstream of the merged (exon + split-read) exontest runs.
#
# Usage: bash scripts/30_merge_downstream.sh
#
# Scores the signature tests (04_signature_pipeline.sh) and the other tools on the
# strand-reconstructed inputs. infer_bipartition_gt.R still selects its lane by
# GT_RULE_STRANDED, so it is exported here; the sweep and the transcript metrics
# read the stranded inputs only.
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
LANE=stranded
export STRANDED=1 GT_RULE_STRANDED=1
CNT_INT=bipartition.merged.stranded.counts
CNT_TSS=bipartition.merged.TSSTTS.stranded.counts
GT=results/gt/bipartition_merged_gt.stranded.txt
SFX=.stranded
# SJAWARE gates the sweep's dpi floor for junction-sourced sides (0.05 vs the
# exonic 0.1). Set to match the shipped add_significant() default, which is
# source-aware as of 2026-09-14. Unset it here to sweep at a uniform floor.
export SJAWARE=1
echo "=== lane: $LANE  (STRANDED=${STRANDED:-unset} SJAWARE=${SJAWARE:-unset}) ==="

# wait for exontest runs of THIS tree only (working directory here); an
# exontest elsewhere on the machine (DICE) must not block the evaluation
sim_exontest () {
  local p
  for p in $(pgrep -f "exec/R.*[e]xontest"); do
    [ "$(readlink /proc/$p/cwd 2>/dev/null)" = "$PWD" ] && return 0
  done
  return 1
}
while sim_exontest; do sleep 60; done

# 1. Count source per event, from the merged COUNTS of this lane. It does not
#    depend on p-values, so a re-run that changed only the tests does not
#    invalidate it. Skipped when present; a re-run of the merge must move it
#    aside, because events are row numbers of the per-gene bipartition files.
echo "=== 1. count source ==="
for t in internal TSSTTS; do
  [ "$t" = internal ] && C=$CNT_INT || C=$CNT_TSS
  CS=results/gt/count_source.$t$SFX.txt
  if [ -s $CS ]; then
    echo "  $CS present, skipping"
  else
    python3 scripts/extract_count_source.py $C $t $CS
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

echo "=== 4. transcript-level metrics ==="
Rscript scripts/transcript_level_metrics.R

echo "=== 5. gene-level tables ==="
python3 scripts/gene_level_table.py

echo "=== Done ($LANE) ==="
