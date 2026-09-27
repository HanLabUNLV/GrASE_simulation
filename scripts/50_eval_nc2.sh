#!/bin/bash
# Evaluate the STRANDED n_choose_2 EBapprox results (run 09-18), replacing the
# 2026-07-13 unstranded evals that the cross panel still reads.
#
# Same evaluator, GT and simulate.rda as scripts/tmp.nc2.sh -- only the test
# directory changes, from n_choose_2.test to n_choose_2.stranded.test.EBapprox.
# Only EBapprox was run on the stranded counts, so EBmap/MLE/wilcoxon are not
# re-evaluated here and their eval dirs remain unstranded.
#
# The previous eval dirs are moved aside to *.pre_stranded rather than
# overwritten, so the old numbers stay recoverable.
#
# Usage: bash scripts/50_eval_nc2.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
EVAL=scripts/evaluate_bipartition_test.R
GT=results/sim_exon_info
SIM=swimdown/simulate/data/simulate.rda
T=n_choose_2.stranded.test.EBapprox

run_one () {   # $1 = file suffix, $2 = eval dir name
  local sfx="$1" out="results/$2"
  local f1="$T/test_n_choose_2.internal_betabinom_EBapprox${sfx}.annotated.txt"
  local f2="$T/test_n_choose_2.TSSTTS_betabinom_EBapprox${sfx}.annotated.txt"
  if [ ! -f "$f1" ] || [ ! -f "$f2" ]; then
    echo "  [skip] $2 -- missing an arm"; return
  fi
  if [ -d "$out" ] && [ ! -d "${out}.pre_stranded" ]; then
    mv -f "$out" "${out}.pre_stranded"
    echo "  backed up $out -> ${out}.pre_stranded"
  fi
  echo "=== $2  $(date +%H:%M:%S) ==="
  Rscript "$EVAL" "$f1,$f2" "$GT" "$out" "$SIM" &> "eval_nc2_stranded_$2.log"
  echo "  -> $(ls "$out" 2>/dev/null | wc -l) files"
}

run_one ""                  eval_n_choose_2_EBapprox
run_one ".mincomb"          eval_n_choose_2_EBapprox_mincomb
run_one ".fisher_combined"  eval_n_choose_2_EBapprox_fishercomb
echo "=== done $(date +%H:%M:%S) ==="
