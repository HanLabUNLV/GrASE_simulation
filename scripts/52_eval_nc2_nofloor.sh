#!/bin/bash
# Evaluate the n_choose_2 NO-FLOOR tests (--min_reads=0, run 09-20) so they can
# be compared against the floored evals already in results/eval_n_choose_2_*.
#
# Same evaluator, GT and simulate.rda as scripts/50_eval_nc2.sh --
# only the test directory changes, so the floor is the single difference.
# Writes to *.nofloor dirs; nothing existing is overwritten.
#
# Usage: bash scripts/52_eval_nc2_nofloor.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
EVAL=scripts/evaluate_bipartition_test.R
GT=results/sim_exon_info
SIM=swimdown/simulate/data/simulate.rda
T=n_choose_2.stranded.test.EBapprox.nofloor

run_one () {   # $1 = file suffix, $2 = eval dir name
  local sfx="$1" out="results/$2"
  local f1="$T/test_n_choose_2.internal_betabinom_EBapprox${sfx}.annotated.txt"
  local f2="$T/test_n_choose_2.TSSTTS_betabinom_EBapprox${sfx}.annotated.txt"
  if [ ! -f "$f1" ] || [ ! -f "$f2" ]; then
    echo "  [skip] $2 -- missing an arm"; return
  fi
  echo "=== $2  $(date +%H:%M:%S) ==="
  Rscript "$EVAL" "$f1,$f2" "$GT" "$out" "$SIM" &> "eval_nc2_nofloor_$2.log"
  echo "  -> $(ls "$out" 2>/dev/null | wc -l) files"
}

run_one ""                  eval_n_choose_2_EBapprox.nofloor
run_one ".mincomb"          eval_n_choose_2_EBapprox_mincomb.nofloor
run_one ".fisher_combined"  eval_n_choose_2_EBapprox_fishercomb.nofloor
echo "=== done $(date +%H:%M:%S) ==="
