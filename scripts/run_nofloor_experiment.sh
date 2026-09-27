#!/bin/bash
# "What if no read floor anywhere in the cross comparison?"
#
# bipartition: FREE -- the pre-filter state was snapshotted as
#   *.annotated.txt.bak_premin_reads by scripts/rerun_stranded_min_reads.sh,
#   so it is evaluated directly, no re-test needed.
# n_choose_2 : needs a real re-test at --min_reads=0. The floor sets
#   p.value <- NA BEFORE adjustment, so the raw p-values of filtered rows are
#   gone from the output and the BH adjustment of every other row changed too;
#   there is no post-hoc shortcut.
# multinomial: already has no floor, so it is unchanged and not rerun.
#
# Usage: bash scripts/run_nofloor_experiment.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
EVAL=scripts/evaluate_bipartition_test.R
GT=results/sim_exon_info
SIM=swimdown/simulate/data/simulate.rda

# --- 1. bipartition, no floor, straight from the snapshots ------------------
B1=bipartition.merged.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt.bak_premin_reads
B2=bipartition.merged.TSSTTS.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt.bak_premin_reads
if [ -f "$B1" ] && [ -f "$B2" ]; then
  echo "=== bipartition no-floor eval $(date +%H:%M:%S) ==="
  Rscript "$EVAL" "$B1,$B2" "$GT" results/eval_bipartition_EBapprox.nofloor "$SIM" \
    &> eval_bip_nofloor.log
  echo "  -> $(ls results/eval_bipartition_EBapprox.nofloor 2>/dev/null | wc -l) files"
else
  echo "  [skip] bipartition snapshots missing"
fi
echo "=== bipartition done $(date +%H:%M:%S) ==="
