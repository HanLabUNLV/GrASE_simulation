#!/bin/bash
# Sweep the effect-size cutoff of the call rule, |delta_pi| against |effect_size|
# (log odds ratio), on the stranded merged EBapprox results, scored under GT_rule
# at the unit level -- same framework as 81_plots_delta_sweep.sh.
#
# Usage: bash scripts/82_plots_effect_sweep.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
Rscript scripts/pr_effect_sweep.R 2>&1 | tail -30
