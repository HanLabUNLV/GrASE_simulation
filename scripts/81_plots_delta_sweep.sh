#!/bin/bash
# Sweep the post-hoc lfc_diff_net threshold (delta) on the stranded merged
# EBapprox results, scored under GT_rule at the unit level -- the same
# framework as the model comparison (pr_curves_models.R) and Figure 2.
#
# The existing 82_plots_posthoc_lfc.sh sweeps the same delta but scores
# under gtIII (results/sim_exon_info, exonic-part level). Both are valid; they
# are not interchangeable, so do not mix their numbers in one paragraph.
#
# Usage: bash scripts/81_plots_delta_sweep.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
Rscript scripts/pr_delta_sweep.R 2>&1 | tail -30
