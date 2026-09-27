#!/bin/bash
# Regenerate the within-bipartition model comparison panels
# (scripts/plots/model_comparison/internalp*.pdf) from the stranded merged
# modelcomp run. phi values above the axis limit are OMITTED, not clamped.
#
# Usage: bash scripts/80_plots_model_comparison.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
Rscript scripts/plot_model_comparison.R
