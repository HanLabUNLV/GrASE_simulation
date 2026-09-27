#!/bin/bash
# Re-run the model comparison PR script now that it stratifies by simulation
# category (Background / DGE / DTE / DTU), assigned per gene from simulate.rda
# exactly as pr_curves_three_levels_gtrule.R does.
cd /mnt/data1/home/mirahan/GrASE_simulation
Rscript scripts/pr_curves_models.R 2>&1 | tail -40
