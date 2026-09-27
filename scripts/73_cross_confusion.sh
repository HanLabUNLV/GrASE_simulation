#!/bin/bash
# Regenerate the cross-comparison panels with the confusion-count figure laid
# out HORIZONTALLY (sim_type panels side by side in one row) instead of
# vertically. Per-panel width is unchanged, so bar width and bar spacing are
# identical to the vertical version -- only the panel arrangement differs.
#
# Usage: bash scripts/73_cross_confusion.sh
set -eu
cd /mnt/data1/home/mirahan/GrASE_simulation
PY=/mnt/data1/home/jaquino/miniconda3/envs/py38/bin/python3

# keep the vertical version for side-by-side comparison
mkdir -p plots/cross.vertical_backup
for f in plots/cross/06_confusion_counts_padj0.01.png \
         plots/cross/06_confusion_counts_padj0.01_restricted.png; do
  [ -f "$f" ] && cp -f "$f" plots/cross.vertical_backup/ || true
done

# --gt-level gtI is stated explicitly rather than left to the script default.
# The level changes ONLY the DTE rows (Background/DGE/DTU are identical at every
# level), and gtI counts all exons of changed transcripts as positive. Relying
# on the default would let a future change to it silently move the DTE numbers.
$PY scripts/visualize_eval.py \
  --results-dir results \
  --out plots \
  --groups cross \
  --gt-level gtI \
  --best-bipartition EBapprox \
  --best-n-choose-2 EBapprox
