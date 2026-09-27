#!/bin/bash
# Cross panel with NO read floor on any structure, for comparison against the
# floored version in plots/cross.
#
# visualize_eval.py resolves eval directories by fixed name, so rather than
# edit the script a shadow results tree is built whose expected names are
# SYMLINKS to the .nofloor dirs. multinomial never had a floor, so its link
# points at the ordinary dir. Nothing is copied and nothing existing is
# modified.
#
# gtI is pinned explicitly, matching scripts/73_cross_confusion.sh.
#
# Usage: bash scripts/71_cross_nofloor.sh
set -u
cd /mnt/data1/home/mirahan/GrASE_simulation
PY=/mnt/data1/home/jaquino/miniconda3/envs/py38/bin/python3
SHADOW=results.nofloor
rm -rf "$SHADOW"; mkdir -p "$SHADOW"
R=$(pwd)/results

link () {  # $1 = expected name, $2 = actual dir
  if [ ! -d "$R/$2" ]; then echo "  MISSING $2"; return 1; fi
  ln -sfn "$R/$2" "$SHADOW/$1"; echo "  $1 -> $2"
}
link eval_bipartition_EBapprox   eval_bipartition_EBapprox.nofloor
link eval_n_choose_2_EBapprox    eval_n_choose_2_EBapprox.nofloor
link eval_multinomial_EBplugin   eval_multinomial_EBplugin
# the cross group also consults the bipartition/nc2 summaries for best-method
# selection; both are specified explicitly below so no other dirs are needed.

$PY scripts/visualize_eval.py \
  --results-dir "$SHADOW" \
  --out plots.nofloor \
  --groups cross \
  --gt-level gtI \
  --best-bipartition EBapprox \
  --best-n-choose-2 EBapprox
