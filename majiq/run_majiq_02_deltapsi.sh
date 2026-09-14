#!/bin/bash
# run_majiq_02_deltapsi.sh
# MAJIQ deltapsi (group1 vs group2) + voila tsv export.
# --show-all is REQUIRED so the tsv includes ALL quantified LSVs (not just
# changing ones); we need the non-changing LSVs to score FP/TN.
# The tsv reports E(dPSI) and P(|dPSI|>=0.20) per junction (MAJIQ default C=0.20).
set -euo pipefail

BASE=$HOME/GrASE_simulation
BUILD=$BASE/majiq/build
OUT=$BASE/majiq/deltapsi
TSV=$BASE/majiq/majiq_deltapsi.tsv
NPROC=20

[ -d "$BUILD" ] || { echo "run run_majiq_01_build.sh first (no build dir)" >&2; exit 1; }

g1=$(for n in $(seq 1 12); do echo "$BUILD/g1_pass2_out_${n}.Aligned.sortedByCoord.out.majiq"; done)
g2=$(for n in $(seq 1 12); do echo "$BUILD/g2_pass2_out_${n}.Aligned.sortedByCoord.out.majiq"; done)

mkdir -p "$OUT"
echo "[deltapsi] running deltapsi group1 vs group2..."
conda run -n majiq2 majiq deltapsi \
  -grp1 $g1 \
  -grp2 $g2 \
  -n group1 group2 \
  -j "$NPROC" --output-type all -o "$OUT"

VOILA=$(ls "$OUT"/*.deltapsi.voila | head -1)
echo "[deltapsi] voila file: $VOILA"
echo "[deltapsi] exporting tsv (--show-all) ..."
conda run -n majiq2 voila tsv "$BUILD/splicegraph.sql" "$VOILA" \
  --show-all -f "$TSV"

echo "[deltapsi] done. tsv: $TSV"
echo "[deltapsi] rows: $(tail -n +2 "$TSV" | wc -l)"
head -1 "$TSV" | tr '\t' '\n' | head -30
