#!/bin/bash
# run_majiq_03_voila_thresholds.sh
# Export voila tsv at TWO effect-size magnitudes (C = 0.20 and C = 0.10) from the
# SAME deltapsi posterior. deltapsi is NOT rerun -- the dPSI posterior is
# threshold-independent; --threshold only sets the P(|dPSI| >= C) column that
# voila reports. --show-all keeps non-changing LSVs so we can score FP/TN.
# Fast (minutes). Run after run_majiq_02_deltapsi.sh has produced the .voila file.
set -euo pipefail

BASE=$HOME/GrASE_simulation
BUILD=$BASE/majiq/build
DPSI=$BASE/majiq/deltapsi
VOILA=$(ls "$DPSI"/*.deltapsi.voila 2>/dev/null | head -1)

[ -n "$VOILA" ] || { echo "no .deltapsi.voila found; run run_majiq_02_deltapsi.sh first" >&2; exit 1; }
echo "[voila] posterior file: $VOILA"

for C in 0.20 0.10; do
  OUT=$BASE/majiq/majiq_deltapsi.thr${C}.tsv
  echo "[voila] exporting tsv at threshold C=$C -> $OUT"
  conda run -n majiq2 voila tsv "$BUILD/splicegraph.sql" "$VOILA" \
    --threshold "$C" --show-all -f "$OUT"
  echo "[voila]   rows: $(tail -n +2 "$OUT" | grep -vc '^#' || true)"
done

echo "[voila] done. Column header:"
head -1 "$BASE/majiq/majiq_deltapsi.thr0.20.tsv" | tr '\t' '\n' | cat -n
