#!/bin/bash
# run_majiq_01_build.sh
# MAJIQ build on the staged BAMs. Writes the build config then runs majiq build.
# Matches the rMATS run: unstranded (strandness=None), de novo junctions on,
# intron retention ON (default -- do NOT pass --disable-ir).
# Long step (~1h+ with -j 20); launch in the background.
set -euo pipefail

BASE=$HOME/GrASE_simulation
GFF3=$BASE/ref/gencode.v28.annotation.gff3
STAGE=$BASE/majiq/bam
CONF=$BASE/majiq/conf.ini
OUT=$BASE/majiq/build
NPROC=20

[ -s "$GFF3" ] || { echo "run run_majiq_00_stage.sh first (no GFF3)" >&2; exit 1; }

# comma-separated experiment names = link basenames WITHOUT .bam
g1=$(for n in $(seq 1 12); do echo -n "g1_pass2_out_${n}.Aligned.sortedByCoord.out,"; done | sed 's/,$//')
g2=$(for n in $(seq 1 12); do echo -n "g2_pass2_out_${n}.Aligned.sortedByCoord.out,"; done | sed 's/,$//')

mkdir -p "$BASE/majiq" "$OUT"
cat > "$CONF" <<EOF
[info]
bamdirs=$STAGE
genome=hg38
strandness=None
[experiments]
group1=$g1
group2=$g2
EOF

echo "[build] config written to $CONF:"
cat "$CONF"
echo "[build] starting majiq build (IR on, de novo on, strandness None)..."
conda run -n majiq2 majiq build "$GFF3" -c "$CONF" -j "$NPROC" -o "$OUT"
echo "[build] done. splicegraph + per-sample .majiq files in $OUT"
ls "$OUT"/*.majiq 2>/dev/null | wc -l
ls "$OUT"/splicegraph.sql 2>/dev/null
