#!/bin/bash
# run_majiq_00_stage.sh
# Stage inputs for the MAJIQ benchmark:
#   1. Convert the gencode v28 GTF to GFF3 (majiq build needs GFF3, not GTF).
#   2. Work around the group1/group2 BAM basename collision: both groups have
#      identical names (pass2_out_N...), but MAJIQ keys experiments by basename
#      across bamdirs. We stage group-prefixed symlinks (g1_/g2_) into one dir.
# Safe/idempotent: only touches files under majiq/ (a dedicated staging area).
set -euo pipefail

BASE=$HOME/GrASE_simulation
GTF=$BASE/ref/gencode.v28.annotation.gtf
GFF3=$BASE/ref/gencode.v28.annotation.gff3
STAGE=$BASE/majiq/bam

# --- 1. GTF -> GFF3 (gffread lives in the rnaseq env) ---
if [ ! -s "$GFF3" ]; then
  echo "[stage] converting GTF -> GFF3 via gffread (rnaseq env)"
  conda run -n rnaseq gffread "$GTF" --keep-genes -o "$GFF3"
else
  echo "[stage] GFF3 already present: $GFF3"
fi
echo "[stage] GFF3 feature types:"
cut -f3 "$GFF3" | grep -v '^#' | sort | uniq -c | sort -rn | head
# NOTE: if 'majiq build' later rejects this GFF3 (it is picky about the
# gene -> mRNA -> exon ID/Parent hierarchy), fall back to the official
# gencode.v28.annotation.gff3 from gencodegenes.org.

# --- 2. Group-prefixed BAM symlinks into one bamdir ---
mkdir -p "$STAGE"
for g in 1 2; do
  SRC=$BASE/STAR/group${g}/bam
  for n in $(seq 1 12); do
    bam=$SRC/pass2_out_${n}.Aligned.sortedByCoord.out.bam
    [ -s "$bam" ]     || { echo "MISSING bam: $bam"  >&2; exit 1; }
    [ -s "$bam.bai" ] || { echo "MISSING .bai: $bam" >&2; exit 1; }
    link=$STAGE/g${g}_pass2_out_${n}.Aligned.sortedByCoord.out.bam
    rm -f "$link" "$link.bai"          # only removes links in the staging dir
    ln -s "$bam"     "$link"
    ln -s "$bam.bai" "$link.bai"
  done
done
echo "[stage] staged $(ls "$STAGE"/*.bam | wc -l) prefixed BAM symlinks in $STAGE"
ls "$STAGE"/*.bam | xargs -n1 basename
