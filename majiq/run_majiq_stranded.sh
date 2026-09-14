#!/bin/bash
# MAJIQ on the strand-reconstructed BAMs -- fairness check for the benchmark.
#
# WHY. MAJIQ was run with strandness=None (correct: the library is unstranded),
# so it assigns reads to genes without strand. It is largely immune to the
# antisense contamination that cost GrASE 223 gene FPs, because a junction is
# identified by donor/acceptor coordinates that must match the gene's own
# annotated splice sites. But it is not FULLY immune: intron retention is
# quantified from raw intronic coverage, which is strand-blind. Measured:
#   MAJIQ_C0.20  33 of 52 gene FPs (63.5%) overlap a manipulated neighbour
#                vs a 15.7% baseline among tested-but-uncalled null genes
#   MAJIQ_C0.10  58 of 135 (43.0%) vs 15.5%
# a 4x enrichment (GrASE's was 6.4x). So MAJIQ carries some of the same artifact
# and the cross-tool comparison is only fair if it gets the same correction.
#
# HOW. Same device as the DEXSeq recount: the read's true strand of origin
# (recovered from the polyester read name) is carried as BAM membership, and the
# annotation is split by strand, so a plus-origin read can only reach a plus
# strand gene. Two independent builds, concatenated at the TSV stage -- a gene
# lives on exactly one strand, so LSV ids cannot collide.
#
# Usage: bash majiq/run_majiq_stranded.sh
set -euo pipefail
BASE=$HOME/GrASE_simulation
GFF3=$BASE/ref/gencode.v28.annotation.gff3
SR=$BASE/strand_recon/bam
NPROC=20

for S in plus minus; do
  SYM=$([ "$S" = plus ] && echo '+' || echo '-')
  G3=$BASE/ref/gencode.v28.annotation.${S}.gff3
  [ -s "$G3" ] || awk -F'\t' -v s="$SYM" '/^#/ || $7==s' "$GFF3" > "$G3"

  STAGE=$BASE/majiq/bam_$S
  mkdir -p "$STAGE"
  for g in 1 2; do
    for n in $(seq 1 12); do
      src=$SR/group${g}_out_${n}.${S}.bam
      dst=$STAGE/g${g}_pass2_out_${n}.Aligned.sortedByCoord.out.bam
      ln -sf "$src" "$dst"; ln -sf "$src.bai" "$dst.bai"
    done
  done

  CONF=$BASE/majiq/conf.${S}.ini
  g1=$(for n in $(seq 1 12); do echo -n "g1_pass2_out_${n}.Aligned.sortedByCoord.out,"; done | sed 's/,$//')
  g2=$(for n in $(seq 1 12); do echo -n "g2_pass2_out_${n}.Aligned.sortedByCoord.out,"; done | sed 's/,$//')
  cat > "$CONF" <<INI
[info]
bamdirs=$STAGE
genome=hg38
strandness=None
[experiments]
group1=$g1
group2=$g2
INI

  BUILD=$BASE/majiq/build_$S
  if [ ! -s "$BUILD/splicegraph.sql" ]; then
    echo "[$S] build..."
    conda run -n majiq2 majiq build "$G3" -c "$CONF" -j "$NPROC" -o "$BUILD"
  fi

  DP=$BASE/majiq/deltapsi_$S
  if [ -z "$(ls "$DP"/*.deltapsi.voila 2>/dev/null || true)" ]; then
    echo "[$S] deltapsi..."
    a=$(for n in $(seq 1 12); do echo "$BUILD/g1_pass2_out_${n}.Aligned.sortedByCoord.out.majiq"; done)
    b=$(for n in $(seq 1 12); do echo "$BUILD/g2_pass2_out_${n}.Aligned.sortedByCoord.out.majiq"; done)
    mkdir -p "$DP"
    conda run -n majiq2 majiq deltapsi -grp1 $a -grp2 $b -n group1 group2 \
        -j "$NPROC" --output-type all -o "$DP"
  fi

  V=$(ls "$DP"/*.deltapsi.voila | head -1)
  for C in 0.20 0.10; do
    conda run -n majiq2 voila tsv "$BUILD/splicegraph.sql" "$V" \
        --threshold "$C" --show-all -f "$BASE/majiq/majiq_deltapsi.stranded.${S}.thr${C}.tsv"
  done
done

# concatenate the two strands (drop the second file's comment/header block)
for C in 0.20 0.10; do
  OUT=$BASE/majiq/majiq_deltapsi.stranded.thr${C}.tsv
  cp "$BASE/majiq/majiq_deltapsi.stranded.plus.thr${C}.tsv" "$OUT"
  grep -v '^#' "$BASE/majiq/majiq_deltapsi.stranded.minus.thr${C}.tsv" | tail -n +2 >> "$OUT"
  echo "[merge] C=$C rows: $(grep -vc '^#' "$OUT")"
done
echo "=== done ==="
