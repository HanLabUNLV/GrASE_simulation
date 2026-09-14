#!/bin/bash
# rMATS on the strand-reconstructed BAMs.
#
# WHY. rMATS quantifies retained introns from intronic coverage, which is
# strand-blind, so it is exposed to the same antisense contamination that cost
# GrASE 223 gene FPs and DEXSeq 497. It has the lowest gene precision of any
# tool in the benchmark (0.683, 378 gene FPs), and is the only one still on
# uncorrected input.
#
# HOW. Same device as the DEXSeq and MAJIQ reruns: the read's true strand of
# origin is carried as BAM membership, the GTF is split by strand, and the two
# halves are processed independently. A gene lives on exactly one strand, so the
# two event catalogs are disjoint.
#
# CAUTION. rMATS numbers its events from 0 in each run, so the plus and minus
# event IDs COLLIDE. The two strands are therefore kept in separate output
# directories and never blindly concatenated. Gene- and transcript-level scoring
# needs only GeneID / coordinates and is unaffected, but UNIT-level scoring joins
# on "event_type:ID" against results.beforegatefix/sim_junction_gt.txt, which was
# built from the ORIGINAL run's IDs -- that ground truth must be rebuilt
# (scripts/infer_rmats_junctions_gt.R) before any unit-level rMATS row is quoted.
#
# Flags mirror the original run exactly (-t paired for prep, -t single for post,
# --novelSS) so the only variable is the input reads. The prep/post -t mismatch
# is a known harmless inconsistency: rerunning post with -t paired gave
# byte-identical output.
#
# Usage: bash rMATS/run_rmats_stranded.sh
set -euxo pipefail
BASE=$HOME/GrASE_simulation
RM=$HOME/miniconda3/envs/rmats.4.2.0/rmats_turbo_v4_2_0
PY=$HOME/miniconda3/envs/rmats.4.2.0/bin/python
GTF=$BASE/ref/gencode.v28.annotation.gtf
SR=$BASE/strand_recon/bam

for S in plus minus; do
  SYM=$([ "$S" = plus ] && echo '+' || echo '-')
  G=$BASE/ref/gencode.v28.annotation.${S}.gtf
  [ -s "$G" ] || awk -F'\t' -v s="$SYM" '/^#/ || $7==s' "$GTF" > "$G"

  for CT in group1 group2; do
    D=$BASE/rMATS/stranded_$S
    mkdir -p "$D/prep/$CT" "$D/prep/tmp_$CT" "$D/tmp_post"
    ls $SR/${CT}_out_*.${S}.bam | sed -n -e ':a;N;$!ba;s/\n/,/g;p' > "$D/prep/${CT}.prep.txt"
    if [ ! -d "$D/prep/tmp_$CT" ] || [ -z "$(ls "$D/prep/tmp_$CT"/*.rmats 2>/dev/null || true)" ]; then
      $PY $RM/rmats.py --b1 "$D/prep/${CT}.prep.txt" --gtf "$G" -t paired \
          --readLength 100 --nthread 16 --od "$D/prep/$CT" \
          --tmp "$D/prep/tmp_$CT" --task prep
      $PY $RM/cp_with_prefix.py prep_${CT}_ "$D/tmp_post/" "$D/prep/tmp_$CT"/*.rmats
    fi
  done

  D=$BASE/rMATS/stranded_$S
  mkdir -p "$D/post"
  $PY $RM/rmats.py --b1 "$D/prep/group1.prep.txt" --b2 "$D/prep/group2.prep.txt" \
      --gtf "$G" -t single --readLength 100 --nthread 16 \
      --od "$D/post/" --tmp "$D/tmp_post/" --task post --novelSS
done

echo "=== event counts per strand ==="
for S in plus minus; do
  for et in SE A3SS A5SS RI; do
    f=$BASE/rMATS/stranded_$S/post/${et}.MATS.JCEC.txt
    [ -s "$f" ] && echo "  $S $et: $(( $(wc -l < "$f") - 1 ))"
  done
done
echo "=== done ==="
