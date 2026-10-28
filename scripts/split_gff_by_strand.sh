#!/bin/bash
# Split the flattened DEXSeq GFF by gene strand. Paired with the origin-split
# BAMs and --stranded no, this is what makes the strand filter work: against the
# FULL gff a read in an antisense overlap sees features of both genes and is
# dropped as _ambiguous, whereas against the matching half it sees only the gene
# it actually came from.
set -eu
BASE=/mnt/data1/home/mirahan/GrASE_simulation
GFF=$BASE/ref/gencode.v28.dexseq.bygene.gff
awk -F'\t' '$7=="+"' "$GFF" > $BASE/ref/gencode.v28.dexseq.bygene.plus.gff
awk -F'\t' '$7=="-"' "$GFF" > $BASE/ref/gencode.v28.dexseq.bygene.minus.gff
wc -l "$GFF" $BASE/ref/gencode.v28.dexseq.bygene.plus.gff $BASE/ref/gencode.v28.dexseq.bygene.minus.gff
