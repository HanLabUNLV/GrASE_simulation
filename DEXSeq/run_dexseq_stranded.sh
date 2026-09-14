#!/bin/bash
# DEXSeq on the strand-reconstructed counts.
#
# The original run used DEXSeq/count_files, made with dexseq_count --stranded
# reverse on an UNSTRANDED library: half of every bin's reads discarded at
# random, and antisense neighbours' reads counted into a gene's bins. DEXSeq had
# the worst gene precision of any tool in the benchmark (0.755, 604 gene FPs),
# largely from that artifact -- the same one that cost GrASE 223 of its 257 gene
# FPs. Re-run here on count_files_stranded so the cross-tool comparison is fair.
#
# New outdir; the original dexseq_group1_group2/ is untouched. The dxd.*.rds
# caches are per-outdir, so nothing stale is picked up.
#
# Usage: bash run_dexseq_stranded.sh
set -eux
cd /mnt/data1/home/mirahan/GrASE_simulation/DEXSeq
Rscript dexseq.R \
    --cntdir ../DEXSeq/count_files_stranded \
    --outdir dexseq_group1_group2_stranded \
    --gff    ../ref/gencode.v28.dexseq.bygene.gff \
    --cell1  group1 --cell2 group2
