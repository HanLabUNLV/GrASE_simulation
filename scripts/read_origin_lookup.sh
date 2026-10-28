#!/bin/bash
# Which transcript did the reads in a genomic window actually come from?
#
# polyester writes the source transcript into the read name
#   >read33488787/ENST00000361624.2;mate1:402-501;mate2:510-609
# and STAR truncates that at the first "/", so the BAM keeps only "read33488787".
# The mapping survives in the simulated FASTA, so the origin is recoverable by
# joining on the read id. This is exact transcript-of-origin, and it is the only
# strand-like handle this simulation has: simulate_reads.R calls
# simulate_experiment() WITHOUT strand_specific, so polyester defaulted to
# unstranded and there is no protocol strand in the data to filter on.
#
# Usage: bash scripts/read_origin_lookup.sh <group> <sample> <chr:start-end>
#   e.g. bash scripts/read_origin_lookup.sh group1 out_1 chr2:169810081-169810344
set -eu
G=$1; S=$2; REGION=$3
BASE=/mnt/data1/home/mirahan/GrASE_simulation
BAM=$BASE/STAR/$G/bam/pass2_$S.Aligned.sortedByCoord.out.bam
FA=$BASE/STAR/$G/${S}_sample_*_1_shuffled.fa.gz

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
samtools view "$BAM" "$REGION" | cut -f1 | sort -u > "$TMP/ids.txt"
echo "reads in $REGION ($G/$S): $(wc -l < "$TMP/ids.txt")"

zcat $FA | awk -v idf="$TMP/ids.txt" '
  BEGIN { while ((getline l < idf) > 0) want[l] = 1 }
  /^>/ {
    name = substr($0, 2)
    slash = index(name, "/")
    id = substr(name, 1, slash - 1)
    if (id in want) {
      rest = substr(name, slash + 1)
      semi = index(rest, ";")
      tx = (semi ? substr(rest, 1, semi - 1) : rest)
      cnt[tx]++
    }
  }
  END { for (t in cnt) print cnt[t], t }' | sort -rn
