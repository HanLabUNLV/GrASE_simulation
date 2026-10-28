#!/bin/bash
# Re-enumerate the simulation bipartitions from graphml.v28 with the corrected
# dedup in grase::bipartition_paths. Among rows sharing (setdiff1, setdiff2) it
# now keeps the ref_ex_part CLOSEST to the distinct parts first, then the
# fewest reference parts -- counting R and L as reference slots -- and on a
# remaining tie keeps the internal row.
#
# The previous order took fewest reference parts first and did not count R/L.
# A skipped exon has two flanking reference parts; an outer bubble anchored at
# R or L had only one. So the skip's internal bubble was discarded and the exon
# was tested as a TERMINAL event -- on the order of 17,800 rMATS SE events in
# this simulation landed on such a terminal survivor.
#
# Writes to bipartition.v28.dedupfix, NOT bipartition.v28. The existing
# benchmark was built on bipartition.v28 and the before/after comparison needs
# it intact. Arguments match the README and DICE's run_bipartition.sh except
# --outdir and --mc_cores. Unlike the original v28 run this also writes
# <gene>.bipartition.all.txt, the table BEFORE the dedup, so a later rule
# change can be evaluated without re-enumerating.
#
# Usage: bash scripts/00_bipartition_v28_dedupfix.sh
set -eu -o pipefail
BASE=/mnt/data1/home/mirahan/GrASE_simulation
R=/mnt/data1/home/mirahan/GrASE/Rpkg/scripts
cd $BASE

# grase resolves via ~/.Renviron to the shared library at
# /mnt/data1/home/jaquino/R/... -- re-install there after any edit to
# Rpkg/R/*.R, or the previous build is used silently:
#   R CMD INSTALL --library=/mnt/data1/home/jaquino/R/x86_64-pc-linux-gnu-library/4.2 \
#                 /mnt/data1/home/mirahan/GrASE/Rpkg
# Refuse to run on a build without the corrected order, rather than produce a
# second copy of the old bipartitions under the new name.
ok=$(Rscript -e 'd <- deparse(grase::bipartition_paths)
  p <- grep(".prox ==", d, fixed = TRUE)[1]; n <- grep("ref_part_cnt ==", d, fixed = TRUE)[1]
  r <- grep("source == \"R\"", d, fixed = TRUE)
  cat(!is.na(p) && !is.na(n) && p < n && length(r) > 0)')
if [ "$ok" != "TRUE" ]; then
  echo "installed grase does not have the corrected dedup -- reinstall first" >&2
  exit 1
fi

OUT=$BASE/bipartition.v28.dedupfix
mkdir -p "$OUT"
echo "start: $(date)"
Rscript $R/bubble_path_split.R \
  --graphdir="$BASE/graphml.v28" \
  --outdir="$OUT" \
  --split=bipartition \
  --mc_cores=48
echo "end: $(date)"
echo "deduped: $(ls "$OUT" | grep -c '\.bipartition\.txt$')  pre-dedup: $(ls "$OUT" | grep -c '\.bipartition\.all\.txt$')"
