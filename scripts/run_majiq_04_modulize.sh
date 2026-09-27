#!/bin/bash
# voila modulize: classify splice-graph modules into MAJIQ's own event types
# (alternate first exon, alternate last exon, cassette, A5/A3, intron
# retention, ...). Runs on the EXISTING build + deltapsi outputs -- no rebuild,
# no realignment. Used to split MAJIQ's transcript recovery by event class
# using MAJIQ's own definitions rather than an ad-hoc positional heuristic.
set -eux
BASE=$HOME/GrASE_simulation
VOILA=$HOME/miniconda3/envs/majiq2/bin/voila
OUT=$BASE/majiq/modulize
mkdir -p $OUT
$VOILA modulize \
  $BASE/majiq/build/splicegraph.sql \
  $BASE/majiq/deltapsi/group1-group2.deltapsi.voila \
  --overwrite --show-all \
  -d $OUT -j 16
echo "MODULIZE COMPLETE"
ls -la $OUT | head -30
