#!/bin/bash
# Structural junction-level ground truth for the STRANDED rMATS runs.
#
# rMATS numbers its events from 0 in every run, so the plus and minus catalogs
# reuse the same IDs and neither matches the original numbering that
# results.beforegatefix/sim_junction_gt.txt is keyed to. The combined table here
# prefixes every event with its strand ("plus:SE:123"), matching the uid the
# sweep builds when STRANDED=1.
#
# Gene sets are disjoint between strands, so the two halves never collide at the
# gene level; the prefix only disambiguates the numeric event IDs.
#
# Usage: bash scripts/run_infer_junctions_gt_stranded.sh [n_shards]   (default 8)
set -euo pipefail
SIM=$HOME/GrASE_simulation
GTF=$SIM/ref/gencode.v28.annotation.gtf
RDA=$SIM/swimdown/simulate/data/simulate.rda
N=${1:-8}

for S in plus minus; do
  RD=$SIM/rMATS/stranded_$S/post
  OD=$SIM/results/sim_junction_gt.stranded_$S
  mkdir -p "$OD"
  echo "[$S] inferring junction GT ($N shards)..."
  for i in $(seq 1 "$N"); do
    Rscript $SIM/scripts/infer_rmats_junctions_gt.R \
      "$RD" "$GTF" "$RDA" "$OD" "$i" "$N" \
      &> "$SIM/results/sim_junction_gt.stranded_$S.log.shard$i" &
  done
  wait
done

COMB=$SIM/results/sim_junction_gt.stranded.txt
first=1
for S in plus minus; do
  OD=$SIM/results/sim_junction_gt.stranded_$S
  files=("$OD"/*.txt)
  [ -e "${files[0]}" ] || { echo "no per-gene files for $S" >&2; exit 1; }
  if [ $first -eq 1 ]; then
    # header + strand column
    printf '%s\tstrand_run\n' "$(head -1 "${files[0]}")" > "$COMB"
    first=0
  fi
  find "$OD" -maxdepth 1 -name '*.txt' -print0 \
    | xargs -0 tail -q -n +2 \
    | awk -v s="$S" 'BEGIN{FS=OFS="\t"} {print $0, s}' >> "$COMB"
done
echo "combined: $COMB ($(( $(wc -l < "$COMB") - 1 )) events)"
