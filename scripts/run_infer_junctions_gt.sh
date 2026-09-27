#!/bin/bash
# Structural junction-level ground truth, written one file per gene.
# Uses fromGTF.*.txt (full annotated universe) + MATS.JCEC.txt (tested flag).
#
# Per-gene files: results/sim_junction_gt/<GeneID>.txt  (resumable, parallel-safe)
# Combined table: results/sim_junction_gt.txt           (built in shell below)
set -euo pipefail

SIM_DIR=~/GrASE_simulation
RMATS_DIR=${SIM_DIR}/rMATS/rmats_post_group1_group2
GTF=${SIM_DIR}/ref/gencode.v28.annotation.gtf
RDA=${SIM_DIR}/swimdown/simulate/data/simulate.rda
OUT_DIR=${SIM_DIR}/results/sim_junction_gt
COMBINED=${SIM_DIR}/results/sim_junction_gt.txt
LOG=${SIM_DIR}/results/sim_junction_gt.log

# Number of parallel shards. Set to 1 for a single serial run.
N_SHARDS=${1:-8}

mkdir -p "${OUT_DIR}"

if [ "${N_SHARDS}" -le 1 ]; then
  Rscript ${SIM_DIR}/scripts/infer_rmats_junctions_gt.R \
    "${RMATS_DIR}" "${GTF}" "${RDA}" "${OUT_DIR}" &> "${LOG}"
else
  # launch N_SHARDS parallel jobs, each handling a 1/N_SHARDS slice of genes
  for i in $(seq 1 "${N_SHARDS}"); do
    Rscript ${SIM_DIR}/scripts/infer_rmats_junctions_gt.R \
      "${RMATS_DIR}" "${GTF}" "${RDA}" "${OUT_DIR}" "${i}" "${N_SHARDS}" \
      &> "${LOG}.shard${i}" &
  done
  wait
  cat "${LOG}".shard* > "${LOG}"
fi

# combine per-gene files into one table: keep one header, drop the rest.
# done in shell because R-side combining over thousands of files is too slow.
# use a glob array (not `ls | head`, which trips SIGPIPE under pipefail) and
# feed files via find+xargs to stay safely under the arg-length limit.
files=("${OUT_DIR}"/*.txt)
head -1 "${files[0]}" > "${COMBINED}"
find "${OUT_DIR}" -maxdepth 1 -name '*.txt' -print0 \
  | xargs -0 tail -q -n +2 >> "${COMBINED}"

echo "Per-gene files: ${OUT_DIR}/"
echo "Combined table: ${COMBINED} ($(wc -l < "${COMBINED}") lines)"
echo "Log: ${LOG}"
