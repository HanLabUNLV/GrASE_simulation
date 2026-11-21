#!/bin/bash
# Bipartitions of graphml.v28 with the partitions read off the paths' features
# (grase::bipartition_paths_by_feature), with the same 15-path cap as
# bipartition_paths, so the output can be compared with bipartition.v28.dedupfix.
# Runs the package from source (devtools::load_all on Rpkg), not the installed
# build. Writes to bipartition.v28.features; resumable: genes with output
# or a .running marker are skipped (delete stale markers after a crash).
#
# Usage: bash scripts/00_bipartition_v28_features.sh [cores, default 60] [max_path, default 15]
set -eu -o pipefail
BASE=/mnt/data1/home/mirahan/GrASE_simulation
PKG=/mnt/data1/home/mirahan/GrASE/Rpkg
CORES=${1:-60}; MAXP=${2:-15}
OUT=$BASE/bipartition.v28.features
mkdir -p "$OUT"
cd $BASE
echo "start: $(date)  cores $CORES  max_path $MAXP"
Rscript - "$BASE/graphml.v28" "$OUT" "$CORES" "$MAXP" "$PKG" <<'RS'
a <- commandArgs(TRUE); graphdir <- a[1]; outdir <- a[2]; cores <- as.integer(a[3])
maxp <- as.numeric(a[4]); pkg <- a[5]
suppressPackageStartupMessages({ devtools::load_all(pkg, quiet = TRUE); library(parallel) })
genes <- sub("\\.graphml$", "", list.files(graphdir, pattern = "\\.graphml$"))
todo <- genes[!file.exists(file.path(outdir, paste0(genes, ".bipartition.txt"))) &
              !file.exists(file.path(outdir, paste0(genes, ".running")))]
cat(length(genes), "graphs,", length(todo), "to split\n")
res <- mclapply(todo, function(gene) {
  mark <- file.path(outdir, paste0(gene, ".running")); file.create(mark)
  tryCatch({
    g <- igraph::read_graph(file.path(graphdir, paste0(gene, ".graphml")), format = "graphml")
    t <- system.time(suppressMessages(bipartition_paths_by_feature(gene, g, outdir, max_path = maxp)))[["elapsed"]]
    unlink(mark); sprintf("%s\t%.1f", gene, t)
  }, error = function(e) {
    cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%F %T"), gene, conditionMessage(e)),
        file = file.path(outdir, "errors.log"), append = TRUE)
    unlink(mark); NA_character_
  })
}, mc.cores = cores, mc.preschedule = FALSE)
res <- unlist(res)
cat("errors:", sum(is.na(res)), "\n")
cat(c("gene\tseconds", res[!is.na(res)]), sep = "\n", file = file.path(outdir, "split_times.tsv"), append = FALSE)
RS
echo "end: $(date)"
echo "deduped: $(ls "$OUT" | grep -c '\.bipartition\.txt$')  pre-dedup: $(ls "$OUT" | grep -c '\.bipartition\.all\.txt$')"
