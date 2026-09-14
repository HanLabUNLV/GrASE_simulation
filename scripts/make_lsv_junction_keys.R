#!/usr/bin/env Rscript
#
# scripts/make_lsv_junction_keys.R
#
# Build the keyed per-junction table for MAJIQ LSVs. This is STRUCTURE, not
# ground truth: it explodes each LSV in the voila tsv into one row per junction,
# keeping lsv_id / gene / junction coordinates / K / MAJIQ's own per-junction
# statistics, plus the transcripts using each junction.
#
# Split out of decision_space_fpr.R, which used to produce it as a side effect.
# Ground-truth generation depending on an evaluation script is how that file
# ended up carrying count-scale columns that were mistaken for truth.
#
# Ground truth is added downstream by scripts/infer_lsv_junction_gt.R.
#
# The `users` column stores the TRANSCRIPT IDENTITIES using each junction, not
# just a count. Evaluators need the identities to build implicated-transcript
# sets for transcript- and gene-level scoring; storing only n_users forced
# pr_curves_three_levels_gtrule.R to rebuild the whole junction->transcript
# mapping from the voila tsv on every run (~30 min), which is exactly the
# duplication splitting this script out was meant to remove.
#
# Usage: Rscript scripts/make_lsv_junction_keys.R
# Out:   results/gt/lsv_junction_keys.txt

suppressPackageStartupMessages({ library(dplyr); library(rtracklayer) })
BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
OUT  <- file.path(BASE, "results/gt"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
load(file.path(BASE, "swimdown/simulate/data/simulate.rda"))
get_st <- function(g) ifelse(g %in% dte.genes, "DTE", ifelse(g %in% dtu.genes, "DTU",
                      ifelse(g %in% dge.genes, "DGE", "Background")))
parse_parts <- function(x) { v <- strsplit(as.character(x), ",")[[1]]
  trimws(v[nchar(trimws(v)) > 0 & trimws(v) != "NA"]) }

cat("Building tx_junc from GTF...\n")
ex <- import(file.path(BASE, "ref/gencode.v28.annotation.gtf"),
             feature.type = "exon", colnames = "transcript_id")
o <- order(ex$transcript_id, start(ex)); tid <- ex$transcript_id[o]
st <- start(ex)[o]; en <- end(ex)[o]; n <- length(tid); sm <- tid[-n] == tid[-1]
tx_junc <- split(paste(en[-n], st[-1], sep = "-")[sm], tid[-n][sm])

## STRANDED=1 enumerates the universe from the origin-split MAJIQ build. This
## matters because the universe is DATA-dependent (it is whatever LSVs MAJIQ
## found), unlike the GT labels, which are design-based. Scoring the stranded
## run against the old build's 50,489 LSVs charges it 55 GT-positive LSVs it
## never tested.
STRANDED <- nzchar(Sys.getenv("STRANDED"))
MAJIQ_TSV <- if (STRANDED) "majiq/majiq_deltapsi.stranded.thr0.20.tsv" else "majiq/majiq_deltapsi.thr0.20.tsv"
KEYS_OUT  <- if (STRANDED) "lsv_junction_keys.stranded.txt" else "lsv_junction_keys.txt"
tsv <- read.table(file.path(BASE, MAJIQ_TSV), header = TRUE,
                  sep = "\t", quote = "", comment.char = "#", stringsAsFactors = FALSE)
tsv$sim_type <- get_st(tsv$gene_id)
cat(sprintf("LSVs in voila tsv: %d\n", nrow(tsv)))

n_mis_p <- n_mis_d <- n_mis_K <- 0L
rows <- vector("list", nrow(tsv))
for (i in seq_len(nrow(tsv))) {
  lj <- parse_parts(gsub(";", ",", tsv$junctions_coords[i]))
  pj <- suppressWarnings(as.numeric(strsplit(tsv$probability_changing[i], "[;,]")[[1]]))
  dj <- suppressWarnings(as.numeric(strsplit(tsv$mean_dpsi_per_lsv_junction[i], "[;,]")[[1]]))
  K <- length(lj); if (K == 0) next
  if (length(pj) != K) { pj <- rep(NA_real_, K); n_mis_p <- n_mis_p + 1L }
  if (length(dj) != K) { dj <- rep(NA_real_, K); n_mis_d <- n_mis_d + 1L }
  if (!is.na(tsv$num_junctions[i]) && tsv$num_junctions[i] != K) n_mis_K <- n_mis_K + 1L
  gtx <- intersect(txdf$TXNAME[txdf$GENEID == tsv$gene_id[i]], names(tx_junc))
  nu <- rep(NA_integer_, K); us <- rep(NA_character_, K)
  if (length(gtx)) {
    ul <- lapply(lj, function(x) gtx[vapply(tx_junc[gtx], function(v) x %in% v, logical(1))])
    nu <- vapply(ul, length, integer(1))
    us <- vapply(ul, function(v) paste(v, collapse = ","), character(1))
  }
  rows[[i]] <- data.frame(lsv_id = tsv$lsv_id[i], gene = tsv$gene_id[i],
                          sim_type = tsv$sim_type[i], lsv_type = tsv$lsv_type[i],
                          junction = lj, j_index = seq_len(K), K = K,
                          prob = pj, dpsi_est = dj, n_users = nu, users = us,
                          sig_p95 = !is.na(pj) & pj >= 0.95,
                          stringsAsFactors = FALSE)
  if (i %% 10000 == 0) cat(sprintf("  ...%d/%d LSVs\n", i, nrow(tsv)))
}
j <- bind_rows(rows[!sapply(rows, is.null)])
cat(sprintf("\nfield-length mismatches vs K: probability_changing %d, mean_dpsi %d, num_junctions %d\n",
            n_mis_p, n_mis_d, n_mis_K))
f <- file.path(OUT, KEYS_OUT)
write.table(j, f, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d junctions from %d LSVs -> %s\n", nrow(j), length(unique(j$lsv_id)), f))
cat("junctions per LSV (K):\n"); print(summary(j$K))
cat("by sim_type:\n"); print(table(j$sim_type))
