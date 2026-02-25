#!/usr/bin/env Rscript

# scripts/infer_diff_exons_gt.R
#
# Infer ground-truth positive/negative exonic bins from simulate.rda + DEXSeq GFF.
#
# Iterates over genes in txdf (the simulation gene/transcript universe).
#
# Two logic branches:
#   Background / DGE  all bins NEGATIVE  (no differential transcript usage)
#   DTU / DTE         bin POSITIVE if it overlaps any transcript whose
#                        fold_change != 1; NEGATIVE otherwise
#
# Usage:
#   Rscript infer_diff_exons_gt.R <gff_dir> <simulate_rda> <out_dir>

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript infer_diff_exons_gt.R <gff_dir> <simulate_rda> <out_dir>\n")
  quit(save = "no", status = 1)
}
gff_dir <- args[1]
sim_rda  <- args[2]
out_dir  <- args[3]

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load simulation data 

cat("Loading simulate.rda...\n")
load(sim_rda)
# Required objects:
#   txdf         data.frame  GENEID, TXNAME, ntx  simulation gene/tx universe
#   fold_changes matrix      rows = transcript IDs, cols = conditions
#   dge.genes / dte.genes / dtu.genes  character vectors

dge_genes <- if (exists("dge.genes")) dge.genes else character(0)
dte_genes <- if (exists("dte.genes")) dte.genes else character(0)
dtu_genes <- if (exists("dtu.genes")) dtu.genes else character(0)

cat(sprintf("  DGE: %d  DTE: %d  DTU: %d\n",
            length(dge_genes), length(dte_genes), length(dtu_genes)))
cat(sprintf("  fold_changes: %d transcripts x %d conditions\n",
            nrow(fold_changes), ncol(fold_changes)))
cat(sprintf("  txdf: %d transcripts across %d genes\n",
            nrow(txdf), length(unique(txdf$GENEID))))

# Changed transcripts: any fold_change != 1.
# Transcripts absent from fold_changes are treated as unchanged.
changed_tx <- rownames(fold_changes)[
  apply(fold_changes, 1, function(x) any(x != 1))
]
cat(sprintf("  Changed transcripts (any FC != 1): %d\n", length(changed_tx)))

# Build gene transcript list from txdf
txs_by_gene <- split(txdf$TXNAME, txdf$GENEID)
all_genes   <- names(txs_by_gene)
cat(sprintf("  Genes to process: %d\n", length(all_genes)))

# Helpers 

# Format fold change values for a set of transcripts as "tx1:v1/v2+tx2:v1/v2"
# (values = one entry per column of fold_changes; "/" separates conditions,
#  "+" separates transcripts)
format_fc <- function(txs, fc_matrix) {
  txs <- txs[txs %in% rownames(fc_matrix)]
  if (length(txs) == 0) return(NA_character_)
  paste(vapply(txs, function(tx) {
    paste0(tx, ":", paste(fc_matrix[tx, ], collapse = "/"))
  }, character(1)), collapse = "+")
}

# GFF parser 

parse_dexseq_gff <- function(gff_path) {
  lines    <- readLines(gff_path, warn = FALSE)
  ep_lines <- lines[grepl("\texonic_part\t", lines, fixed = TRUE)]
  if (length(ep_lines) == 0) return(NULL)
  attr_field <- sub("^([^\t]*\t){8}", "", ep_lines)
  bin_num    <- sub('.*exonic_part_number "([^"]+)".*', "\\1", attr_field)
  tx_raw     <- sub('.*transcripts "([^"]+)".*',        "\\1", attr_field)
  data.frame(
    bin_id  = paste0("E", bin_num),
    tx_str  = tx_raw,
    tx_list = I(strsplit(tx_raw, "+", fixed = TRUE)),
    stringsAsFactors = FALSE
  )
}

# Process genes 

cat("\nProcessing genes...\n")
n_written <- 0L
n_skipped <- 0L

for (gene_name in all_genes) {

  sim_type <- dplyr::case_when(
    gene_name %in% dge_genes ~ "DGE",
    gene_name %in% dte_genes ~ "DTE",
    gene_name %in% dtu_genes ~ "DTU",
    TRUE                      ~ "Background"
  )

  gff_path <- file.path(gff_dir, paste0(gene_name, ".dexseq.gff"))
  if (!file.exists(gff_path)) { n_skipped <- n_skipped + 1L; next }

  bins <- parse_dexseq_gff(gff_path)
  if (is.null(bins) || nrow(bins) == 0) { n_skipped <- n_skipped + 1L; next }

  if (sim_type %in% c("Background", "DGE")) {

    # All bins negative: no differential transcript usage.
    # For DGE, every transcript in the bin carries the gene-level FC; report it.
    fc_vals <- vapply(bins$tx_list, function(txs) {
      format_fc(txs, fold_changes)
    }, character(1))

    out_df <- data.frame(
      gene        = gene_name,
      sim_type    = sim_type,
      exonic_part = bins$bin_id,
      fold_change = fc_vals,
      group       = "negative",
      transcripts = bins$tx_str,
      changed_tx  = NA_character_,
      alt_tx      = NA_character_,
      source      = NA_character_,
      sink        = NA_character_,
      stringsAsFactors = FALSE
    )

  } else {

    # DTU / DTE: bin positive if it overlaps any changed transcript for this gene
    gene_changed <- intersect(txs_by_gene[[gene_name]], changed_tx)

    rows <- lapply(seq_len(nrow(bins)), function(i) {
      hit    <- intersect(bins$tx_list[[i]], gene_changed)
      alt    <- setdiff(gene_changed, bins$tx_list[[i]])
      # DTE: positive if any changed tx is in the bin
      # DTU: positive only if some changed tx are in the bin AND some are outside
      #      (bin must separate changed transcripts, not contain all of them)
      is_pos <- if (sim_type == "DTU") length(hit) > 0 && length(alt) > 0
                else                   length(hit) > 0
      # Report FC for any changed transcripts present in the bin (hit),
      # regardless of whether the bin is positive or negative.
      fc_str <- format_fc(hit, fold_changes)
      data.frame(
        gene        = gene_name,
        sim_type    = sim_type,
        exonic_part = bins$bin_id[i],
        fold_change = fc_str,
        group       = if (is_pos) "changed" else "negative",
        transcripts = bins$tx_str[i],
        changed_tx  = if (is_pos) paste(hit, collapse = "+") else NA_character_,
        alt_tx      = if (is_pos && length(alt) > 0) paste(alt, collapse = "+") else NA_character_,
        source      = NA_character_,
        sink        = NA_character_,
        stringsAsFactors = FALSE
      )
    })
    out_df <- bind_rows(rows)

  }

  write.table(out_df,
              file.path(out_dir, paste0(gene_name, ".exonic_parts_fc.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  n_written <- n_written + 1L
}

cat(sprintf("Done. Written %d  |  Skipped %d (no GFF)   %s\n",
            n_written, n_skipped, out_dir))
