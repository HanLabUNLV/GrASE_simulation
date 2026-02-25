#!/usr/bin/env Rscript

# scripts/evaluate_tools_vs_gt.R
#
# Compare rMATS and Saturn exonic-part detection accuracy against simulation
# ground truth.  Uses the same GT as evaluate_bipartition_test.R so results
# are directly comparable to GrASE.
#
# rMATS:
#   Event files : <rmats_dir>/{SE,A3SS,A5SS,RI,MXE}.MATS.JCEC.txt
#   Mapping     : <map_dir>/combined.fromGTF.{SE,A3SS,A5SS,RI}.txt
#                 columns: ID (rMATS event ID), GeneID, DexseqFragment
#                          (comma-sep exonic parts of the alternative exon)
#   Significance: FDR column
#
# Saturn:
#   File        : <saturn_file>
#                 ids column format: "ENSG...:Exxx"  (gene:exonic_part)
#   Significance: empirical_FDR column (falls back to regular_FDR)
#
# Ground truth (gt_dir, from infer_diff_exons_simulation.R):
#   gene.exonic_parts_fc.txt   columns: gene, exonic_part, group
#                                       group in {changed, constant, negative}
#
# Usage:
#   Rscript evaluate_tools_vs_gt.R \
#     <gt_dir> <rmats_dir> <map_dir> <saturn_file> <out_dir> <simulate_rda>
#
# Outputs (written to out_dir):
#   {rmats,saturn}_per_gene_padj<thr>.txt            per-gene metrics (full GT)
#   {rmats,saturn}_restricted_per_gene_padj<thr>.txt per-gene metrics (GT restricted to each tool's testable exons)
#   {rmats,saturn}_summary_by_simtype.txt            aggregated by sim_type x threshold (full GT)
#   {rmats,saturn}_restricted_summary_by_simtype.txt same, restricted to each tool's testable exons
#   comparison_summary.txt                           side-by-side comparison (rMATS, rMATS_restricted, Saturn, Saturn_restricted)

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  cat(paste(
    "Usage: Rscript evaluate_tools_vs_gt.R",
    "  <gt_dir>       directory with gene.exonic_parts_fc.txt files",
    "  <rmats_dir>    directory with *.MATS.JCEC.txt files",
    "  <map_dir>      directory with combined.fromGTF.*.txt mapping files",
    "  <saturn_file>  Saturn results file",
    "  <out_dir>      output directory",
    "  <simulate_rda> path to simulate.rda (for DGE/DTE/DTU labels)",
    sep = "\n"
  ), "\n")
  quit(save = "no", status = 1)
}

gt_dir      <- args[1]
rmats_dir   <- args[2]
map_dir     <- args[3]
saturn_file <- args[4]
out_dir     <- args[5]
sim_rda     <- args[6]

padj_thresholds <- c(0.01, 0.05, 0.1, 0.2)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# helpers 
f1_safe <- function(p, r) {
  if (is.na(p) || is.na(r) || (p + r) == 0) NA_real_
  else 2 * p * r / (p + r)
}

# Compute TP/FP/FN/TN given detected exonic part set and GT labels data frame.
# testable:      if non-NULL (restricted eval), restricts the evaluation universe.
# restrict_pos:  if TRUE (restricted eval), restrict both gt_pos and gt_neg to
#                testable (coverage-restricted recall; FN excludes the coverage
#                gap). Full eval passes no testable at all (universe = all GT).
expart_metrics <- function(det_exparts, gt_gene_df,
                           testable = NULL, restrict_pos = TRUE) {
  gt_pos <- unique(gt_gene_df$exonic_part[gt_gene_df$group %in% c("changed", "constant")])
  gt_neg <- unique(gt_gene_df$exonic_part[gt_gene_df$group == "negative"])
  if (!is.null(testable)) {
    if (restrict_pos) gt_pos <- intersect(gt_pos, testable)
    # Always restrict gt_neg: only count TN for exons the tool explicitly tested
    gt_neg <- intersect(gt_neg, testable)
  }
  det    <- unique(det_exparts[!is.na(det_exparts) & nchar(det_exparts) > 0])
  TP <- length(intersect(det, gt_pos))
  FP <- length(intersect(det, gt_neg))
  FN <- length(setdiff(gt_pos, det))
  TN <- length(setdiff(gt_neg, det))
  prec   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  list(
    n_gt_pos = length(gt_pos), n_gt_neg = length(gt_neg), n_detected = length(det),
    TP = TP, FP = FP, FN = FN, TN = TN,
    precision = round(prec,           4),
    recall    = round(recall,         4),
    f1        = round(f1_safe(prec, recall), 4)
  )
}

build_summary <- function(per_gene_df) {
  by_type <- per_gene_df %>%
    group_by(sim_type, padj_thr) %>%
    summarise(
      n_genes          = n(),
      total_TP         = sum(TP,  na.rm = TRUE),
      total_FP         = sum(FP,  na.rm = TRUE),
      total_FN         = sum(FN,  na.rm = TRUE),
      total_TN         = sum(TN,  na.rm = TRUE),
      micro_precision  = round(sum(TP) / max(sum(TP) + sum(FP), 1), 4),
      micro_recall     = round(sum(TP) / max(sum(TP) + sum(FN), 1), 4),
      micro_f1         = round(f1_safe(
                           sum(TP) / max(sum(TP) + sum(FP), 1),
                           sum(TP) / max(sum(TP) + sum(FN), 1)), 4),
      macro_precision  = round(mean(precision, na.rm = TRUE), 4),
      macro_recall     = round(mean(recall,    na.rm = TRUE), 4),
      macro_f1         = round(mean(f1,        na.rm = TRUE), 4),
      .groups = "drop"
    )
  all_row <- by_type %>%
    group_by(padj_thr) %>%
    summarise(
      sim_type        = "ALL",
      n_genes         = sum(n_genes),
      total_TP        = sum(total_TP),
      total_FP        = sum(total_FP),
      total_FN        = sum(total_FN),
      total_TN        = sum(total_TN),
      micro_precision = round(total_TP / max(total_TP + total_FP, 1), 4),
      micro_recall    = round(total_TP / max(total_TP + total_FN, 1), 4),
      micro_f1        = round(f1_safe(micro_precision, micro_recall), 4),
      macro_precision = NA_real_,
      macro_recall    = NA_real_,
      macro_f1        = NA_real_,
      .groups = "drop"
    )
  bind_rows(by_type, all_row) %>% arrange(sim_type, padj_thr)
}

# load sim gene type labels
 
cat("Loading simulate.rda...\n")
dge_genes <- dte_genes <- dtu_genes <- character(0)
if (file.exists(sim_rda)) {
  load(sim_rda)
  if (exists("dge.genes")) dge_genes <- dge.genes
  if (exists("dte.genes")) dte_genes <- dte.genes
  if (exists("dtu.genes")) dtu_genes <- dtu.genes
  cat(sprintf("  DGE: %d  DTE: %d  DTU: %d genes\n",
              length(dge_genes), length(dte_genes), length(dtu_genes)))
} else {
  warning("simulate.rda not found; sim_type will be 'Unknown'")
}

get_sim_type <- function(gene) {
  dplyr::case_when(
    gene %in% dge_genes ~ "DGE",
    gene %in% dte_genes ~ "DTE",
    gene %in% dtu_genes ~ "DTU",
    TRUE ~ "Background"
  )
}

# load ground truth 

cat("Loading ground truth...\n")
gt_fc_files      <- list.files(gt_dir, pattern = "\\.exonic_parts_fc\\.txt$",
                                full.names = TRUE)
gt_genes_avail   <- sub("\\.exonic_parts_fc\\.txt$", "", basename(gt_fc_files))

gt_all <- lapply(gt_fc_files, function(f) {
  tryCatch(
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               quote = "", comment.char = "",
               colClasses = c(exonic_part = "character",
                              fold_change = "character",
                              source = "character", sink = "character")),
    error = function(e) NULL
  )
})
gt_all <- bind_rows(gt_all[!sapply(gt_all, is.null)])
gt_all$sim_type <- get_sim_type(gt_all$gene)
gt_by_gene <- split(gt_all, gt_all$gene)
cat(sprintf("  %d exonic part records across %d GT genes\n",
            nrow(gt_all), length(gt_by_gene)))

# build rMATS exonic-part call table 

cat("\nBuilding rMATS exonic-part call table...\n")

# Event types that have mapping files (MXE omitted no mapping available)
rmats_etypes <- c("SE", "A3SS", "A5SS", "RI")

rmats_rows <- list()

for (etype in rmats_etypes) {
  map_path   <- file.path(map_dir, paste0("combined.fromGTF.", etype, ".txt"))
  rmats_path <- file.path(rmats_dir, paste0(etype, ".MATS.JCEC.txt"))

  if (!file.exists(map_path)) {
    message(sprintf("  Skipping %s: mapping file not found", etype))
    next
  }
  if (!file.exists(rmats_path)) {
    message(sprintf("  Skipping %s: rMATS result file not found", etype))
    next
  }

  map_df <- read.table(map_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")
  # Only rows with a real exonic-part mapping
  map_df <- map_df[!is.na(map_df$DexseqFragment) &
                     map_df$DexseqFragment != "NA" &
                     map_df$DexseqFragment != "", ]
  if (nrow(map_df) == 0) {
    message(sprintf("  Skipping %s: no exonic-part mappings in file", etype))
    next
  }

  rmats_df <- read.table(rmats_path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "")
  # Deduplicate the duplicate ID column (rMATS has ID twice)
  names(rmats_df) <- make.unique(names(rmats_df))
  rmats_df$GeneID <- gsub('"', '', rmats_df$GeneID)  # strip surrounding quotes

  # Keep only FDR and ID for joining
  rmats_sig <- rmats_df[, c("ID", "FDR")]
  rmats_sig$ID <- as.character(rmats_sig$ID)
  map_df$ID    <- as.character(map_df$ID)

  merged <- merge(map_df[, c("ID", "GeneID", "DexseqFragment")],
                  rmats_sig, by = "ID", all.x = FALSE)
  if (nrow(merged) == 0) next

  # Expand comma-separated DexseqFragment into one row per exonic part
  expanded <- lapply(seq_len(nrow(merged)), function(i) {
    parts <- trimws(unlist(strsplit(merged$DexseqFragment[i], ",")))
    data.frame(
      gene        = merged$GeneID[i],
      exonic_part = parts,
      FDR         = merged$FDR[i],
      event_type  = etype,
      stringsAsFactors = FALSE
    )
  })
  rmats_rows[[etype]] <- bind_rows(expanded)
  message(sprintf("  %s: %d mapped exonic-part rows", etype, nrow(rmats_rows[[etype]])))
}

rmats_exparts <- bind_rows(rmats_rows)
# If the same exonic part is covered by multiple events, keep the minimum FDR
rmats_exparts <- rmats_exparts %>%
  group_by(gene, exonic_part) %>%
  summarise(FDR = min(FDR, na.rm = TRUE), .groups = "drop")

cat(sprintf("  rMATS: %d unique gene x exonic_part entries\n", nrow(rmats_exparts)))
rmats_by_gene <- split(rmats_exparts, rmats_exparts$gene)

# Build per-gene set of all exonic parts rMATS can test (regardless of significance)
# Used for coverage-restricted evaluation so GT denominator is limited to rMATS-testable exons
rmats_testable_by_gene <- lapply(
  split(rmats_exparts$exonic_part, rmats_exparts$gene),
  unique
)
cat(sprintf("  rMATS: %d genes with testable exonic parts\n", length(rmats_testable_by_gene)))

#  build Saturn exonic-part call table 

cat("Loading Saturn results...\n")
saturn_df <- read.table(saturn_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, quote = "", comment.char = "")

# ids format: "ENSG00000000419.12:E001"
saturn_df$gene        <- sub(":.*$", "", saturn_df$ids)
saturn_df$exonic_part <- sub("^.*:", "", saturn_df$ids)

# Use empirical_FDR; fall back to regular_FDR if absent or all-NA
fdr_col <- if ("empirical_FDR" %in% names(saturn_df) &&
               any(!is.na(saturn_df$empirical_FDR))) {
  "empirical_FDR"
} else {
  message("  Using regular_FDR (empirical_FDR missing or all-NA)")
  "regular_FDR"
}
saturn_df$FDR <- saturn_df[[fdr_col]]
cat(sprintf("  Saturn: %d exonic-part rows (FDR column: %s)\n",
            nrow(saturn_df), fdr_col))
saturn_by_gene <- split(saturn_df[, c("gene", "exonic_part", "FDR")],
                        saturn_df$gene)

# Build per-gene set of all exonic parts Saturn tested (all rows in the file,
# regardless of significance; Saturn pre-filters low-count/single-exon genes
# and exons before running tests, so every row in the file was actually tested)
saturn_testable_by_gene <- lapply(
  split(saturn_df$exonic_part, saturn_df$gene),
  unique
)
cat(sprintf("  Saturn: %d genes with testable exonic parts\n", length(saturn_testable_by_gene)))

# Also build a regular_FDR variant for Saturn (same testable universe)
if ("regular_FDR" %in% names(saturn_df)) {
  saturn_regfdr_df <- saturn_df
  saturn_regfdr_df$FDR <- saturn_regfdr_df$regular_FDR
  saturn_regfdr_by_gene <- split(
    saturn_regfdr_df[, c("gene", "exonic_part", "FDR")],
    saturn_regfdr_df$gene
  )
  cat(sprintf("  Saturn regular_FDR: %d genes\n",
              length(saturn_regfdr_by_gene)))
} else {
  saturn_regfdr_by_gene <- NULL
  cat("  Saturn regular_FDR: column not found, skipping\n")
}

#  evaluate both tools across padj thresholds 

eval_genes <- unique(gt_all$gene)

evaluate_tool <- function(tool_by_gene, tool_name,
                          testable_by_gene = NULL, restrict_pos = TRUE) {
  cat(sprintf("\nEvaluating %s...\n", tool_name))
  rows_all <- list()

  for (thr in padj_thresholds) {
    rows_thr <- lapply(eval_genes, function(gene) {
      gt_df <- gt_by_gene[[gene]]
      if (is.null(gt_df)) return(NULL)

      testable <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]] else NULL
      # Skip genes where the tool has no testable exonic parts
      if (!is.null(testable) && length(testable) == 0) return(NULL)

      tool_gene <- tool_by_gene[[gene]]
      if (!is.null(tool_gene)) {
        det_exparts <- tool_gene$exonic_part[
          !is.na(tool_gene$FDR) & tool_gene$FDR < thr
        ]
      } else {
        det_exparts <- character(0)
      }

      m <- expart_metrics(det_exparts, gt_df,
                          testable = testable, restrict_pos = restrict_pos)
      data.frame(
        gene       = gene,
        sim_type   = gt_df$sim_type[1],
        tool       = tool_name,
        padj_thr   = thr,
        n_gt_pos   = m$n_gt_pos,
        n_gt_neg   = m$n_gt_neg,
        n_detected = m$n_detected,
        TP         = m$TP, FP = m$FP, FN = m$FN, TN = m$TN,
        precision  = m$precision,
        recall     = m$recall,
        f1         = m$f1,
        stringsAsFactors = FALSE
      )
    })
    rows_all <- c(rows_all, rows_thr)
  }
  bind_rows(rows_all[!sapply(rows_all, is.null)])
}

# Full evaluation (framework i - total space): universe = all GT exons.
# TN = all GT-negative exons not in each tool's significant set (including
# untested ones; non-tested = implicit negative call).
# No testable restriction; FN includes the coverage gap.
rmats_per_gene  <- evaluate_tool(rmats_by_gene,  "rMATS")
saturn_per_gene <- evaluate_tool(saturn_by_gene, "Saturn")

# Coverage-restricted: both GT-positive and GT-negative restricted to each
# tool's testable universe. FN = tested but not significant (no coverage gap).
rmats_restricted_per_gene <- evaluate_tool(
  rmats_by_gene, "rMATS_restricted",
  testable_by_gene = rmats_testable_by_gene,
  restrict_pos = TRUE
)
saturn_restricted_per_gene <- evaluate_tool(
  saturn_by_gene, "Saturn_restricted",
  testable_by_gene = saturn_testable_by_gene,
  restrict_pos = TRUE
)

# Saturn regular_FDR evaluations
if (!is.null(saturn_regfdr_by_gene)) {
  saturn_regfdr_per_gene <- evaluate_tool(
    saturn_regfdr_by_gene, "Saturn_regularFDR")
  saturn_regfdr_restricted_per_gene <- evaluate_tool(
    saturn_regfdr_by_gene, "Saturn_regularFDR_restricted",
    testable_by_gene = saturn_testable_by_gene,
    restrict_pos = TRUE
  )
}

#  write per-gene results 

for (thr in padj_thresholds) {
  write.table(
    rmats_per_gene[rmats_per_gene$padj_thr == thr, ],
    file.path(out_dir, sprintf("rmats_per_gene_padj%.2f.txt", thr)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    saturn_per_gene[saturn_per_gene$padj_thr == thr, ],
    file.path(out_dir, sprintf("saturn_per_gene_padj%.2f.txt", thr)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    rmats_restricted_per_gene[rmats_restricted_per_gene$padj_thr == thr, ],
    file.path(out_dir, sprintf("rmats_restricted_per_gene_padj%.2f.txt", thr)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    saturn_restricted_per_gene[saturn_restricted_per_gene$padj_thr == thr, ],
    file.path(out_dir, sprintf("saturn_restricted_per_gene_padj%.2f.txt", thr)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  if (!is.null(saturn_regfdr_by_gene)) {
    write.table(
      saturn_regfdr_per_gene[saturn_regfdr_per_gene$padj_thr == thr, ],
      file.path(out_dir, sprintf("saturn_regularFDR_per_gene_padj%.2f.txt", thr)),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    write.table(
      saturn_regfdr_restricted_per_gene[
        saturn_regfdr_restricted_per_gene$padj_thr == thr, ],
      file.path(out_dir,
                sprintf("saturn_regularFDR_restricted_per_gene_padj%.2f.txt", thr)),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
}

#  summaries 

rmats_summary              <- build_summary(rmats_per_gene)
saturn_summary             <- build_summary(saturn_per_gene)
rmats_restricted_summary   <- build_summary(rmats_restricted_per_gene)
saturn_restricted_summary  <- build_summary(saturn_restricted_per_gene)

write.table(rmats_summary,
            file.path(out_dir, "rmats_summary_by_simtype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(saturn_summary,
            file.path(out_dir, "saturn_summary_by_simtype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rmats_restricted_summary,
            file.path(out_dir, "rmats_restricted_summary_by_simtype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(saturn_restricted_summary,
            file.path(out_dir, "saturn_restricted_summary_by_simtype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(saturn_regfdr_by_gene)) {
  saturn_regfdr_summary            <- build_summary(saturn_regfdr_per_gene)
  saturn_regfdr_restricted_summary <- build_summary(saturn_regfdr_restricted_per_gene)
  write.table(saturn_regfdr_summary,
              file.path(out_dir, "saturn_regularFDR_summary_by_simtype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(saturn_regfdr_restricted_summary,
              file.path(out_dir, "saturn_regularFDR_restricted_summary_by_simtype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# print results 

cat("\n=== rMATS Summary (full GT) ===\n")
print(as.data.frame(rmats_summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

cat("\n=== rMATS Summary (restricted to rMATS-testable exonic parts) ===\n")
print(as.data.frame(rmats_restricted_summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

cat("\n=== Saturn Summary (full GT) ===\n")
print(as.data.frame(saturn_summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

cat("\n=== Saturn Summary (restricted to Saturn-testable exonic parts) ===\n")
print(as.data.frame(saturn_restricted_summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

if (!is.null(saturn_regfdr_by_gene)) {
  cat("\n=== Saturn Summary regular_FDR (full GT) ===\n")
  print(as.data.frame(saturn_regfdr_summary %>%
    select(sim_type, padj_thr, n_genes,
           micro_precision, micro_recall, micro_f1,
           total_TP, total_FP, total_FN)))

  cat("\n=== Saturn Summary regular_FDR (restricted) ===\n")
  print(as.data.frame(saturn_regfdr_restricted_summary %>%
    select(sim_type, padj_thr, n_genes,
           micro_precision, micro_recall, micro_f1,
           total_TP, total_FP, total_FN)))
}

cat(sprintf("\nResults written to: %s\n", out_dir))
