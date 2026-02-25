#!/usr/bin/env Rscript

# scripts/evaluate_bipartition_test.R
#
# Evaluates GrASE bipartition statistical test results against simulation
# ground truth.
#
# Input test results file columns (tab-separated):
#   gene, event, LRT, p.value, model, phi, effect_size, padj,
#   source, sink, ref_ex_part, setdiff1, setdiff2,
#   transcripts1, transcripts2, path1, path2,
#   ref_mean, diff1_mean, diff2_mean, diff_mean, which, setdiff
#
# Ground truth (gt_dir, from infer_diff_exons_gt.R):
#   gene.exonic_parts_fc.txt   exonic_part, fold_change, group, transcripts,
#                               changed_tx, alt_tx, source, sink
#
# Exonic part detection:
#   Detected positive = exonic parts in `setdiff` column with padj < threshold
#   GT positive       = exonic parts labelled "changed" or "constant"
#   GT negative       = exonic parts labelled "negative"
#   Reports TP/FP/FN/TN, precision, recall, F1 at multiple padj thresholds
#
#   Two evaluations:
#     Full (total space)   : universe = all GT genes; untested = implicit neg
#     Restricted           : universe = only GrASE-testable exons (setdiff)
#
# Usage:
#   Rscript evaluate_bipartition_test.R <test_results> <gt_dir> <out_dir> <simulate_rda>
#
#   <test_results> may be a single file path or a comma-separated list of paths.
#   When multiple files are given they are combined: the minimum
#   padj across files is used per (gene, setdiff).
#
# Outputs:
#   grase_per_gene_padj<thr>.txt                  per-gene metrics (full GT)
#   grase_restricted_per_gene_padj<thr>.txt       per-gene metrics (restricted)
#   grase_summary_by_simtype.txt                  P/R/F1 by sim_type (full GT)
#   grase_restricted_summary_by_simtype.txt       same, restricted

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
n_cores <- 30

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript evaluate_bipartition_test.R <test_results> <gt_dir> <out_dir> <simulate_rda>\n")
  quit(save = "no", status = 1)
}

test_files  <- trimws(unlist(strsplit(args[1], ",")))
gt_dir      <- args[2]
out_dir     <- args[3]
sim_rda     <- args[4]

padj_thresholds <- c(0.01, 0.05, 0.1, 0.2)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#  helpers 

parse_exparts <- function(x) {
  if (is.na(x) || trimws(x) == "" || x == "NA") return(character(0))
  trimws(unlist(strsplit(x, ",")))
}

f1_safe <- function(p, r) {
  if (is.na(p) || is.na(r) || (p + r) == 0) NA_real_
  else 2 * p * r / (p + r)
}

#  load test results 

cat(sprintf("loading %d test file(s)...\n", length(test_files)))
tests_list <- lapply(test_files, function(f) {
  cat(sprintf("  %s\n", basename(f)))
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", comment.char = "", na.strings = c("NA", ""),
             colClasses = c(source = "character", sink = "character"))
})
tests <- bind_rows(tests_list)
cat(sprintf("  %d total rows, %d unique genes across all files\n",
            nrow(tests), length(unique(tests$gene))))

# when multiple files are combined, take the minimum padj per (gene, setdiff)
if (length(test_files) > 1) {
  tests_eval <- tests %>%
    group_by(gene, setdiff) %>%
    summarise(padj = min(padj, na.rm = TRUE), .groups = "drop")
} else {
  tests_eval <- tests %>% select(gene, setdiff, padj)
}

# build per-gene set of exonic parts grase explicitly tested (setdiff only)
# setdiff may be a comma-separated list of exonic parts; split before deduplicating
grase_testable_by_gene <- lapply(
  split(tests_eval$setdiff, tests_eval$gene),
  function(x) unique(unlist(lapply(x[!is.na(x) & nchar(x) > 0], parse_exparts)))
)
cat(sprintf("  grase: %d genes with testable exonic parts (setdiff only)\n",
            length(grase_testable_by_gene)))

#  load simulation gene type labels 

cat("loading simulate.rda...\n")
dge_genes <- dte_genes <- dtu_genes <- character(0)
if (file.exists(sim_rda)) {
  load(sim_rda)
  if (exists("dge.genes")) dge_genes <- dge.genes
  if (exists("dte.genes")) dte_genes <- dte.genes
  if (exists("dtu.genes")) dtu_genes <- dtu.genes
  cat(sprintf("  DGE: %d  DTE: %d  DTU: %d genes\n",
              length(dge_genes), length(dte_genes), length(dtu_genes)))
} else {
  warning("simulate.rda not found; sim_type will be 'Unknown' for all genes")
}

get_sim_type <- function(gene) {
  if (gene %in% dge_genes) return("DGE")
  if (gene %in% dte_genes) return("DTE")
  if (gene %in% dtu_genes) return("DTU")
  return("Background")
}

#  load ground truth 

cat("loading ground truth exonic_parts_fc files...\n")
gt_fc_files <- list.files(gt_dir, pattern = "\\.exonic_parts_fc\\.txt$",
                           full.names = TRUE)
gt_genes_available <- sub("\\.exonic_parts_fc\\.txt$", "",
                          basename(gt_fc_files))

test_genes <- unique(tests$gene)
eval_genes <- gt_genes_available   # all gt genes (total-space universe)
cat(sprintf("  gt genes: %d  |  test genes: %d  |  gt-only (untested): %d\n",
            length(gt_genes_available), length(test_genes),
            length(setdiff(gt_genes_available, test_genes))))

gt_all <- mclapply(gt_fc_files, function(f) {
  tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "",
                      colClasses = c(exonic_part = "character",
                                     fold_change = "character",
                                     source      = "character",
                                     sink        = "character")),
           error = function(e) NULL)
}, mc.cores = n_cores)
gt_all <- bind_rows(gt_all[!sapply(gt_all, is.null)])
cat(sprintf("  loaded %d exonic part records across %d genes\n",
            nrow(gt_all), length(unique(gt_all$gene))))

gt_all$sim_type <- sapply(gt_all$gene, get_sim_type)

# pre-split gt_all by gene for o(1) lookup 
gt_by_gene <- split(gt_all, gt_all$gene)

# exonic part detection 

cat("\nRunning exonic part detection...\n")

add_all_row <- function(summary_df) {
  all_rows <- summary_df %>%
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
  bind_rows(summary_df, all_rows) %>% arrange(sim_type, padj_thr)
}

eval_exonic_parts <- function(label, testable_by_gene = NULL, restrict_pos = TRUE,
                       gt_col = "group") {
  cat(sprintf("  [%s]\n", label))
  all_thresholds <- list()

  for (thr in padj_thresholds) {
    cat(sprintf("    padj < %.2f\n", thr))

    sig_rows    <- tests_eval[!is.na(tests_eval$padj) & tests_eval$padj < thr, ]
    det_by_gene <- split(sig_rows$setdiff, sig_rows$gene)

    rows <- mclapply(eval_genes, mc.cores = n_cores, function(gene) {
      gt_gene <- gt_by_gene[[gene]]
      if (is.null(gt_gene) || nrow(gt_gene) == 0) return(NULL)

      testable <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]] else NULL
      if (!is.null(testable_by_gene) &&
          (is.null(testable) || length(testable) == 0)) return(NULL)

      # Use gt_col if present; fall back to "group" for old GT files
      col <- if (gt_col %in% names(gt_gene)) gt_col else "group"
      gt_pos <- unique(gt_gene$exonic_part[
        gt_gene[[col]] %in% c("changed", "constant")])
      gt_neg <- unique(gt_gene$exonic_part[gt_gene[[col]] == "negative"])
      if (!is.null(testable)) {
        if (restrict_pos) gt_pos <- intersect(gt_pos, testable)
        gt_neg <- intersect(gt_neg, testable)
      }

      sim_type    <- gt_gene$sim_type[1]
      det_exparts <- unique(unlist(lapply(det_by_gene[[gene]], parse_exparts)))

      TP <- length(intersect(det_exparts, gt_pos))
      FP <- length(intersect(det_exparts, gt_neg))
      FN <- length(setdiff(gt_pos, det_exparts))
      TN <- length(setdiff(gt_neg, det_exparts))

      precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
      recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
      f1        <- f1_safe(precision, recall)

      data.frame(
        gene       = gene,
        sim_type   = sim_type,
        padj_thr   = thr,
        n_gt_pos   = length(gt_pos),
        n_gt_neg   = length(gt_neg),
        n_detected = length(det_exparts),
        TP = TP, FP = FP, FN = FN, TN = TN,
        precision  = round(precision, 4),
        recall     = round(recall,    4),
        f1         = round(f1,        4),
        stringsAsFactors = FALSE
      )
    })
    grase_df <- bind_rows(rows[!sapply(rows, is.null)])

    all_thresholds[[as.character(thr)]] <- grase_df

    write.table(grase_df,
                file.path(out_dir, sprintf("%s_per_gene_padj%.2f.txt",
                                          label, thr)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  summary_rows <- lapply(names(all_thresholds), function(thr_str) {
    df  <- all_thresholds[[thr_str]]
    thr <- as.numeric(thr_str)
    df %>%
      group_by(sim_type) %>%
      summarise(
        padj_thr        = thr,
        n_genes         = n(),
        total_TP        = sum(TP,  na.rm = TRUE),
        total_FP        = sum(FP,  na.rm = TRUE),
        total_FN        = sum(FN,  na.rm = TRUE),
        total_TN        = sum(TN,  na.rm = TRUE),
        micro_precision = round(sum(TP) / max(sum(TP) + sum(FP), 1), 4),
        micro_recall    = round(sum(TP) / max(sum(TP) + sum(FN), 1), 4),
        micro_f1        = round(f1_safe(
                            sum(TP) / max(sum(TP) + sum(FP), 1),
                            sum(TP) / max(sum(TP) + sum(FN), 1)), 4),
        macro_precision = round(mean(precision, na.rm = TRUE), 4),
        macro_recall    = round(mean(recall,    na.rm = TRUE), 4),
        macro_f1        = round(mean(f1,        na.rm = TRUE), 4),
        .groups = "drop"
      )
  })
  summary1 <- bind_rows(summary_rows) %>% arrange(sim_type, padj_thr)
  summary1  <- add_all_row(summary1)

  write.table(summary1,
              file.path(out_dir, sprintf("%s_summary_by_simtype.txt", label)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  list(thresholds = all_thresholds, summary = summary1)
}

# Full evaluation (total space): universe = all GT exons
eval_full  <- eval_exonic_parts("grase")
# Coverage-restricted: universe = GrASE-testable exons (setdiff only)
eval_restr <- eval_exonic_parts("grase_restricted",
                        testable_by_gene = grase_testable_by_gene,
                        restrict_pos = TRUE)


cat("\n=== Summary (full GT, micro metrics) ===\n")
print(as.data.frame(eval_full$summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

cat("\n=== Summary (restricted to GrASE-testable exons) ===\n")
print(as.data.frame(eval_restr$summary %>%
  select(sim_type, padj_thr, n_genes,
         micro_precision, micro_recall, micro_f1,
         total_TP, total_FP, total_FN)))

cat(sprintf("\nResults written to: %s\n", out_dir))
