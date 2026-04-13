#!/usr/bin/env Rscript
#
# scripts/evaluate_dte_bubble.R
#
# For DTE genes, evaluates GrASE detection at the bubble level.
# A bubble is "detected" if either _s1 or _s2 side has padj < threshold,
# regardless of which side the changed transcript's unique exons are on.
#
# Also identifies "wrong-side" detections: bubble detected, but the
# significant side does not contain any GT-positive (changed) exons --
# i.e. GrASE called the other side of the bipartition.
#
# Outputs (in <out_dir>):
#   dte_bubble_per_gene.txt   -- per-gene bubble-level + exon-level metrics
#   dte_bubble_summary.txt    -- aggregated summary across all DTE genes
#
# Usage:
#   Rscript evaluate_dte_bubble.R <test_results> <gt_dir> <out_dir> <simulate_rda>

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
n_cores <- 30

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript evaluate_dte_bubble.R <test_results> <gt_dir> <out_dir> <simulate_rda> [<delta>]\n")
  quit(save = "no", status = 1)
}

test_files <- trimws(unlist(strsplit(args[1], ",")))
gt_dir     <- args[2]
out_dir    <- args[3]
sim_rda    <- args[4]
delta      <- if (length(args) >= 5) as.numeric(args[5]) else 0

cat(sprintf("lfc_diff_net delta = %g\n", delta))

padj_thresholds <- c(0.01, 0.05, 0.1, 0.2)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

parse_exparts <- function(x) {
  if (is.na(x) || trimws(x) == "" || x == "NA") return(character(0))
  trimws(unlist(strsplit(x, ",")))
}

# load test results -----------------------------------------------------

cat(sprintf("loading %d test file(s)...\n", length(test_files)))
tests_list <- lapply(test_files, function(f) {
  cat(sprintf("  %s\n", basename(f)))
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", comment.char = "", na.strings = c("NA", ""),
             colClasses = c(source = "character", sink = "character",
                            comparison = "character",
                            setdiff1 = "character", setdiff2 = "character"))
})
tests <- bind_rows(tests_list)
cat(sprintf("  %d rows, %d genes\n", nrow(tests), length(unique(tests$gene))))

# derive comparison-relevant setdiff: each row is one comparison direction,
# so use setdiff1 for diff1_vs_ref and setdiff2 for diff2_vs_ref.
if ("comparison" %in% names(tests) && "setdiff1" %in% names(tests)) {
  tests <- tests %>% mutate(
    setdiff = ifelse(comparison == "diff1_vs_ref", setdiff1, setdiff2)
  )
}

# keep all comparison rows; each row is evaluated independently.
tests_ev <- tests %>%
  select(gene, event, comparison, setdiff, padj, any_of("lfc_diff_net")) %>%
  mutate(event_base = sub("_s[12]$", "", event))

# load simulation type labels -------------------------------------------

cat("loading simulate.rda...\n")
dge_genes <- dte_genes <- dtu_genes <- character(0)
if (file.exists(sim_rda)) {
  load(sim_rda)
  if (exists("dge.genes")) dge_genes <- dge.genes
  if (exists("dte.genes")) dte_genes <- dte.genes
  if (exists("dtu.genes")) dtu_genes <- dtu.genes
  cat(sprintf("  DGE: %d  DTE: %d  DTU: %d\n",
              length(dge_genes), length(dte_genes), length(dtu_genes)))
} else {
  warning("simulate.rda not found")
}

# load ground truth -----------------------------------------------------

cat("loading GT...\n")
gt_files <- list.files(gt_dir, pattern = "\\.exonic_parts_fc\\.txt$",
                        full.names = TRUE)
gt_genes  <- sub("\\.exonic_parts_fc\\.txt$", "", basename(gt_files))

gt_all <- bind_rows(mclapply(gt_files, function(f) {
  tryCatch(
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               quote = "", comment.char = "",
               colClasses = c(exonic_part  = "character",
                              fold_change  = "character",
                              source       = "character",
                              sink         = "character")),
    error = function(e) NULL)
}, mc.cores = n_cores))

gt_all$sim_type <- sapply(gt_all$gene, function(g) {
  if (g %in% dge_genes) "DGE"
  else if (g %in% dte_genes) "DTE"
  else if (g %in% dtu_genes) "DTU"
  else "Background"
})
gt_by_gene      <- split(gt_all, gt_all$gene)
tests_by_gene   <- split(tests_ev, tests_ev$gene)

dte_genes_in_gt <- gt_genes[gt_genes %in% dte_genes]
cat(sprintf("  DTE genes in GT: %d\n", length(dte_genes_in_gt)))

# bubble-level evaluation -----------------------------------------------

eval_one_threshold <- function(thr, use_lfc = FALSE) {
  cat(sprintf("  padj < %.2f%s\n", thr, if (use_lfc) "  [+lfc_diff_net]" else ""))

  lfc_ok_ev <- if (use_lfc && "lfc_diff_net" %in% names(tests_ev))
                 (is.na(tests_ev$lfc_diff_net) | tests_ev$lfc_diff_net > delta)
               else TRUE
  sig_ev <- tests_ev %>%
    filter(!is.na(padj) & padj < thr & lfc_ok_ev) %>%
    select(gene, event, event_base) %>%
    distinct()
  sig_by_gene <- split(sig_ev, sig_ev$gene)

  rows <- mclapply(dte_genes_in_gt, mc.cores = n_cores, function(gene) {
    tryCatch({
      gt_gene    <- gt_by_gene[[gene]]
      gene_tests <- tests_by_gene[[gene]]

      if (is.null(gt_gene) || nrow(gt_gene) == 0) return(NULL)

      gt_pos_all <- unique(gt_gene$exonic_part[gt_gene$group == "changed"])

      if (is.null(gene_tests) || nrow(gene_tests) == 0) {
        return(data.frame(
          gene = gene, padj_thr = thr,
          n_gt_pos_bubbles = NA_integer_,
          bubble_TP = 0L, bubble_FN = NA_integer_,
          wrong_side = NA_integer_,
          bubble_recall = NA_real_, pct_wrong_side = NA_real_,
          exon_TP = 0L, exon_FP = 0L, exon_FN = NA_integer_,
          exon_recall = NA_real_, exon_precision = NA_real_,
          stringsAsFactors = FALSE))
      }

      testable_exons    <- unique(unlist(lapply(gene_tests$setdiff, parse_exparts)))
      gt_pos_restricted <- intersect(gt_pos_all, testable_exons)
      gt_neg_restricted <- setdiff(testable_exons, gt_pos_all)

      if (length(gt_pos_restricted) == 0) return(NULL)

      gene_tests$setdiff_exons <- lapply(gene_tests$setdiff, parse_exparts)
      gene_tests$has_gt_pos <- sapply(gene_tests$setdiff_exons,
                                      function(ex) any(ex %in% gt_pos_restricted))

      pos_event_bases <- unique(gene_tests$event_base[gene_tests$has_gt_pos])
      n_pos_bubbles   <- length(pos_event_bases)

      gene_sig  <- sig_by_gene[[gene]]
      sig_bases <- if (is.null(gene_sig)) character(0) else unique(gene_sig$event_base)

      detected <- intersect(pos_event_bases, sig_bases)
      missed   <- setdiff(pos_event_bases, sig_bases)
      bubble_TP <- length(detected)
      bubble_FN <- length(missed)

      wrong_side <- 0L
      for (eb in detected) {
        lfc_ok_gt <- if (use_lfc && "lfc_diff_net" %in% names(gene_tests))
                       (is.na(gene_tests$lfc_diff_net) | gene_tests$lfc_diff_net > delta)
                     else TRUE
        sig_rows <- gene_tests[gene_tests$event_base == eb &
                                 !is.na(gene_tests$padj) &
                                 gene_tests$padj < thr & lfc_ok_gt, ]
        if (!any(sig_rows$has_gt_pos)) wrong_side <- wrong_side + 1L
      }

      lfc_ok_all <- if (use_lfc && "lfc_diff_net" %in% names(gene_tests))
                      (is.na(gene_tests$lfc_diff_net) | gene_tests$lfc_diff_net > delta)
                    else TRUE
      sig_rows_all <- gene_tests[!is.na(gene_tests$padj) & gene_tests$padj < thr &
                                   lfc_ok_all, ]
      det_exparts  <- unique(unlist(sig_rows_all$setdiff_exons))
      exon_TP <- length(intersect(det_exparts, gt_pos_restricted))
      exon_FP <- length(intersect(det_exparts, gt_neg_restricted))
      exon_FN <- length(setdiff(gt_pos_restricted, det_exparts))

      bubble_recall  <- bubble_TP / max(bubble_TP + bubble_FN, 1)
      pct_wrong_side <- if (bubble_TP > 0) wrong_side / bubble_TP else NA_real_
      exon_recall    <- if ((exon_TP + exon_FN) > 0) exon_TP / (exon_TP + exon_FN) else NA_real_
      exon_prec      <- if ((exon_TP + exon_FP) > 0) exon_TP / (exon_TP + exon_FP) else NA_real_

      data.frame(
        gene = gene, padj_thr = thr,
        n_gt_pos_bubbles = n_pos_bubbles,
        bubble_TP = bubble_TP, bubble_FN = bubble_FN,
        wrong_side = wrong_side,
        bubble_recall = round(bubble_recall, 4),
        pct_wrong_side = round(pct_wrong_side, 4),
        exon_TP = exon_TP, exon_FP = exon_FP, exon_FN = exon_FN,
        exon_recall = round(exon_recall, 4),
        exon_precision = round(exon_prec, 4),
        stringsAsFactors = FALSE)
    }, error = function(e) NULL)
  })

  per_gene <- bind_rows(Filter(is.data.frame, rows))

  n_pos  <- sum(per_gene$n_gt_pos_bubbles, na.rm = TRUE)
  bTP    <- sum(per_gene$bubble_TP,  na.rm = TRUE)
  bFN    <- sum(per_gene$bubble_FN,  na.rm = TRUE)
  bWrong <- sum(per_gene$wrong_side, na.rm = TRUE)
  eTP    <- sum(per_gene$exon_TP, na.rm = TRUE)
  eFP    <- sum(per_gene$exon_FP, na.rm = TRUE)
  eFN    <- sum(per_gene$exon_FN, na.rm = TRUE)

  summary <- data.frame(
    padj_thr         = thr,
    n_dte_genes      = nrow(per_gene),
    n_gt_pos_bubbles = n_pos,
    bubble_TP        = bTP,
    bubble_FN        = bFN,
    bubble_recall    = round(bTP / max(bTP + bFN, 1), 4),
    wrong_side       = bWrong,
    pct_wrong_side   = round(100 * bWrong / max(bTP, 1), 2),
    exon_TP          = eTP, exon_FP = eFP, exon_FN = eFN,
    exon_recall      = round(eTP / max(eTP + eFN, 1), 4),
    exon_precision   = round(eTP / max(eTP + eFP, 1), 4),
    stringsAsFactors = FALSE)

  cat(sprintf("    bubble recall: %.3f  (%d/%d)\n",
              summary$bubble_recall, bTP, bTP + bFN))
  cat(sprintf("    wrong-side:    %d / %d detected (%.1f%%)\n",
              bWrong, bTP, summary$pct_wrong_side))
  cat(sprintf("    exon recall (restricted): %.3f  precision: %.3f\n",
              summary$exon_recall, summary$exon_precision))

  list(per_gene = per_gene, summary = summary)
}

write_results <- function(results, per_gene_file, summary_file) {
  per_gene_all <- bind_rows(lapply(results, `[[`, "per_gene"))
  summary_all  <- bind_rows(lapply(results, `[[`, "summary"))
  write.table(per_gene_all, per_gene_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_all,  summary_file,  sep = "\t", quote = FALSE, row.names = FALSE)
  summary_all
}

cat("\nRunning bubble-level DTE evaluation (baseline)...\n")
res_base <- lapply(padj_thresholds, eval_one_threshold, use_lfc = FALSE)
summary_base <- write_results(res_base,
  file.path(out_dir, "dte_bubble_per_gene.txt"),
  file.path(out_dir, "dte_bubble_summary.txt"))

cat("\nRunning bubble-level DTE evaluation (lfc_diff_net filtered)...\n")
res_lfc <- lapply(padj_thresholds, eval_one_threshold, use_lfc = TRUE)
summary_lfc <- write_results(res_lfc,
  file.path(out_dir, "dte_bubble_per_gene_lfc_filtered.txt"),
  file.path(out_dir, "dte_bubble_summary_lfc_filtered.txt"))

cat(sprintf("\nResults written to: %s\n", out_dir))
cat("\n=== Summary (baseline) ===\n")
print(summary_base)
cat("\n=== Summary (lfc_diff_net filtered) ===\n")
print(summary_lfc)
