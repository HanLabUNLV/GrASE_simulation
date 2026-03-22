#!/usr/bin/env Rscript
#
# scripts/evaluate_dte_fp_source.R
#
# For DTE genes, decomposes false positives by source:
#
#   "paired_FP"   -- FP exons from the OTHER side of a bubble whose correct
#                    side (GT-positive setdiff) was ALSO called significant.
#                    These are the "other-side" FPs: GrASE found the real
#                    signal but also called the mirror side.
#
#   "wrong_side_FP" -- FP exons from a bubble where the correct side was NOT
#                    called significant. GrASE called the wrong side only.
#
#   "unrelated_FP" -- FP exons from bubbles that have NO GT-positive side at
#                    all (neither side's setdiff contains a GT-positive exon).
#                    These are entirely spurious calls.
#
# Outputs (in <out_dir>):
#   dte_fp_source_summary.txt    -- aggregate FP breakdown per threshold
#   dte_fp_source_per_gene.txt   -- per-gene FP breakdown
#
# Usage:
#   Rscript evaluate_dte_fp_source.R <test_results> <gt_dir> <out_dir> <simulate_rda>

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
n_cores <- 30

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: Rscript evaluate_dte_fp_source.R <test_results> <gt_dir> <out_dir> <simulate_rda>\n")
  quit(save = "no", status = 1)
}

test_files <- trimws(unlist(strsplit(args[1], ",")))
gt_dir     <- args[2]
out_dir    <- args[3]
sim_rda    <- args[4]

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
             colClasses = c(source = "character", sink = "character"))
})
tests <- bind_rows(tests_list)

# min padj per (gene, event) across files; keep setdiff
tests_ev <- tests %>%
  group_by(gene, event, setdiff) %>%
  summarise(padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
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

gt_by_gene    <- split(gt_all, gt_all$gene)
tests_by_gene <- split(tests_ev, tests_ev$gene)

dte_genes_in_gt <- gt_genes[gt_genes %in% dte_genes]
cat(sprintf("  DTE genes in GT: %d\n", length(dte_genes_in_gt)))

# per-threshold FP decomposition ----------------------------------------

eval_one_threshold <- function(thr) {
  cat(sprintf("  padj < %.2f\n", thr))

  rows <- mclapply(dte_genes_in_gt, mc.cores = n_cores, function(gene) {
  tryCatch({

    gt_gene    <- gt_by_gene[[gene]]
    gene_tests <- tests_by_gene[[gene]]

    if (is.null(gt_gene) || nrow(gt_gene) == 0) return(NULL)
    if (is.null(gene_tests) || nrow(gene_tests) == 0) return(NULL)

    gt_pos_all <- unique(gt_gene$exonic_part[gt_gene$group == "changed"])

    # restrict to testable (setdiff) exons
    testable_exons    <- unique(unlist(lapply(gene_tests$setdiff, parse_exparts)))
    gt_pos_restricted <- intersect(gt_pos_all, testable_exons)
    if (length(gt_pos_restricted) == 0) return(NULL)

    # annotate each test row: parse setdiff, classify GT
    gene_tests$setdiff_exons <- lapply(gene_tests$setdiff, parse_exparts)
    gene_tests$has_gt_pos <- sapply(gene_tests$setdiff_exons,
                                    function(ex) any(ex %in% gt_pos_restricted))

    # only significant rows
    sig <- gene_tests[!is.na(gene_tests$padj) & gene_tests$padj < thr, ]
    if (nrow(sig) == 0) return(NULL)

    # for each event_base: is there any significant GT-positive side?
    base_has_correct_sig <- sig %>%
      group_by(event_base) %>%
      summarise(correct_side_sig = any(has_gt_pos), .groups = "drop")

    # classify each significant row
    sig <- sig %>%
      left_join(base_has_correct_sig, by = "event_base") %>%
      mutate(row_class = case_when(
        has_gt_pos                        ~ "correct_side",
        correct_side_sig & !has_gt_pos   ~ "paired_FP",
        !correct_side_sig & !has_gt_pos  ~ "wrong_or_unrelated_FP",
        TRUE                             ~ "other"
      ))

    # for wrong_or_unrelated: split into wrong_side (pos bubble missed) vs unrelated
    # a bubble is a "positive bubble" if it has at least one row with has_gt_pos
    pos_bases <- unique(gene_tests$event_base[gene_tests$has_gt_pos])

    sig <- sig %>%
      mutate(row_class = ifelse(
        row_class == "wrong_or_unrelated_FP" & event_base %in% pos_bases,
        "wrong_side_FP",   # positive bubble but only wrong side called
        row_class
      )) %>%
      mutate(row_class = ifelse(
        row_class == "wrong_or_unrelated_FP",
        "unrelated_FP",    # event_base has no GT-positive side at all
        row_class
      ))

    # count FP exons per class
    # TP exons: from correct_side rows
    correct_rows  <- sig[sig$row_class == "correct_side", ]
    paired_rows   <- sig[sig$row_class == "paired_FP", ]
    wrong_rows    <- sig[sig$row_class == "wrong_side_FP", ]
    unrel_rows    <- sig[sig$row_class == "unrelated_FP", ]

    exon_TP         <- length(unique(intersect(
      unlist(correct_rows$setdiff_exons), gt_pos_restricted)))
    exon_FP_paired  <- length(unique(intersect(
      unlist(paired_rows$setdiff_exons),  testable_exons)))
    exon_FP_wrong   <- length(unique(intersect(
      unlist(wrong_rows$setdiff_exons),   testable_exons)))
    exon_FP_unrel   <- length(unique(intersect(
      unlist(unrel_rows$setdiff_exons),   testable_exons)))

    data.frame(
      gene            = gene,
      padj_thr        = thr,
      exon_TP         = exon_TP,
      exon_FP_paired  = exon_FP_paired,   # other side of a TP bubble
      exon_FP_wrong   = exon_FP_wrong,    # wrong side only (no TP side called)
      exon_FP_unrel   = exon_FP_unrel,    # unrelated bubble
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
  })

  per_gene <- bind_rows(Filter(is.data.frame, rows))

  total_TP      <- sum(per_gene$exon_TP,         na.rm = TRUE)
  total_paired  <- sum(per_gene$exon_FP_paired,  na.rm = TRUE)
  total_wrong   <- sum(per_gene$exon_FP_wrong,   na.rm = TRUE)
  total_unrel   <- sum(per_gene$exon_FP_unrel,   na.rm = TRUE)
  total_FP      <- total_paired + total_wrong + total_unrel

  summary <- data.frame(
    padj_thr       = thr,
    exon_TP        = total_TP,
    exon_FP        = total_FP,
    FP_paired      = total_paired,
    FP_wrong_side  = total_wrong,
    FP_unrelated   = total_unrel,
    pct_paired     = round(100 * total_paired / max(total_FP, 1), 1),
    pct_wrong_side = round(100 * total_wrong  / max(total_FP, 1), 1),
    pct_unrelated  = round(100 * total_unrel  / max(total_FP, 1), 1),
    precision      = round(total_TP / max(total_TP + total_FP, 1), 4),
    stringsAsFactors = FALSE
  )

  cat(sprintf("    TP: %d  FP: %d  (precision: %.3f)\n",
              total_TP, total_FP, summary$precision))
  cat(sprintf("    FP breakdown -- paired: %d (%.1f%%)  wrong_side: %d (%.1f%%)  unrelated: %d (%.1f%%)\n",
              total_paired,  summary$pct_paired,
              total_wrong,   summary$pct_wrong_side,
              total_unrel,   summary$pct_unrelated))

  list(per_gene = per_gene, summary = summary)
}

cat("\nRunning FP source decomposition...\n")
results <- lapply(padj_thresholds, eval_one_threshold)

per_gene_all <- bind_rows(lapply(results, `[[`, "per_gene"))
summary_all  <- bind_rows(lapply(results, `[[`, "summary"))

write.table(per_gene_all,
            file.path(out_dir, "dte_fp_source_per_gene.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_all,
            file.path(out_dir, "dte_fp_source_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nResults written to: %s\n", out_dir))
cat("\n=== Summary ===\n")
print(summary_all)
