#!/usr/bin/env Rscript
#
# scripts/evaluate_posthoc_lfc_filter.R
#
# Evaluates the post-hoc LFC filter by comparing two detection criteria on
# existing bipartition test results:
#
#   baseline : padj < threshold
#   filtered : padj < threshold & lfc_diff_net > delta
#
# lfc_diff_net = abs(lfc_diff) - abs(lfc_ref), pre-computed by exontest.R and
# present in the annotated files.  Positive values indicate the diff exon
# changes more than the ref exon (genuine DTU signal); negative values suggest
# the ref exon is the driver (denominator-effect FP).
#
# Two evaluation universes (both run, saved to separate files):
#   full       : universe = all GT exons
#   restricted : universe = GrASE-testable exons only (setdiff column)
#
# Outputs (in <out_dir>):
#   precision_recall_by_simtype.{full,restricted}.txt
#   dte_fp_breakdown.{full,restricted}.txt
#   tp_attrition.{full,restricted}.txt
#   roc_data.{full,restricted}.txt
#
# Usage:
#   Rscript evaluate_posthoc_lfc_filter.R \
#     <test_files> <gt_dir> <out_dir> <simulate_rda> [<delta>]
#
#   test_files  : comma-separated .annotated.txt result files
#   gt_dir      : directory with ground-truth .exonic_parts_fc.txt files
#   out_dir     : output directory
#   simulate_rda: path to simulate.rda with dge.genes / dte.genes / dtu.genes
#   delta       : lfc_diff_net threshold (default 0); events with
#                 lfc_diff_net <= delta are filtered out

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(grase))
n_cores <- 30

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat(paste0(
    "Usage: Rscript evaluate_posthoc_lfc_filter.R",
    " <test_files> <gt_dir> <out_dir> <simulate_rda>",
    " [<delta>]\n"))
  quit(save = "no", status = 1)
}

test_files  <- trimws(unlist(strsplit(args[1], ",")))
gt_dir      <- args[2]
out_dir     <- args[3]
sim_rda     <- args[4]
delta       <- if (length(args) >= 5) as.numeric(args[5]) else 0

cat(sprintf("delta = %g\n", delta))

padj_thresholds <- c(0.01, 0.05, 0.10, 0.20)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

parse_exparts <- function(x) {
  if (is.na(x) || trimws(x) == "" || x == "NA") return(character(0))
  trimws(unlist(strsplit(x, ",")))
}

f1_safe <- function(p, r) {
  if (is.na(p) || is.na(r) || (p + r) == 0) NA_real_
  else 2 * p * r / (p + r)
}

# --- load test results -------------------------------------------------------

cat(sprintf("loading %d test file(s)...\n", length(test_files)))
tests_list <- lapply(test_files, function(tf) {
  cat(sprintf("  %s\n", basename(tf)))
  t <- read.table(tf, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  quote = "", comment.char = "", na.strings = c("NA", ""),
                  colClasses = c(source = "character", sink = "character",
                                 comparison = "character",
                                 setdiff1 = "character", setdiff2 = "character"))
  t %>% mutate(setdiff = ifelse(comparison == "diff1_vs_ref", setdiff1, setdiff2))
})
tests <- bind_rows(tests_list)
cat(sprintf("  %d total rows, %d unique genes\n",
            nrow(tests), length(unique(tests$gene))))
n_filtered <- sum(!is.na(tests$lfc_diff_net) & tests$lfc_diff_net <= delta)
n_na       <- sum(is.na(tests$lfc_diff_net))
cat(sprintf(
  "  lfc_diff_net <= %g: %d rows would be filtered, %d NA (out of %d)\n",
  delta, n_filtered, n_na, nrow(tests)))

# --- min padj per (gene, setdiff) across all files ---------------------------
tests_eval <- tests %>%
  group_by(gene, setdiff) %>%
  slice(which.min(replace(padj, is.na(padj), Inf))) %>%
  ungroup()

# GrASE-testable exonic parts per gene (from setdiff column)
grase_testable_by_gene <- lapply(
  split(tests_eval$setdiff, tests_eval$gene),
  function(x) unique(unlist(lapply(x[!is.na(x) & nchar(x) > 0], parse_exparts)))
)
cat(sprintf("  %d genes with GrASE-testable exonic parts\n",
            length(grase_testable_by_gene)))

tests_eval_by_gene <- split(tests_eval, tests_eval$gene)

# for FP breakdown
tests$event_base <- sub("_s[12]$", "", tests$event)
tests_ev <- tests %>%
  group_by(gene, event, setdiff) %>%
  slice(which.min(replace(padj, is.na(padj), Inf))) %>%
  ungroup()
tests_ev$event_base <- sub("_s[12]$", "", tests_ev$event)
tests_ev_by_gene <- split(tests_ev, tests_ev$gene)

# --- load simulation gene type labels ----------------------------------------

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

get_sim_type <- function(gene) {
  if (gene %in% dge_genes) return("DGE")
  if (gene %in% dte_genes) return("DTE")
  if (gene %in% dtu_genes) return("DTU")
  return("Background")
}

# --- load ground truth -------------------------------------------------------

cat("loading ground truth...\n")
gt_files <- list.files(gt_dir, pattern = "\\.exonic_parts_fc\\.txt$",
                       full.names = TRUE)
gt_genes  <- sub("\\.exonic_parts_fc\\.txt$", "", basename(gt_files))

gt_all <- bind_rows(mclapply(gt_files, function(f) {
  tryCatch(
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               quote = "", comment.char = "",
               colClasses = c(exonic_part = "character",
                              fold_change = "character",
                              source      = "character",
                              sink        = "character")),
    error = function(e) NULL)
}, mc.cores = n_cores))

gt_all$sim_type <- sapply(gt_all$gene, get_sim_type)
gt_by_gene      <- split(gt_all, gt_all$gene)
eval_genes_full <- unique(gt_all$gene)
# restricted: only genes that GrASE produced test results for
eval_genes_restr <- intersect(eval_genes_full, unique(tests_eval$gene))
cat(sprintf("  GT: %d genes (%d GrASE-tested, %d GT-only)\n",
            length(eval_genes_full),
            length(eval_genes_restr),
            length(setdiff(eval_genes_full, eval_genes_restr))))

# Returns gt_pos and gt_neg for one gene.
# When testable_by_gene is not NULL (restricted mode), intersects with the
# GrASE-testable setdiff exons and applies a DTE correction: a GT positive
# exon is only kept if at least one test event has a ref_ex_part that is NOT
# in the GT changed set, i.e. the ratio test is meaningful for that exon.
get_gt_sets <- function(gene, gt_gene, testable) {
  gt_changed <- unique(gt_gene$exonic_part[gt_gene$group == "changed"])
  gt_neg_all <- unique(gt_gene$exonic_part[gt_gene$group == "negative"])

  if (is.null(testable)) {
    return(list(gt_pos = gt_changed, gt_neg = gt_neg_all))
  }

  gt_pos <- intersect(gt_changed, testable)
  gt_neg <- intersect(gt_neg_all, testable)

  if (get_sim_type(gene) == "DTE" && length(gt_pos) > 0) {
    gene_rows <- tests_eval_by_gene[[gene]]
    if (!is.null(gene_rows) && nrow(gene_rows) > 0) {
      gt_pos <- Filter(function(e) {
        in_setdiff <- sapply(gene_rows$setdiff,
                             function(s) e %in% parse_exparts(s))
        ev_rows <- gene_rows[in_setdiff, ]
        if (nrow(ev_rows) == 0) return(FALSE)
        # keep e if any event has ref_ex_part entirely outside gt_changed
        any(sapply(ev_rows$ref_ex_part, function(r) {
          !any(parse_exparts(r) %in% gt_changed)
        }))
      }, gt_pos)
    }
  }

  list(gt_pos = gt_pos, gt_neg = gt_neg)
}

# =============================================================================
# Output A: precision/recall by sim_type
# =============================================================================

cat("\n--- Output A: precision/recall by sim_type ---\n")

eval_pr_one <- function(thr, lfc_filter, testable_by_gene) {
  if (lfc_filter) {
    sig_rows <- tests_eval[!is.na(tests_eval$padj) &
                             tests_eval$padj < thr &
                             !is.na(tests_eval$lfc_diff_net) &
                             tests_eval$lfc_diff_net > delta, ]
  } else {
    sig_rows <- tests_eval[!is.na(tests_eval$padj) & tests_eval$padj < thr, ]
  }
  det_by_gene <- split(sig_rows$setdiff, sig_rows$gene)
  eval_genes  <- if (is.null(testable_by_gene)) eval_genes_full
                 else eval_genes_restr

  rows <- mclapply(eval_genes, mc.cores = n_cores, function(gene) {
    gt_gene  <- gt_by_gene[[gene]]
    testable <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]] else NULL
    if (is.null(gt_gene)) return(NULL)
    if (!is.null(testable_by_gene) &&
        (is.null(testable) || length(testable) == 0)) return(NULL)

    sets   <- get_gt_sets(gene, gt_gene, testable)
    gt_pos <- sets$gt_pos
    gt_neg <- sets$gt_neg
    if (length(gt_pos) == 0 && length(gt_neg) == 0) return(NULL)

    det_exparts <- unique(unlist(lapply(det_by_gene[[gene]], parse_exparts)))
    TP <- length(intersect(det_exparts, gt_pos))
    FP <- length(intersect(det_exparts, gt_neg))
    FN <- length(setdiff(gt_pos, det_exparts))

    data.frame(gene = gene, sim_type = get_sim_type(gene),
               padj_thr = thr, TP = TP, FP = FP, FN = FN,
               stringsAsFactors = FALSE)
  })
  bind_rows(rows[!sapply(rows, is.null)])
}

run_pr_output <- function(label, testable_by_gene) {
  cat(sprintf("  [%s]\n", label))
  pr_rows <- list()
  for (thr in padj_thresholds) {
    cat(sprintf("    padj < %.2f\n", thr))
    b <- eval_pr_one(thr, lfc_filter = FALSE, testable_by_gene)
    b$criterion <- "baseline"
    f <- eval_pr_one(thr, lfc_filter = TRUE,  testable_by_gene)
    f$criterion <- "filtered"
    pr_rows[[length(pr_rows) + 1]] <- b
    pr_rows[[length(pr_rows) + 1]] <- f
  }
  pr_all <- bind_rows(pr_rows)
  pr_summary <- pr_all %>%
    group_by(sim_type, padj_thr, criterion) %>%
    summarise(
      n_genes         = n(),
      total_TP        = sum(TP),
      total_FP        = sum(FP),
      total_FN        = sum(FN),
      micro_precision = round(sum(TP) / max(sum(TP) + sum(FP), 1), 4),
      micro_recall    = round(sum(TP) / max(sum(TP) + sum(FN), 1), 4),
      micro_f1        = round(f1_safe(
                          sum(TP) / max(sum(TP) + sum(FP), 1),
                          sum(TP) / max(sum(TP) + sum(FN), 1)), 4),
      .groups = "drop"
    ) %>%
    arrange(sim_type, padj_thr, criterion)
  outfile <- file.path(out_dir,
    sprintf("precision_recall_by_simtype.%s.txt", label))
  write.table(pr_summary, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n=== Precision/recall [%s] ===\n", label))
  print(as.data.frame(pr_summary %>%
    select(sim_type, padj_thr, criterion,
           micro_precision, micro_recall, micro_f1,
           total_TP, total_FP, total_FN)))
}

run_pr_output("full",       testable_by_gene = NULL)
run_pr_output("restricted", testable_by_gene = grase_testable_by_gene)

# =============================================================================
# Output B: FP source breakdown for DTE genes
# =============================================================================

cat("\n--- Output B: DTE FP source breakdown ---\n")

dte_genes_in_gt <- gt_genes[gt_genes %in% dte_genes]
cat(sprintf("  DTE genes in GT: %d\n", length(dte_genes_in_gt)))

eval_fp_one <- function(thr, lfc_filter, testable_by_gene) {
  rows <- mclapply(dte_genes_in_gt, mc.cores = n_cores, function(gene) {
    tryCatch({
      gt_gene    <- gt_by_gene[[gene]]
      gene_tests <- tests_ev_by_gene[[gene]]
      if (is.null(gt_gene)    || nrow(gt_gene)    == 0) return(NULL)
      if (is.null(gene_tests) || nrow(gene_tests) == 0) return(NULL)

      all_testable <- unique(unlist(lapply(gene_tests$setdiff, parse_exparts)))
      testable     <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]]
                      else all_testable
      if (!is.null(testable_by_gene) &&
          (is.null(testable) || length(testable) == 0)) return(NULL)

      sets        <- get_gt_sets(gene, gt_gene, testable)
      gt_pos_eval <- sets$gt_pos
      if (length(gt_pos_eval) == 0) return(NULL)

      gene_tests$setdiff_exons <- lapply(gene_tests$setdiff, parse_exparts)
      gene_tests$has_gt_pos    <- sapply(gene_tests$setdiff_exons,
                                         function(ex) any(ex %in% gt_pos_eval))

      if (lfc_filter) {
        sig <- gene_tests[!is.na(gene_tests$padj) & gene_tests$padj < thr &
                            !is.na(gene_tests$lfc_diff_net) &
                            gene_tests$lfc_diff_net > delta, ]
      } else {
        sig <- gene_tests[!is.na(gene_tests$padj) & gene_tests$padj < thr, ]
      }
      if (nrow(sig) == 0) return(NULL)

      base_has_correct_sig <- sig %>%
        group_by(event_base) %>%
        summarise(correct_side_sig = any(has_gt_pos), .groups = "drop")

      sig <- sig %>%
        left_join(base_has_correct_sig, by = "event_base") %>%
        mutate(row_class = case_when(
          has_gt_pos                        ~ "correct_side",
          correct_side_sig & !has_gt_pos   ~ "paired_FP",
          !correct_side_sig & !has_gt_pos  ~ "wrong_or_unrelated_FP",
          TRUE                             ~ "other"
        ))

      pos_bases <- unique(gene_tests$event_base[gene_tests$has_gt_pos])
      sig <- sig %>%
        mutate(row_class = ifelse(
          row_class == "wrong_or_unrelated_FP" & event_base %in% pos_bases,
          "wrong_side_FP", row_class)) %>%
        mutate(row_class = ifelse(
          row_class == "wrong_or_unrelated_FP",
          "unrelated_FP", row_class))

      correct_rows <- sig[sig$row_class == "correct_side", ]
      paired_rows  <- sig[sig$row_class == "paired_FP", ]
      wrong_rows   <- sig[sig$row_class == "wrong_side_FP", ]
      unrel_rows   <- sig[sig$row_class == "unrelated_FP", ]

      data.frame(
        gene           = gene,
        padj_thr       = thr,
        criterion      = if (lfc_filter) "filtered" else "baseline",
        exon_TP        = length(unique(intersect(
                           unlist(correct_rows$setdiff_exons), gt_pos_eval))),
        exon_FP_paired = length(unique(intersect(
                           unlist(paired_rows$setdiff_exons),  testable))),
        exon_FP_wrong  = length(unique(intersect(
                           unlist(wrong_rows$setdiff_exons),   testable))),
        exon_FP_unrel  = length(unique(intersect(
                           unlist(unrel_rows$setdiff_exons),   testable))),
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
  })
  bind_rows(Filter(is.data.frame, rows))
}

run_fp_output <- function(label, testable_by_gene) {
  cat(sprintf("  [%s]\n", label))
  fp_rows <- list()
  for (thr in padj_thresholds) {
    cat(sprintf("    padj < %.2f\n", thr))
    fp_rows[[length(fp_rows) + 1]] <-
      eval_fp_one(thr, lfc_filter = FALSE, testable_by_gene)
    fp_rows[[length(fp_rows) + 1]] <-
      eval_fp_one(thr, lfc_filter = TRUE,  testable_by_gene)
  }
  fp_all <- bind_rows(fp_rows)
  fp_summary <- fp_all %>%
    group_by(padj_thr, criterion) %>%
    summarise(
      exon_TP        = sum(exon_TP),
      FP_paired      = sum(exon_FP_paired),
      FP_wrong_side  = sum(exon_FP_wrong),
      FP_unrelated   = sum(exon_FP_unrel),
      .groups = "drop"
    ) %>%
    mutate(
      exon_FP        = FP_paired + FP_wrong_side + FP_unrelated,
      pct_paired     = round(100 * FP_paired     / pmax(exon_FP, 1), 1),
      pct_wrong_side = round(100 * FP_wrong_side / pmax(exon_FP, 1), 1),
      pct_unrelated  = round(100 * FP_unrelated  / pmax(exon_FP, 1), 1),
      precision      = round(exon_TP / pmax(exon_TP + exon_FP, 1), 4)
    ) %>%
    select(padj_thr, criterion, exon_TP, exon_FP,
           FP_paired, FP_wrong_side, FP_unrelated,
           pct_paired, pct_wrong_side, pct_unrelated, precision) %>%
    arrange(padj_thr, criterion)
  outfile <- file.path(out_dir,
    sprintf("dte_fp_breakdown.%s.txt", label))
  write.table(fp_summary, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n=== DTE FP breakdown [%s] ===\n", label))
  print(as.data.frame(fp_summary))
}

run_fp_output("full",       testable_by_gene = NULL)
run_fp_output("restricted", testable_by_gene = grase_testable_by_gene)

# =============================================================================
# Output C: TP attrition
# =============================================================================

cat("\n--- Output C: TP attrition ---\n")

run_attrition_output <- function(label, testable_by_gene) {
  cat(sprintf("  [%s]\n", label))
  eval_genes <- if (is.null(testable_by_gene)) eval_genes_full
                else eval_genes_restr

  attrition_gene_rows <- mclapply(eval_genes, mc.cores = n_cores, function(gene) {
    gt_gene  <- gt_by_gene[[gene]]
    testable <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]] else NULL
    if (is.null(gt_gene)) return(NULL)
    if (!is.null(testable_by_gene) &&
        (is.null(testable) || length(testable) == 0)) return(NULL)

    sets   <- get_gt_sets(gene, gt_gene, testable)
    gt_pos <- sets$gt_pos
    if (length(gt_pos) == 0) return(NULL)

    gene_rows <- tests_eval_by_gene[[gene]]
    if (is.null(gene_rows) || nrow(gene_rows) == 0) return(NULL)
    stype <- get_sim_type(gene)

    gene_rows$setdiff_exons <- lapply(gene_rows$setdiff, parse_exparts)
    gene_rows$is_tp <- sapply(gene_rows$setdiff_exons,
                              function(ex) any(ex %in% gt_pos))

    lapply(padj_thresholds, function(thr) {
      sig <- gene_rows[!is.na(gene_rows$padj) & gene_rows$padj < thr, ]
      if (nrow(sig) == 0) return(NULL)
      tp_rows <- sig[sig$is_tp, ]
      n_tp_removed <- sum(
        is.na(tp_rows$lfc_diff_net) | tp_rows$lfc_diff_net <= delta)
      data.frame(
        gene             = gene,
        sim_type         = stype,
        padj_thr         = thr,
        n_TP_baseline    = nrow(tp_rows),
        n_TP_lfc_removed = n_tp_removed,
        stringsAsFactors = FALSE
      )
    })
  })

  attrition_all <- bind_rows(unlist(attrition_gene_rows, recursive = FALSE))
  attrition_summary <- bind_rows(
    attrition_all %>%
      group_by(sim_type, padj_thr) %>%
      summarise(
        n_TP_baseline    = sum(n_TP_baseline),
        n_TP_lfc_removed = sum(n_TP_lfc_removed),
        .groups = "drop"
      ),
    attrition_all %>%
      group_by(padj_thr) %>%
      summarise(
        sim_type         = "ALL",
        n_TP_baseline    = sum(n_TP_baseline),
        n_TP_lfc_removed = sum(n_TP_lfc_removed),
        .groups = "drop"
      )
  ) %>%
    mutate(pct_TP_lost =
             round(100 * n_TP_lfc_removed / pmax(n_TP_baseline, 1), 2)) %>%
    arrange(sim_type, padj_thr)

  outfile <- file.path(out_dir, sprintf("tp_attrition.%s.txt", label))
  write.table(attrition_summary, outfile,
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n=== TP attrition [%s] ===\n", label))
  print(as.data.frame(attrition_summary))
}

run_attrition_output("full",       testable_by_gene = NULL)
run_attrition_output("restricted", testable_by_gene = grase_testable_by_gene)

# =============================================================================
# Output D: ROC data -- sweep delta on lfc_diff_net > delta
# =============================================================================

cat("\n--- Output D: ROC sweep (padj < 0.01) ---\n")

delta_values <- c(-Inf, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0)
roc_padj <- 0.01

eval_pr_delta <- function(d, thr, testable_by_gene) {
  if (is.infinite(d) && d < 0) {
    sig_rows <- tests_eval[!is.na(tests_eval$padj) & tests_eval$padj < thr, ]
  } else {
    sig_rows <- tests_eval[!is.na(tests_eval$padj) & tests_eval$padj < thr &
                             !is.na(tests_eval$lfc_diff_net) &
                             tests_eval$lfc_diff_net > d, ]
  }
  det_by_gene <- split(sig_rows$setdiff, sig_rows$gene)
  eval_genes  <- if (is.null(testable_by_gene)) eval_genes_full
                 else eval_genes_restr

  rows <- mclapply(eval_genes, mc.cores = n_cores, function(gene) {
    gt_gene  <- gt_by_gene[[gene]]
    testable <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]] else NULL
    if (is.null(gt_gene)) return(NULL)
    if (!is.null(testable_by_gene) &&
        (is.null(testable) || length(testable) == 0)) return(NULL)

    sets   <- get_gt_sets(gene, gt_gene, testable)
    gt_pos <- sets$gt_pos
    gt_neg <- sets$gt_neg
    if (length(gt_pos) == 0 && length(gt_neg) == 0) return(NULL)
    det_exparts <- unique(unlist(lapply(det_by_gene[[gene]], parse_exparts)))
    TP <- length(intersect(det_exparts, gt_pos))
    FP <- length(intersect(det_exparts, gt_neg))
    FN <- length(setdiff(gt_pos, det_exparts))
    data.frame(gene = gene, sim_type = get_sim_type(gene),
               delta = d, padj_thr = thr,
               TP = TP, FP = FP, FN = FN,
               stringsAsFactors = FALSE)
  })
  bind_rows(rows[!sapply(rows, is.null)])
}

run_roc_output <- function(label, testable_by_gene) {
  cat(sprintf("  [%s]\n", label))
  roc_rows <- list()
  for (d in delta_values) {
    cat(sprintf("    delta = %s\n",
                ifelse(is.infinite(d), "-Inf", sprintf("%.2f", d))))
    roc_rows[[length(roc_rows) + 1]] <-
      eval_pr_delta(d, roc_padj, testable_by_gene)
  }
  roc_all <- bind_rows(roc_rows)
  roc_summary <- roc_all %>%
    group_by(sim_type, delta, padj_thr) %>%
    summarise(
      n_genes         = n(),
      total_TP        = sum(TP),
      total_FP        = sum(FP),
      total_FN        = sum(FN),
      micro_precision = round(sum(TP) / max(sum(TP) + sum(FP), 1), 4),
      micro_recall    = round(sum(TP) / max(sum(TP) + sum(FN), 1), 4),
      micro_f1        = round(f1_safe(
                          sum(TP) / max(sum(TP) + sum(FN), 1),
                          sum(TP) / max(sum(TP) + sum(FP), 1)), 4),
      .groups = "drop"
    ) %>%
    arrange(sim_type, padj_thr, delta)
  outfile <- file.path(out_dir, sprintf("roc_data.%s.txt", label))
  write.table(roc_summary, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("\n=== ROC summary [%s] (padj < 0.01, DTE + DTU) ===\n", label))
  print(as.data.frame(roc_summary %>%
    filter(sim_type %in% c("DTE", "DTU")) %>%
    select(sim_type, delta, micro_precision, micro_recall, micro_f1,
           total_TP, total_FP)))
}

run_roc_output("full",       testable_by_gene = NULL)
run_roc_output("restricted", testable_by_gene = grase_testable_by_gene)

cat(sprintf("\nResults written to: %s\n", out_dir))
