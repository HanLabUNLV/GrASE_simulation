#!/usr/bin/env Rscript

# scripts/evaluate_bipartition_sj_test.R
#
# SJ-specific fork of evaluate_bipartition_test.R.
# Differences from the exon-based version:
# 1. Only the d1d2 comparison (diff2_vs_diff1) is used. d1r and d2r are
#    dropped: junction "ref" counts are noisier than exon body counts and the
#    lfc_diff_net filter does not sufficiently suppress DTE confounding for
#    junction tests.
# 2. d1d2 is credited to the union of setdiff1 and setdiff2 (not setdiff2
#    alone) because path assignment is arbitrary.
# 3. lfc_diff_net filter is not applied to d1d2 (direction-agnostic).
#
# Evaluates GrASE SJ bipartition statistical test results against simulation
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
#   GT positive       = exonic parts labelled "changed"
#   GT negative       = exonic parts labelled "negative"
#   Reports TP/FP/FN/TN, precision, recall, F1 at multiple padj thresholds
#
#   Two evaluations:
#     Full (total space)   : universe = all GT genes; untested = implicit neg
#     Restricted           : universe = only GrASE-testable exons (setdiff)
#
# Usage:
#   Rscript evaluate_bipartition_test.R <test_results> <gt_dir> <out_dir> <simulate_rda> [<delta>]
#
#   <test_results> may be a single file path or a comma-separated list of paths.
#   When multiple files are given they are combined: the row with minimum
#   padj is kept per (gene, setdiff).
#   <delta>        lfc_diff_net threshold (default 0); significant calls also
#                  require lfc_diff_net > delta (NA lfc_diff_net rows pass).
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
  cat("Usage: Rscript evaluate_bipartition_test.R <test_results> <gt_dir> <out_dir> <simulate_rda> [<delta>]\n")
  quit(save = "no", status = 1)
}

test_files  <- trimws(unlist(strsplit(args[1], ",")))
gt_dir      <- args[2]
out_dir     <- args[3]
sim_rda     <- args[4]
delta       <- if (length(args) >= 5) as.numeric(args[5]) else 0

cat(sprintf("lfc_diff_net delta = %g\n", delta))

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

# Return GT-positive exonic parts for one gene's GT data frame.
# For DTE genes, dte_gt_level controls stringency:
#   "i"   -- all exons of changed transcripts (group == "changed")
#   "ii"  -- minus constitutive exons (dte_gt_level %in% c("shared","unique"))
#   "iii" -- unique to changed transcripts only (dte_gt_level == "unique")
# For non-DTE genes the level is ignored; group == "changed" is always used.
# gt_neg is always group == "negative" at every level -- exons removed from
# gt_pos at levels ii/iii (constitutive, shared) are excluded from the
# evaluation entirely (not counted as FP when detected).
get_gt_pos <- function(gt_gene, dte_gt_level = "i") {
  is_dte <- !is.null(gt_gene$sim_type) && any(gt_gene$sim_type == "DTE")
  if (!is_dte || dte_gt_level == "i") {
    return(unique(gt_gene$exonic_part[gt_gene$group == "changed"]))
  }
  if (!"dte_gt_level" %in% names(gt_gene))
    stop("dte_gt_level column missing -- rerun infer_diff_exons_gt.R")
  lvl <- gt_gene$dte_gt_level
  keep <- if (dte_gt_level == "ii") lvl %in% c("shared", "unique")
          else                      lvl == "unique"   # "iii"
  unique(gt_gene$exonic_part[!is.na(keep) & keep])
}

#  load test results 

cat(sprintf("loading %d test file(s)...\n", length(test_files)))
tests_list <- lapply(test_files, function(f) {
  cat(sprintf("  %s\n", basename(f)))
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", comment.char = "", na.strings = c("NA", ""),
             colClasses = c(source = "character", sink = "character",
                            comparison = "character",
                            setdiff1 = "character", setdiff2 = "character",
                            setdiff_union = "character"))
})
tests <- bind_rows(tests_list)
cat(sprintf("  %d total rows, %d unique genes across all files\n",
            nrow(tests), length(unique(tests$gene))))

# SJ pipeline uses d1d2 only; drop d1r/d2r rows.
tests <- tests[tests$comparison == "diff2_vs_diff1" | is.na(tests$comparison), ]
cat(sprintf("  %d rows after filtering to d1d2 only\n", nrow(tests)))

# Each row is one comparison (d1d2 only for SJ).
tests_eval <- tests %>%
  select(gene, event, comparison, padj, lfc_diff_net, setdiff1, setdiff2,
         any_of("setdiff_union"))

# build per-gene set of exonic parts grase explicitly tested (setdiff only)
# setdiff may be a comma-separated list of exonic parts; split before deduplicating
grase_testable_by_gene <- lapply(
  split(tests_eval %>% select(setdiff1, setdiff2), tests_eval$gene),
  function(df) {
    all_setdiffs <- c(df$setdiff1, df$setdiff2)
    unique(unlist(lapply(all_setdiffs[!is.na(all_setdiffs) & nchar(all_setdiffs) > 0], parse_exparts)))
  }
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

# per-event lfc_diff_net vs TP/FP labeling

cat("\nWriting per-event lfc_diff_net labels...\n")
per_event_rows <- mclapply(seq_len(nrow(tests_eval)), mc.cores = n_cores, function(i) {
  row <- tests_eval[i, ]
  gene <- row$gene
  event <- row$event # Get event ID
  gt_gene <- gt_by_gene[[gene]]
  sim_type <- get_sim_type(gene)

  # Resolve the exonic parts to check for overlaps (SJ-specific rules):
  # - combination files: use setdiff_union.
  # - diff1_vs_ref  -> setdiff1 only.
  # - diff2_vs_ref  -> setdiff2 only.
  # - diff2_vs_diff1 -> union of setdiff1 and setdiff2.
  if (!is.null(row$setdiff_union) && !is.na(row$setdiff_union) && nchar(row$setdiff_union) > 0) {
    combined_setdiff <- row$setdiff_union
    event_setdiffs   <- row$setdiff_union
  } else if (!is.null(row$comparison) && !is.na(row$comparison) &&
             row$comparison == "diff2_vs_diff1") {
    raw <- c(row$setdiff1, row$setdiff2)
    raw <- raw[!is.na(raw) & nchar(raw) > 0]
    combined_setdiff <- if (length(raw) > 0) paste(raw, collapse = ",") else NA_character_
    event_setdiffs   <- raw
  } else {
    sd <- if (!is.null(row$comparison) && !is.na(row$comparison) &&
              row$comparison == "diff1_vs_ref") row$setdiff1 else row$setdiff2
    combined_setdiff <- if (!is.na(sd) && nchar(sd) > 0) sd else NA_character_
    event_setdiffs   <- if (!is.na(sd) && nchar(sd) > 0) sd else character(0)
  }
  ex <- unique(unlist(lapply(
    event_setdiffs[!is.na(event_setdiffs) & nchar(event_setdiffs) > 0], parse_exparts)))

  if (is.null(gt_gene) || nrow(gt_gene) == 0) {
    return(data.frame(gene = gene, event = event, comparison = row$comparison,
                      setdiff = combined_setdiff,
                      padj = row$padj, lfc_diff_net = row$lfc_diff_net,
                      sim_type = sim_type,
                      overlaps_gt_pos = NA, overlaps_gt_neg = NA,
                      stringsAsFactors = FALSE))
  }

  gt_pos <- get_gt_pos(gt_gene, dte_gt_level = "i")   # level i for per-event labeling
  gt_neg <- unique(gt_gene$exonic_part[gt_gene$group == "negative"])

  data.frame(gene = gene, event = event, comparison = row$comparison,
             setdiff = combined_setdiff,
             padj = row$padj, lfc_diff_net = row$lfc_diff_net,
             sim_type = sim_type,
             overlaps_gt_pos = any(ex %in% gt_pos),
             overlaps_gt_neg = any(ex %in% gt_neg),
             stringsAsFactors = FALSE)
})
per_event_df <- bind_rows(per_event_rows)
write.table(per_event_df,
            file.path(out_dir, "per_event_lfc_diff_net.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  wrote %d rows to per_event_lfc_diff_net.txt\n", nrow(per_event_df)))

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
                              dte_gt_level = "i") {
  cat(sprintf("  [%s  GT level %s]\n", label, dte_gt_level))
  all_thresholds <- list()

  for (thr in padj_thresholds) {
    cat(sprintf("    padj < %.2f\n", thr))

    # d1d2 is direction-agnostic (path assignment is arbitrary), so lfc_diff_net
    # filter does not apply to diff2_vs_diff1 rows.
    sig_rows <- tests_eval[!is.na(tests_eval$padj) & tests_eval$padj < thr &
                             (tests_eval$comparison == "diff2_vs_diff1" |
                              is.na(tests_eval$lfc_diff_net) |
                              tests_eval$lfc_diff_net > delta), ]
    
    # Detected exonic parts (SJ-specific rules):
    # - combination files (setdiff_union present): use the precomputed union.
    # - diff1_vs_ref  -> setdiff1
    # - diff2_vs_ref  -> setdiff2
    # - diff2_vs_diff1 (d1d2) -> union of setdiff1 and setdiff2: the test
    #   detects differential bipartition usage regardless of which path carries
    #   the distinct exon.
    has_union <- "setdiff_union" %in% names(sig_rows)
    det_by_gene <- lapply(
      split(sig_rows %>% select(comparison, setdiff1, setdiff2,
                                any_of("setdiff_union")), sig_rows$gene),
      function(df) {
        parts <- if (has_union) {
          df$setdiff_union
        } else {
          vapply(seq_len(nrow(df)), function(j) {
            comp <- df$comparison[j]
            sd1  <- df$setdiff1[j]
            sd2  <- df$setdiff2[j]
            if (comp == "diff1_vs_ref")
              return(if (!is.na(sd1) && nchar(sd1) > 0) sd1 else NA_character_)
            if (comp == "diff2_vs_ref")
              return(if (!is.na(sd2) && nchar(sd2) > 0) sd2 else NA_character_)
            # diff2_vs_diff1: union of both sides
            ex <- unique(c(parse_exparts(sd1), parse_exparts(sd2)))
            if (length(ex) == 0L) return(NA_character_)
            paste(ex, collapse = ",")
          }, character(1L))
        }
        unique(unlist(lapply(parts[!is.na(parts) & nchar(parts) > 0], parse_exparts)))
      }
    )

    rows <- mclapply(eval_genes, mc.cores = n_cores, function(gene) {
      gt_gene <- gt_by_gene[[gene]]
      if (is.null(gt_gene) || nrow(gt_gene) == 0) return(NULL)

      testable <- if (!is.null(testable_by_gene)) testable_by_gene[[gene]] else NULL
      if (!is.null(testable_by_gene) &&
          (is.null(testable) || length(testable) == 0)) return(NULL)

      gt_pos <- get_gt_pos(gt_gene, dte_gt_level)
      gt_neg <- unique(gt_gene$exonic_part[gt_gene$group == "negative"])
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

for (lv in c("i", "ii", "iii")) {
  lv_tag <- paste0("gt", toupper(lv))
  # Full evaluation (total space): universe = all GT exons
  eval_full  <- eval_exonic_parts(sprintf("grase_%s",            lv_tag),
                                  dte_gt_level = lv)
  # Coverage-restricted: universe = GrASE-testable exons (setdiff only)
  eval_restr <- eval_exonic_parts(sprintf("grase_restricted_%s", lv_tag),
                                  testable_by_gene = grase_testable_by_gene,
                                  restrict_pos     = TRUE,
                                  dte_gt_level     = lv)

  cat(sprintf("\n=== Summary (full GT level %s, micro metrics) ===\n", lv))
  print(as.data.frame(eval_full$summary %>%
    select(sim_type, padj_thr, n_genes,
           micro_precision, micro_recall, micro_f1,
           total_TP, total_FP, total_FN)))

  cat(sprintf("\n=== Summary (restricted, GT level %s) ===\n", lv))
  print(as.data.frame(eval_restr$summary %>%
    select(sim_type, padj_thr, n_genes,
           micro_precision, micro_recall, micro_f1,
           total_TP, total_FP, total_FN)))
}

cat(sprintf("\nResults written to: %s\n", out_dir))
