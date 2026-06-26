#!/usr/bin/env Rscript

# scripts/evaluate_tools_vs_gt.R
#
# Evaluate rMATS (junction-level GT) and Saturn (exon-bin GT) against
# simulation ground truth.
#
# rMATS: evaluated using junction-level GT from infer_junctions_gt.R
#   Event files : <rmats_dir>/{SE,A3SS,A5SS,RI}.MATS.JCEC.txt
#   Junction GT : <junction_gt_file> (sim_junction_gt.txt)
#   Significance: FDR column
#
# Saturn:
#   File        : <saturn_file>
#                 ids column format: "ENSG...:Exxx"  (gene:exonic_part)
#   Significance: empirical_FDR column (falls back to regular_FDR)
#
# Usage:
#   Rscript evaluate_tools_vs_gt.R \
#     <gt_dir> <rmats_dir> <map_dir> <saturn_file> <out_dir> <simulate_rda> \
#     [junction_gt_file]
#   map_dir is unused but kept for positional compatibility.
#   junction_gt_file defaults to <results_dir>/sim_junction_gt.txt.
#
# Outputs (written to out_dir):
#   rmats_junctiongt_per_gene_padj<thr>.txt              per-gene (full junction GT)
#   rmats_junctiongt_restricted_per_gene_padj<thr>.txt   per-gene (rMATS-tested events only)
#   rmats_junctiongt_summary_by_simtype.txt              aggregated (full)
#   rmats_junctiongt_restricted_summary_by_simtype.txt   aggregated (restricted)
#   {saturn}_per_gene_padj<thr>.txt                      per-gene (full exon-bin GT)
#   {saturn}_restricted_per_gene_padj<thr>.txt           per-gene (Saturn-testable only)
#   {saturn}_summary_by_simtype.txt                      aggregated (full)
#   {saturn}_restricted_summary_by_simtype.txt           aggregated (restricted)

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

gt_dir          <- args[1]
rmats_dir       <- args[2]
map_dir         <- args[3]
saturn_file     <- args[4]
out_dir         <- args[5]
sim_rda         <- args[6]
junction_gt_file <- if (length(args) >= 7) args[7] else
  file.path(dirname(args[5]), "sim_junction_gt.txt")

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
  gt_pos <- unique(gt_gene_df$exonic_part[gt_gene_df$group == "changed"])
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

#  evaluate tools across padj thresholds

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
      # (also skips genes absent from testable_by_gene, where testable is NULL)
      if (!is.null(testable_by_gene) &&
          (is.null(testable) || length(testable) == 0)) return(NULL)

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

# junction-level GT evaluation for rMATS (runs before Saturn)

rmats_junctiongt_per_gene        <- NULL
rmats_junctiongt_summary         <- NULL

if (file.exists(junction_gt_file)) {
  cat(sprintf("\nLoading junction GT: %s\n", junction_gt_file))
  jgt <- read.table(junction_gt_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
  # exclude untestable events: either side empty -> PSI constant/undefined,
  # cannot be scored. Often our exact-exon matching missed an isoform rMATS
  # sees in reads, so these are untestable, not true negatives.
  jgt <- jgt[!(jgt$n_inc_tx == 0 | jgt$n_skip_tx == 0), ]
  cat(sprintf("  %d testable events (%d GT-positive, %d GT-negative)\n",
              nrow(jgt), sum(jgt$gt_positive), sum(!jgt$gt_positive)))

  # helper: mean of comma-separated numeric string (NA if all NA)
  csv_mean <- function(x) {
    v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  }

  # load raw rMATS FDR + read count/PSI filter columns per event
  rmats_fdr_rows <- list()
  for (etype in c("SE", "A3SS", "A5SS", "RI")) {
    f <- file.path(rmats_dir, paste0(etype, ".MATS.JCEC.txt"))
    if (!file.exists(f)) next
    df <- read.table(f, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "")
    names(df) <- make.unique(names(df))
    df$GeneID <- gsub('"', '', df$GeneID)

    # avg total reads per group: mean(IJC + SJC) across replicates
    ijc1 <- sapply(df$IJC_SAMPLE_1, csv_mean)
    sjc1 <- sapply(df$SJC_SAMPLE_1, csv_mean)
    ijc2 <- sapply(df$IJC_SAMPLE_2, csv_mean)
    sjc2 <- sapply(df$SJC_SAMPLE_2, csv_mean)
    avg_count1 <- ijc1 + sjc1
    avg_count2 <- ijc2 + sjc2

    # avg IncLevel (PSI) per group
    psi1 <- sapply(df$IncLevel1, csv_mean)
    psi2 <- sapply(df$IncLevel2, csv_mean)

    # SIGNIFICANCE filter (forces non-significant if failed): avg count >= 10 in
    # both groups AND avg PSI in [0.05, 0.95] in both groups.
    pass_filter <- !is.na(avg_count1) & avg_count1 >= 10 &
                   !is.na(avg_count2) & avg_count2 >= 10 &
                   !is.na(psi1) & psi1 >= 0.05 & psi1 <= 0.95 &
                   !is.na(psi2) & psi2 >= 0.05 & psi2 <= 0.95

    # RESTRICTED-UNIVERSE membership: count >= 10 and a PSI was computed. A
    # computed PSI means rMATS tested the junction, so it belongs in the tested
    # universe even if PSI is outside [0.05, 0.95]. Such events stay in the
    # restricted denominator but are still forced non-significant ("Negative")
    # by the pass_filter above -> they count as FN (GT-pos) or TN (GT-neg),
    # rather than being dropped (which would inflate recall).
    in_restricted <- !is.na(avg_count1) & avg_count1 >= 10 &
                     !is.na(avg_count2) & avg_count2 >= 10 &
                     !is.na(psi1) & !is.na(psi2)

    rmats_fdr_rows[[etype]] <- data.frame(
      event_type       = etype,
      ID               = as.character(df$ID),
      GeneID           = df$GeneID,
      FDR              = df$FDR,
      rmats_pass_filter = pass_filter,
      rmats_in_restricted = in_restricted,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  JCEC %s: %d events, %d pass sig filter, %d in restricted universe\n",
                etype, nrow(df), sum(pass_filter), sum(in_restricted)))
  }
  rmats_fdr <- bind_rows(rmats_fdr_rows)

  # join rMATS FDR with junction GT by event_type + ID
  jgt$ID <- as.character(jgt$ID)
  rmats_fdr$ID <- as.character(rmats_fdr$ID)
  jgt_merged <- merge(jgt,
                      rmats_fdr[, c("event_type", "ID", "FDR",
                                    "rmats_pass_filter", "rmats_in_restricted")],
                      by = c("event_type", "ID"), all.x = TRUE)
  # flag tested events BEFORE imputing FDR (NA means rMATS did not test this event)
  jgt_merged$rmats_tested <- !is.na(jgt_merged$FDR)
  jgt_merged$rmats_pass_filter[is.na(jgt_merged$rmats_pass_filter)] <- FALSE
  jgt_merged$rmats_in_restricted[is.na(jgt_merged$rmats_in_restricted)] <- FALSE
  # untested or filter-failing events cannot be significant
  jgt_merged$FDR[is.na(jgt_merged$FDR)] <- 1
  jgt_merged$FDR[!jgt_merged$rmats_pass_filter] <- 1

  jgt_by_gene    <- split(jgt_merged, jgt_merged$GeneID)
  jgt_eval_genes <- unique(jgt_merged$GeneID)

  # restricted=FALSE: full universe (all annotated events;
  #                   untested or filter-failing = non-significant)
  # restricted=TRUE:  restricted universe (only tested+filter-passing events)
  eval_jgt_universe <- function(tool_label, restricted) {
    rows_all <- list()
    for (thr in padj_thresholds) {
      rows_thr <- lapply(jgt_eval_genes, function(gene) {
        ev <- jgt_by_gene[[gene]]
        if (is.null(ev) || nrow(ev) == 0) return(NULL)
        if (restricted) {
          ev <- ev[ev$rmats_tested & ev$rmats_in_restricted, ]
          if (nrow(ev) == 0) return(NULL)
        }
        sig <- !is.na(ev$FDR) & ev$FDR < thr
        TP <- sum( sig &  ev$gt_positive, na.rm = TRUE)
        FP <- sum( sig & !ev$gt_positive, na.rm = TRUE)
        FN <- sum(!sig &  ev$gt_positive, na.rm = TRUE)
        TN <- sum(!sig & !ev$gt_positive, na.rm = TRUE)
        prec   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
        recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
        data.frame(
          gene = gene, sim_type = ev$sim_type[1], tool = tool_label,
          padj_thr = thr, n_gt_pos = sum(ev$gt_positive),
          n_gt_neg = sum(!ev$gt_positive), n_detected = sum(sig),
          TP = TP, FP = FP, FN = FN, TN = TN,
          precision = round(prec, 4), recall = round(recall, 4),
          f1 = round(f1_safe(prec, recall), 4),
          stringsAsFactors = FALSE
        )
      })
      rows_all <- c(rows_all, rows_thr)
    }
    bind_rows(rows_all[!sapply(rows_all, is.null)])
  }

  rmats_junctiongt_per_gene            <- eval_jgt_universe("rMATS_junctionGT",            restricted = FALSE)
  rmats_junctiongt_restricted_per_gene <- eval_jgt_universe("rMATS_junctionGT_restricted",  restricted = TRUE)
  rmats_junctiongt_summary             <- build_summary(rmats_junctiongt_per_gene)
  rmats_junctiongt_restricted_summary  <- build_summary(rmats_junctiongt_restricted_per_gene)

  for (thr in padj_thresholds) {
    write.table(
      rmats_junctiongt_per_gene[rmats_junctiongt_per_gene$padj_thr == thr, ],
      file.path(out_dir, sprintf("rmats_junctiongt_per_gene_padj%.2f.txt", thr)),
      sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(
      rmats_junctiongt_restricted_per_gene[
        rmats_junctiongt_restricted_per_gene$padj_thr == thr, ],
      file.path(out_dir,
                sprintf("rmats_junctiongt_restricted_per_gene_padj%.2f.txt", thr)),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  write.table(rmats_junctiongt_summary,
              file.path(out_dir, "rmats_junctiongt_summary_by_simtype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(rmats_junctiongt_restricted_summary,
              file.path(out_dir, "rmats_junctiongt_restricted_summary_by_simtype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # --- rMATS breakdown by event type (SE/A3SS/A5SS/RI + ALL) ----------------
  # Aggregates event-level confusion (not per-gene) directly from jgt_merged,
  # which already excludes GT-untestable events and has FDR imputed to 1 for
  # untested / filter-failing events. Pooled over all sim types, so FP includes
  # Background/DGE false positives.
  rmats_etype_summary <- function(restricted) {
    rows <- list()
    for (thr in padj_thresholds) {
      d <- jgt_merged
      if (restricted) d <- d[d$rmats_tested & d$rmats_in_restricted, ]
      for (et in c("SE", "A3SS", "A5SS", "RI", "ALL")) {
        dd  <- if (et == "ALL") d else d[d$event_type == et, ]
        if (nrow(dd) == 0) next
        sig <- !is.na(dd$FDR) & dd$FDR < thr
        TP  <- sum( sig &  dd$gt_positive, na.rm = TRUE)
        FP  <- sum( sig & !dd$gt_positive, na.rm = TRUE)
        FN  <- sum(!sig &  dd$gt_positive, na.rm = TRUE)
        TN  <- sum(!sig & !dd$gt_positive, na.rm = TRUE)
        prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
        rec  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
        rows[[length(rows) + 1]] <- data.frame(
          event_type = et, padj_thr = thr, n_events = nrow(dd),
          total_TP = TP, total_FP = FP, total_FN = FN, total_TN = TN,
          micro_precision = round(prec, 4), micro_recall = round(rec, 4),
          micro_f1 = round(f1_safe(prec, rec), 4),
          stringsAsFactors = FALSE)
      }
    }
    bind_rows(rows)
  }
  rmats_etype_summary_full       <- rmats_etype_summary(FALSE)
  rmats_etype_summary_restricted <- rmats_etype_summary(TRUE)
  write.table(rmats_etype_summary_full,
              file.path(out_dir, "rmats_junctiongt_summary_by_eventtype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(rmats_etype_summary_restricted,
              file.path(out_dir, "rmats_junctiongt_restricted_summary_by_eventtype.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\n=== rMATS Summary (junction GT, full universe + filter) ===\n")
  print(as.data.frame(rmats_junctiongt_summary %>%
    select(sim_type, padj_thr, n_genes,
           micro_precision, micro_recall, micro_f1,
           total_TP, total_FP, total_FN)))
  cat("\n=== rMATS Summary (junction GT, restricted to tested+filtered events) ===\n")
  print(as.data.frame(rmats_junctiongt_restricted_summary %>%
    select(sim_type, padj_thr, n_genes,
           micro_precision, micro_recall, micro_f1,
           total_TP, total_FP, total_FN)))
  cat("\n=== rMATS by event type (restricted, all padj) ===\n")
  print(as.data.frame(rmats_etype_summary_restricted))
} else {
  cat(sprintf("\n[skip] junction GT file not found: %s\n", junction_gt_file))
}

# Saturn evaluation (after junction GT)

saturn_per_gene <- evaluate_tool(saturn_by_gene, "Saturn")
saturn_restricted_per_gene <- evaluate_tool(
  saturn_by_gene, "Saturn_restricted",
  testable_by_gene = saturn_testable_by_gene,
  restrict_pos = TRUE
)

#  write Saturn per-gene results

for (thr in padj_thresholds) {
  write.table(
    saturn_per_gene[saturn_per_gene$padj_thr == thr, ],
    file.path(out_dir, sprintf("saturn_per_gene_padj%.2f.txt", thr)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    saturn_restricted_per_gene[saturn_restricted_per_gene$padj_thr == thr, ],
    file.path(out_dir, sprintf("saturn_restricted_per_gene_padj%.2f.txt", thr)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

#  Saturn summaries and print

saturn_summary            <- build_summary(saturn_per_gene)
saturn_restricted_summary <- build_summary(saturn_restricted_per_gene)

write.table(saturn_summary,
            file.path(out_dir, "saturn_summary_by_simtype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(saturn_restricted_summary,
            file.path(out_dir, "saturn_restricted_summary_by_simtype.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

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

cat(sprintf("\nResults written to: %s\n", out_dir))
