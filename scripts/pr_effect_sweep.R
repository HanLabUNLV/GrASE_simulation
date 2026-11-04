#!/usr/bin/env Rscript
#
# scripts/pr_effect_sweep.R
#
# Sweep the effect-size cutoff of the call rule, comparing four gates: raw |delta_pi|,
# per-base |delta_pi| (length-normalized, grase::annotate_pi_perbase; raw where a side has
# no length, i.e. use_perbase = TRUE), raw AND per-base (use_perbase = "both"), and
# |effect_size| (the log odds ratio of the fitted model), on the stranded merged
# EBapprox results, scored under GT_rule at the unit level -- the same inputs,
# ground truth, universes and two-stage padj as pr_delta_sweep.R.
#
# A unit is called at cutoff c iff grase::is_significant() passes it with the two-stage
# padj < 0.01, lfc_diff_net > 0 and that one measure >= c (the other effect gate off, no
# read floor). One cutoff is applied to every side; the pipeline's separate min_dpi_sj for
# junction-sourced sides is not modelled. The shipped rule (is_significant() defaults:
# both delta_pi scales >= 0.05, |effect_size| >= 0.5, d_support >= 10) is scored alongside
# as a reference point.
#
# Usage: Rscript scripts/pr_effect_sweep.R
#   outputs: results/eval_metricfair/effect_sweep.GT_rule.txt
#            results/eval_metricfair/effect_sweep.GT_rule.table.txt
#            plots/effect_sweep.GT_rule.{full,restricted}.pdf
suppressMessages({library(dplyr); library(grase)})

BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
OUT  <- file.path(BASE, "results/eval_metricfair")
GT   <- file.path(BASE, "results/gt/bipartition_merged_gt.stranded.txt")
DIRS <- file.path(BASE, c("bipartition.merged.stranded.test.EBapprox",
                          "bipartition.merged.TSSTTS.stranded.test.EBapprox"))
MODEL <- "betabinom_EBapprox"
GFF   <- file.path(BASE, "dexseq.gff")   # per-gene DEXSeq GFFs of the simulation annotation
READ_LEN <- 100                           # 2x100 reads; effective exonic-part lengths of per-base pi
PADJ  <- 0.01
DPI_CUTS <- c(0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3)
CUTS  <- list(delta_pi = DPI_CUTS, delta_pi_perbase = DPI_CUTS, delta_pi_both = DPI_CUTS,
              effect_size = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3))

## Two-stage rule, identical to pr_delta_sweep.R / pr_curves_three_levels_gtrule.R.
two_stage <- function(pv, gene, padj_gene, src) {
  out <- rep(NA_real_, length(pv))
  for (s in unique(src)) {
    i  <- which(src == s)
    w  <- ave(pv[i], gene[i], FUN = function(x) p.adjust(x, method = "BH"))
    gp <- padj_gene[i][!duplicated(gene[i])]
    gp <- sort(gp[!is.na(gp)])
    g  <- unique(gp)
    cnt <- findInterval(g, gp)
    x  <- w * length(gp)
    j  <- findInterval(x, g * cnt, left.open = TRUE) + 1
    above <- ifelse(j <= length(g), g[pmin(j, length(g))], Inf)
    below <- ifelse(j > 1, x / cnt[pmax(j - 1, 1)], Inf)
    out[i] <- pmax(padj_gene[i], pmin(above, below))
  }
  out
}

load(file.path(BASE, "swimdown/simulate/data/simulate.rda"))
get_st <- function(g) ifelse(g %in% dte.genes, "DTE", ifelse(g %in% dtu.genes, "DTU",
                     ifelse(g %in% dge.genes, "DGE", "Background")))

gt <- read.table(GT, header = TRUE, sep = "\t", quote = "", comment.char = "",
                 stringsAsFactors = FALSE)
gt$uid <- paste(gt$src, gt$gene, gt$event, gt$comparison, sep = "|")
pos <- setNames(gt$gt_positive %in% c(TRUE, "TRUE"), gt$uid)

fs <- file.path(DIRS, sprintf("test_bipartition.merged_%s.annotated.txt", MODEL))
rd <- lapply(fs, function(f) read.table(f, header = TRUE, sep = "\t", quote = "",
        comment.char = "", stringsAsFactors = FALSE)[, c("gene","event","comparison",
        "padj","padj_gene","p.value","lfc_diff_net","delta_pi","effect_size","d_support",
        "pi_trt", intersect(c("pi_ctrl", "pi_ref"), names(read.table(f, header = TRUE, sep = "\t", nrows = 1, quote = "", comment.char = ""))),
        "ref_ex_part","setdiff1","setdiff2")])
d <- bind_rows(rd)
if ("pi_ref" %in% names(d)) d <- rename(d, pi_ctrl = pi_ref)   # results written before the rename
d <- annotate_pi_perbase(d, GFF, READ_LEN)
cat(sprintf("per-base delta_pi resolved for %d of %d rows\n", sum(!is.na(d$delta_pi_perbase)), nrow(d)))
## the gate of each measure as an is_significant() call at cutoff ct
GATE <- list(
  delta_pi         = function(ct) is_significant(d, PADJ, 0, ct, ct, 0, padj = d$stat, use_perbase = FALSE, min_effect = 0),
  delta_pi_perbase = function(ct) is_significant(d, PADJ, 0, ct, ct, 0, padj = d$stat, use_perbase = TRUE, min_effect = 0),
  delta_pi_both    = function(ct) is_significant(d, PADJ, 0, ct, ct, 0, padj = d$stat, use_perbase = "both", min_effect = 0),
  effect_size      = function(ct) is_significant(d, PADJ, 0, 0, 0, 0, padj = d$stat, min_effect = ct))
d$src <- rep(c("internal","TSSTTS"), times = vapply(rd, nrow, integer(1)))
d$uid <- paste(d$src, d$gene, d$event, d$comparison, sep = "|")
d$stat <- two_stage(d$p.value, d$gene, d$padj_gene, d$src)
d$testable <- !is.na(d$stat)
d$stat[is.na(d$stat)] <- 1
d$truth <- pos[d$uid]; d$truth[is.na(d$truth)] <- FALSE
d$cat <- get_st(d$gene)
d$pipeline <- is_significant(d, PADJ, 0, 0.05, 0.05, 10, padj = d$stat)   # the shipped defaults
cat(sprintf("%d units, %d GT-positive, %d testable\n", nrow(d), sum(d$truth), sum(d$testable)))

f1 <- function(p, r) ifelse(p + r > 0, 2 * p * r / (p + r), 0)
CATS <- c("ALL", "DTE", "DTU", "DGE", "Background")
score <- function(sig, keep, uni, measure, cut) {
  bind_rows(lapply(CATS, function(cc) {
    k <- if (cc == "ALL") which(keep) else which(d$cat == cc & keep)
    tp <- sum(sig[k] & d$truth[k]); fp <- sum(sig[k] & !d$truth[k]); np <- sum(d$truth[k])
    p <- tp / max(tp + fp, 1); r <- tp / max(np, 1)
    data.frame(universe = uni, category = cc, measure = measure, cutoff = cut,
               TP = tp, FP = fp, n_pos = np, precision = p, recall = r, F1 = f1(p, r))
  }))
}
rows <- list()
for (uni in c("full", "restricted")) {
  keep <- if (uni == "full") rep(TRUE, nrow(d)) else d$testable
  for (m in names(CUTS)) for (ct in CUTS[[m]])
    rows[[length(rows) + 1]] <- score(keep & GATE[[m]](ct), keep, uni, m, ct)
  rows[[length(rows) + 1]] <- score(keep & d$pipeline, keep, uni, "pipeline_significant", NA)
}
tab <- bind_rows(rows)
write.table(tab, file.path(OUT, "effect_sweep.GT_rule.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

con <- file(file.path(OUT, "effect_sweep.GT_rule.table.txt"), "w")
for (uni in c("full", "restricted")) for (cc in CATS) {
  x <- tab[tab$universe == uni & tab$category == cc, ]
  writeLines(sprintf("\n=== %s | %s universe | GT_rule unit level | padj %.2f ===", cc, uni, PADJ), con)
  writeLines(sprintf("%-22s %7s %7s %7s %9s %9s %8s", "measure", "cutoff", "TP", "FP", "precision", "recall", "F1"), con)
  for (i in seq_len(nrow(x)))
    writeLines(sprintf("%-22s %7s %7d %7d %9.4f %9.4f %8.4f", x$measure[i],
                       ifelse(is.na(x$cutoff[i]), "-", format(x$cutoff[i])),
                       x$TP[i], x$FP[i], x$precision[i], x$recall[i], x$F1[i]), con)
}
close(con)

## Figures: one per universe; panels ALL, DTE, DTU. Each measure is a PR trajectory
## over its cutoffs; the pipeline's own calls are the open square.
COL_M <- c(delta_pi = "#1f78b4", delta_pi_perbase = "#33a02c", delta_pi_both = "#6a3d9a", effect_size = "#e66101")
for (uni in c("full", "restricted")) {
  tu <- tab[tab$universe == uni, ]
  pp <- tu[tu$category %in% c("ALL", "DTE", "DTU"), ]
  ## a cutoff that calls nothing has precision 0 by convention; keep it off the axis range
  XL <- c(0, max(pp$recall)); YL <- range(pp$precision[pp$TP > 0])
  pdf(file.path(BASE, sprintf("plots/effect_sweep.GT_rule.%s.pdf", uni)), width = 10.8, height = 4.3)
  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2.4, 0))
  for (cc in c("ALL", "DTE", "DTU")) {
    plot(NA, xlim = XL, ylim = YL, xlab = "recall", ylab = "precision", main = cc)
    for (m in names(CUTS)) {
      x <- tu[tu$category == cc & tu$measure == m, ]
      lines(x$recall, x$precision, col = COL_M[[m]], lwd = 2)
      points(x$recall, x$precision, col = COL_M[[m]], pch = 19, cex = 0.7)
      text(x$recall, x$precision, labels = format(x$cutoff), pos = 4, cex = 0.55, col = COL_M[[m]])
    }
    b <- tu[tu$category == cc & tu$measure == "pipeline_significant", ]
    points(b$recall, b$precision, pch = 0, cex = 1.5, lwd = 2)
    if (cc == "ALL")
      legend("bottomleft", c("raw |delta_pi|", "per-base |delta_pi|", "raw and per-base |delta_pi|",
                             "|effect_size|", "pipeline calls"),
             pch = c(19, 19, 19, 19, 0), col = c(COL_M, "black"), lty = c(1, 1, 1, 1, NA), lwd = 2, bty = "n", cex = 0.7)
  }
  mtext(sprintf("Effect-size cutoff sweep, EBapprox, stranded merged (GT_rule, unit level, padj %.2f, %s universe)", PADJ, uni),
        side = 3, line = 0.6, outer = TRUE, cex = 0.85, font = 2)
  invisible(dev.off())
}
cat("wrote:\n  results/eval_metricfair/effect_sweep.GT_rule.txt\n",
    " results/eval_metricfair/effect_sweep.GT_rule.table.txt\n",
    " plots/effect_sweep.GT_rule.{full,restricted}.pdf\n", sep = "")
