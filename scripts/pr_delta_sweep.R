#!/usr/bin/env Rscript
#
# scripts/pr_curves_models.R
#
# PR curves for the four statistical models, UNIT LEVEL ONLY.
#
# Deliberately separate from pr_curves_three_levels_gtrule.R. That script
# compares TOOLS with genuinely different units, which is why it carries
# per-tool ground truths (GrASE -> GT_rule_bipartition, MAJIQ -> GT_rule_lsv,
# DEXSeq -> part-level gtI), a restricted/full universe split, and three
# levels. None of that applies here: the four models share ONE unit, ONE
# implicated-transcript set and ONE FDR architecture, so transcript- and
# gene-level curves would differ between them only through which sides get
# called -- the unit level shows that directly and without the indirection.
#
# Ground truth is read from the EXTERNAL GT_rule table written by
# infer_bipartition_gt.R; nothing is re-derived here.
#   results/gt/bipartition_merged_gt.stranded.txt   (gt_positive column)
#
# Everything is held constant except the model: same merged stranded counts,
# and within an arm the same raw phi (EBapprox and EBmap differ only in how it
# is moderated; MLE and wilcoxon ignore it). Across arms phi necessarily
# differs -- the arms enumerate different events.
#
# Usage: Rscript scripts/pr_curves_models.R
#   inputs : bipartition.merged{,.TSSTTS}.stranded.test.modelcomp/
#            test_bipartition.merged_<model>.annotated.txt
#   outputs: results/eval_metricfair/pr_models.GT_rule.txt
#            plots/pr_models.GT_rule.png
suppressMessages({library(dplyr)})

BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
OUT  <- file.path(BASE, "results/eval_metricfair")
GT   <- file.path(BASE, "results/gt/bipartition_merged_gt.stranded.txt")
DIRS <- file.path(BASE, c("bipartition.merged.stranded.test.modelcomp",
                          "bipartition.merged.TSSTTS.stranded.test.modelcomp"))
MODELS <- c(EBapprox = "betabinom_EBapprox")
COL <- c(EBapprox="#1f78b4", EBmap="#33a02c", MLE="#e31a1c", wilcoxon="#6a3d9a")
THR <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2)
DELTAS <- c(-Inf, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2)
PADJ   <- 0.01

## Two-stage rule, identical to the tool sweep: within-gene BH RECOMPUTED from
## raw p-values, then max with the gene-level padj. Must not use the stored
## `padj`, which pvalueAdjustment.R populates only for genes clearing a
## hardcoded alpha at stage 1 and leaves NA elsewhere -- combining that with
## padj_gene screens genes twice and truncates the sweep above that alpha.
two_stage <- function(pv, gene, padj_gene) {
  within <- ave(pv, gene, FUN = function(x) p.adjust(x, method = "BH"))
  pmax(within, padj_gene)
}

## Simulation category per GENE, assigned exactly as pr_curves_three_levels_gtrule.R
## does (from simulate.rda, not from the GT table's sim_type -- that column labels
## only DTE and DTU, so Background and DGE would be invisible).
load(file.path(BASE, "swimdown/simulate/data/simulate.rda"))
get_st <- function(g) ifelse(g %in% dte.genes, "DTE", ifelse(g %in% dtu.genes, "DTU",
                     ifelse(g %in% dge.genes, "DGE", "Background")))

gt <- read.table(GT, header = TRUE, sep = "\t", quote = "", comment.char = "",
                 stringsAsFactors = FALSE)
gt$uid <- paste(gt$src, gt$gene, gt$event, gt$comparison, sep = "|")
pos <- setNames(gt$gt_positive %in% c(TRUE, "TRUE"), gt$uid)
## gt$gene is version-tagged (ENSG...\.NN), exactly as simulate.rda stores
## dte/dtu/dge.genes -- match on it directly, do not strip the version.
## Category is a property of the GENE, so derive it from the uid's gene field
## (uid = src|gene|event|comparison) rather than from the GT table, which lists
## only DTE/DTU genes.
gene_of_uid <- function(u) sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", u)
cat(sprintf("GT_rule table: %d units (%d positive)\n", length(pos), sum(pos)))

## Universe. The models do NOT all emit a row per event: betabinom_MLE returns
## nothing where its fit fails (22,753 rows against 29,097 for the others on the
## internal arm), so a per-model denominator would shrink MLE's universe and
## INFLATE its recall -- the same artefact the tool sweep documents for the
## restricted universe. Both views are therefore reported:
##   full       common universe = every GT unit any model emitted; a model that
##              produced no row for a unit simply does not call it. This is the
##              comparable view and the one to read across models.
##   restricted per-model universe = only the units that model emitted. A
##              within-model view; must NOT be read across models.
mdl_data <- list()
for (nm in names(MODELS)) {
  fs <- file.path(DIRS, sprintf("test_bipartition.merged_%s.annotated.txt", MODELS[[nm]]))
  if (!all(file.exists(fs))) { cat(sprintf("  %-9s missing an arm, skipping\n", nm)); next }
  n  <- vapply(fs, function(f) nrow(read.table(f, header=TRUE, sep="\t", quote="",
               comment.char="", stringsAsFactors=FALSE)), integer(1))
  d  <- bind_rows(lapply(fs, function(f) read.table(f, header=TRUE, sep="\t", quote="",
        comment.char="", stringsAsFactors=FALSE)[, c("gene","event","comparison",
        "padj","padj_gene","p.value","lfc_diff_net")]))
  d$src <- rep(c("internal","TSSTTS"), times = n)
  d$uid <- paste(d$src, d$gene, d$event, d$comparison, sep = "|")
  d$stat <- two_stage(d$p.value, d$gene, d$padj_gene)
  ## Testability must be recorded BEFORE stat is coerced to 1: a unit whose
  ## p-value is NA (read floor, independent filtering) is a non-call, not a
  ## unit the model could have called. The restricted universe keeps only the
  ## testable ones, so recall is not charged for units GrASE can never reach.
  d$testable <- !is.na(d$stat)
  d$stat[is.na(d$stat)] <- 1
  ## A unit absent from the GT table is a NEGATIVE, not a unit to drop -- this
  ## is the sweep's rule (pr_curves_three_levels_gtrule.R, "Any unit not in the
  ## GT table is a negative"). Dropping them would delete every Background and
  ## DGE unit, i.e. exactly the null units whose calls are the FPs we report.
  d$truth <- pos[d$uid]
  d$truth[is.na(d$truth)] <- FALSE
  mdl_data[[nm]] <- d
  cat(sprintf("  %-9s %6d units emitted, %5d GT-positive\n", nm, nrow(d), sum(d$truth)))
}

universe <- Reduce(union, lapply(mdl_data, `[[`, "uid"))
u_truth  <- pos[universe]; u_truth[is.na(u_truth)] <- FALSE
cat_by_uid <- setNames(get_st(gene_of_uid(universe)), universe)
cat(sprintf("common universe: %d units, %d positive (union across models)\n",
            length(universe), sum(u_truth)))


## --- Delta sweep on the post-hoc lfc_diff_net filter -----------------------
## Same scoring path as pr_curves_models.R (GT_rule, unit level, stranded
## merged, nested-BH two-stage), model fixed to EBapprox, sweeping the
## lfc_diff_net threshold. delta = -Inf is the unfiltered BASELINE.
##
## Two universes, matching evaluate_posthoc_lfc_filter.R:
##   full       every emitted unit; units the model cannot test stay in the
##              recall denominator as non-calls.
##   restricted GrASE-testable units only (a computable p-value), so recall is
##              not charged for units removed by the read floor or by
##              independent filtering.
f1 <- function(p, r) ifelse(p + r > 0, 2 * p * r / (p + r), 0)
d   <- mdl_data[["EBapprox"]]
stopifnot(!is.null(d))

lfc  <- setNames(d$lfc_diff_net, d$uid)[universe]
st   <- setNames(rep(1, length(universe)), universe);        st[d$uid]   <- d$stat
tst  <- setNames(rep(FALSE, length(universe)), universe);    tst[d$uid]  <- d$testable
ucat <- cat_by_uid[universe]
cat(sprintf("universe %d units; GrASE-testable %d (%.1f%%)\n",
            length(universe), sum(tst), 100 * sum(tst) / length(universe)))

CATS <- c("ALL", "DTE", "DTU", "DGE", "Background")
rows <- list()
for (uni in c("full", "restricted")) {
  keep_u <- if (uni == "full") rep(TRUE, length(universe)) else tst
  for (dl in DELTAS) {
    ## a unit passes only if significant AND clears the filter; an NA
    ## lfc_diff_net never clears it (matches the main pipeline).
    sig <- keep_u & (st < PADJ) & !is.na(lfc) & lfc > dl
    for (cc in CATS) {
      k <- if (cc == "ALL") which(keep_u) else which(ucat == cc & keep_u)
      if (!length(k)) next
      tp <- sum(sig[k] & u_truth[k]); fp <- sum(sig[k] & !u_truth[k])
      np <- sum(u_truth[k])
      rows[[length(rows)+1]] <- data.frame(
        universe = uni, category = cc, delta = dl, padj_thr = PADJ,
        TP = tp, FP = fp, n_pos = np,
        precision = tp/max(tp+fp,1), recall = tp/max(np,1),
        F1 = f1(tp/max(tp+fp,1), tp/max(np,1)))
    }
  }
}
tab <- bind_rows(rows)
write.table(tab, file.path(OUT, "delta_sweep.GT_rule.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)

con <- file(file.path(OUT, "delta_sweep.GT_rule.table.txt"), "w")
for (uni in c("full", "restricted")) for (cc in CATS) {
  x <- tab[tab$universe == uni & tab$category == cc, ]
  if (!nrow(x)) next
  writeLines(sprintf("\n=== %s | %s universe | GT_rule unit level | padj %.2f ===",
                     cc, uni, PADJ), con)
  writeLines(sprintf("%8s %8s %8s %9s %9s %8s", "delta", "TP", "FP",
                     "precision", "recall", "F1"), con)
  for (i in seq_len(nrow(x)))
    writeLines(sprintf("%8s %8d %8d %9.4f %9.4f %8.4f",
                       ifelse(is.infinite(x$delta[i]), "none", sprintf("%.1f", x$delta[i])),
                       x$TP[i], x$FP[i], x$precision[i], x$recall[i], x$F1[i]), con)
}
close(con)

## --- Figures: one per universe ---------------------------------------------
## Panel 1 ALL, panel 2 DTE and DTU TOGETHER (colour-coded, shared axes with
## panel 1), panel 3 null-gene FP counts by delta.
COL_CAT <- c(ALL = "#1f78b4", DTE = "#e66101", DTU = "#1b7837")
PANEL_CATS <- c("ALL", "DTE", "DTU")

draw_traj <- function(x, bl, col, lab_delta) {
  lines(x$recall, x$precision, col = col, lwd = 2)
  points(x$recall, x$precision, col = col, pch = 19, cex = 0.7)
  if (lab_delta) text(x$recall, x$precision, labels = sprintf("%.1f", x$delta),
                      pos = 4, cex = 0.55, col = col)
  if (nrow(bl)) {
    d0 <- x[x$delta == 0, ]
    if (nrow(d0)) segments(bl$recall, bl$precision, d0$recall, d0$precision,
                           col = col, lty = 2)
    points(bl$recall, bl$precision, pch = 0, cex = 1.5, col = col, lwd = 2)
  }
}

for (uni in c("full", "restricted")) {
  tu <- tab[tab$universe == uni, ]
  ## Common axes across both PR panels, over every plotted point incl. baseline.
  pp <- tu[tu$category %in% PANEL_CATS, ]
  XL <- range(pp$recall); YL <- range(pp$precision)

  ## H includes the outer top margin reserved by oma below for the overall
  ## title. Without oma, mtext(outer=TRUE) has no strip to draw in and a
  ## negative `line` pushes the title down onto the panel headings.
  W <- 3.6 * 3; H <- 4.3
for (dev_i in 1:1) {
    if (dev_i == 1)
      pdf(file.path(BASE, sprintf("plots/delta_sweep.GT_rule.%s.pdf", uni)),
          width = W, height = H)
    else
      png(file.path(BASE, sprintf("plots/delta_sweep.GT_rule.%s.png", uni)),
          width = W, height = H, units = "in", res = 300)
    par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2.4, 0))

    ## panel 1 -- ALL
    plot(NA, xlim = XL, ylim = YL, xlab = "recall", ylab = "precision", main = "ALL")
    draw_traj(tu[tu$category == "ALL" &  is.finite(tu$delta), ],
              tu[tu$category == "ALL" & !is.finite(tu$delta), ], COL_CAT[["ALL"]], TRUE)
    legend("bottomleft", c("lfc_diff_net delta", "no filter (baseline)"),
           pch = c(19, 0), col = COL_CAT[["ALL"]], lty = c(1, NA), lwd = 2,
           bty = "n", cex = 0.7)

    ## panel 2 -- DTE and DTU together
    plot(NA, xlim = XL, ylim = YL, xlab = "recall", ylab = "precision",
         main = "DTE and DTU")
    for (cc in c("DTE", "DTU"))
      draw_traj(tu[tu$category == cc &  is.finite(tu$delta), ],
                tu[tu$category == cc & !is.finite(tu$delta), ], COL_CAT[[cc]], TRUE)
    legend("bottomleft", c("DTE", "DTU", "no filter (baseline)"),
           pch = c(19, 19, 0), col = c(COL_CAT[["DTE"]], COL_CAT[["DTU"]], "black"),
           lty = c(1, 1, NA), lwd = 2, bty = "n", cex = 0.7)

    ## panel 3 -- null-gene false positives by delta, baseline leftmost
    nb  <- tu[tu$category %in% c("Background", "DGE"), ]
    dls <- sort(unique(nb$delta))                  # -Inf sorts first = baseline
    m   <- sapply(dls, function(dd)
            c(Background = sum(nb$FP[nb$delta == dd & nb$category == "Background"]),
              DGE        = sum(nb$FP[nb$delta == dd & nb$category == "DGE"])))
    colnames(m) <- ifelse(is.finite(dls), sprintf("%.1f", dls), "none")
    bp <- barplot(m, col = c("#bdbdbd", "#636363"), las = 2, cex.names = 0.7,
                  ylab = "null-gene false positives",
                  main = "Background + DGE (leftmost = no filter)",
                  ylim = c(0, max(colSums(m)) * 1.25))
    text(bp, colSums(m), labels = colSums(m), pos = 3, cex = 0.7, xpd = NA)
    legend("topright", c("Background", "DGE"), fill = c("#bdbdbd", "#636363"),
           bty = "n", cex = 0.8)

    mtext(sprintf("lfc_diff_net delta sweep, EBapprox, stranded merged (GT_rule, unit level, %s universe)", uni),
          side = 3, line = 0.6, outer = TRUE, cex = 0.85, font = 2)
    invisible(dev.off())
  }
}

cat("wrote:\n  results/eval_metricfair/delta_sweep.GT_rule.txt\n",
    " results/eval_metricfair/delta_sweep.GT_rule.table.txt\n",
    " plots/delta_sweep.GT_rule.{full,restricted}.{pdf,png}\n", sep="")
