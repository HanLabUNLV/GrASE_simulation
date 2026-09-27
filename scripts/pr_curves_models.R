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
MODELS <- c(EBapprox = "betabinom_EBapprox", EBmap = "betabinom_EBmap",
            MLE = "betabinom_MLE", wilcoxon = "wilcoxon")
COL <- c(EBapprox="#1f78b4", EBmap="#33a02c", MLE="#e31a1c", wilcoxon="#6a3d9a")
THR <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2)
DELTA <- 0

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
  ok <- !is.na(d$lfc_diff_net) & d$lfc_diff_net > DELTA
  d$stat[is.na(d$stat) | !ok] <- 1
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

rows <- list()
for (nm in names(mdl_data)) {
  d <- mdl_data[[nm]]
  st <- setNames(rep(1, length(universe)), universe)   # absent -> never called
  st[d$uid] <- d$stat
  ucat <- cat_by_uid[universe]
  for (t in THR) {
    sig <- st < t
    tp <- sum(sig & u_truth); fp <- sum(sig & !u_truth)
    rows[[length(rows)+1]] <- data.frame(model = nm, universe = "full", category = "ALL",
      thr = t, TP = tp, FP = fp, n_pos = sum(u_truth),
      precision = tp/max(tp+fp,1), recall = tp/max(sum(u_truth),1))
    ## per simulation category. Background and DGE genes carry no GT positive,
    ## so their rows are pure false-positive counts -- which is what the results
    ## text quotes as "null FPs".
    for (cc in c("DTU","DTE","DGE","Background")) {
      k <- which(ucat == cc)
      if (!length(k)) next
      tpc <- sum(sig[k] & u_truth[k]); fpc <- sum(sig[k] & !u_truth[k]); npc <- sum(u_truth[k])
      rows[[length(rows)+1]] <- data.frame(model = nm, universe = "full", category = cc,
        thr = t, TP = tpc, FP = fpc, n_pos = npc,
        precision = tpc/max(tpc+fpc,1), recall = tpc/max(npc,1))
    }
    sig2 <- d$stat < t
    tp2 <- sum(sig2 & d$truth); fp2 <- sum(sig2 & !d$truth)
    rows[[length(rows)+1]] <- data.frame(model = nm, universe = "restricted", category = "ALL",
      thr = t, TP = tp2, FP = fp2, n_pos = sum(d$truth),
      precision = tp2/max(tp2+fp2,1), recall = tp2/max(sum(d$truth),1))
  }
}
tab <- bind_rows(rows)
if (!nrow(tab)) stop("no models scored -- are the modelcomp arms present?")
write.table(tab, file.path(OUT, "pr_models.GT_rule.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)

PANELS <- c("ALL", "DTE", "DTU")

## --- Figure: PR per category (a) + null-gene FP counts (b) ------------------
for (uni in c("full", "restricted")) {
  ## Only the ALL panel is defined on the restricted universe: the per-category
  ## rows are scored on the common universe only (see the universe note above).
  cats <- if (uni == "full") PANELS else "ALL"
  ## Vector PDF is the figure of record (these are manuscript panels and get
  ## scaled for print); the PNG is a convenience copy at 300 dpi. The former
  ## single PNG at res=130 gave ~2.9 in per panel, which is too coarse to
  ## enlarge. Dimensions are in INCHES here, not pixels.
  ## H includes the oma strip reserved below for the overall title.
  W <- 3.4 * (length(cats) + 1); H <- 4.2
  for (dev_i in 1:2) {
    if (dev_i == 1)
      pdf(file.path(BASE, sprintf("plots/pr_models.GT_rule.%s.pdf", uni)),
          width = W, height = H)
    else
      png(file.path(BASE, sprintf("plots/pr_models.GT_rule.%s.png", uni)),
          width = W, height = H, units = "in", res = 300)
  par(mfrow = c(1, length(cats) + 1), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2.4, 0))

  ## Common axes across the category panels, for the same reason as the delta
  ## sweep: per-panel ranges make equal screen positions mean unequal values.
  pp <- tab[tab$universe == uni & tab$category %in% cats, ]
  XL <- range(pp$recall); YL <- range(pp$precision)
  for (cc in cats) {
    pf <- tab[tab$universe == uni & tab$category == cc, ]
    plot(NA, xlim = XL, ylim = YL,
         xlab = "recall", ylab = "precision", main = cc)
    for (nm in unique(pf$model)) {
      x <- pf[pf$model == nm, ]; x <- x[order(x$recall), ]
      lines(x$recall, x$precision, col = COL[[nm]], lwd = 2)
      points(x$recall, x$precision, col = COL[[nm]], pch = 19, cex = 0.7)
    }
    if (cc == cats[1])
      legend("bottomleft", names(COL), col = unname(COL), lwd = 2, bty = "n", cex = 0.8)
  }

  ## (b) false positives on null genes at the operating threshold. Background
  ## and DGE genes carry no GT positive, so every call there is an FP.
  nb <- tab[tab$universe == uni & tab$category %in% c("Background", "DGE") &
            tab$thr == 0.01, ]
  if (nrow(nb)) {
    m  <- sapply(names(COL), function(nm)
            c(Background = sum(nb$FP[nb$model == nm & nb$category == "Background"]),
              DGE        = sum(nb$FP[nb$model == nm & nb$category == "DGE"])))
    bp <- barplot(m, beside = FALSE, col = c("#bdbdbd", "#636363"),
                  ylab = "false positives on null genes", main = "Background + DGE (padj 0.01)",
                  las = 2, cex.names = 0.8, ylim = c(0, max(colSums(m)) * 1.25))
    text(bp, colSums(m), labels = colSums(m), pos = 3, cex = 0.8, xpd = NA)
    legend("topleft", c("Background", "DGE"), fill = c("#bdbdbd", "#636363"),
           bty = "n", cex = 0.8)
  } else plot.new()

  mtext(sprintf("Statistical models, merged stranded unit (GT_rule, unit level, %s universe)", uni),
        side = 3, line = 0.6, outer = TRUE, cex = 0.85, font = 2)
  invisible(dev.off())
  }
}

## --- Tables ----------------------------------------------------------------
f1 <- function(p, r) ifelse(p + r > 0, 2 * p * r / (p + r), 0)
tab$F1 <- f1(tab$precision, tab$recall)
write.table(tab, file.path(OUT, "pr_models.GT_rule.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## Operating-point table at padj 0.01, one block per category, both universes.
con <- file(file.path(OUT, "pr_models.GT_rule.table.txt"), "w")
for (uni in c("full", "restricted")) {
  for (cc in c("ALL", "DTE", "DTU", "Background", "DGE")) {
    x <- tab[tab$universe == uni & tab$category == cc & tab$thr == 0.01, ]
    if (!nrow(x)) next
    writeLines(sprintf("\n=== %s | %s universe | padj 0.01 ===", cc, uni), con)
    writeLines(sprintf("%-9s %7s %7s %8s %9s %9s %7s", "model", "TP", "FP",
                       "n_pos", "precision", "recall", "F1"), con)
    for (i in seq_len(nrow(x)))
      writeLines(sprintf("%-9s %7d %7d %8d %9.4f %9.4f %7.4f", x$model[i], x$TP[i],
                         x$FP[i], x$n_pos[i], x$precision[i], x$recall[i], x$F1[i]), con)
  }
}
close(con)

cat("wrote:\n  results/eval_metricfair/pr_models.GT_rule.txt (full sweep)\n",
    " results/eval_metricfair/pr_models.GT_rule.table.txt (padj 0.01 by category)\n",
    " plots/pr_models.GT_rule.full.png\n  plots/pr_models.GT_rule.restricted.png\n", sep="")
