#!/usr/bin/env Rscript
#
# Partial ROC laid out as DTE and DTU ROWS, one figure per universe -- the ROC
# counterpart to scripts/pr_dte_dtu_by_universe.R.
#
# PARTIAL, for the same reason scripts/plot_roc_partial.R is: the threshold grid
# stops well short of calling everything and structural non-coverage caps TPR
# far below 1, so the curves do not reach (1,1) and AUC is not computable.
# Axes are zoomed to the operating quadrant -- a 0-1 FPR axis would compress
# every curve onto the y-axis, since the classes are imbalanced roughly 1:35 at
# gene level and observed FPR never exceeds 0.22.
#
# The x-limit is shared between the DTE and DTU panels of a COLUMN so the two
# rows are directly comparable; it is not shared across columns, because the
# levels have different negative pools and a common limit would flatten three
# of the four.
#
# CAVEAT drawn on the figure: at GENE level DTE and DTU share the same negative
# pool (the null genes), so FP, TN and FPR are identical between the two rows by
# construction and only TPR differs there. The transcript and unit levels do not
# have this property.
#
# Reads the table the sweep already wrote; recomputes nothing.
# Usage: Rscript scripts/roc_dte_dtu_by_universe.R
B <- "/mnt/data1/home/mirahan/GrASE_simulation"
source(file.path(B, "scripts/plot_palette.R"))     # COL, disp, plab
tab <- read.table(file.path(B, "results/eval_metricfair/pr_three_levels.GT_rule.txt"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

LEVELS <- c("unit","transcript_strict","transcript_tolerant","gene")
LVLAB  <- c(unit = "Unit", transcript_strict = "transcript_strict",
            transcript_tolerant = "transcript_tolerant", gene = "Gene")
PLOT_TOOLS <- setdiff(names(COL),
                      c("GrASE", "GrASE_dpi0.1", "GrASE_dpi0.2",
                        "GrASE_internal", "GrASE_merged_internal"))

d <- subset(tab, category %in% c("DTE","DTU") & level %in% LEVELS &
                 !is.na(n_neg) & n_neg > 0 & !is.na(n_pos) & n_pos > 0 &
                 tool %in% PLOT_TOOLS)
d$FPR <- d$FP / d$n_neg
d$TPR <- d$TP / d$n_pos
## Match the PR figures (pr_curves_three_levels_gtrule.R): curves start at padj
## 0.001 (1e-4 trimmed; MAJIQ keeps its probability grid), and the UNIT panel is
## a within-GrASE view -- other tools' units are scored against different GTs.
## The written table keeps every row.
d <- d[!grepl("^MAJIQ", d$tool) & d$thr >= 0.001 | grepl("^MAJIQ", d$tool), ]
d <- d[d$level != "unit" | grepl("^GrASE", d$tool), ]

for (uni in c("full","restricted")) {
  du <- subset(d, universe == uni)
  ## shared x-limit per column (level), across both category rows
  xmax <- sapply(LEVELS, function(lv) {
    v <- du$FPR[du$level == lv]; if (!length(v)) 1 else max(v) * 1.05 })
  ymax <- sapply(LEVELS, function(lv) {
    v <- du$TPR[du$level == lv]; if (!length(v)) 1 else min(1, max(v) * 1.08) })

  W <- 4.0 * length(LEVELS); H <- 8.2
for (dev_i in 1:1) {
    f <- file.path(B, sprintf("plots/roc_dte_dtu.%s.%s", uni,
                              if (dev_i == 1) "pdf" else "png"))
    if (dev_i == 1) pdf(f, width = W, height = H)
    else png(f, width = W, height = H, units = "in", res = 300)
    par(mfrow = c(2, length(LEVELS)), mar = c(4.2, 4.4, 3.2, 1),
        oma = c(4.6, 1.2, 3.0, 0))

    for (categ in c("DTE","DTU")) {
      for (li in seq_along(LEVELS)) {
        lv <- LEVELS[li]
        x <- subset(du, level == lv & category == categ)
        plot(NA, xlim = c(0, xmax[li]), ylim = c(0, ymax[li]),
             xlab = "False positive rate", ylab = "True positive rate",
             main = sprintf("%s -- %s", categ,
                            ifelse(is.na(LVLAB[lv]), lv, LVLAB[lv])))
        ## chance line, for reference only; at these axis limits it is nearly flat
        abline(0, 1, col = "grey80", lty = 3)
        if (!nrow(x)) { text(xmax[li]/2, ymax[li]/2, "no data", col = "grey60"); next }
        for (tl in unique(x$tool)) {
          z <- x[x$tool == tl, ]; z <- z[order(z$FPR, z$TPR), ]
          lines(z$FPR, z$TPR, col = COL[[tl]], lwd = 2)
          points(z$FPR, z$TPR, col = COL[[tl]], pch = 19, cex = 0.55)
        }
      }
    }

    par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    tl_in <- intersect(PLOT_TOOLS, unique(du$tool))
    legend("bottom", legend = plab(tl_in), col = unname(COL[tl_in]),
           lwd = 2, pch = 19, pt.cex = 0.6, ncol = 5, bty = "n",
           cex = 0.78, inset = c(0, 0.042))
    mtext(sprintf("Partial ROC by evaluation level, %s universe (GT_rule). Axes zoomed to the operating quadrant; AUC not computable.", uni),
          side = 3, line = 1.0, outer = TRUE, cex = 0.92, font = 2)
    mtext("At GENE level DTE and DTU share the same negative pool, so FPR is identical between rows by construction and only TPR differs.",
          side = 1, line = 2.6, outer = TRUE, cex = 0.72, col = "#555555")
    invisible(dev.off())
  }
  cat(sprintf("wrote plots/roc_dte_dtu.%s.{pdf,png}\n", uni))
}
