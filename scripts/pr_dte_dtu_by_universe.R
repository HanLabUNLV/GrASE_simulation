#!/usr/bin/env Rscript
#
# PR curves laid out as DTE and DTU ROWS, one figure per universe.
#
# The sweep's own figures (plots/pr_three_levels.GT_rule_{DTE,DTU}.png) put the
# two universes in rows and split the categories across separate FILES, which
# makes the DTE-vs-DTU contrast -- the comparison the results text actually
# makes -- impossible to see without flipping between images. This transposes
# that: rows are DTE and DTU, columns are the four evaluation levels, and the
# universe becomes the file split.
#
# Reads the table the sweep already wrote, so it does NOT recompute anything
# and cannot drift from the published numbers. Regenerate the table first with
#   STRANDED=1 Rscript scripts/pr_curves_three_levels_gtrule.R
# if the underlying results have changed.
#
# Usage: Rscript scripts/pr_dte_dtu_by_universe.R
BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
TAB  <- file.path(BASE, "results/eval_metricfair/pr_three_levels.GT_rule.txt")
stopifnot(file.exists(TAB))
tab <- read.table(TAB, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

LEVELS <- c("unit","transcript_strict","transcript_tolerant","gene")

## Styling copied verbatim from pr_curves_three_levels_gtrule.R so these panels
## are colour-consistent with the sweep's own figures.
COL <- c("GrASE_BH"="#fec44f",
         GrASE="#fb6a4a", "GrASE_merged_all"="#fb6a4a",
         "GrASE_dpi0.1"="#de2d26", "GrASE_merged_dpi0.1"="#de2d26",
         "GrASE_dpi0.2"="#a50f15", "GrASE_merged_dpi0.2"="#a50f15",
         "GrASE_internal"="#67000d", "GrASE_merged_internal"="#67000d",
         DEXSeq="#3182bd",
         MAJIQ_C0.10="#a1d99b", MAJIQ_C0.20="#238b45",
         rMATS="#bcbddc",
         "rMATS_dpsi0.1"="#756bb1", "rMATS_dpsi0.2"="#3f007d")
PLOT_TOOLS <- setdiff(names(COL),
                      c("GrASE", "GrASE_dpi0.1", "GrASE_dpi0.2",
                        "GrASE_internal", "GrASE_merged_internal"))
disp <- function(x) {
  x <- sub("^GrASE_BH$", "GrASE (plain BH)", x)
  sub("^GrASE_merged_", "GrASE_", sub("^GrASE_merged_all$", "GrASE_dpi0", x))
}
LVLAB <- c(unit = "Unit", transcript_strict = "transcript_strict",
           transcript_tolerant = "transcript_tolerant", gene = "Gene")

## Two figure families in the same style: the DTE/DTU contrast as two rows,
## and ALL on its own as a single row. Same columns, same colours, same axes.
FAMS <- list(list(tag = "pr_dte_dtu", cats = c("DTE", "DTU")),
             list(tag = "pr_all",     cats = "ALL"))

for (fam in FAMS) {
 CATS <- fam$cats
 for (uni in c("full", "restricted")) {
  W <- 4.0 * length(LEVELS); H <- if (length(CATS) == 2) 8.2 else 5.0
  for (dev_i in 1:2) {
    f <- file.path(BASE, sprintf("plots/%s.%s.%s", fam$tag, uni,
                                 if (dev_i == 1) "pdf" else "png"))
    if (dev_i == 1) pdf(f, width = W, height = H)
    else png(f, width = W, height = H, units = "in", res = 300)
    ## rows = category, columns = level. oma reserves the title strip at the
    ## top and the legend strip at the bottom.
    par(mfrow = c(length(CATS), length(LEVELS)), mar = c(4.2, 4.4, 3.2, 1),
        oma = c(if (length(CATS) == 2) 4.2 else 5.6, 1.2, 3.0, 0))

    for (categ in CATS) {
      for (lv in LEVELS) {
        d <- tab[tab$level == lv & tab$universe == uni & tab$category == categ, ]
        d <- d[d$tool %in% PLOT_TOOLS, ]
        ## Match pr_curves_three_levels_gtrule.R: start at padj 0.001 (1e-4 trimmed,
        ## MAJIQ keeps its probability grid); UNIT is a within-GrASE view.
        d <- d[!grepl("^MAJIQ", d$tool) & d$thr >= 0.001 | grepl("^MAJIQ", d$tool), ]
        if (lv == "unit") d <- d[grepl("^GrASE", d$tool), ]
        plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Recall",
             ylab = "Precision", main = sprintf("%s -- %s", categ,
             ifelse(is.na(LVLAB[lv]), lv, LVLAB[lv])))
        if (!nrow(d)) { text(0.5, 0.5, "no data", col = "grey60"); next }
        for (tl in unique(d$tool)) {
          x <- d[d$tool == tl, ]; x <- x[order(x$recall), ]
          lines(x$recall, x$precision, col = COL[[tl]], lwd = 2)
          points(x$recall, x$precision, col = COL[[tl]], pch = 19, cex = 0.55)
        }
      }
    }

    ## one shared legend for all eight panels
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    tl_in <- intersect(PLOT_TOOLS, unique(tab$tool))
    legend("bottom", legend = disp(tl_in), col = unname(COL[tl_in]),
           lwd = 2, pch = 19, pt.cex = 0.6, horiz = FALSE, ncol = 5,
           bty = "n", cex = 0.78,
           inset = c(0, if (length(CATS) == 2) 0.005 else 0.03))
    mtext(sprintf("Precision-recall by evaluation level, %s universe (GT_rule)", uni),
          side = 3, line = 1.0, outer = TRUE, cex = 1.0, font = 2)
    invisible(dev.off())
  }
  cat(sprintf("wrote plots/%s.%s.{pdf,png}\n", fam$tag, uni))
 }
}
