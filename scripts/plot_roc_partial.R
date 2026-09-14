#!/usr/bin/env Rscript
# Partial ROC, transcript and gene space, including GrASE under BOTH FDR
# architectures.
#
# exontest.R ships with --padj_method nested_BH (gene screen, then units inside
# surviving genes), while DEXSeq applies plain BH across bins. That difference
# alone shapes the curves: under nested_BH only 1,345 null units are reachable
# at any threshold, against 209,981 under plain BH. GrASE_BH here is plain BH
# applied post hoc to the SAME raw p-values, so the two curves differ only in
# the multiple-testing procedure. At padj 0.01 they are nearly identical (132 vs
# 127 null calls); they diverge only as the threshold loosens.
#
# Curves start at padj 0.01 (MAJIQ prob 0.99) and loosen. PARTIAL ROC: the grid
# stops at 0.2 and structural non-coverage caps TPR well below 1, so AUC is not
# computable and the axes are zoomed to the operating quadrant.
#
# GrASE_BH is read from pr_three_levels.GT_rule.txt like every other curve.
# It used to be recomputed here from the annotated test files (p.adjust on the
# raw p-values, implicated sets re-parsed, TPR/FPR re-derived against
# simulate.rda). That second implementation agreed with the table to all
# printed digits at every threshold in both spaces, so it was redundant -- and
# a standing drift risk, since a change to the sweep's call rule or universe
# would have moved the table without moving this script. One source now.
#
# Usage: Rscript scripts/plot_roc_partial.R
# Out:   plots/roc_partial.GT_rule.png
B <- "/mnt/data1/home/mirahan/GrASE_simulation"
source(file.path(B, "scripts/plot_palette.R"))
tab <- read.table(file.path(B,"results/eval_metricfair/pr_three_levels.GT_rule.txt"),
                  header=TRUE, sep="\t", stringsAsFactors=FALSE)

# GrASE_BH first: it is a GrASE variant, so it belongs with the other GrASE
# entries in the legend rather than orphaned at the end. intersect() below
# preserves this order.
TOOLS <- c("GrASE_BH","GrASE_merged_all","GrASE_merged_dpi0.1","GrASE_merged_dpi0.2",
           "DEXSeq","MAJIQ_C0.10","MAJIQ_C0.20","rMATS","rMATS_dpsi0.1","rMATS_dpsi0.2")
# COL, disp(), plab() come from plot_palette.R (sourced above).
get <- function(lv){ d <- tab[tab$level==lv & tab$universe=="full" & tab$category=="ALL", ]
  ## MAJIQ keeps its full grid including 0.99 (a posterior probability, not an
  ## FDR -- 0.95 vs 0.99 is a routine choice, unlike padj 1e-3/1e-4).
  d <- d[!grepl("^MAJIQ",d$tool) & d$thr>=0.01 | grepl("^MAJIQ",d$tool), ]
  d$TPR <- d$TP/d$n_pos; d$FPR <- d$FP/d$n_neg; d[d$tool %in% TOOLS, ] }

png(file.path(B,"plots/roc_partial.GT_rule.png"), width=1520, height=780, res=140)
par(mfrow=c(1,2), mar=c(4.6,4.8,3.2,1), oma=c(2.0,0,1.6,0))
for (cfg in list(c("transcript_strict","transcript space"), c("gene","gene space"))) {
  d <- get(cfg[1]); tl <- intersect(TOOLS, d$tool)
  plot(NA, xlim=c(0,max(d$FPR)*1.06), ylim=c(0,max(d$TPR)*1.08), xaxs="i", yaxs="i",
       xlab="FPR", ylab="TPR", main=cfg[2])
  grid(col="grey92", lty=1)
  for (t in tl){ z <- d[d$tool==t,]; z <- z[order(z$FPR),]
    lines(z$FPR,z$TPR,col=COL[t],lwd=2.2); points(z$FPR,z$TPR,col=COL[t],pch=16,cex=0.85) }
  legend("bottomright", legend=plab(tl), col=COL[tl], lwd=2, bty="n", cex=0.66)
}
mtext("Partial ROC: TPR vs FPR, curves start at padj 0.01 (MAJIQ prob 0.99) and loosen. UP-AND-LEFT is better.",
      outer=TRUE, cex=0.85, font=2)
mtext("GrASE curves use the merged exon+SJ unit with nested_BH (exontest.R default); GrASE (plain BH) is plain BH on the same raw p-values, matching DEXSeq's FDR architecture",
      side=1, outer=TRUE, cex=0.62, line=0.3, font=3)
invisible(dev.off())
cat("wrote plots/roc_partial.GT_rule.png\n\n")
for (lv in c("transcript_strict","gene")) {
  d <- get(lv)
  cat(sprintf("=== %s ===\n", lv))
  for (t in intersect(TOOLS, d$tool)) { z <- d[d$tool==t,]; z <- z[order(z$FPR),]
    cat(sprintf("  %-22s FPR %.5f -> %.5f   TPR %.3f -> %.3f\n", plab(t),
        min(z$FPR), max(z$FPR), min(z$TPR), max(z$TPR))) }
}
