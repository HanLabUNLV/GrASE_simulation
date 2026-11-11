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
# Curves start at padj 0.001 (MAJIQ prob 0.99) and loosen. PARTIAL ROC: the grid
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
# Both universes are drawn. full divides every tool by the same n_pos/n_neg;
# restricted divides each tool by the units it can reach, so TP and FP are the
# same but both rates rise, most for the tools that reach least (MAJIQ, rMATS),
# and FPR is no longer on a common scale across tools.
#
# Usage: Rscript scripts/plot_roc_partial.R
# Out:   plots/roc_partial.GT_rule.{tolerant,strict}.pdf             (full)
#        plots/roc_partial.GT_rule.{tolerant,strict}.restricted.pdf  (restricted)
#        plots/roc_partial.GT_rule.all3{,.restricted}.pdf  strict, tolerant and gene together
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
get <- function(lv, U="full"){ d <- tab[tab$level==lv & tab$universe==U & tab$category=="ALL", ]
  ## MAJIQ keeps its full grid including 0.99 (a posterior probability, not an
  ## FDR -- 0.95 vs 0.99 is a routine choice, unlike padj 1e-3/1e-4).
  d <- d[!grepl("^MAJIQ",d$tool) & d$thr>=0.001 | grepl("^MAJIQ",d$tool), ]
  d$TPR <- d$TP/d$n_pos; d$FPR <- d$FP/d$n_neg; d[d$tool %in% TOOLS, ] }

## Vector PDF is the figure of record; the PNG is a 300 dpi convenience copy.
## Dimensions in INCHES, so the panel size is explicit rather than an artifact
## of pixel count divided by res.
## Both transcript FP rules are emitted. They share TP and n_pos, so TPR is
## IDENTICAL between them and only FPR moves -- by about an order of magnitude.
## strict   charges every co-travelling transcript, so its FPR is dominated by
##          implicated-set SIZE (call resolution).
## tolerant charges only transcripts implicated exclusively by calls that hit
##          nothing, so its FPR reflects decision ACCURACY.
## Use tolerant to compare tools on decision quality; strict is the resolution
## penalty shown alongside. The gene panel is identical in both files (gene
## level has no strict/tolerant split).
W <- 10.9; H <- 5.6
for (U in c("full","restricted")) {
SUF <- if (U == "full") "" else ".restricted"
for (TXLV in c("transcript_tolerant","transcript_strict")) {
TAG <- sub("^transcript_","",TXLV)
for (dev_i in 1:1) {
if (dev_i == 1) pdf(file.path(B,sprintf("plots/roc_partial.GT_rule.%s%s.pdf",TAG,SUF)), width=W, height=H)
else png(file.path(B,sprintf("plots/roc_partial.GT_rule.%s%s.png",TAG,SUF)), width=W, height=H,
         units="in", res=300)
par(mfrow=c(1,2), mar=c(4.6,4.8,3.2,1), oma=c(2.0,0,1.6,0))
for (cfg in list(c(TXLV, sprintf("transcript space (%s FP rule)",TAG)), c("gene","gene space"))) {
  d <- get(cfg[1], U); tl <- intersect(TOOLS, d$tool)
  plot(NA, xlim=c(0,max(d$FPR)*1.06), ylim=c(0,max(d$TPR)*1.08), xaxs="i", yaxs="i",
       xlab="FPR", ylab="TPR", main=cfg[2])
  grid(col="grey92", lty=1)
  for (t in tl){ z <- d[d$tool==t,]; z <- z[order(z$FPR),]
    lines(z$FPR,z$TPR,col=COL[t],lwd=2.2); points(z$FPR,z$TPR,col=COL[t],pch=16,cex=0.85) }
  legend("bottomright", legend=plab(tl), col=COL[tl], lwd=2, bty="n", cex=0.66)
}
mtext(sprintf("Partial ROC (%s universe): TPR vs FPR, curves start at padj 0.001 (MAJIQ prob 0.99) and loosen. UP-AND-LEFT is better.", U),
      outer=TRUE, cex=0.85, font=2)
mtext("GrASE curves use the merged exon+SJ unit with nested_BH (exontest.R default); GrASE (plain BH) is plain BH on the same raw p-values, matching DEXSeq's FDR architecture",
      side=1, outer=TRUE, cex=0.62, line=0.3, font=3)
invisible(dev.off())
}
cat(sprintf("wrote plots/roc_partial.GT_rule.%s%s.pdf\n", TAG, SUF))
}
}

## All three spaces side by side, one figure per universe.
for (U in c("full","restricted")) {
SUF <- if (U == "full") "" else ".restricted"
pdf(file.path(B,sprintf("plots/roc_partial.GT_rule.all3%s.pdf",SUF)), width=15.6, height=H)
par(mfrow=c(1,3), mar=c(4.6,4.8,3.2,1), oma=c(2.0,0,1.6,0))
for (cfg in list(c("transcript_strict","transcript space (strict FP rule)"),
                 c("transcript_tolerant","transcript space (tolerant FP rule)"),
                 c("gene","gene space"))) {
  d <- get(cfg[1], U); tl <- intersect(TOOLS, d$tool)
  plot(NA, xlim=c(0,max(d$FPR)*1.06), ylim=c(0,max(d$TPR)*1.08), xaxs="i", yaxs="i",
       xlab="FPR", ylab="TPR", main=cfg[2])
  grid(col="grey92", lty=1)
  for (t in tl){ z <- d[d$tool==t,]; z <- z[order(z$FPR),]
    lines(z$FPR,z$TPR,col=COL[t],lwd=2.2); points(z$FPR,z$TPR,col=COL[t],pch=16,cex=0.85) }
  legend("bottomright", legend=plab(tl), col=COL[tl], lwd=2, bty="n", cex=0.66)
}
mtext(sprintf("Partial ROC (%s universe): TPR vs FPR, curves start at padj 0.001 (MAJIQ prob 0.99) and loosen. UP-AND-LEFT is better.", U),
      outer=TRUE, cex=0.85, font=2)
mtext("GrASE curves use the merged exon+SJ unit with nested_BH (exontest.R default); GrASE (plain BH) is plain BH on the same raw p-values, matching DEXSeq's FDR architecture",
      side=1, outer=TRUE, cex=0.62, line=0.3, font=3)
invisible(dev.off())
cat(sprintf("wrote plots/roc_partial.GT_rule.all3%s.pdf\n", SUF))
}
cat("\n")
for (U in c("full","restricted")) for (lv in c("transcript_tolerant","transcript_strict","gene")) {
  d <- get(lv, U)
  cat(sprintf("=== %s, %s ===\n", lv, U))
  for (t in intersect(TOOLS, d$tool)) { z <- d[d$tool==t,]; z <- z[order(z$FPR),]
    cat(sprintf("  %-22s FPR %.5f -> %.5f   TPR %.3f -> %.3f\n", plab(t),
        min(z$FPR), max(z$FPR), min(z$TPR), max(z$TPR))) }
}
