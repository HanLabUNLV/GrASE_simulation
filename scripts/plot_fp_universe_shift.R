#!/usr/bin/env Rscript
# Where each tool's false positives LIVE: signal genes vs null genes.
#
# Replaces the null-gene ROC (plot_fpr_recall_curve.R, deleted 2026-09-11).
# Because TPR is identical between the full-universe and null-gene ROCs (same
# TP, same n_pos), the only thing that differed between those two figures was
# ONE number per tool -- the FPR. Drawing that as a second two-axis ROC made the
# plots look alike: the shared TPR fixes every y position, and a linear FPR axis
# dominated by rMATS/DEXSeq squashes the region where the reordering actually
# happens. Nine of ten tools change FPR rank between the two universes and the
# full/null ratio spans 2.4x to 108x, so the signal is real -- it just was not
# visible as two ROCs.
#
# Left  : FPR on null genes vs FPR on the full universe, log-log. Diagonals are
#         constant ratios. A tool low-and-right converts almost all of its false
#         positives into co-travellers inside genes that genuinely changed; a
#         tool near the top calls genes with no splicing change at all.
# Right : the same thing as one interpretable number -- what share of a tool's
#         false positives sit in genes that have no splicing change by design.
#
# Usage: Rscript scripts/plot_fp_universe_shift.R
# Out:   plots/fp_universe_shift.GT_rule.png
B <- "/mnt/data1/home/mirahan/GrASE_simulation"
source(file.path(B, "scripts/plot_palette.R"))
tab <- read.table(file.path(B,"results/eval_metricfair/pr_three_levels.GT_rule.txt"),
                  header=TRUE, sep="\t", stringsAsFactors=FALSE)
d <- tab[tab$level=="transcript_strict" & tab$universe=="full", ]
TOOLS <- c("GrASE_BH","GrASE_merged_all","GrASE_merged_dpi0.1","GrASE_merged_dpi0.2",
           "DEXSeq","MAJIQ_C0.10","MAJIQ_C0.20","rMATS","rMATS_dpsi0.1","rMATS_dpsi0.2")
# Structural negative class: all transcripts of Background and DGE genes. Its
# size must not depend on whether a tool happened to make calls in a category --
# taking n_neg from the tool's own rows dropped 12,193 transcripts for any tool
# with no DGE calls (GrASE_dpi0.1 read 0.00030 instead of 0.00028).
NNULL <- sum(vapply(c("Background","DGE"), function(k) max(d$n_neg[d$category==k]), numeric(1)))
op <- function(t) if (grepl("^MAJIQ", t)) 0.95 else 0.01

z <- do.call(rbind, lapply(TOOLS, function(t) {
  a <- d[d$tool==t & d$category=="ALL" & d$thr==op(t), ]
  if (!nrow(a)) return(NULL)
  fpn <- sum(d$FP[d$tool==t & d$thr==op(t) & d$category %in% c("Background","DGE")])
  data.frame(tool=t, fpr_full=a$FP[1]/a$n_neg[1], fpr_null=fpn/NNULL,
             share=100*fpn/a$FP[1], stringsAsFactors=FALSE)
}))

png(file.path(B,"plots/fp_universe_shift.GT_rule.png"), width=1560, height=740, res=140)
par(mfrow=c(1,2), mar=c(4.6,4.8,3.0,1), oma=c(1.8,0,1.6,0))

plot(z$fpr_full, z$fpr_null, log="xy", pch=16, cex=1.5, col=COL[z$tool],
     xlab="FPR, full universe (all non-manipulated transcripts)",
     ylab="FPR, null genes only (Background + DGE)",
     main="where the false positives live",
     xlim=range(z$fpr_full)*c(0.75,1.6), ylim=range(z$fpr_null)*c(0.6,2.2))
for (r in c(2,10,50,100)) {
  xs <- range(z$fpr_full)*c(0.5,2)
  lines(xs, xs/r, col="grey85", lty=2)
  text(xs[2], xs[2]/r, sprintf("%dx", r), col="grey55", cex=0.6, pos=4, xpd=NA)
}
points(z$fpr_full, z$fpr_null, pch=16, cex=1.5, col=COL[z$tool])
text(z$fpr_full, z$fpr_null, plab(z$tool), pos=3, cex=0.58, col=COL[z$tool])

o <- order(z$share)
bp <- barplot(z$share[o], horiz=TRUE, col=COL[z$tool[o]], border=NA, las=1,
              xlab="% of a tool's false positives that are in null genes",
              main="share of FPs with no splicing signal at all",
              xlim=c(0, max(z$share)*1.35), names.arg=rep("", nrow(z)))
text(z$share[o], bp, sprintf(" %.1f%%  %s", z$share[o], plab(z$tool[o])),
     pos=4, cex=0.62, xpd=NA)
mtext("Same TPR in both universes, so the ONLY difference between the two ROCs is this shift. Lower is better in both panels.",
      outer=TRUE, cex=0.80, font=2)
invisible(dev.off())
cat("wrote plots/fp_universe_shift.GT_rule.png\n\n")
print(z[order(z$share), ], row.names=FALSE, digits=3)
