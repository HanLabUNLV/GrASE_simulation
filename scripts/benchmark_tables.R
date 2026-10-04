#!/usr/bin/env Rscript
# All four benchmark tables, from one row order and one set of labels.
#
#   native_unit_table.txt        each tool on its OWN unit (WITHIN-TOOL ONLY)
#   structural_reach_table.txt   what each universe can address (shared currency)
#   transcript_level_table.txt   transcript currency (cross-tool valid)
#   gene_level_table_fmt.txt     gene currency, restricted and full universes
#
# Kept in one script so the tables cannot drift apart in row order or naming.
#
# LABELS. The ungated GrASE row is GrASE_dpi0, never plain "GrASE": exontest.R
# ships --min_dpi 0.1, so its `significant` column means padj < 0.01 AND
# |delta_pi| >= 0.1. Labelling the ungated run "GrASE" would present something
# 65% larger than the shipped default (1,012 vs 615 calls on the merged stranded
# internal file). Naming all three by their gate makes the axis explicit.
#
# The sweep and the transcript metrics use different internal tool names for the
# same configuration, so both maps are given below.
#
# Usage: Rscript scripts/benchmark_tables.R
# In:   results/eval_metricfair/pr_three_levels.GT_rule.txt   (unit, gene)
#       results/eval_metricfair/pr_calls_cache.GT_rule.rds    (genes tested)
#       results/eval_metricfair/transcript_level_metrics.txt  (transcript)
suppressPackageStartupMessages(library(dplyr))
BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
OUT  <- file.path(BASE, "results/eval_metricfair")
tab  <- read.table(file.path(OUT,"pr_three_levels.GT_rule.txt"), header=TRUE,
                   sep="\t", stringsAsFactors=FALSE)
txm  <- read.table(file.path(OUT,"transcript_level_metrics.txt"), header=TRUE,
                   sep="\t", stringsAsFactors=FALSE)
cl   <- readRDS(file.path(OUT,"pr_calls_cache.GT_rule.rds"))
ngene <- vapply(cl, function(z) length(unique(z$gene)), integer(1))

LAB   <- c("GrASE (plain BH)","GrASE_dpi0","GrASE_dpi0.1","GrASE_dpi0.2","GrASE_internal",
           "MAJIQ_C0.20","MAJIQ_C0.10","DEXSeq","rMATS","rMATS_dpsi0.1","rMATS_dpsi0.2")
SWEEP <- c("GrASE_BH","GrASE_merged_all","GrASE_merged_dpi0.1","GrASE_merged_dpi0.2",
           "GrASE_merged_internal","MAJIQ_C0.20","MAJIQ_C0.10","DEXSeq",
           "rMATS","rMATS_dpsi0.1","rMATS_dpsi0.2")
## Every tool here is on strand-corrected input. rMATS used to be the exception
## (its rows were the uncorrected run) -- strand correction costs it zero TPs but
## removes ~1.3% of its strict FPs, so leaving it uncorrected mixed conventions
## within one table for no reason.
TXNM  <- c("GrASE_BH","GrASE_merged_stranded","GrASE_merged_stranded_dpi0.1",
           "GrASE_merged_stranded_dpi0.2","GrASE_merged_stranded_internal",
           "MAJIQ_C0.20_stranded","MAJIQ_C0.10_stranded","DEXSeq_stranded",
           "rMATS_stranded","rMATS_stranded_dpsi0.1","rMATS_stranded_dpsi0.2")
thr_of <- function(t) if (grepl("^MAJIQ", t)) 0.95 else 0.01
cm  <- function(x) formatC(x, format="d", big.mark=",")
md  <- function(hdr, lines) { cat("|", paste(hdr, collapse=" | "), "|\n")
  cat("|", paste(rep("---", length(hdr)), collapse="|"), "|\n"); cat(lines, sep="") }

pick <- function(tool, lv, uni) {
  r <- tab[tab$tool==tool & tab$level==lv & tab$universe==uni &
           tab$category=="ALL" & tab$thr==thr_of(tool), ]
  if (nrow(r)) r[1, ] else NULL
}

## --- 1. native unit --------------------------------------------------------
rows <- list(); out <- character(0)
for (i in seq_along(SWEEP)) {
  r <- pick(SWEEP[i], "unit", "restricted"); if (is.null(r)) next
  U <- r$TP + r$FP + r$FN + r$TN
  rows[[length(rows)+1]] <- data.frame(tool=LAB[i], Universe=U,
    Genes=if (SWEEP[i] %in% names(ngene)) ngene[[SWEEP[i]]] else NA_integer_,
    TP=r$TP, FP=r$FP, FN=r$FN, TN=r$TN, n_pos=r$n_pos, n_neg=r$n_neg,
    Precision=round(r$precision,3), Recall=round(r$recall,3), stringsAsFactors=FALSE)
  out <- c(out, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %.3f | %.3f |\n",
    LAB[i], cm(U), cm(rows[[length(rows)]]$Genes), cm(r$TP), cm(r$FP), cm(r$FN),
    cm(r$TN), cm(r$n_pos), cm(r$n_neg), r$precision, r$recall))
}
d1 <- bind_rows(rows)
write.table(d1, file.path(OUT,"native_unit_table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n### native unit (WITHIN-TOOL ONLY: different units, different ground truths)\n")
md(c("tool","Universe","Genes","TP","FP","FN","TN","n_pos","n_neg","Precision","Recall"), out)

## --- 2. structural reach ---------------------------------------------------
REACH <- c(GrASE_merged_all="GrASE", MAJIQ_C0.20="MAJIQ", DEXSeq="DEXSeq", rMATS="rMATS")
NG <- 3002; NT <- 4503; rr <- list(); out <- character(0)
for (t in names(REACH)) {
  g <- pick(t,"gene","restricted"); x <- pick(t,"transcript_strict","restricted")
  if (is.null(g) || is.null(x)) next
  rr[[length(rr)+1]] <- data.frame(tool=REACH[[t]], genes_reachable=g$n_pos, genes_total=NG,
    pct_genes=round(100*g$n_pos/NG,1), tx_reachable=x$n_pos, tx_total=NT,
    pct_tx=round(100*x$n_pos/NT,1), stringsAsFactors=FALSE)
  out <- c(out, sprintf("| %s | %s / %s (%.1f%%) | %s / %s (%.1f%%) |\n", REACH[[t]],
    cm(g$n_pos), cm(NG), 100*g$n_pos/NG, cm(x$n_pos), cm(NT), 100*x$n_pos/NT))
}
write.table(bind_rows(rr), file.path(OUT,"structural_reach_table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n### structural reach (shared currency: genes and transcripts)\n")
md(c("tool","signal genes reachable","manipulated transcripts reachable"), out)

## --- 3. transcript level ---------------------------------------------------
rows <- list(); out <- character(0)
for (i in seq_along(TXNM)) {
  r <- txm[txm$tool==TXNM[i] & txm$category=="ALL", ]; if (!nrow(r)) next
  r <- r[1, ]
  rows[[length(rows)+1]] <- data.frame(tool=LAB[i], TP=r$TP, FP_strict=r$FP_strict,
    FP_tolerant=r$FP_tolerant, transcript_recall=r$tx_recall,
    precision_strict=r$prec_strict, precision_tolerant=r$prec_tolerant,
    n_calls=r$n_sig_calls, n_calls_hit=r$calls_hit,
    n_transcripts_per_call=r$sharp_median, stringsAsFactors=FALSE)
  out <- c(out, sprintf("| %s | %s | %s | %s | %.4f | %.4f | %.4f | %s | %s | %s |\n",
    LAB[i], cm(r$TP), cm(r$FP_strict), cm(r$FP_tolerant), r$tx_recall,
    r$prec_strict, r$prec_tolerant, cm(r$n_sig_calls), cm(r$calls_hit), r$sharp_median))
}
write.table(bind_rows(rows), file.path(OUT,"transcript_level_table.txt"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n### transcript level (shared currency; strict and tolerant bracket unit granularity)\n")
md(c("tool","TP","FP strict","FP tolerant","transcript recall","precision strict",
     "precision tolerant","n_calls","n_calls_hit","n_tx per call"), out)

## --- 4. gene level, both universes -----------------------------------------
rows <- list()
for (uni in c("restricted","full")) {
  out <- character(0)
  for (i in seq_along(SWEEP)) {
    r <- pick(SWEEP[i], "gene", uni); if (is.null(r)) next
    U <- r$TP + r$FP + r$FN + r$TN
    rows[[length(rows)+1]] <- data.frame(universe=uni, tool=LAB[i], Universe=U,
      TP=r$TP, FP=r$FP, FN=r$FN, TN=r$TN, n_pos=r$n_pos, n_neg=r$n_neg,
      Precision=round(r$precision,3), Recall=round(r$recall,3), stringsAsFactors=FALSE)
    out <- c(out, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %.3f | %.3f |\n",
      LAB[i], cm(U), cm(r$TP), cm(r$FP), cm(r$FN), cm(r$TN), cm(r$n_pos), cm(r$n_neg),
      r$precision, r$recall))
  }
  cat(sprintf("\n### gene level, %s universe\n", uni))
  md(c("tool","Universe","TP","FP","FN","TN","n_pos","n_neg","Precision","Recall"), out)
}
write.table(bind_rows(rows), file.path(OUT,"gene_level_table_fmt.txt"), sep="\t", quote=FALSE, row.names=FALSE)
## --- 5. bubble-class attribution (Table 6) ---------------------------------
## Computed by transcript_level_metrics.R, which owns the manipulated-transcript
## truth and the internal/TSSTTS src tag; formatted here with the other tables.
## The Total row reconciles with the transcript table: recovered(internal) +
## recovered(TSSTTS) - both = the GrASE_dpi0 row's TP.
attf <- file.path(OUT, "tsstts_attribution.txt")
if (file.exists(attf)) {
  at <- read.table(attf, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                   check.names=FALSE)
  out <- character(0)
  for (i in seq_len(nrow(at)))
    out <- c(out, sprintf("| %s | %s | %s | %s |\n", at$class[i],
             cm(at$ALL[i]), cm(at$DTE[i]), cm(at$DTU[i])))
  cat("\n### bubble class recovering each manipulated transcript\n")
  md(c(" ","ALL","DTE","DTU"), out)
  tot <- at$ALL[at$class == "Total"]
  pc  <- round(100 * at$ALL / tot)
  cat(sprintf("\nin-text: %s (%d%%) TSS/TTS only, %s (%d%%) internal only, %s (%d%%) both, of %s recovered\n",
      cm(at$ALL[at$class=="TSS/TTS only"]),  pc[at$class=="TSS/TTS only"],
      cm(at$ALL[at$class=="internal only"]), pc[at$class=="internal only"],
      cm(at$ALL[at$class=="both"]),          pc[at$class=="both"], cm(tot)))
} else {
  cat("\nNOTE: tsstts_attribution.txt absent -- run transcript_level_metrics.R\n")
}

cat(sprintf("\nwrote 5 tables to %s\n", OUT))
