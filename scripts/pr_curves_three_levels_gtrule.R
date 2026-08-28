#!/usr/bin/env Rscript
#
# scripts/pr_curves_three_levels.R
#
# Precision-recall curves at THREE scoring currencies:
#   1. UNIT   -- each tool on its own native unit, scored by its own GT_rule
#                where one exists:
#                  GrASE  bipartition test  GT_rule_bipartition
#                  MAJIQ  LSV               GT_rule_lsv
#                  DEXSeq exonic bin        part-level gtI  (no GT_rule built)
#                  rMATS  event             structural junction GT (unchanged)
#                WITHIN-TOOL ONLY: per-tool references make the unit panel
#                non-comparable across tools. Use the TRANSCRIPT and GENE panels
#                for cross-tool reading.
#                NOTE the FULL-universe recall denominator is dropped for the
#                unit level here: GT_rule can only label units the tool actually
#                quantified, so a full-universe positive count does not exist for
#                MAJIQ under it. Unit recall below is RESTRICTED-universe.
#   LEVEL NAMES: transcript_strict and transcript_tolerant are the same recall
#   (manipulated transcripts implicated by >=1 call) with two precisions --
#   strict charges every co-travelling transcript, tolerant charges only calls
#   that hit nothing. Formerly "transcript" and "call".
#
#   2. TRANSCRIPT_STRICT -- common currency. A transcript is CALLED if it is in the
#                implicated set of >= 1 significant call; positive if it is a
#                manipulated transcript (iso.dtu / iso.dte). Negatives are all
#                other transcripts of tested genes (incl. Background/DGE).
#   3. GENE   -- common currency. A gene is CALLED if it has >= 1 significant
#                call; positive if DTE or DTU; negative if Background or DGE.
#
# Every tool is swept over its own statistic (GrASE/DEXSeq padj, MAJIQ
# P(|dPSI|>=C), rMATS FDR) -- points are NOT threshold-matched across tools;
# each curve is that tool's own operating characteristic.
#
# Usage:  Rscript scripts/pr_curves_three_levels_gtrule.R
# Output: results/eval_metricfair/pr_three_levels.GT_rule.txt
#         plots/pr_three_levels.GT_rule.*.png

suppressPackageStartupMessages({ library(dplyr); library(rtracklayer) })
BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
OUT  <- file.path(BASE, "results/eval_metricfair")
CACHE <- file.path(OUT, "pr_calls_cache.GT_rule.rds")
load(file.path(BASE, "swimdown/simulate/data/simulate.rda"))
get_st <- function(g) ifelse(g %in% dte.genes, "DTE", ifelse(g %in% dtu.genes, "DTU",
                     ifelse(g %in% dge.genes, "DGE", "Background")))
parse_l <- function(x, sep = "[,+]") { v <- strsplit(x, sep)[[1]]; trimws(v[nchar(trimws(v)) > 0 & trimws(v) != "NA"]) }
manip_tx <- unique(c(names(which(iso.dtu)), names(which(iso.dte))))
sig_genes <- unique(c(dte.genes, dtu.genes))
gene_tx <- split(txdf$TXNAME, txdf$GENEID)   # needed on both the cached and uncached path

## --- GT tables (cheap reads; needed on both cached and uncached paths so the
##     sweep can be subset by simulation category) ----------------------------
ep <- readRDS(file.path(OUT, "ep_gt_cache.rds"))
ep$key <- paste(ep$gene, ep$exonic_part, sep = ":")
pos_key <- new.env(); for (k in ep$key[ep$gt_pos]) assign(k, TRUE, envir = pos_key)
neg_key <- new.env(); for (k in ep$key[ep$gt_neg]) assign(k, TRUE, envir = neg_key)
lsvgt   <- read.table(file.path(BASE, "results/sim_lsv_gt.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lsv_pos <- setNames(lsvgt$gt_positive, lsvgt$lsv_id)
lsv_full<- read.table(file.path(BASE, "results/sim_lsv_gt.full.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
## --- GT_rule unit truth -----------------------------------------------------
bpgt <- read.table(file.path(BASE, "results/gt/bipartition_gt.txt"), header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
bpgt$uid <- paste(bpgt$src, bpgt$gene, bpgt$event, bpgt$comparison, sep = "|")
bp_pos <- setNames(bpgt$GT_rule_bipartition, bpgt$uid)
n_pos_bp_cat <- table(bpgt$sim_type[bpgt$GT_rule_bipartition])
jngt <- read.table(file.path(BASE, "results/gt/lsv_junction_gt.txt"), header = TRUE,
                   sep = "\t", quote = "", stringsAsFactors = FALSE)
jngt$uid <- paste(jngt$lsv_id, jngt$j_index, sep = "|")
jn_pos <- setNames(jngt$GT_rule_junction, jngt$uid)
n_pos_jn_cat <- table(jngt$sim_type[jngt$GT_rule_junction])
lsvr <- read.table(file.path(BASE, "results/gt/lsv_gt.txt"), header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
lsv_pos_rule <- setNames(lsvr$GT_rule_lsv, lsvr$lsv_id)
n_pos_lsvrule_cat <- table(lsvr$sim_type[lsvr$GT_rule_lsv])
jgt_all <- read.table(file.path(BASE, "results.beforegatefix/sim_junction_gt.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
jgt_all <- jgt_all[!(jgt_all$n_inc_tx == 0 | jgt_all$n_skip_tx == 0), ]
jgt_all$uid <- paste(jgt_all$event_type, jgt_all$ID, sep = ":")
rm_pos  <- setNames(jgt_all$gt_positive, jgt_all$uid)
# per-category FULL positive counts
cat_of_gene <- function(g) get_st(g)
n_pos_parts_cat <- table(ep$sim_type[ep$gt_pos])
n_pos_lsv_cat   <- table(lsv_full$sim_type[lsv_full$gt_positive])
n_pos_evt_cat   <- table(jgt_all$sim_type[jgt_all$gt_positive])

## ===========================================================================
## Build per-tool CALL tables once: (gene, stat, unit_id, imp = transcripts,
## unit_pos = is the unit GT-positive). stat is oriented so that SMALLER = more
## significant for padj/FDR tools and is negated for MAJIQ probabilities.
## ===========================================================================
if (file.exists(CACHE)) { calls <- readRDS(CACHE) } else {
calls <- list()

## --- GrASE (exonic parts; unit GT from the cached exonic-part GT) ----------
ep <- readRDS(file.path(OUT, "ep_gt_cache.rds"))
ep$key <- paste(ep$gene, ep$exonic_part, sep = ":")
pos_key <- new.env(); for (k in ep$key[ep$gt_pos]) assign(k, TRUE, envir = pos_key)
neg_key <- new.env(); for (k in ep$key[ep$gt_neg]) assign(k, TRUE, envir = neg_key)

grase_files <- c("bipartition.test.fulldesign/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
                 "bipartition.test.fulldesign/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
grase_n <- vapply(grase_files, function(f) nrow(read.table(file.path(BASE,f), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)), integer(1))
gr <- bind_rows(lapply(grase_files,
  function(f) read.table(file.path(BASE, f), header = TRUE, sep = "\t", quote = "",
    comment.char = "", stringsAsFactors = FALSE)[, c("gene","event","comparison","padj","padj_gene",
    "lfc_diff_net","delta_pi","setdiff1","setdiff2","transcripts1","transcripts2")]))
gr$src <- rep(c("internal","TSSTTS"), times = grase_n)
## UNIVERSE = EVERY enumerated bipartition test (339,684). A test whose gene was
## screened out at nested-BH stage 1 went through stage 1 and was rejected; that
## is a call, not an absence of testing. Both the stage-1 screen (padj == NA) and
## lfc_diff_net > 0 are therefore CALL conditions, applied via `stat` below.
## This matches how the other tools are treated -- MAJIQ's universe is every
## quantified junction whether called or not, DEXSeq's every bin with a padj.
## Filtering either condition into the universe removes GT positives from
## recall's denominator and inflates GrASE recall relative to those tools.
gr$sd  <- ifelse(grepl("diff1", gr$comparison), as.character(gr$setdiff1), as.character(gr$setdiff2))
gr$txs <- ifelse(grepl("diff1", gr$comparison), gr$transcripts1, gr$transcripts2)
gr$uid <- paste(gr$src, gr$gene, gr$event, gr$comparison, sep = "|")
gr$imp <- lapply(gr$txs, parse_l, sep = ",")
lfcok  <- !is.na(gr$lfc_diff_net) & gr$lfc_diff_net > 0
st <- gr$padj
st[is.na(st)] <- 1                    # gene screened out at stage 1: a CALL
st[!lfcok]    <- 1                    # lfc is a CALL condition
calls$GrASE <- list(
  gene = gr$gene, stat = st, gene_stat = gr$padj_gene,
  units = as.list(gr$uid), imp = gr$imp, lower_better = TRUE)
cat(sprintf("  GrASE universe %d tests; callable (padj computed AND lfc>0): %d\n",
            nrow(gr), sum(!is.na(gr$padj) & lfcok)))
## GrASE restricted to INTERNAL bubbles. MAJIQ cannot test a TSS/TTS boundary
## at all, so the strict like-for-like against a junction uses this entry --
## pooling TSS/TTS compares GrASE's extra reach against MAJIQ's absence rather
## than comparing the units.
## NOTE this is a UNIVERSE restriction, not a call-side filter, so it SUBSETS.
## The dpi variants above must NOT subset (they keep every tested unit and make
## the filtered ones un-callable) because a call-side filter does not change
## what was tested. Getting these two backwards puts recall out by ~12x.
ki <- gr$src == "internal"
sti <- gr$padj[ki]; sti[is.na(sti)] <- 1; sti[!lfcok[ki]] <- 1
calls$GrASE_internal <- list(
  gene = gr$gene[ki], stat = sti, gene_stat = gr$padj_gene[ki],
  units = as.list(gr$uid[ki]), imp = gr$imp[ki], lower_better = TRUE)
cat(sprintf("  GrASE_internal: %d of %d tests (universe restricted)\n", sum(ki), nrow(gr)))

## Call-side effect-size variants, symmetric to MAJIQ's built-in C:
##   GrASE_dpi0.1 <-> MAJIQ_C0.10      GrASE_dpi0.2 <-> MAJIQ_C0.20
## Without these the curves compare MAJIQ WITH an effect-size threshold against
## GrASE WITHOUT one, which is not a like-for-like operating characteristic.
for (dp in c(0.1, 0.2)) {
  k <- lfcok & !is.na(gr$padj) & !is.na(gr$delta_pi) & abs(gr$delta_pi) >= dp
  ## Keep EVERY tested unit in the call list and make filtered-out tests
  ## un-callable (padj set to 1, lower_better). Subsetting gr instead would
  ## shrink `all_units`, and recall's denominator (upos_restr) is computed from
  ## it -- so an aggressive filter would appear to RAISE recall, which is
  ## impossible for a subset of the same calls.
  st <- gr$padj; st[is.na(st)] <- 1; st[!k] <- 1
  calls[[sprintf("GrASE_dpi%.1f", dp)]] <- list(
    gene = gr$gene, stat = st, gene_stat = gr$padj_gene,
    units = as.list(gr$uid),
    imp = gr$imp, lower_better = TRUE)
  cat(sprintf("  GrASE_dpi%.1f: %d of %d tests callable (universe unchanged)\n",
              dp, sum(k), nrow(gr)))
}

## --- DEXSeq (bins) ---------------------------------------------------------
dex <- read.table(file.path(BASE, "DEXSeq/dexseq_group1_group2/all.group1_group2.dxd_filteredbyCountMultiExon.txt"),
                  header = FALSE, skip = 1, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
dp <- suppressWarnings(as.numeric(dex[[8]])); keep <- !is.na(dp)
dgene <- dex[[2]][keep]; dexon <- dex[[3]][keep]; dp <- dp[keep]
# bin -> transcripts (only needed for signal genes + nulls we might call)
tx_of_bin <- new.env()
for (g in unique(dgene)) {
  f <- sprintf("%s/results/sim_exon_info/%s.exonic_parts_fc.txt", BASE, g)
  d <- tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                           quote = "", comment.char = ""), error = function(e) NULL)
  if (is.null(d)) next
  for (i in seq_len(nrow(d)))
    assign(paste(g, d$exonic_part[i], sep = ":"), parse_l(d$transcripts[i], "\\+"), envir = tx_of_bin)
}
calls$DEXSeq <- list(gene = dgene, stat = dp, gene_stat = dp,
  units = as.list(paste(dgene, dexon, sep = ":")),
  imp = lapply(paste(dgene, dexon, sep = ":"), function(k)
    if (exists(k, envir = tx_of_bin)) get(k, envir = tx_of_bin) else character(0)),
  lower_better = TRUE)

## --- persisted junction keys: structure + per-junction transcript users -----
## make_lsv_junction_keys.R stores the transcript IDENTITIES per junction, so
## both MAJIQ blocks below JOIN instead of rebuilding the junction->transcript
## map from the voila tsv (which cost ~20 min per run and duplicated work the
## GT pipeline already does).
jk <- read.table(file.path(BASE, "results/gt/lsv_junction_keys.txt"), header = TRUE,
                 sep = "\t", quote = "", stringsAsFactors = FALSE)
stopifnot("users" %in% names(jk))
jk$uid <- paste(jk$lsv_id, jk$j_index, sep = "|")
jk$imp <- lapply(jk$users, function(x)
  if (is.na(x) || !nzchar(x)) character(0) else strsplit(x, ",", fixed = TRUE)[[1]])
cat(sprintf("  junction keys: %d junctions, %d with >=1 user\n",
            nrow(jk), sum(lengths(jk$imp) > 0)))

## per-junction probabilities at the other C setting (cheap: no transcript work)
prob_at <- function(Cv) {
  tv <- read.table(file.path(BASE, sprintf("majiq/majiq_deltapsi.thr%s.tsv", Cv)),
                   header = TRUE, sep = "\t", quote = "", comment.char = "#",
                   stringsAsFactors = FALSE)
  u <- unlist(lapply(seq_len(nrow(tv)), function(i) {
    pj <- suppressWarnings(as.numeric(strsplit(tv$probability_changing[i], "[;,]")[[1]]))
    if (!length(pj)) return(NULL)
    setNames(pj, paste(tv$lsv_id[i], seq_along(pj), sep = "|"))
  }))
  u
}
PROB <- list("0.20" = setNames(jk$prob, jk$uid), "0.10" = prob_at("0.10"))

## --- MAJIQ (LSVs; per-LSV max prob, transcripts of the max-prob junction) --
exons_gr <- import(file.path(BASE, "ref/gencode.v28.annotation.gtf"),
                   feature.type = "exon", colnames = c("transcript_id", "gene_id"))
{ o <- order(exons_gr$transcript_id, start(exons_gr))
  tid <- exons_gr$transcript_id[o]; s_ <- start(exons_gr)[o]; e_ <- end(exons_gr)[o]
  n <- length(tid); same <- tid[-n] == tid[-1]
  tx_junc <- split(paste(e_[-n], s_[-1], sep = "-")[same], tid[-n][same]) }
for (Cv in c("0.20", "0.10")) {
  pv <- PROB[[Cv]][jk$uid]; pv[is.na(pv)] <- 0
  d  <- data.frame(lsv_id = jk$lsv_id, gene = jk$gene, p = as.numeric(pv),
                   idx = seq_len(nrow(jk)), stringsAsFactors = FALSE)
  # per LSV: the max probability and the row index achieving it
  o    <- order(d$lsv_id, -d$p)
  dd   <- d[o, ]
  keep <- !duplicated(dd$lsv_id)
  top  <- dd[keep, ]
  calls[[paste0("MAJIQ_C", Cv)]] <- list(
    gene = top$gene, stat = top$p, gene_stat = top$p,
    units = as.list(top$lsv_id), imp = jk$imp[top$idx], lower_better = FALSE)
  cat(sprintf("  MAJIQ_C%s: %d LSVs (joined)\n", Cv, nrow(top)))
}

## --- MAJIQ at JUNCTION level -----------------------------------------------
## A bipartition and a junction are the same KIND of object: one alternative
## against a local reference, one df each. An LSV is a family of K-1 such tests,
## so bipartition-vs-LSV compares a test to a family. These entries put the
## junction in the sweep so the like-for-like comparison is a curve, not a
## hand-matched pair of operating points.
for (Cv in c("0.20", "0.10")) {
  pv <- PROB[[Cv]][jk$uid]; pv[is.na(pv)] <- 0
  calls[[paste0("MAJIQjunc_C", Cv)]] <- list(
    gene = jk$gene, stat = as.numeric(pv), gene_stat = as.numeric(pv),
    units = as.list(jk$uid), imp = jk$imp, lower_better = FALSE)
  cat(sprintf("  MAJIQjunc_C%s: %d junctions (joined)\n", Cv, nrow(jk)))
}

## --- rMATS (events; inc U skip transcripts) --------------------------------
csv_mean <- function(x) { v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
                          if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE) }
jgt <- read.table(file.path(BASE, "results.beforegatefix/sim_junction_gt.txt"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
jgt <- jgt[!(jgt$n_inc_tx == 0 | jgt$n_skip_tx == 0), ]
jgt$uid <- paste(jgt$event_type, jgt$ID, sep = ":")
rm_pos <- setNames(jgt$gt_positive, jgt$uid)
rg <- list(); rstat <- c(); rgene <- c(); rimp <- list(); ruid <- c()
for (et in c("SE","A3SS","A5SS","RI")) {
  d <- read.table(file.path(BASE, "rMATS/rmats_post_group1_group2", paste0(et, ".MATS.JCEC.txt")),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  names(d) <- make.unique(names(d)); d$GeneID <- gsub('"', "", d$GeneID)
  c1 <- sapply(d$IJC_SAMPLE_1, csv_mean) + sapply(d$SJC_SAMPLE_1, csv_mean)
  c2 <- sapply(d$IJC_SAMPLE_2, csv_mean) + sapply(d$SJC_SAMPLE_2, csv_mean)
  p1 <- sapply(d$IncLevel1, csv_mean); p2 <- sapply(d$IncLevel2, csv_mean)
  testable <- !is.na(c1) & c1 >= 10 & !is.na(c2) & c2 >= 10 &
              !is.na(p1) & p1 >= 0.05 & p1 <= 0.95 & !is.na(p2) & p2 >= 0.05 & p2 <= 0.95
  d <- d[testable, ]
  if (nrow(d) == 0) next
  for (i in seq_len(nrow(d))) {
    g <- d$GeneID[i]; gtx <- intersect(gene_tx[[g]], names(tx_junc))
    js <- if (et == "SE")
        c(paste(d$upstreamEE[i], d$exonStart_0base[i]+1, sep="-"),
          paste(d$exonEnd[i], d$downstreamES[i]+1, sep="-"),
          paste(d$upstreamEE[i], d$downstreamES[i]+1, sep="-"))
      else if (et %in% c("A3SS","A5SS"))
        c(paste(d$flankingEE[i], d$longExonStart_0base[i]+1, sep="-"),
          paste(d$flankingEE[i], d$shortES[i]+1, sep="-"),
          paste(d$longExonEnd[i], d$flankingES[i]+1, sep="-"),
          paste(d$shortEE[i], d$flankingES[i]+1, sep="-"))
      else paste(d$upstreamEE[i], d$downstreamES[i]+1, sep="-")
    rimp[[length(rimp)+1]] <- if (length(gtx)) gtx[vapply(tx_junc[gtx], function(v) any(js %in% v), logical(1))] else character(0)
    rgene <- c(rgene, g); rstat <- c(rstat, d$FDR[i]); ruid <- c(ruid, paste(et, d$ID[i], sep=":"))
  }
}
calls$rMATS <- list(gene = rgene, stat = rstat, gene_stat = rstat,
                    units = as.list(ruid), imp = rimp, lower_better = TRUE)
attr(calls, "pos_key") <- pos_key; attr(calls, "neg_key") <- neg_key
attr(calls, "lsv_pos") <- lsv_pos; attr(calls, "rm_pos") <- rm_pos
saveRDS(calls, CACHE)
}
# (GT tables and pos/neg keys are built above, outside the cache, so that the
#  sweep can subset by simulation category.)

## ===========================================================================
## Sweep each tool at three levels
## ===========================================================================
GRIDS <- list(GrASE = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_dpi0.1 = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_dpi0.2 = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_internal = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              MAJIQjunc_C0.20 = c(0.99,0.95,0.9,0.8,0.5,0.3,0.1),
              MAJIQjunc_C0.10 = c(0.99,0.95,0.9,0.8,0.5,0.3,0.1),
              DEXSeq = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              MAJIQ_C0.20 = c(0.99,0.95,0.9,0.8,0.5,0.3,0.1),
              MAJIQ_C0.10 = c(0.99,0.95,0.9,0.8,0.5,0.3,0.1),
              rMATS = c(1e-4,1e-3,0.01,0.05,0.1,0.2))
# FULL-universe positive counts per unit type (denominator = every GT-positive
# unit that exists, tested or not). RESTRICTED = positives among units the tool
# actually tested. Precision is identical in both (same TP/FP); only recall moves.
n_pos_parts_full <- sum(vapply(ls(pos_key), function(k) TRUE, logical(1)))
lsv_full <- read.table(file.path(BASE, "results/sim_lsv_gt.full.txt"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
n_pos_lsv_full <- sum(lsv_full$gt_positive)
n_pos_evt_full <- sum(rm_pos)                    # all annotated testable events
n_pos_gene_full <- length(sig_genes)
n_pos_tx_full   <- length(manip_tx)

rows <- list()
for (tool in names(calls)) {
 cl0 <- calls[[tool]]
 st0 <- get_st(cl0$gene)
 for (categ in c("ALL", "DTE", "DTU", "Background", "DGE")) {
  # subset the tool's calls to genes of this category. The DTE and DTU panels
  # include BACKGROUND genes as the NEGATIVE class (positives = DTE/DTU genes) so
  # gene precision is not 1 by construction -- a Background-gene call is a real
  # gene FP. This makes DTE/DTU proper signal-vs-null detection problems sharing
  # the Background null; they therefore no longer partition (Background is shared),
  # but the pooled ALL panel is unchanged. DGE stays its own (null) panel.
  keep <- if (categ == "ALL") rep(TRUE, length(cl0$gene))
          else if (categ == "DTE") st0 %in% c("DTE", "Background")
          else if (categ == "DTU") st0 %in% c("DTU", "Background")
          else st0 == categ
  if (!any(keep)) next
  cl <- list(gene = cl0$gene[keep], stat = cl0$stat[keep],
             units = cl0$units[keep], imp = cl0$imp[keep], lower_better = cl0$lower_better)
  cat_genes <- switch(categ, ALL = sig_genes, DTE = dte.genes, DTU = dtu.genes,
                      Background = character(0), DGE = character(0))
  is_null_cat <- categ %in% c("Background", "DGE")
  cat_tx <- if (categ == "ALL") manip_tx else
            intersect(manip_tx, unique(unlist(gene_tx[cat_genes])))
  all_units <- unique(unlist(cl$units))
  # restricted denominators (positives among units this tool tested in this category)
  upos_restr <- if (grepl("^MAJIQjunc", tool)) sum(jn_pos[intersect(all_units, names(jn_pos))], na.rm = TRUE)
                else if (grepl("^GrASE", tool)) sum(bp_pos[intersect(all_units, names(bp_pos))], na.rm = TRUE)
                else if (grepl("^MAJIQ", tool)) sum(lsv_pos_rule[intersect(all_units, names(lsv_pos_rule))], na.rm = TRUE)
                else if (tool == "rMATS")  sum(rm_pos[intersect(all_units, names(rm_pos))], na.rm = TRUE)
                else sum(vapply(all_units, function(k) exists(k, envir = pos_key), logical(1)))
  pick <- function(tabl, full_all) if (categ == "ALL") full_all else
            as.numeric(tabl[categ]) %||% 0
  `%||%` <- function(a, b) if (length(a) == 0 || is.na(a)) b else a
  ## GT_rule cannot label units the tool never quantified, so for GrASE and
  ## MAJIQ the full-universe unit denominator is set to the restricted one.
  upos_full  <- if (grepl("^GrASE", tool)) upos_restr
                else if (grepl("^MAJIQ", tool)) upos_restr
                else if (tool == "rMATS")  (if (categ=="ALL") n_pos_evt_full   else as.numeric(n_pos_evt_cat[categ]))
                else                       (if (categ=="ALL") n_pos_parts_full else as.numeric(n_pos_parts_cat[categ]))
  tested_genes <- unique(cl$gene)
  ## negative-class sizes, for TN. Each level has its own universe:
  ##   unit       units the tool tested (category-local, matching its FP rule)
  ##   transcript transcripts of the genes the tool tested
  ##   gene       genes the tool tested
  ## transcript_tolerant counts DECISIONS, so a "true negative decision" is not
  ## a defined object -- its TN is NA rather than a fabricated number.
  ## SCOPE MUST MATCH THE FP RULE. All four levels now partition by category,
  ## so every universe is the category subset (cl).
  n_units_tested  <- length(unique(unlist(cl$units)))
  n_tx_tested0    <- length(unique(unlist(gene_tx[tested_genes])))
  n_genes_tested0 <- length(tested_genes)
  gpos_restr <- length(intersect(tested_genes, cat_genes))
  tpos_restr <- length(intersect(unique(unlist(gene_tx[tested_genes])), cat_tx))
  n_pos_gene_c <- length(cat_genes); n_pos_tx_c <- length(cat_tx)
  for (thr in GRIDS[[tool]]) {
    sig <- if (cl$lower_better) cl$stat < thr else cl$stat >= thr
    if (!any(sig)) next
    ## --- UNIT ---
    su <- unique(unlist(cl$units[sig]))
    if (grepl("^MAJIQjunc", tool)) {
      v <- jn_pos[su]; v[is.na(v)] <- FALSE
      tp <- sum(v); fp <- sum(!v)
    } else if (grepl("^GrASE", tool)) {
      v <- bp_pos[su]; v[is.na(v)] <- FALSE      # null-gene tests are negatives
      tp <- sum(v); fp <- sum(!v)
    } else if (grepl("^MAJIQ", tool)) {
      v <- lsv_pos_rule[su]; tp <- sum(v, na.rm = TRUE); fp <- sum(!v, na.rm = TRUE)
    } else if (tool == "rMATS") {
      v <- rm_pos[su];  tp <- sum(v, na.rm = TRUE); fp <- sum(!v, na.rm = TRUE)
    } else {
      tp <- sum(vapply(su, function(k) exists(k, envir = pos_key), logical(1)))
      fp <- sum(vapply(su, function(k) exists(k, envir = neg_key), logical(1)))
    }
    ## --- TRANSCRIPT / GENE ---
    ## Each call is scored against the positives of its panel. In the DTE and DTU
    ## panels the negative class is BACKGROUND (see the `keep` block above), so a
    ## Background-gene call is a real gene FP and gene precision is NOT 1 by
    ## construction. (Previously the panels were pure DTE / pure DTU, which left
    ## every called gene a positive -> FP = 0 -> gene precision = 1, an artifact;
    ## sharing the Background null fixes it, at the cost of DTE/DTU no longer
    ## partitioning with Background. The pooled ALL panel is unchanged.)
    sig0 <- if (cl0$lower_better) cl0$stat < thr else cl0$stat >= thr
    stx  <- unique(unlist(cl$imp[sig]));  stx0 <- unique(unlist(cl0$imp[sig0]))
    sg   <- unique(cl$gene[sig]);         sg0  <- unique(cl0$gene[sig0])
    if (is_null_cat) {           # null categories: every call here is false
      ttp <- 0; tfp <- length(stx)
      gtp <- 0; gfp <- length(sg)
    } else {
      ttp <- length(intersect(stx, cat_tx));  tfp <- length(setdiff(stx, manip_tx))
      gtp <- length(intersect(sg, cat_genes)); gfp <- length(setdiff(sg, sig_genes))
    }
    ## --- transcript_tolerant: SAME transcript currency as strict; ONLY the FP
    ## rule differs. TP and recall are shared with strict (TP = ttp = manipulated
    ## transcripts implicated). A co-traveller is charged only if it rides a call
    ## that hit NO manipulated tx -- i.e. transcripts implicated exclusively by
    ## signal-missing calls. Bystanders of a genuine hit are forgiven, so
    ## tol_fp <= tfp and precision_tolerant >= precision_strict at the SAME recall.
    hit_c     <- vapply(cl$imp[sig], function(v) any(v %in% cat_tx), logical(1))
    impl_hit  <- unique(unlist(cl$imp[sig][ hit_c]))
    impl_miss <- unique(unlist(cl$imp[sig][!hit_c]))
    tol_fp    <- if (is_null_cat) length(stx)
                 else length(setdiff(setdiff(impl_miss, impl_hit), manip_tx))
    for (uni in c("restricted", "full")) {
      ## TN / n_neg exist only in the RESTRICTED universe. The full universe
      ## enlarges the positive count with units the tool never tested, and says
      ## nothing about how many negatives lie outside the tested set, so the
      ## confusion table cannot close there. NA rather than an invented number.
      isr <- uni == "restricted"
      du <- if (uni == "restricted") upos_restr else upos_full
      dg <- if (uni == "restricted") gpos_restr else n_pos_gene_c
      dt <- if (uni == "restricted") tpos_restr else n_pos_tx_c
      rows[[length(rows)+1]] <- data.frame(category=categ, level="unit", universe=uni, tool=tool, thr=thr,
        TP=tp, FP=fp, FN=max(du-tp,0), TN=if (isr) max(n_units_tested-tp-fp-max(du-tp,0),0) else NA_integer_,
        n_pos=du, n_neg=if (isr) max(n_units_tested-du,0) else NA_integer_,
        precision=tp/max(tp+fp,1), recall=tp/max(du,1))
      rows[[length(rows)+1]] <- data.frame(category=categ, level="transcript_strict", universe=uni, tool=tool, thr=thr,
        TP=ttp, FP=tfp, FN=max(dt-ttp,0), TN=if (isr) max(n_tx_tested0-ttp-tfp-max(dt-ttp,0),0) else NA_integer_,
        n_pos=dt, n_neg=if (isr) max(n_tx_tested0-dt,0) else NA_integer_,
        precision=ttp/max(ttp+tfp,1), recall=ttp/max(dt,1))
      rows[[length(rows)+1]] <- data.frame(category=categ, level="gene", universe=uni, tool=tool, thr=thr,
        TP=gtp, FP=gfp, FN=max(dg-gtp,0), TN=if (isr) max(n_genes_tested0-gtp-gfp-max(dg-gtp,0),0) else NA_integer_,
        n_pos=dg, n_neg=if (isr) max(n_genes_tested0-dg,0) else NA_integer_,
        precision=gtp/max(gtp+gfp,1), recall=gtp/max(dg,1))
      # transcript_tolerant: same transcript TP / recall / n_pos as strict; only
      # FP differs (tol_fp <= tfp -- co-travellers of a hit forgiven). Closed table
      # in transcript currency, defined in the FULL universe (TN/n_neg restricted-
      # only, as for the other levels).
      rows[[length(rows)+1]] <- data.frame(category=categ, level="transcript_tolerant", universe=uni, tool=tool, thr=thr,
        TP=ttp, FP=tol_fp, FN=max(dt-ttp,0),
        TN=if (isr) max(n_tx_tested0-ttp-tol_fp-max(dt-ttp,0),0) else NA_integer_,
        n_pos=dt, n_neg=if (isr) max(n_tx_tested0-dt,0) else NA_integer_,
        precision=ttp/max(ttp+tol_fp,1), recall=ttp/max(dt,1))
    }
  }
 }
}
tab <- bind_rows(rows)
tab$precision <- round(tab$precision, 4); tab$recall <- round(tab$recall, 4)
write.table(tab, file.path(OUT, "pr_three_levels.GT_rule.txt"), sep="\t", quote=FALSE, row.names=FALSE)
cat("=== PR operating points by category (GrASE*/DEXSeq/rMATS 0.01, MAJIQ 0.95) ===\n")
print(as.data.frame(tab %>% filter((grepl("^GrASE|^DEXSeq|^rMATS", tool) & thr == 0.01) |
                                   (grepl("MAJIQ", tool) & thr == 0.95)) %>%
      arrange(category, level, universe, tool) %>%
      select(category, level, universe, tool, TP, FP, n_pos, precision, recall)), row.names = FALSE)

## ===========================================================================
## Plots: one file per category; 2 rows (restricted/full) x 3 cols (levels)
## ===========================================================================
LEVELS <- c("unit","transcript_strict","transcript_tolerant","gene")
for (categ in c("ALL", "DTE", "DTU", "Background", "DGE")) {
if (!(categ %in% c("Background","DGE"))) {   # PR curves need positives
fn <- if (categ == "ALL") "plots/pr_three_levels.GT_rule.png" else sprintf("plots/pr_three_levels.GT_rule_%s.png", categ)
png(file.path(BASE, fn), width = 1900, height = 1000, res = 130)
par(mfrow = c(2,4), mar = c(4.2,4.2,3,1), oma = c(2.3,0,1.6,0))
COL <- c(GrASE="#ff7f0e", "GrASE_dpi0.1"="#ff7f0e", "GrASE_dpi0.2"="#ff7f0e",
         DEXSeq="#9467bd", MAJIQ_C0.20="#66c2a5", MAJIQ_C0.10="#66c2a5",
         "GrASE_internal"="#d62728",
         "MAJIQjunc_C0.20"="#1b7837", "MAJIQjunc_C0.10"="#1b7837",
         rMATS="#4c72b0")
# GrASE variants share a hue and separate by line type; same for the two MAJIQ C
# settings, so the effect-size knob reads as a variant rather than a new tool.
LTY <- c(GrASE=1, "GrASE_dpi0.1"=2, "GrASE_dpi0.2"=3,
         "GrASE_internal"=1,
         DEXSeq=1, MAJIQ_C0.20=1, MAJIQ_C0.10=2,
         "MAJIQjunc_C0.20"=1, "MAJIQjunc_C0.10"=2, rMATS=1)
cap <- function(x) paste0(toupper(substring(x,1,1)), substring(x,2))
## Panel labels name what each precision actually measures. The `level` column
## in the written table keeps its original values (unit/transcript/call/gene) so
## downstream joins are unaffected -- only the display changes.
##   transcript precision == 1 - (co-traveller rate among called transcripts).
##     It is a function of implicated-set SIZE, not of whether calls were right:
##     rMATS 0.132 with median |imp| 5, GrASE 0.270 with median 2.
##   call precision == fraction of DECISIONS landing on >=1 manipulated
##     transcript -- accuracy, tolerant of co-travellers.
##   Their RECALL is identical by construction (both: manipulated transcripts
##   implicated by >=1 call), so only the precisions differ.
LVLAB <- c(unit = "Unit", transcript_strict = "transcript_strict",
           transcript_tolerant = "transcript_tolerant", gene = "Gene")
lvlab <- function(x) ifelse(is.na(LVLAB[x]), cap(x), LVLAB[x])
for (uni in c("restricted","full")) {
  for (lv in LEVELS) {
    d <- tab[tab$level == lv & tab$universe == uni & tab$category == categ, ]
    ## UNIT is a WITHIN-TOOL panel: only GrASE-family variants share the
    ## bipartition unit and GT_rule_bipartition, so only they may be overlaid.
    ## Other tools' units (bins / LSVs / junctions / events) are scored against
    ## different GTs and are NOT comparable here -- cross-tool reading lives in the
    ## transcript/gene panels only.
    if (lv == "unit") d <- d[grepl("^GrASE", d$tool), ]
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="Recall", ylab="Precision",
         main=sprintf("%s (%s)%s", lvlab(lv), uni, if (lv == "unit") "  within-GrASE" else ""))
    abline(0,1,col="grey80",lty=3)
    for (tl in unique(d$tool)) {
      dd <- d[d$tool == tl, ]; dd <- dd[order(dd$recall), ]
      lines(dd$recall, dd$precision, col=COL[tl], lty=LTY[tl], lwd=2)
      points(dd$recall, dd$precision, col=COL[tl], pch=16, cex=0.8)
    }
    if (lv == "unit" && uni == "restricted") {          # unit legend: GrASE variants only
      gk <- grep("^GrASE", names(COL), value = TRUE)
      legend("bottomleft", legend=gk, col=COL[gk], lty=LTY[gk], lwd=2, bty="n", cex=0.7)
    }
    if (lv == "transcript_strict" && uni == "restricted")  # cross-tool legend on a fair panel
      legend("bottomleft", legend=names(COL), col=COL, lty=LTY, lwd=2, bty="n", cex=0.7)
  }
}
mtext(sprintf("Simulation category: %s", categ), outer = TRUE, cex = 0.9, font = 2)
mtext(paste0("unit = tested unit (within-GrASE only) | transcript_strict = charges every co-traveller | ",
             "transcript_tolerant = forgives co-travellers riding a hit (same recall) | ",
             "DTE/DTU panels use Background as the negative class; ALL is pooled"),
      side = 1, outer = TRUE, cex = 0.6, line = 0.15)
if (categ %in% c("DTE","DTU"))
  mtext(paste0("positives = ", categ, " genes; negatives = Background genes ",
               "(a Background call is a gene FP), so gene precision is informative here"),
        side = 1, outer = TRUE, cex = 0.6, line = 0.95, font = 3)
dev.off()
cat(sprintf("  wrote %s\n", fn))
}

## --- stacked TP / FN / FP bars at the operating points --------------------
fn2 <- sprintf("plots/confusion_three_levels.GT_rule%s.png",
               if (categ == "ALL") "" else paste0("_", categ))
png(file.path(BASE, fn2), width = 1900, height = 1000, res = 130)
par(mfrow = c(2,4), mar = c(4,8.5,3,2.5), oma = c(2.3,0,1.6,0))
op <- tab[tab$category == categ &
          ((grepl("^GrASE|^DEXSeq|^rMATS", tab$tool) & tab$thr == 0.01) |
           (grepl("MAJIQ", tab$tool) & tab$thr == 0.95)), ]
BARCOL <- c(TP = "#1d9e75", FN = "#888780", FP = "#e24b4a")
for (uni in c("restricted","full")) {
  for (lv in LEVELS) {
    d <- op[op$level == lv & op$universe == uni, ]
    ord <- names(COL)
    if (lv == "unit") ord <- grep("^GrASE", ord, value = TRUE)  # unit: within-GrASE only
    d <- d[match(ord, d$tool), ]; d <- d[!is.na(d$tool), ]
    m <- rbind(TP = d$TP, FN = pmax(d$n_pos - d$TP, 0), FP = d$FP)
    colnames(m) <- d$tool
    bp <- barplot(m, horiz = TRUE, las = 1, col = BARCOL, border = "white",
                  cex.names = 0.7, xlab = "count",
                  main = sprintf("%s (%s)%s", lvlab(lv), uni, if (lv == "unit") "  within-GrASE" else ""))
    tot <- colSums(m)
    text(tot, bp, labels = format(tot, big.mark = ","), pos = 4, cex = 0.6, xpd = NA)
    if (lv == "unit" && uni == "restricted")
      legend("bottomright", legend = c("TP","FN","FP"), fill = BARCOL,
             bty = "n", cex = 0.75, border = NA)
  }
}
mtext(sprintf("Confusion counts at operating points -- category: %s", categ),
      outer = TRUE, cex = 0.9, font = 2)
mtext(paste0("unit = tested unit (within-GrASE only) | transcript_strict = charges every co-traveller | ",
             "transcript_tolerant = forgives co-travellers riding a hit (same recall) | ",
             "DTE/DTU panels use Background as the negative class; ALL is pooled"),
      side = 1, outer = TRUE, cex = 0.6, line = 0.15)
if (categ %in% c("DTE","DTU"))
  mtext(paste0("positives = ", categ, " genes; negatives = Background genes ",
               "(a Background call is a gene FP), so gene precision is informative here"),
        side = 1, outer = TRUE, cex = 0.6, line = 0.95, font = 3)
dev.off()
cat(sprintf("  wrote %s\n", fn2))
}
cat(sprintf("\nTable: %s\n", file.path(OUT, "pr_three_levels.GT_rule.txt")))
