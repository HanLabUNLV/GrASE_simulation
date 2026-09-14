#!/usr/bin/env Rscript
#
# scripts/transcript_level_metrics.R
#
# Transcript-level metrics: LOCALIZATION. The simulation truth is native at
# this level (iso.dtu pairs, iso.dte singletons). None of the tools call
# transcripts directly, but every call implicates a transcript set:
#   GrASE : transcripts of the tested distinct path group (transcripts1/2)
#   DEXSeq: transcripts containing the significant bin (GT 'transcripts' col)
#   MAJIQ : transcripts using ANY junction of the significant LSV (whole LSV)
#   rMATS : whole-event calls cannot disambiguate a transcript; localization
#           precision only, via the junction GT's n_manip_inc/skip columns.
#
# Metrics, reported PER CATEGORY (DTE, DTU, Background, DGE, ALL):
#  1. Transcript recall (SHARED by both precisions): a manipulated transcript is
#     RECOVERED if some significant call's implicated set contains it.
#     Denominator = ALL manipulated transcripts of the category (no floor).
#  2. Two precisions bracketing the same recall:
#       prec_strict   -- among implicated transcripts, fraction manipulated
#                        (charges every co-traveller; sensitive to |imp| size)
#       prec_tolerant -- among calls, fraction hitting >= 1 manipulated tx
#                        (charges only calls that miss the signal entirely)
#  3. Sharpness: implicated-set size per call (median [IQR]) -- explains the
#     strict/tolerant gap.
# Null-gene (Background/DGE) significant calls ARE kept and charged as FP in both
# precisions and in the pooled ALL row (FULL-universe / null-FP convention);
# categories partition: DTE + DTU + Background + DGE = ALL.
#
# Significance rules as elsewhere: GrASE padj < 0.01 & lfc_diff_net > 0;
# DEXSeq padj < 0.01; MAJIQ per-junction P(|dPSI|>=0.20) >= 0.95;
# rMATS FDR < 0.01 + count/PSI filter.

suppressPackageStartupMessages({ library(dplyr); library(rtracklayer) })
BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
OUT_DIR <- file.path(BASE, "results/eval_metricfair")
load(file.path(BASE, "swimdown/simulate/data/simulate.rda"))
counts1 <- sim.counts.mat[, 1]; counts2 <- sim.counts.mat[, 2]
if (is.null(names(counts1))) { names(counts1) <- txdf$TXNAME; names(counts2) <- txdf$TXNAME }
FLOOR <- 0.1
strip_ver <- function(x) sub("\\.\\d+$", "", x)

## --- manipulated transcripts + their true proportion shift ------------------
manip_df <- bind_rows(
  data.frame(tx = names(which(iso.dtu)), sim_type = "DTU", stringsAsFactors = FALSE),
  data.frame(tx = names(which(iso.dte)), sim_type = "DTE", stringsAsFactors = FALSE))
manip_df$gene <- txdf$GENEID[match(manip_df$tx, txdf$TXNAME)]
# transcripts absent from the quantified counts matrix are unexpressed: count 0
cts <- function(v, tx) { x <- v[tx]; x[is.na(x)] <- 0; x }
gene_tot1 <- tapply(cts(counts1, txdf$TXNAME), txdf$GENEID, sum)
gene_tot2 <- tapply(cts(counts2, txdf$TXNAME), txdf$GENEID, sum)
manip_df$true_dprop <- abs(cts(counts2, manip_df$tx) / gene_tot2[manip_df$gene] -
                           cts(counts1, manip_df$tx) / gene_tot1[manip_df$gene])
manip_df$true_dprop[!is.finite(manip_df$true_dprop)] <- 0
# NO detectability floor at this level: the transcript IS the manipulated
# object (iso.dtu/iso.dte is the literal design); the realized proportion shift
# is an intended part of the difficulty distribution. It is kept only as a
# STRATIFICATION covariate (recall by shift magnitude), never as a GT filter.
manip_det <- manip_df
manip_det$shift_bin <- cut(manip_det$true_dprop,
                           breaks = c(-Inf, 0.05, 0.1, 0.2, Inf),
                           labels = c("<0.05", "0.05-0.1", "0.1-0.2", ">=0.2"))
cat(sprintf("Manipulated transcripts (all, no floor): %d DTU, %d DTE\n",
    sum(manip_df$sim_type == "DTU"), sum(manip_df$sim_type == "DTE")))
manip_by_gene <- split(manip_df$tx, manip_df$gene)   # all manip tx per gene (for precision)
manip_all     <- manip_df$tx                         # every manipulated transcript
get_st <- function(g) ifelse(g %in% dte.genes, "DTE", ifelse(g %in% dtu.genes, "DTU",
                     ifelse(g %in% dge.genes, "DGE", "Background")))
parse_list <- function(x, sep = "[,+]") { v <- strsplit(x, sep)[[1]]; trimws(v[nchar(trimws(v)) > 0 & trimws(v) != "NA"]) }

## metric helpers: calls = data.frame(gene, sim_type, imp) over ALL genes
## (signal + null). Reports per category (DTE, DTU, Background, DGE, ALL), all in
## TRANSCRIPT currency:
##   tx_recall     -- manip transcripts recovered by >= 1 sig call (SHARED)
##   prec_strict   -- TP/(TP + every implicated co-traveller)
##   prec_tolerant -- TP/(TP + co-travellers riding a call that hit NO manip tx;
##                    bystanders of a genuine hit forgiven)
##   sharp_*       -- implicated-set size per call (median [IQR])
## TP = manip transcripts implicated is shared by both precisions and by recall;
## strict and tolerant differ ONLY in the FP rule. Null-gene calls are kept and
## charged as FP.
score_tool <- function(calls, tool) {
  hit_call <- mapply(function(g, imp) any(imp %in% manip_by_gene[[g]]), calls$gene, calls$imp)
  sizes    <- vapply(calls$imp, length, integer(1))
  # transcript recall: manip tx recovered by a sig call on ITS OWN gene
  sigc <- calls[calls$sim_type %in% c("DTE", "DTU"), ]
  imp_by_gene <- split(sigc$imp, sigc$gene)
  rec <- manip_det
  rec$hit <- mapply(function(t, g) {
    s <- imp_by_gene[[g]]; if (is.null(s)) FALSE else any(vapply(s, function(v) t %in% v, logical(1)))
  }, rec$tx, rec$gene)

  block <- function(idx, cat_label, rec_sub) {
    imp <- calls$imp[idx]; nc <- length(imp); h <- hit_call[idx]
    implic    <- unique(unlist(imp))
    impl_hit  <- unique(unlist(imp[h]))
    impl_miss <- unique(unlist(imp[!h]))
    tp        <- length(intersect(implic, manip_all))
    fp_strict <- length(setdiff(implic, manip_all))
    fp_tol    <- length(setdiff(setdiff(impl_miss, impl_hit), manip_all))
    row <- data.frame(tool = tool, category = cat_label,
      n_manip_tx     = if (!is.null(rec_sub)) nrow(rec_sub) else NA_integer_,
      n_tx_recovered = if (!is.null(rec_sub)) sum(rec_sub$hit) else NA_integer_,
      tx_recall      = if (!is.null(rec_sub)) round(mean(rec_sub$hit), 4) else NA_real_,
      n_sig_calls = nc,
      ## calls_hit -- calls whose implicated set contains >=1 manipulated
      ## transcript (call currency, NOT transcript currency). calls_hit /
      ## n_sig_calls separates "did the call land" from "how many bystanders
      ## came along" (sharp_*), the two factors that jointly set prec_tolerant.
      ## NOTE this criterion is generous to WIDE implicated sets -- more
      ## transcripts implicated is more chances to contain a manipulated one --
      ## so read it alongside sharp_median, not alone.
      calls_hit = sum(h), hit_rate = round(sum(h) / max(nc, 1), 4),
      n_impl_tx = length(implic), n_impl_manip = tp,
      ## explicit TP/FP in TRANSCRIPT currency. TP is shared by both precisions
      ## (= n_impl_manip); the two FP columns are what differ:
      ##   FP_strict   -- every implicated non-manipulated transcript
      ##   FP_tolerant -- only those implicated EXCLUSIVELY by calls that hit
      ##                  no manipulated transcript (co-travellers of a genuine
      ##                  hit forgiven). FP_tolerant <= FP_strict always.
      ## NOTE these are per-TRANSCRIPT, so a wrong call is charged once per
      ## transcript it implicates -- hence far larger than a per-call FP count.
      TP = tp, FP_strict = fp_strict, FP_tolerant = fp_tol,
      prec_strict   = round(tp / max(tp + fp_strict, 1), 4),
      prec_tolerant = round(tp / max(tp + fp_tol, 1), 4),
      sharp_median = if (nc > 0) as.numeric(median(sizes[idx])) else NA_real_,
      sharp_q1     = if (nc > 0) as.numeric(quantile(sizes[idx], .25)) else NA_real_,
      sharp_q3     = if (nc > 0) as.numeric(quantile(sizes[idx], .75)) else NA_real_,
      stringsAsFactors = FALSE)
    # recall stratified by realized proportion shift (descriptive, NOT a filter)
    if (!is.null(rec_sub)) for (b in levels(manip_det$shift_bin)) {
      rb <- rec_sub[!is.na(rec_sub$shift_bin) & rec_sub$shift_bin == b, ]
      row[[paste0("recall_", b)]] <- round(if (nrow(rb) > 0) mean(rb$hit) else NA_real_, 4)
      row[[paste0("n_", b)]] <- nrow(rb)
    }
    row
  }

  out <- list()
  # DTE/DTU panels use Background as the NEGATIVE class (positives = DTE/DTU manip
  # transcripts; Background-gene calls are FP), matching pr_curves. Background/DGE
  # also keep their own null-panel rows; ALL is the pooled view.
  for (st in c("DTE", "DTU"))
    out[[st]] <- block(which(calls$sim_type %in% c(st, "Background")), st, rec[rec$sim_type == st, ])
  for (st in c("Background", "DGE")) {
    idx <- which(calls$sim_type == st)
    if (length(idx)) out[[st]] <- block(idx, st, NULL)
  }
  out[["ALL"]] <- block(seq_len(nrow(calls)), "ALL", rec)          # pooled incl. null FP
  bind_rows(out)
}

## TWO-STAGE THRESHOLD (see pr_curves_three_levels_gtrule.R for the full note).
## exontest.R's nested_BH screens genes at a FIXED alpha = 0.05 and then runs
## within-gene BH inside survivors, so the shipped padj already encodes a gene
## screen at 0.05 regardless of the unit threshold applied here. Filtering
## padj < 0.01 therefore accepts units from genes admitted at a LOOSER screen
## than the unit cut. Both stages are recoverable from the raw p-values, and the
## correct rule at level t is padj_gene < t AND padj_within < t, encoded as
## pmax(). On the internal file this is 972 significant units at 0.01 rather
## than 1,012.
two_stage <- function(pv, gene, padj_gene) {
  within <- ave(pv, gene, FUN = function(x) p.adjust(x, method = "BH"))
  pmax(within, padj_gene)
}

## --- GrASE ------------------------------------------------------------------
cat("GrASE...\n")
grase_src <- c(internal = "bipartition.test.fulldesign/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
               TSSTTS   = "bipartition.test.fulldesign/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
g_n <- vapply(grase_src, function(f) nrow(read.table(file.path(BASE, f), header = TRUE, sep = "\t",
              quote = "", comment.char = "", stringsAsFactors = FALSE)), integer(1))
g <- bind_rows(lapply(grase_src,
  function(f) read.table(file.path(BASE, f), header = TRUE, sep = "\t", quote = "",
                         comment.char = "", stringsAsFactors = FALSE)[,
    c("gene", "comparison", "padj", "padj_gene", "p.value", "lfc_diff_net", "delta_pi",
      "transcripts1", "transcripts2")]))
g$src <- rep(names(grase_src), times = g_n)
## base call condition: BOTH nested_BH stages below 0.01, AND lfc_diff_net > 0
g$padj2 <- two_stage(g$p.value, g$gene, g$padj_gene)
g <- g[!is.na(g$padj2) & g$padj2 < 0.01 & !is.na(g$lfc_diff_net) & g$lfc_diff_net > 0, ]
g$imp <- lapply(ifelse(grepl("diff1", g$comparison), g$transcripts1, g$transcripts2), parse_list)
g$sim_type <- get_st(g$gene)
grase_tab <- score_tool(g[, c("gene", "sim_type", "imp")], "GrASE")
## Call-side effect-size variants, symmetric to MAJIQ's built-in C
## (GrASE_dpi0.1 <-> MAJIQ_C0.10, GrASE_dpi0.2 <-> MAJIQ_C0.20). Without these
## the table compares MAJIQ WITH an effect-size threshold against GrASE WITHOUT
## one, which is not a like-for-like operating point.
for (dp in c(0.1, 0.2)) {
  gk <- g[!is.na(g$delta_pi) & abs(g$delta_pi) >= dp, ]
  grase_tab <- bind_rows(grase_tab,
    score_tool(gk[, c("gene", "sim_type", "imp")], sprintf("GrASE_dpi%.1f", dp)))
}
## GrASE on the STRAND-RECONSTRUCTED counts, exonic unit, no junction merge.
## Same bipartitions, same model, same call rule -- only the counts differ:
## dexseq_count --stranded reverse on an unstranded library (half depth,
## antisense contamination) vs strand recovered from read origin.
## See ../stranded_reconstruction.md.
str_src <- c(internal = "bipartition.internal.stranded.test.EBapprox/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
             TSSTTS   = "bipartition.TSSTTS.stranded.test.EBapprox/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
if (all(file.exists(file.path(BASE, str_src)))) {
  sn <- vapply(str_src, function(f) nrow(read.table(file.path(BASE, f), header = TRUE,
        sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)), integer(1))
  sg <- bind_rows(lapply(str_src,
    function(f) read.table(file.path(BASE, f), header = TRUE, sep = "\t", quote = "",
      comment.char = "", stringsAsFactors = FALSE)[,
      c("gene","comparison","padj","padj_gene","p.value","lfc_diff_net","delta_pi",
        "transcripts1","transcripts2")]))
  sg$src <- rep(names(str_src), times = sn)
  sg$padj2 <- two_stage(sg$p.value, sg$gene, sg$padj_gene)
  sg <- sg[!is.na(sg$padj2) & sg$padj2 < 0.01 &
           !is.na(sg$lfc_diff_net) & sg$lfc_diff_net > 0, ]
  sg$imp <- lapply(ifelse(grepl("diff1", sg$comparison), sg$transcripts1, sg$transcripts2), parse_list)
  sg$sim_type <- get_st(sg$gene)
  K <- c("gene","sim_type","imp")
  grase_tab <- bind_rows(grase_tab, score_tool(sg[, K], "GrASE_stranded"))
  for (dp in c(0.1, 0.2)) {
    sk <- sg[!is.na(sg$delta_pi) & abs(sg$delta_pi) >= dp, ]
    grase_tab <- bind_rows(grase_tab,
      score_tool(sk[, K], sprintf("GrASE_stranded_dpi%.1f", dp)))
  }
  grase_tab <- bind_rows(grase_tab,
    score_tool(sg[sg$src == "internal", K], "GrASE_stranded_internal"))
}

## internal-only bubbles: the strict like-for-like against a junction, since
## MAJIQ cannot test a TSS/TTS boundary at all.
gi <- g[g$src == "internal", ]
grase_tab <- bind_rows(grase_tab,
  score_tool(gi[, c("gene", "sim_type", "imp")], "GrASE_internal"))

## GrASE with junction-merged bipartitions, BOTH types. A side with no
## exclusive exonic part emits no test in the exonic pipeline; where that side
## owns an exclusive intron, the merged counts substitute its split-read count.
## The rule is the same for internal and TSS/TTS (merge_exon_sj_counts.R fires
## on is.na(setdiff) & !is.na(intron_distinct)). It was first applied only to
## internal on the assumption that a TSS/TTS empty side is always a nested
## terminal path with no unique junction; that holds for part of them, but a
## large share are alternative donors/acceptors riding inside a terminal bubble
## and DO own an exclusive junction (scan: scripts/tsstts_empty_scan.py).
merged_src <- c(
  internal = "bipartition.merged.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt",
  TSSTTS   = "bipartition.merged.TSSTTS.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt")
read_merged <- function(f) {
  p <- file.path(BASE, f)
  if (!file.exists(p)) return(NULL)
  d <- read.table(p, header = TRUE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE)
  d$padj2 <- two_stage(d$p.value, d$gene, d$padj_gene)
  d <- d[!is.na(d$padj2) & d$padj2 < 0.01 &
         !is.na(d$lfc_diff_net) & d$lfc_diff_net > 0, ]
  d$imp <- lapply(ifelse(grepl("diff1", d$comparison), d$transcripts1, d$transcripts2), parse_list)
  d$sim_type <- get_st(d$gene)
  d$src <- if (grepl("TSSTTS", f)) "TSSTTS" else "internal"
  d[, c("gene", "sim_type", "src", "delta_pi", "imp")]
}
mg_int <- read_merged(merged_src[["internal"]])
mg_tss <- read_merged(merged_src[["TSSTTS"]])
KEEP <- c("gene", "sim_type", "imp")
## internal merged + exonic TSS/TTS -- the earlier operating point, kept so the
## TSS/TTS extension can be read off as a separate increment.
grase_tab <- bind_rows(grase_tab, score_tool(
  bind_rows(mg_int[, KEEP], g[g$src == "TSSTTS", KEEP]), "GrASE_merged"))
if (!is.null(mg_tss)) {
  mg <- bind_rows(mg_int, mg_tss)
  grase_tab <- bind_rows(grase_tab, score_tool(mg[, KEEP], "GrASE_merged_all"))
  ## Effect-size variants on the merged unit, so every matched triple has a
  ## merged member:  GrASE_merged_dpi0.1 <-> MAJIQ_C0.10 <-> rMATS_dpsi0.1, etc.
  for (dp in c(0.1, 0.2)) {
    mk <- mg[!is.na(mg$delta_pi) & abs(mg$delta_pi) >= dp, ]
    grase_tab <- bind_rows(grase_tab,
      score_tool(mk[, KEEP], sprintf("GrASE_merged_dpi%.1f", dp)))
  }
  ## internal-only merged: the like-for-like against a junction tool, since
  ## MAJIQ cannot test a TSS/TTS boundary at all.
  grase_tab <- bind_rows(grase_tab,
    score_tool(mg_int[, KEEP], "GrASE_merged_internal"))
}
## same merged unit, on the STRAND-RECONSTRUCTED counts
mstr <- c(internal = "bipartition.merged.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt",
          TSSTTS   = "bipartition.merged.TSSTTS.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt")
if (all(file.exists(file.path(BASE, mstr)))) {
  ms <- bind_rows(lapply(names(mstr), function(k) read_merged(mstr[[k]])))
  KEEP <- c("gene","sim_type","imp")
  grase_tab <- bind_rows(grase_tab, score_tool(ms[, KEEP], "GrASE_merged_stranded"))
  for (dp in c(0.1, 0.2)) {
    mk <- ms[!is.na(ms$delta_pi) & abs(ms$delta_pi) >= dp, ]
    grase_tab <- bind_rows(grase_tab,
      score_tool(mk[, KEEP], sprintf("GrASE_merged_stranded_dpi%.1f", dp)))
  }
  ## internal-only merged+stranded: the strict like-for-like against a junction
  ## tool, since MAJIQ cannot test a TSS/TTS boundary at all.
  grase_tab <- bind_rows(grase_tab,
    score_tool(ms[ms$src == "internal", KEEP], "GrASE_merged_stranded_internal"))

  ## Same merged unit and the same raw p-values, but PLAIN BH instead of
  ## exontest.R's nested_BH gene screen -- matching DEXSeq's FDR architecture.
  ## Under nested_BH only 1,345 null-gene units are reachable at any threshold;
  ## under plain BH 209,981 are. The two agree closely at padj 0.01 and diverge
  ## as it loosens, so carrying both makes the architecture explicit.
  bh_raw <- bind_rows(lapply(names(mstr), function(k) {
    d <- read.table(file.path(BASE, mstr[[k]]), header = TRUE, sep = "\t",
                    quote = "", comment.char = "", stringsAsFactors = FALSE)
    ## source/sink are node ids: integer in the internal file, character ("R"/"L")
    ## in TSS/TTS, so bind_rows needs them coerced before stacking.
    d$source <- as.character(d$source); d$sink <- as.character(d$sink)
    d$src <- k; d }))
  bh_raw$padj_BH <- p.adjust(bh_raw$p.value, method = "BH")
  bh <- bh_raw[!is.na(bh_raw$padj_BH) & bh_raw$padj_BH < 0.01 &
               !is.na(bh_raw$lfc_diff_net) & bh_raw$lfc_diff_net > 0, ]
  bh$imp <- lapply(ifelse(grepl("diff1", bh$comparison), bh$transcripts1, bh$transcripts2), parse_list)
  bh$sim_type <- get_st(bh$gene)
  grase_tab <- bind_rows(grase_tab, score_tool(bh[, KEEP], "GrASE_BH"))
}

## --- DEXSeq -----------------------------------------------------------------
cat("DEXSeq...\n")
## Scored for both count sources: the original (dexseq_count --stranded reverse
## on an unstranded library) and the strand-reconstructed rerun.
dexseq_tab <- NULL
for (dxcfg in list(c("DEXSeq",          "DEXSeq/dexseq_group1_group2/all.group1_group2.dxd_filteredbyCountMultiExon.txt"),
                   c("DEXSeq_stranded", "DEXSeq/dexseq_group1_group2_stranded/all.group1_group2.dxd_filteredbyCountMultiExon.txt"))) {
dxlab <- dxcfg[1]; dxpath <- dxcfg[2]
if (!file.exists(file.path(BASE, dxpath))) next
dex <- read.table(file.path(BASE, dxpath),
                  header = FALSE, skip = 1, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
dexp <- suppressWarnings(as.numeric(dex[[8]]))
dd <- data.frame(gene = dex[[2]], exon = dex[[3]], stringsAsFactors = FALSE)[!is.na(dexp) & dexp < 0.01, ]
dd$sim_type <- get_st(dd$gene)
# null-gene (Background/DGE) calls are KEPT and charged as FP downstream
# bin -> transcripts from the GT files of those genes
tx_of_bin <- new.env()
for (gene in unique(dd$gene)) {
  f <- sprintf("%s/results/sim_exon_info/%s.exonic_parts_fc.txt", BASE, gene)
  d <- tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                           quote = "", comment.char = ""), error = function(e) NULL)
  if (is.null(d)) next
  for (i in seq_len(nrow(d)))
    assign(paste(gene, d$exonic_part[i], sep = ":"), parse_list(d$transcripts[i], "\\+"), envir = tx_of_bin)
}
dd$imp <- lapply(paste(dd$gene, dd$exon, sep = ":"), function(k)
  if (exists(k, envir = tx_of_bin)) get(k, envir = tx_of_bin) else character(0))
dexseq_tab <- bind_rows(dexseq_tab, score_tool(dd[, c("gene", "sim_type", "imp")], dxlab))
}

## --- MAJIQ (LSV is the unit; whole-LSV implication) -------------------------
cat("MAJIQ...\n")
exons_gr <- import(file.path(BASE, "ref/gencode.v28.annotation.gtf"),
                   feature.type = "exon", colnames = c("transcript_id", "gene_id"))
{ o <- order(exons_gr$transcript_id, start(exons_gr))
  tid <- exons_gr$transcript_id[o]; st_ <- start(exons_gr)[o]; en_ <- end(exons_gr)[o]
  n <- length(tid); same <- tid[-n] == tid[-1]
  jstr <- paste(en_[-n], st_[-1], sep = "-")[same]
  tx_junc <- split(jstr, tid[-n][same]) }
majiq_tab <- NULL
## Scored for both BAM sources: the original unsplit BAMs and the
## origin-split ones. MAJIQ is largely but not wholly immune to the antisense
## contamination -- junctions are protected by splice-site matching, intron
## retention is quantified from strand-blind intronic coverage.
majiq_cfg <- list(c("", "majiq/majiq_deltapsi.thr%s.tsv"),
                  c("_stranded", "majiq/majiq_deltapsi.stranded.thr%s.tsv"))
for (mcfg in majiq_cfg) {
mlab <- mcfg[1]; mpat <- mcfg[2]
for (Cv in c("0.20", "0.10")) {
if (!file.exists(file.path(BASE, sprintf(mpat, Cv)))) next
tsv <- read.table(file.path(BASE, sprintf(mpat, Cv)), header = TRUE,
                  sep = "\t", quote = "", comment.char = "#", stringsAsFactors = FALSE)
tsv$sim_type <- get_st(tsv$gene_id)
# null-gene LSVs KEPT -- a significant LSV in a null gene is an FP call.
# LSV is the unit: significant iff its MAX junction clears C; when significant it
# implicates the WHOLE LSV (union of ALL its junctions' users), regardless of
# per-junction significance -- matching pr_curves MAJIQ_C and rMATS inc+skip.
mrows <- list()
for (i in seq_len(nrow(tsv))) {
  pj <- suppressWarnings(as.numeric(strsplit(tsv$probability_changing[i], "[;,]")[[1]]))
  lj <- parse_list(gsub(";", ",", tsv$junctions_coords[i]))
  K <- min(length(pj), length(lj)); if (K == 0) next
  mx <- suppressWarnings(max(pj[1:K], na.rm = TRUE))
  if (!is.finite(mx) || mx < 0.95) next
  gtx <- intersect(txdf$TXNAME[txdf$GENEID == tsv$gene_id[i]], names(tx_junc))
  users <- unique(unlist(lapply(lj[1:K], function(jc)
    gtx[vapply(tx_junc[gtx], function(v) jc %in% v, logical(1))])))
  mrows[[length(mrows) + 1]] <- list(gene = tsv$gene_id[i], sim_type = tsv$sim_type[i], imp = users)
}
mj <- data.frame(gene = vapply(mrows, `[[`, "", "gene"),
                 sim_type = vapply(mrows, `[[`, "", "sim_type"), stringsAsFactors = FALSE)
mj$imp <- lapply(mrows, `[[`, "imp")
majiq_tab <- bind_rows(majiq_tab, score_tool(mj, paste0("MAJIQ_C", Cv, mlab)))
}
}

## --- rMATS: implicated set = inc_tx UNION skip_tx of the significant event --
## (symmetric with MAJIQ: mirrored LSV significance implicates both sides'
## users, so whole-event inc+skip implication is the same logic.)
cat("rMATS (inc+skip implicated sets)...\n")
# per-transcript exon spans (for RI inclusion side)
tx_ex <- split(data.frame(s = start(exons_gr), e = end(exons_gr)), exons_gr$transcript_id)
add_ev <- function(gene, jset, span = NULL) {
  gtx <- intersect(txdf$TXNAME[txdf$GENEID == gene], names(tx_junc))
  imp <- gtx[vapply(tx_junc[gtx], function(v) any(jset %in% v), logical(1))]
  if (!is.null(span)) {   # RI: add transcripts with an exon spanning the intron
    gtx2 <- intersect(txdf$TXNAME[txdf$GENEID == gene], names(tx_ex))
    sp <- gtx2[vapply(tx_ex[gtx2], function(d) any(d$s <= span[1] & d$e >= span[2]), logical(1))]
    imp <- union(imp, sp)
  }
  imp
}
csv_mean <- function(x) { v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
                          if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE) }
## Reported twice: the documented per-group PSI/count filter ON (as rMATS's
## docs recommend) and OFF. The filter is a RAW-PSI boundary guard, not an
## effect-size threshold, so it is not the analogue of GrASE_dpi / MAJIQ_C --
## both rows are given so its cost is visible.
rmats_tab <- NULL
## variants: documented gate ON (base), gate OFF, and call-side effect-size
## thresholds on rMATS' own |IncLevelDifference| -- the true analogue of
## GrASE_dpi / MAJIQ_C (the PSI gate is a raw-PSI boundary guard, not one).
## Scored on BOTH count sources, matching every other tool in this table.
## The stranded rMATS run is split by strand (rMATS renumbers events from 0 in
## each half-run, so the halves are separate post dirs and are simply pooled
## here -- a gene contributes to one strand only). Directory convention mirrors
## rm_dirs in pr_curves_three_levels_gtrule.R so the two scripts agree.
RM_UNSTR <- c(all   = "rMATS/rmats_post_group1_group2")
RM_STR   <- c(plus  = "rMATS/stranded_plus/post",
              minus = "rMATS/stranded_minus/post")
rm_variants <- list(list(nm = "rMATS",                  gate = TRUE, dp = 0,   dirs = RM_UNSTR),
                    list(nm = "rMATS_dpsi0.1",          gate = TRUE, dp = 0.1, dirs = RM_UNSTR),
                    list(nm = "rMATS_dpsi0.2",          gate = TRUE, dp = 0.2, dirs = RM_UNSTR),
                    list(nm = "rMATS_stranded",         gate = TRUE, dp = 0,   dirs = RM_STR),
                    list(nm = "rMATS_stranded_dpsi0.1", gate = TRUE, dp = 0.1, dirs = RM_STR),
                    list(nm = "rMATS_stranded_dpsi0.2", gate = TRUE, dp = 0.2, dirs = RM_STR))
for (VAR in rm_variants) {
GATE <- VAR$gate
if (!any(file.exists(file.path(BASE, VAR$dirs)))) {
  cat(sprintf("  %s: no post dir, skipped\n", VAR$nm)); next
}
rmats_imp_rows <- list()
for (rmk in names(VAR$dirs)) {
for (et in c("SE", "A3SS", "A5SS", "RI")) {
  rmf <- file.path(BASE, VAR$dirs[[rmk]], paste0(et, ".MATS.JCEC.txt"))
  if (!file.exists(rmf)) next
  d <- read.table(rmf, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  names(d) <- make.unique(names(d)); d$GeneID <- gsub('"', "", d$GeneID)
  c1 <- sapply(d$IJC_SAMPLE_1, csv_mean) + sapply(d$SJC_SAMPLE_1, csv_mean)
  c2 <- sapply(d$IJC_SAMPLE_2, csv_mean) + sapply(d$SJC_SAMPLE_2, csv_mean)
  p1 <- sapply(d$IncLevel1, csv_mean); p2 <- sapply(d$IncLevel2, csv_mean)
  ok <- !is.na(d$FDR) & d$FDR < 0.01
  if (GATE) ok <- ok & !is.na(c1) & c1 >= 10 & !is.na(c2) & c2 >= 10 &
        !is.na(p1) & p1 >= 0.05 & p1 <= 0.95 & !is.na(p2) & p2 >= 0.05 & p2 <= 0.95
  if (VAR$dp > 0) { dps <- suppressWarnings(as.numeric(d$IncLevelDifference))
                    dps[is.na(dps)] <- 0; ok <- ok & abs(dps) >= VAR$dp }
  d <- d[ok, ]                       # null-gene events KEPT -- charged as FP
  if (nrow(d) == 0) next
  for (i in seq_len(nrow(d))) {
    gene <- d$GeneID[i]
    if (et == "SE") {
      jset <- c(paste(d$upstreamEE[i], d$exonStart_0base[i] + 1, sep = "-"),
                paste(d$exonEnd[i],    d$downstreamES[i] + 1,    sep = "-"),
                paste(d$upstreamEE[i], d$downstreamES[i] + 1,    sep = "-"))
      imp <- add_ev(gene, jset)
    } else if (et %in% c("A3SS", "A5SS")) {
      # implicate via all four possible flanking<->long/short junction strings;
      # genomic left-right orientation covers both strands
      jset <- c(paste(d$flankingEE[i], d$longExonStart_0base[i] + 1, sep = "-"),
                paste(d$flankingEE[i], d$shortES[i] + 1,             sep = "-"),
                paste(d$longExonEnd[i], d$flankingES[i] + 1,         sep = "-"),
                paste(d$shortEE[i],     d$flankingES[i] + 1,         sep = "-"))
      imp <- add_ev(gene, jset)
    } else {  # RI: skip junction + transcripts whose exon spans the intron
      jset <- paste(d$upstreamEE[i], d$downstreamES[i] + 1, sep = "-")
      imp <- add_ev(gene, jset, span = c(d$upstreamEE[i], d$downstreamES[i] + 1))
    }
    rmats_imp_rows[[length(rmats_imp_rows) + 1]] <-
      list(gene = gene, sim_type = get_st(gene), imp = imp)
  }
}
}
if (!length(rmats_imp_rows)) { cat(sprintf("  %s: no calls\n", VAR$nm)); next }
rm_df <- data.frame(gene = vapply(rmats_imp_rows, `[[`, "", "gene"),
                    sim_type = vapply(rmats_imp_rows, `[[`, "", "sim_type"),
                    stringsAsFactors = FALSE)
rm_df$imp <- lapply(rmats_imp_rows, `[[`, "imp")
rmats_tab <- bind_rows(rmats_tab, score_tool(rm_df, VAR$nm))
}

## --- report ------------------------------------------------------------------
tab <- bind_rows(grase_tab, dexseq_tab, majiq_tab, rmats_tab)
tab$category <- factor(tab$category, levels = c("DTE", "DTU", "Background", "DGE", "ALL"))
tab <- tab %>% arrange(category, tool)
write.table(tab, file.path(OUT_DIR, "transcript_level_metrics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== Transcript-level metrics (localization) ===\n")
print(as.data.frame(tab), row.names = FALSE)
cat(sprintf("\nWritten: %s\n", file.path(OUT_DIR, "transcript_level_metrics.txt")))
cat("NOTE: all four tools scored symmetrically in TRANSCRIPT currency. TP = manip\n",
    "transcripts implicated; recall is shared; prec_strict charges every co-\n",
    "traveller, prec_tolerant only co-travellers riding a call that hit nothing.\n",
    "DTE/DTU panels use Background as the negative class (Background calls are\n",
    "FP); Background/DGE also kept as null panels; ALL is pooled. Implication:\n",
    "GrASE tested distinct path group;",
    "group; DEXSeq bin membership; MAJIQ/rMATS junction users. Recall over ALL\n",
    "manip transcripts (no floor); shift bins are descriptive strata.\n")
