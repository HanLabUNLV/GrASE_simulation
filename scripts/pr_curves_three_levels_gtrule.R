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
bg_genes <- setdiff(names(gene_tx), c(dte.genes, dtu.genes, dge.genes))  # Background genes (full-universe negatives)

## --- GT tables (cheap reads; needed on both cached and uncached paths so the
##     sweep can be subset by simulation category) ----------------------------
ep <- readRDS(file.path(OUT, "ep_gt_cache.rds"))
ep$key <- paste(ep$gene, ep$exonic_part, sep = ":")
pos_key <- new.env(); for (k in ep$key[ep$gt_pos]) assign(k, TRUE, envir = pos_key)
neg_key <- new.env(); for (k in ep$key[ep$gt_neg]) assign(k, TRUE, envir = neg_key)
lsvgt   <- read.table(file.path(BASE, "results/sim_lsv_gt.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
lsv_pos <- setNames(lsvgt$gt_positive, lsvgt$lsv_id)
lsv_full<- read.table(file.path(BASE, "results/sim_lsv_gt.full.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
## --- strand-reconstruction switch -------------------------------------------
## STRANDED=1 scores every tool on the strand-corrected inputs: exonic counts
## recounted with the read's true strand of origin (see stranded_reconstruction.md),
## MAJIQ and rMATS rerun on the origin-split BAMs, DEXSeq rerun on the new bins.
## Unset reproduces the original run exactly.
STRANDED <- nzchar(Sys.getenv("STRANDED"))
if (STRANDED) cat("*** STRANDED inputs ***\n")
## rMATS renumbers its events from 0 in every run, so the stranded plus/minus
## IDs collide with each other and with the original numbering that
## results.beforegatefix/sim_junction_gt.txt is keyed to. Event IDs are
## strand-prefixed below to keep them unique, which means they no longer join to
## that GT -- rMATS UNIT rows are therefore emitted as NA under STRANDED rather
## than silently scored against stale truth. Transcript and gene levels are
## unaffected (they key on GeneID and coordinates). Rebuild with
## scripts/infer_rmats_junctions_gt.R to restore unit-level rMATS.
## Once scripts/run_infer_junctions_gt_stranded.sh has rebuilt the junction GT
## against the new IDs, it is picked up automatically and unit rows become real.
RMATS_STR_GT <- file.path(BASE, "results/sim_junction_gt.stranded.txt")
RMATS_UNIT_NA <- STRANDED && !file.exists(RMATS_STR_GT)
if (STRANDED && !RMATS_UNIT_NA) cat("  rMATS unit GT: using stranded rebuild\n")

## --- GT_rule unit truth -----------------------------------------------------
bpgt <- read.table(file.path(BASE, if (STRANDED) "results/gt/bipartition_gt.stranded.txt"
                                   else "results/gt/bipartition_gt.txt"), header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
bpgt$uid <- paste(bpgt$src, bpgt$gene, bpgt$event, bpgt$comparison, sep = "|")
bp_pos <- setNames(bpgt$GT_rule_bipartition, bpgt$uid)
n_pos_bp_cat <- table(bpgt$sim_type[bpgt$GT_rule_bipartition])
## Merged (exon + split-read) unit truth, same rule, built by
## scripts/infer_bipartition_gt.R --merged. bipartition_gt.txt cannot label a side
## whose exonic distinct set is empty -- it drops those rows -- so scoring
## GrASE_merged_all against bp_pos would silently call every junction-sourced
## unit a negative (v[is.na(v)] <- FALSE below). Each GrASE variant is scored
## against the truth built on the unit it actually tested.
BPMGT <- file.path(BASE, if (STRANDED) "results/gt/bipartition_merged_gt.stranded.txt"
                         else "results/gt/bipartition_merged_gt.txt")
bp_pos_merged <- NULL
if (file.exists(BPMGT)) {
  bpm <- read.table(BPMGT, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  bpm$uid <- paste(bpm$src, bpm$gene, bpm$event, bpm$comparison, sep = "|")
  bp_pos_merged <- setNames(bpm$GT_rule_bipartition, bpm$uid)
  cat(sprintf("merged unit GT: %d tests (%d exon-sourced, %d sj-sourced)\n",
              nrow(bpm), sum(bpm$D_source == "exon"), sum(bpm$D_source == "sj")))
}
## unit truth map for a GrASE variant
bp_map <- function(tool) {
  if (grepl("^GrASE_merged", tool) || identical(tool, "GrASE_BH")) {
    if (is.null(bp_pos_merged))
      stop("GrASE_merged_all needs results/gt/bipartition_merged_gt.txt -- ",
           "run scripts/infer_bipartition_gt.R --merged first")
    return(bp_pos_merged)
  }
  bp_pos
}
jngt <- read.table(file.path(BASE, "results/gt/lsv_junction_gt.txt"), header = TRUE,
                   sep = "\t", quote = "", stringsAsFactors = FALSE)
jngt$uid <- paste(jngt$lsv_id, jngt$j_index, sep = "|")
jn_pos <- setNames(jngt$GT_rule_junction, jngt$uid)
n_pos_jn_cat <- table(jngt$sim_type[jngt$GT_rule_junction])
## LSV universe and labels from the origin-split build when STRANDED: the
## universe is data-dependent (whatever LSVs MAJIQ found), so reusing the old
## build's 50,489 charges the stranded run 55 GT-positive LSVs it never tested.
lsvr <- read.table(file.path(BASE, if (STRANDED) "results/gt/lsv_gt.stranded.txt"
                                   else "results/gt/lsv_gt.txt"), header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
lsv_pos_rule <- setNames(lsvr$GT_rule_lsv, lsvr$lsv_id)
n_pos_lsvrule_cat <- table(lsvr$sim_type[lsvr$GT_rule_lsv])
## This is the AUTHORITATIVE rm_pos. It must live OUTSIDE the call-cache gate
## below: the rMATS block re-derives it, but that block is skipped whenever the
## cache is reused, which silently left rm_pos keyed on the ORIGINAL run's
## unprefixed ids ("SE:123") while calls$rMATS$units held strand-prefixed ones
## ("plus:SE:39"). Every lookup then returned NA and all rMATS unit rows came
## out as TP 0 / FP 0 / n_pos 0 -- but only on a cached rerun, which is why the
## first stranded sweep looked fine.
jgt_all <- read.table(if (STRANDED && !RMATS_UNIT_NA) RMATS_STR_GT else
                      file.path(BASE, "results.beforegatefix/sim_junction_gt.txt"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
jgt_all <- jgt_all[!(jgt_all$n_inc_tx == 0 | jgt_all$n_skip_tx == 0), ]
jgt_all$uid <- if (STRANDED && !RMATS_UNIT_NA)
  paste(jgt_all$strand_run, jgt_all$event_type, jgt_all$ID, sep = ":") else
  paste(jgt_all$event_type, jgt_all$ID, sep = ":")
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

grase_files <- if (STRANDED)
  c("bipartition.internal.stranded.test.EBapprox/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
    "bipartition.TSSTTS.stranded.test.EBapprox/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt") else
  c("bipartition.test.fulldesign/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
    "bipartition.test.fulldesign/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
grase_n <- vapply(grase_files, function(f) nrow(read.table(file.path(BASE,f), header=TRUE, sep="\t",
                  quote="", comment.char="", stringsAsFactors=FALSE)), integer(1))
gr <- bind_rows(lapply(grase_files,
  function(f) read.table(file.path(BASE, f), header = TRUE, sep = "\t", quote = "",
    comment.char = "", stringsAsFactors = FALSE)[, c("gene","event","comparison","padj","padj_gene",
    "p.value","lfc_diff_net","delta_pi","setdiff1","setdiff2","transcripts1","transcripts2")]))
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
## TWO-STAGE THRESHOLD. exontest.R's nested_BH runs a gene screen at a FIXED
## alpha = 0.05 and then within-gene BH inside the surviving genes, so the
## shipped `padj` is NA for every unit of a screened-out gene. Sweeping `padj`
## alone therefore relaxes only stage 2 while the gene screen stays pinned at
## 0.05: the gene set never changes, which both flattens the curve and makes
## every threshold except 0.05 incoherent (at 0.01 we were accepting units from
## genes admitted at a LOOSER screen than the unit threshold -- anti-conservative
## by ~4%; above 0.05 no new genes could ever enter).
## nested_bh() computes padj_gene as p.adjust(simes_p, "BH") independently of
## alpha, and stage 2 as within-gene BH independently of which other genes
## passed -- so both are recoverable for EVERY gene from the raw p-values, and
## the correct threshold-t rule is padj_gene < t AND padj_within < t. Encoded as
## pmax(), since max(a,b) < t iff a < t and b < t. Verified: padj_within
## reproduces the shipped padj on all 2,989 units that have one (max abs diff
## 1.5e-15), and the two rules agree exactly at t = 0.05.
## lfc_diff_net gate. exontest.R applies `lfc_diff_net > delta` when building its
## `significant` column, with --delta defaulting to 0; the sweep must apply the
## same call condition or it would credit GrASE with calls the tool would not
## make. This was hardcoded as `> 0` in two places -- named here so it tracks
## --delta rather than silently assuming the default.
##
## The comparison is SIGNED and that is deliberate -- do NOT wrap it in abs().
## lfc_diff_net is defined as abs(lfc_diff) - abs(lfc_ref), so the magnitudes
## are already inside the quantity: it is positive when the DISTINCT set moved
## more than the shared reference did, negative when the reference moved more.
## `> delta` therefore selects genuine differential path usage and excludes
## gene-level changes that move distinct and reference together. abs() would
## readmit exactly the reference-dominated events the filter exists to drop.
DELTA <- 0

two_stage <- function(pv, gene, padj_gene) {
  within <- ave(pv, gene, FUN = function(x) p.adjust(x, method = "BH"))
  pmax(within, padj_gene)
}
gr$padj2 <- two_stage(gr$p.value, gr$gene, gr$padj_gene)
lfcok  <- !is.na(gr$lfc_diff_net) & gr$lfc_diff_net > DELTA
st <- gr$padj2
st[is.na(st)] <- 1                    # no statistic: a CALL condition, not a gap
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
sti <- gr$padj2[ki]; sti[is.na(sti)] <- 1; sti[!lfcok[ki]] <- 1
calls$GrASE_internal <- list(
  gene = gr$gene[ki], stat = sti, gene_stat = gr$padj_gene[ki],
  units = as.list(gr$uid[ki]), imp = gr$imp[ki], lower_better = TRUE)
cat(sprintf("  GrASE_internal: %d of %d tests (universe restricted)\n", sum(ki), nrow(gr)))

## Call-side effect-size variants, symmetric to MAJIQ's built-in C:
##   GrASE_dpi0.1 <-> MAJIQ_C0.10      GrASE_dpi0.2 <-> MAJIQ_C0.20
## Without these the curves compare MAJIQ WITH an effect-size threshold against
## GrASE WITHOUT one, which is not a like-for-like operating characteristic.
for (dp in c(0.1, 0.2)) {
  k <- lfcok & !is.na(gr$padj2) & !is.na(gr$delta_pi) & abs(gr$delta_pi) >= dp
  ## Keep EVERY tested unit in the call list and make filtered-out tests
  ## un-callable (padj set to 1, lower_better). Subsetting gr instead would
  ## shrink `all_units`, and recall's denominator (upos_restr) is computed from
  ## it -- so an aggressive filter would appear to RAISE recall, which is
  ## impossible for a subset of the same calls.
  st <- gr$padj2; st[is.na(st)] <- 1; st[!k] <- 1
  calls[[sprintf("GrASE_dpi%.1f", dp)]] <- list(
    gene = gr$gene, stat = st, gene_stat = gr$padj_gene,
    units = as.list(gr$uid),
    imp = gr$imp, lower_better = TRUE)
  cat(sprintf("  GrASE_dpi%.1f: %d of %d tests callable (universe unchanged)\n",
              dp, sum(k), nrow(gr)))
}

## --- GrASE_merged_all (exon+SJ merged counts, internal AND TSS/TTS) --------
## A bipartition side with no exclusive exonic part emits no exonic test. Where
## that side owns an exclusive intron, merge_exon_sj_counts.R substitutes its
## split-read count, so the side becomes testable. Same rule for both types
## (is.na(setdiff) & !is.na(intron_distinct)).
## Unit truth for this tool comes from results/gt/bipartition_merged_gt.txt
## (same GT_rule, junction-sourced sides included), selected by bp_map() above.
## Its exon-sourced rows are identical to bipartition_gt.txt, so the only
## difference from GrASE is the added junction-sourced units.
merged_files <- if (STRANDED) c(
  "bipartition.merged.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt",
  "bipartition.merged.TSSTTS.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt") else c(
  "bipartition.merged.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt",
  "bipartition.merged.TSSTTS.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt")
if (all(file.exists(file.path(BASE, merged_files)))) {
  mgn <- vapply(merged_files, function(f) nrow(read.table(file.path(BASE,f), header=TRUE, sep="\t",
                quote="", comment.char="", stringsAsFactors=FALSE)), integer(1))
  mg <- bind_rows(lapply(merged_files,
    function(f) read.table(file.path(BASE, f), header = TRUE, sep = "\t", quote = "",
      comment.char = "", stringsAsFactors = FALSE)[, c("gene","event","comparison","padj","padj_gene",
      "p.value","lfc_diff_net","delta_pi","setdiff1","setdiff2","transcripts1","transcripts2")]))
  mg$src <- rep(c("internal","TSSTTS"), times = mgn)
  mg$txs <- ifelse(grepl("diff1", mg$comparison), mg$transcripts1, mg$transcripts2)
  mg$uid <- paste(mg$src, mg$gene, mg$event, mg$comparison, sep = "|")
  mg$imp <- lapply(mg$txs, parse_l, sep = ",")
  mg$padj2 <- two_stage(mg$p.value, mg$gene, mg$padj_gene)
  mlfcok <- !is.na(mg$lfc_diff_net) & mg$lfc_diff_net > DELTA
  mst <- mg$padj2; mst[is.na(mst)] <- 1; mst[!mlfcok] <- 1
  calls$GrASE_merged_all <- list(
    gene = mg$gene, stat = mst, gene_stat = mg$padj_gene,
    units = as.list(mg$uid), imp = mg$imp, lower_better = TRUE)
  cat(sprintf("  GrASE_merged_all universe %d tests; callable: %d\n",
              nrow(mg), sum(!is.na(mg$padj) & mlfcok)))
  ## Effect-size variants on the merged unit, matching the exonic ones:
  ##   GrASE_merged_dpi0.1 <-> MAJIQ_C0.10 <-> rMATS_dpsi0.1
  ## Call-side, so the universe is UNCHANGED (see the GrASE_dpi note above --
  ## subsetting here would shrink all_units and make recall rise under a
  ## stricter filter).
  for (dp in c(0.1, 0.2)) {
    mk <- mlfcok & !is.na(mg$padj2) & !is.na(mg$delta_pi) & abs(mg$delta_pi) >= dp
    mstk <- mg$padj2; mstk[is.na(mstk)] <- 1; mstk[!mk] <- 1
    calls[[sprintf("GrASE_merged_dpi%.1f", dp)]] <- list(
      gene = mg$gene, stat = mstk, gene_stat = mg$padj_gene,
      units = as.list(mg$uid), imp = mg$imp, lower_better = TRUE)
    cat(sprintf("  GrASE_merged_dpi%.1f: %d of %d tests callable (universe unchanged)\n",
                dp, sum(mk), nrow(mg)))
  }
  ## internal-only merged bubbles: UNIVERSE restriction (as GrASE_internal is),
  ## the strict like-for-like against a junction tool.
  mki <- mg$src == "internal"
  msti <- mg$padj2[mki]; msti[is.na(msti)] <- 1; msti[!mlfcok[mki]] <- 1
  calls$GrASE_merged_internal <- list(
    gene = mg$gene[mki], stat = msti, gene_stat = mg$padj_gene[mki],
    units = as.list(mg$uid[mki]), imp = mg$imp[mki], lower_better = TRUE)
  cat(sprintf("  GrASE_merged_internal: %d of %d tests (universe restricted)\n",
              sum(mki), nrow(mg)))
  ## Same merged unit and the same raw p-values, but PLAIN BH instead of
  ## exontest.R's nested_BH gene screen -- matching DEXSeq's FDR architecture.
  ## Under nested_BH only 1,345 null-gene units are reachable at any threshold
  ## (99.4% get padj = NA at stage 1); under plain BH 209,981 are. The two agree
  ## closely at padj 0.01 (132 vs 127 null calls) and diverge as it loosens, so
  ## the curve shapes are partly an FDR-architecture artifact rather than a
  ## property of the test. Carrying both makes that explicit.
  mbh <- p.adjust(mg$p.value, method = "BH")
  mbh[!mlfcok] <- 1
  calls$GrASE_BH <- list(
    gene = mg$gene, stat = mbh, gene_stat = mg$padj_gene,
    units = as.list(mg$uid), imp = mg$imp, lower_better = TRUE)
  cat(sprintf("  GrASE_BH (plain BH): %d of %d tests callable\n",
              sum(!is.na(mbh) & mbh < 0.01), nrow(mg)))
} else {
  cat("  GrASE_merged_all: merged test files absent, skipping\n")
}

## --- DEXSeq (bins) ---------------------------------------------------------
dex <- read.table(file.path(BASE, if (STRANDED)
    "DEXSeq/dexseq_group1_group2_stranded/all.group1_group2.dxd_filteredbyCountMultiExon.txt" else
    "DEXSeq/dexseq_group1_group2/all.group1_group2.dxd_filteredbyCountMultiExon.txt"),
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
jk <- read.table(file.path(BASE, if (STRANDED) "results/gt/lsv_junction_keys.stranded.txt"
                                 else "results/gt/lsv_junction_keys.txt"), header = TRUE,
                 sep = "\t", quote = "", stringsAsFactors = FALSE)
stopifnot("users" %in% names(jk))
jk$uid <- paste(jk$lsv_id, jk$j_index, sep = "|")
jk$imp <- lapply(jk$users, function(x)
  if (is.na(x) || !nzchar(x)) character(0) else strsplit(x, ",", fixed = TRUE)[[1]])
cat(sprintf("  junction keys: %d junctions, %d with >=1 user\n",
            nrow(jk), sum(lengths(jk$imp) > 0)))

## per-junction probabilities at the other C setting (cheap: no transcript work)
prob_at <- function(Cv) {
  tv <- read.table(file.path(BASE, sprintf(if (STRANDED)
                     "majiq/majiq_deltapsi.stranded.thr%s.tsv" else
                     "majiq/majiq_deltapsi.thr%s.tsv", Cv)),
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

## --- MAJIQ (LSVs; per-LSV max prob; LSV is the unit so it implicates the WHOLE
## LSV -- the union of ALL its junctions' users, regardless of per-junction
## significance -- analogous to rMATS implicating inc+skip. (PSI sums to 1 across
## the LSV, so the non-max junctions are the other side of the same switch.) The
## per-junction fine view is MAJIQjunc below.) --------------------------------
exons_gr <- import(file.path(BASE, "ref/gencode.v28.annotation.gtf"),
                   feature.type = "exon", colnames = c("transcript_id", "gene_id"))
{ o <- order(exons_gr$transcript_id, start(exons_gr))
  tid <- exons_gr$transcript_id[o]; s_ <- start(exons_gr)[o]; e_ <- end(exons_gr)[o]
  n <- length(tid); same <- tid[-n] == tid[-1]
  tx_junc <- split(paste(e_[-n], s_[-1], sep = "-")[same], tid[-n][same]) }
# whole-LSV implicated set: union of every junction's users, keyed by lsv_id
imp_by_lsv <- tapply(seq_len(nrow(jk)), jk$lsv_id,
                     function(ix) unique(unlist(jk$imp[ix])))
for (Cv in c("0.20", "0.10")) {
  pv <- PROB[[Cv]][jk$uid]; pv[is.na(pv)] <- 0
  d  <- data.frame(lsv_id = jk$lsv_id, gene = jk$gene, p = as.numeric(pv),
                   idx = seq_len(nrow(jk)), stringsAsFactors = FALSE)
  # per LSV: the max probability (the LSV's significance statistic)
  o    <- order(d$lsv_id, -d$p)
  dd   <- d[o, ]
  keep <- !duplicated(dd$lsv_id)
  top  <- dd[keep, ]
  calls[[paste0("MAJIQ_C", Cv)]] <- list(
    gene = top$gene, stat = top$p, gene_stat = top$p,
    units = as.list(top$lsv_id), imp = unname(imp_by_lsv[top$lsv_id]),
    lower_better = FALSE)
  cat(sprintf("  MAJIQ_C%s: %d LSVs (whole-LSV implication)\n", Cv, nrow(top)))
}


## --- rMATS (events; inc U skip transcripts) --------------------------------
csv_mean <- function(x) { v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
                          if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE) }
jgt <- read.table(if (STRANDED && !RMATS_UNIT_NA) RMATS_STR_GT else
                  file.path(BASE, "results.beforegatefix/sim_junction_gt.txt"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
jgt <- jgt[!(jgt$n_inc_tx == 0 | jgt$n_skip_tx == 0), ]
## uid must match how the rMATS block keys its events below: strand-prefixed
## under STRANDED, because rMATS renumbers from 0 in each half-run.
jgt$uid <- if (STRANDED && !RMATS_UNIT_NA)
  paste(jgt$strand_run, jgt$event_type, jgt$ID, sep = ":") else
  paste(jgt$event_type, jgt$ID, sep = ":")
## NOT reassigning rm_pos here: it is set once above, outside the cache gate, so
## that a cached rerun cannot end up with a stale key convention.
rg <- list(); rstat <- c(); rgene <- c(); rimp <- list(); ruid <- c()
rstat_nogate <- c()   # same events, FDR unmodified (PSI gate OFF)
rdpsi <- c()          # rMATS' own |dPSI| (IncLevelDifference), for effect-size variants
## per-transcript exon spans: an RI event's INCLUSION isoform is defined by the
## ABSENCE of splicing, so it has no junction to match. Junction-only
## implication therefore misses every intron-retaining transcript. Mirrors
## transcript_level_metrics.R's add_ev(span=) so the two scripts agree.
tx_ex <- split(data.frame(s = start(exons_gr), e = end(exons_gr)), exons_gr$transcript_id)
rm_dirs <- if (STRANDED) c(plus = "rMATS/stranded_plus/post", minus = "rMATS/stranded_minus/post")
           else c(all = "rMATS/rmats_post_group1_group2")
for (rmk in names(rm_dirs)) {
for (et in c("SE","A3SS","A5SS","RI")) {
  rmf <- file.path(BASE, rm_dirs[[rmk]], paste0(et, ".MATS.JCEC.txt"))
  if (!file.exists(rmf)) next
  d <- read.table(rmf, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  names(d) <- make.unique(names(d)); d$GeneID <- gsub('"', "", d$GeneID)
  c1 <- sapply(d$IJC_SAMPLE_1, csv_mean) + sapply(d$SJC_SAMPLE_1, csv_mean)
  c2 <- sapply(d$IJC_SAMPLE_2, csv_mean) + sapply(d$SJC_SAMPLE_2, csv_mean)
  p1 <- sapply(d$IncLevel1, csv_mean); p2 <- sapply(d$IncLevel2, csv_mean)
  ## rMATS PSI/count filter is a CALL-SIDE condition (decided project
  ## convention; see ver12.md and the Table 2 eval, which force FDR to 1 and
  ## KEEP the event). Do NOT subset the universe here: dropping the filtered
  ## events removes GT positives from recall's denominator and inflates rMATS
  ## recall (0.425 -> 0.690 at padj 0.01). This mirrors how the GrASE block
  ## above treats screened-out / lfc-failed tests -- kept, made un-callable.
  testable <- !is.na(c1) & c1 >= 10 & !is.na(c2) & c2 >= 10 &
              !is.na(p1) & p1 >= 0.05 & p1 <= 0.95 & !is.na(p2) & p2 >= 0.05 & p2 <= 0.95
  fdr_call <- d$FDR; fdr_call[!testable] <- 1
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
    imp_i <- if (length(gtx)) gtx[vapply(tx_junc[gtx], function(v) any(js %in% v), logical(1))] else character(0)
    if (et == "RI") {   # add transcripts whose exon spans the retained intron
      gtx2 <- intersect(gene_tx[[g]], names(tx_ex))
      sp <- gtx2[vapply(tx_ex[gtx2], function(z)
              any(z$s <= d$upstreamEE[i] & z$e >= d$downstreamES[i] + 1), logical(1))]
      imp_i <- union(imp_i, sp)
    }
    rimp[[length(rimp)+1]] <- imp_i
    ## strand-prefixed so the two half-runs cannot collide (rMATS numbers from 0)
    rgene <- c(rgene, g); rstat <- c(rstat, fdr_call[i])
    ruid <- c(ruid, if (STRANDED) paste(rmk, et, d$ID[i], sep=":") else paste(et, d$ID[i], sep=":"))
    rstat_nogate <- c(rstat_nogate, d$FDR[i])
    rdpsi <- c(rdpsi, suppressWarnings(as.numeric(d$IncLevelDifference[i])))
  }
}
}
calls$rMATS <- list(gene = rgene, stat = rstat, gene_stat = rstat,
                    units = as.list(ruid), imp = rimp, lower_better = TRUE)
## Call-side effect-size variants using rMATS' OWN dPSI estimate
## (IncLevelDifference), the true analogue of GrASE_dpi / MAJIQ_C:
##   rMATS_dpsi0.1 <-> GrASE_dpi0.1 <-> MAJIQ_C0.10
##   rMATS_dpsi0.2 <-> GrASE_dpi0.2 <-> MAJIQ_C0.20
## Built on the GATED base (the documented configuration); like the GrASE dpi
## variants these do NOT subset -- below-threshold events stay in the universe
## and are made un-callable.
rdpsi[is.na(rdpsi)] <- 0
for (dp in c(0.1, 0.2)) {
  st_dp <- rstat; st_dp[abs(rdpsi) < dp] <- 1
  calls[[sprintf("rMATS_dpsi%.1f", dp)]] <- list(
    gene = rgene, stat = st_dp, gene_stat = st_dp,
    units = as.list(ruid), imp = rimp, lower_better = TRUE)
}
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
              GrASE_merged_all = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              "GrASE_merged_dpi0.1" = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              "GrASE_merged_dpi0.2" = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_merged_internal = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_BH = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_dpi0.1 = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_dpi0.2 = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              GrASE_internal = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              DEXSeq = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              MAJIQ_C0.20 = c(0.99,0.95,0.9,0.8,0.7,0.5),
              MAJIQ_C0.10 = c(0.99,0.95,0.9,0.8,0.7,0.5),
              rMATS = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              "rMATS_dpsi0.1" = c(1e-4,1e-3,0.01,0.05,0.1,0.2),
              "rMATS_dpsi0.2" = c(1e-4,1e-3,0.01,0.05,0.1,0.2))
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
  ## full-universe negative populations (per panel): every enumerable gene / tx /
  ## unit in the panel's scope. DTE/DTU panels include Background as the negatives.
  panel_genes <- switch(categ, ALL = names(gene_tx),
                        DTE = union(dte.genes, bg_genes), DTU = union(dtu.genes, bg_genes),
                        Background = bg_genes, DGE = dge.genes)
  panel_st    <- switch(categ, ALL = c("DTE","DTU","Background","DGE"),
                        DTE = c("DTE","Background"), DTU = c("DTU","Background"),
                        Background = "Background", DGE = "DGE")
  n_genes_full_c <- length(panel_genes)
  n_tx_full_c    <- length(unique(unlist(gene_tx[panel_genes])))
  n_parts_full_c <- sum(ep$sim_type %in% panel_st & (ep$gt_pos | ep$gt_neg))  # DEXSeq bins
  n_evt_full_c   <- sum(jgt_all$sim_type %in% panel_st)                        # rMATS events
  all_units <- unique(unlist(cl$units))
  # restricted denominators (positives among units this tool tested in this category)
  upos_restr <- if (grepl("^MAJIQjunc", tool)) sum(jn_pos[intersect(all_units, names(jn_pos))], na.rm = TRUE)
                else if (grepl("^GrASE", tool)) { bpm_ <- bp_map(tool)
                     sum(bpm_[intersect(all_units, names(bpm_))], na.rm = TRUE) }
                else if (grepl("^MAJIQ", tool)) sum(lsv_pos_rule[intersect(all_units, names(lsv_pos_rule))], na.rm = TRUE)
                else if (grepl("^rMATS", tool))  sum(rm_pos[intersect(all_units, names(rm_pos))], na.rm = TRUE)
                else sum(vapply(all_units, function(k) exists(k, envir = pos_key), logical(1)))
  pick <- function(tabl, full_all) if (categ == "ALL") full_all else
            as.numeric(tabl[categ]) %||% 0
  `%||%` <- function(a, b) if (length(a) == 0 || is.na(a)) b else a
  ## GT_rule cannot label units the tool never quantified, so for GrASE and
  ## MAJIQ the full-universe unit denominator is set to the restricted one.
  upos_full  <- if (grepl("^GrASE", tool)) upos_restr
                else if (grepl("^MAJIQ", tool)) upos_restr
                else if (grepl("^rMATS", tool))  (if (categ=="ALL") n_pos_evt_full   else as.numeric(n_pos_evt_cat[categ]))
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
  ## whether each TESTED unit's implicated set hits >= 1 manip tx -- independent of
  ## the threshold, so compute ONCE per category (hoisted out of the thr loop;
  ## recomputing this vapply over GrASE's ~339k units per threshold was the cost).
  hit_unit <- vapply(cl$imp, function(v) any(v %in% cat_tx), logical(1))
  for (thr in GRIDS[[tool]]) {
    sig <- if (cl$lower_better) cl$stat < thr else cl$stat >= thr
    if (!any(sig)) next
    ## --- UNIT ---
    su <- unique(unlist(cl$units[sig]))
    if (grepl("^MAJIQjunc", tool)) {
      v <- jn_pos[su]; v[is.na(v)] <- FALSE
      tp <- sum(v); fp <- sum(!v)
    } else if (grepl("^GrASE", tool)) {
      v <- bp_map(tool)[su]; v[is.na(v)] <- FALSE  # null-gene tests are negatives
      tp <- sum(v); fp <- sum(!v)
    } else if (grepl("^MAJIQ", tool)) {
      v <- lsv_pos_rule[su]; tp <- sum(v, na.rm = TRUE); fp <- sum(!v, na.rm = TRUE)
    } else if (grepl("^rMATS", tool)) {
      if (RMATS_UNIT_NA) { tp <- NA_integer_; fp <- NA_integer_ } else {
        v <- rm_pos[su];  tp <- sum(v, na.rm = TRUE); fp <- sum(!v, na.rm = TRUE) }
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
    hit_c     <- hit_unit[sig]
    impl_hit  <- unique(unlist(cl$imp[sig][ hit_c]))
    impl_miss <- unique(unlist(cl$imp[sig][!hit_c]))
    tol_fp    <- if (is_null_cat) length(stx)
                 else length(setdiff(setdiff(impl_miss, impl_hit), manip_tx))
    ## call currency (per-call/per-unit): a call is TP if its implicated set hits
    ## >= 1 manipulated tx, FP if it hits none; recall over tested units bearing
    ## signal. A wrong call is 1 FP no matter how many transcripts it drags --
    ## contrast transcript_tolerant, which charges those dragged transcripts. This
    ## is the most coarse-unit-tolerant view (see doc: rMATS call vs transcript).
    ctp <- sum( sig & hit_unit); cfp <- sum( sig & !hit_unit)
    cfn <- sum(!sig & hit_unit); ctn <- sum(!sig & !hit_unit)
    for (uni in c("restricted", "full")) {
      ## TN / n_neg are counted over the enumerable NEGATIVE population of the
      ## panel: restricted = the units/tx/genes the tool tested; full = every
      ## enumerable one in the panel's scope (all genes, all their transcripts, all
      ## enumerated units). Unit universe by tool: GrASE = all enumerated
      ## bipartitions (already held in cl$units), DEXSeq = all bins, rMATS = all
      ## events, MAJIQ = quantified LSVs/junctions (no all-possible catalog, so its
      ## full == restricted).
      isr <- uni == "restricted"
      du <- if (isr) upos_restr else upos_full
      dg <- if (isr) gpos_restr else n_pos_gene_c
      dt <- if (isr) tpos_restr else n_pos_tx_c
      N_gene <- if (isr) n_genes_tested0 else n_genes_full_c
      N_tx   <- if (isr) n_tx_tested0    else n_tx_full_c
      N_unit <- if (isr) n_units_tested
                else if (grepl("^GrASE", tool)) n_units_tested   # cl$units holds all bipartitions
                else if (tool == "DEXSeq") n_parts_full_c
                else if (grepl("^rMATS", tool))  n_evt_full_c
                else n_units_tested                              # MAJIQ: quantified universe
      rows[[length(rows)+1]] <- data.frame(category=categ, level="unit", universe=uni, tool=tool, thr=thr,
        TP=tp, FP=fp, FN=max(du-tp,0), TN=max(N_unit-tp-fp-max(du-tp,0),0),
        n_pos=du, n_neg=max(N_unit-du,0),
        precision=tp/max(tp+fp,1), recall=tp/max(du,1))
      rows[[length(rows)+1]] <- data.frame(category=categ, level="transcript_strict", universe=uni, tool=tool, thr=thr,
        TP=ttp, FP=tfp, FN=max(dt-ttp,0), TN=max(N_tx-ttp-tfp-max(dt-ttp,0),0),
        n_pos=dt, n_neg=max(N_tx-dt,0),
        precision=ttp/max(ttp+tfp,1), recall=ttp/max(dt,1))
      rows[[length(rows)+1]] <- data.frame(category=categ, level="gene", universe=uni, tool=tool, thr=thr,
        TP=gtp, FP=gfp, FN=max(dg-gtp,0), TN=max(N_gene-gtp-gfp-max(dg-gtp,0),0),
        n_pos=dg, n_neg=max(N_gene-dg,0),
        precision=gtp/max(gtp+gfp,1), recall=gtp/max(dg,1))
      # transcript_tolerant: same transcript TP / recall / n_pos as strict; only FP
      # differs (tol_fp <= tfp -- co-travellers of a hit forgiven). Transcript
      # currency, so TN closes in both universes just like strict.
      rows[[length(rows)+1]] <- data.frame(category=categ, level="transcript_tolerant", universe=uni, tool=tool, thr=thr,
        TP=ttp, FP=tol_fp, FN=max(dt-ttp,0), TN=max(N_tx-ttp-tol_fp-max(dt-ttp,0),0),
        n_pos=dt, n_neg=max(N_tx-dt,0),
        precision=ttp/max(ttp+tol_fp,1), recall=ttp/max(dt,1))
      # call currency: TP/FP/FN/TN counted in CALLS (one per tested unit). recall
      # and the closed table are restricted-only (untested units have no call).
      rows[[length(rows)+1]] <- data.frame(category=categ, level="call", universe=uni, tool=tool, thr=thr,
        TP=ctp, FP=cfp,
        FN=if (isr) cfn else NA_integer_, TN=if (isr) ctn else NA_integer_,
        n_pos=if (isr) ctp+cfn else NA_integer_, n_neg=if (isr) cfp+ctn else NA_integer_,
        precision=ctp/max(ctp+cfp,1), recall=if (isr) ctp/max(ctp+cfn,1) else NA_real_)
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

## --- plot styling + tool selection ------------------------------------------
## Effect size is shown by SHADE (darker = stronger threshold); all lines solid:
##   GrASE (none) -> GrASE_dpi0.1 -> GrASE_dpi0.2   light -> dark orange
##   MAJIQ_C0.10  -> MAJIQ_C0.20                     light -> dark teal
## Distinct HUE per tool; effect size = SHADE within the hue (darker = stronger):
##   GrASE = red, DEXSeq = orange, MAJIQ = green, rMATS = navy/purple.
## A merged variant keeps its exonic twin's colour, so the same shade means the
## same operating point in either family. The PLOTS show only the MERGED GrASE
## family: the merged unit is a strict superset of the exonic one (a junction
## count is substituted only where the exonic pipeline had nothing to count), so
## the exonic curves are a subset view and drawing both doubled the red curves
## for no added information. Both families stay in the TABLE.
## GrASE_BH first so it heads the legend, then the nested_BH GrASE family.
## DEXSeq moved to blue to free the orange range for GrASE_BH.
COL <- c("GrASE_BH"="#fec44f",
         GrASE="#fb6a4a", "GrASE_merged_all"="#fb6a4a",
         "GrASE_dpi0.1"="#de2d26", "GrASE_merged_dpi0.1"="#de2d26",
         "GrASE_dpi0.2"="#a50f15", "GrASE_merged_dpi0.2"="#a50f15",
         "GrASE_internal"="#67000d", "GrASE_merged_internal"="#67000d",
         DEXSeq="#3182bd",
         MAJIQ_C0.10="#a1d99b", MAJIQ_C0.20="#238b45",
         ## rMATS purples were #3f007d / #6a51a3 / #54278f -- three near-identical
         ## dark purples, and the ladder ran the WRONG WAY (the no-threshold
         ## variant was the darkest). Widened to a proper light->dark ramp so
         ## the shade convention holds here as it does for GrASE and MAJIQ.
         rMATS="#bcbddc",
         "rMATS_dpsi0.1"="#756bb1", "rMATS_dpsi0.2"="#3f007d")
LTY <- setNames(rep(1L, length(COL)), names(COL))   # all solid; effect size = shade
## EDIT this to choose which curves are drawn (set to names(COL) for all).
## Dropped by default: every EXONIC GrASE variant (superseded by its merged
## twin) and GrASE_merged_internal (a universe-restricted view, not comparable
## to the pooled curves).
PLOT_TOOLS <- setdiff(names(COL),
                      c("GrASE", "GrASE_dpi0.1", "GrASE_dpi0.2",
                        "GrASE_internal", "GrASE_merged_internal"))
## Display labels only. The merged variants are the ONLY GrASE family drawn, so
## the "_merged" tag carries no information in the figure and just lengthens
## every label. Names in the written table are untouched -- figure and table are
## joined through the tool column, not the label.
##   GrASE_merged_all -> GrASE   GrASE_merged_dpi0.1 -> GrASE_dpi0.1   etc.
disp <- function(x) {
  ## ungated is GrASE_dpi0, not "GrASE": exontest.R ships --min_dpi 0.1, so its
  ## `significant` column is padj<0.01 AND |delta_pi|>=0.1. Labelling the ungated
  ## run "GrASE" would present something 65% larger than the shipped default
  ## (1,012 vs 615 calls on the merged stranded internal file).
  x <- sub("^GrASE_BH$", "GrASE (plain BH)", x)
  sub("^GrASE_merged_", "GrASE_", sub("^GrASE_merged_all$", "GrASE_dpi0", x))
}

for (categ in c("ALL", "DTE", "DTU", "Background", "DGE")) {
if (!(categ %in% c("Background","DGE"))) {   # PR curves need positives
fn <- if (categ == "ALL") "plots/pr_three_levels.GT_rule.png" else sprintf("plots/pr_three_levels.GT_rule_%s.png", categ)
png(file.path(BASE, fn), width = 1900, height = 1000, res = 130)
par(mfrow = c(2,4), mar = c(4.2,4.2,3,1), oma = c(3.1,0,1.6,0))
# COL, LTY, PLOT_TOOLS defined once above (effect size = shade; all solid).
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
    d <- d[d$tool %in% PLOT_TOOLS, ]                     # user-selected tools only
    ## Curves start at the standard operating point and loosen from there.
    ## padj < 0.01 is where these tools are reported; 1e-3 and 1e-4 are not
    ## operating points anyone uses, and plotting them is misleading -- on the
    ## null-space FPR view DEXSeq looked dominant purely on a padj 1e-4 point
    ## (22 null FP) while being 2-8x worse than GrASE at every usable threshold.
    ## The rows stay in the written table; only the curves are trimmed.
    ## MAJIQ's statistic runs the other way (probability floor, larger =
    ## stricter), so its standard point is 0.95 and it loosens downward.
    ## MAJIQ keeps 0.99: its statistic is a posterior probability, not an FDR,
    ## and 0.95 vs 0.99 is a routine user choice in the MAJIQ literature --
    ## unlike padj 1e-3/1e-4, which nobody reports. The extra point also lets
    ## MAJIQ_C0.20 reach FP 24 in gene space, which finally overlaps
    ## GrASE_dpi0.1's range (21-22) and makes a matched comparison possible.
    d <- d[!grepl("^MAJIQ", d$tool) & d$thr >= 0.01 | grepl("^MAJIQ", d$tool), ]
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
      gk <- intersect(grep("^GrASE", names(COL), value = TRUE), PLOT_TOOLS)
      legend("bottomleft", legend=disp(gk), col=COL[gk], lty=LTY[gk], lwd=2, bty="n", cex=0.7)
    }
    if (lv == "transcript_strict" && uni == "restricted") {  # cross-tool legend on a fair panel
      lk <- intersect(names(COL), PLOT_TOOLS)
      legend("bottomleft", legend=disp(lk), col=COL[lk], lty=LTY[lk], lwd=2, bty="n", cex=0.7)
    }
  }
}
mtext(sprintf("Simulation category: %s", categ), outer = TRUE, cex = 0.9, font = 2)
mtext(paste0("unit = tested unit (within-GrASE only) | transcript_strict = charges every co-traveller | ",
             "transcript_tolerant = forgives co-travellers riding a hit (same recall) | ",
             "DTE/DTU panels use Background as the negative class; ALL is pooled"),
      side = 1, outer = TRUE, cex = 0.6, line = 0.15)
## Row identity matters for how these panels may be read. The restricted row
## gives each tool its OWN negative class (only what its event definitions can
## reach: 113,616 transcripts for GrASE vs 85,822 for rMATS), so those
## denominators are tool-specific and the row is a WITHIN-TOOL view. The full
## row shares one negative class (198,524) across every tool, and is the
## cross-tool comparison.
mtext(paste0("TOP ROW restricted universe: each tool's own reachable negatives ",
             "(tool-specific denominator -- WITHIN-TOOL only)   |   ",
             "BOTTOM ROW full universe: one shared negative class for all tools ",
             "-- this is the cross-tool comparison"),
      side = 1, outer = TRUE, cex = 0.6, line = 0.85, font = 2)
if (categ %in% c("DTE","DTU"))
  mtext(paste0("positives = ", categ, " genes; negatives = Background genes ",
               "(a Background call is a gene FP), so gene precision is informative here"),
        side = 1, outer = TRUE, cex = 0.6, line = 1.55, font = 3)
dev.off()
cat(sprintf("  wrote %s\n", fn))
}

## --- stacked TP / FN / FP bars at the operating points --------------------
fn2 <- sprintf("plots/confusion_three_levels.GT_rule%s.png",
               if (categ == "ALL") "" else paste0("_", categ))
png(file.path(BASE, fn2), width = 1900, height = 1000, res = 130)
par(mfrow = c(2,4), mar = c(4,8.5,3,2.5), oma = c(3.1,0,1.6,0))
op <- tab[tab$category == categ &
          ((grepl("^GrASE|^DEXSeq|^rMATS", tab$tool) & tab$thr == 0.01) |
           (grepl("MAJIQ", tab$tool) & tab$thr == 0.95)), ]
BARCOL <- c(TP = "#1d9e75", FN = "#888780", FP = "#e24b4a")
for (uni in c("restricted","full")) {
  for (lv in LEVELS) {
    d <- op[op$level == lv & op$universe == uni, ]
    ord <- intersect(names(COL), PLOT_TOOLS)                     # user-selected tools only
    if (lv == "unit") ord <- grep("^GrASE", ord, value = TRUE)  # unit: within-GrASE only
    d <- d[match(ord, d$tool), ]; d <- d[!is.na(d$tool), ]
    m <- rbind(TP = d$TP, FN = pmax(d$n_pos - d$TP, 0), FP = d$FP)
    colnames(m) <- disp(d$tool)
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
mtext(paste0("TOP ROW restricted universe (tool-specific negatives, WITHIN-TOOL only)",
             "   |   BOTTOM ROW full universe (shared negatives, cross-tool)"),
      side = 1, outer = TRUE, cex = 0.6, line = 0.85, font = 2)
if (categ %in% c("DTE","DTU"))
  mtext(paste0("positives = ", categ, " genes; negatives = Background genes ",
               "(a Background call is a gene FP), so gene precision is informative here"),
        side = 1, outer = TRUE, cex = 0.6, line = 1.55, font = 3)
dev.off()
cat(sprintf("  wrote %s\n", fn2))
}
cat(sprintf("\nTable: %s\n", file.path(OUT, "pr_three_levels.GT_rule.txt")))
