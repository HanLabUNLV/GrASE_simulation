#!/usr/bin/env Rscript
#
# scripts/infer_junctions_gt.R
#
# Derive junction-level ground truth for rMATS evaluation, using the
# threshold-free STRUCTURAL ground truth from the simulation.
#
# Uses fromGTF.*.txt files as the full event universe (all annotated events,
# regardless of read support). Each event is assigned GT status from which
# transcripts the simulation actually manipulated (iso.dtu / iso.dte), matched
# to the inclusion vs exclusion path of the event via GTF coordinate matching.
# Events also present in *.MATS.JCEC.txt are flagged rmats_tested=TRUE.
#
# GT logic (structural, no PSI threshold):
#   DTU gene: GT-positive iff the manipulated transcripts are NOT all on the
#             same side, where each manip tx is assigned to exactly one of
#             inc / skip / neither. (Spanning >1 side -> PSI shifts.)
#   DTE gene: GT-positive iff the single manipulated tx is on either side.
#   DGE/Background: always GT-negative.
#
# PSI / dPSI are retained as DIAGNOSTIC columns only (used by check_rmats_fps.R)
# and no longer drive gt_positive.
#
# This supports two evaluation universes:
#   full       : all annotated events (untested GT-positives count as FN)
#   restricted : only rmats_tested events (fair comparison for tested events)
#
# Output: one file per gene <out_dir>/<GeneID>.txt (resumable, parallel-safe).
#   columns: event_type, ID, GeneID, chr, strand, coord,
#            n_inc_tx, n_skip_tx, n_manip_inc, n_manip_skip,
#            PSI1, PSI2, dPSI (diagnostic), sim_type,
#            gt_positive, rmats_tested
#   Combine the per-gene files into one table with shell (head + tail);
#   see run_infer_junctions_gt.sh. R-side combining is too slow.
#
# Resume: a present <GeneID>.txt means that gene is done and is skipped.
#         Writes are atomic (.tmp then rename) so partial files never appear.
#
# Usage:
#   Rscript infer_junctions_gt.R <rmats_dir> <gtf> <simulate_rda> <out_dir>
#                                [shard_idx n_shards]
#   With shard_idx (1-based) and n_shards, only genes in that shard are
#   processed -- launch n_shards jobs in parallel, then combine in shell.

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  rmats_dir   <- "/mnt/data1/home/mirahan/GrASE_simulation/rMATS/rmats_post_group1_group2"
  gtf_file    <- "/mnt/data1/home/mirahan/GrASE_simulation/ref/gencode.v28.annotation.gtf"
  sim_rda     <- "/mnt/data1/home/mirahan/GrASE_simulation/swimdown/simulate/data/simulate.rda"
  out_dir     <- "/mnt/data1/home/mirahan/GrASE_simulation/results/sim_junction_gt"
  shard_idx   <- NA_integer_
  n_shards    <- NA_integer_
} else {
  rmats_dir  <- args[1]
  gtf_file   <- args[2]
  sim_rda    <- args[3]
  out_dir    <- args[4]  # per-gene output directory
  # optional sharding for parallel runs: process only genes in this shard
  shard_idx  <- if (length(args) >= 6) as.integer(args[5]) else NA_integer_
  n_shards   <- if (length(args) >= 6) as.integer(args[6]) else NA_integer_
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sharded <- !is.na(shard_idx) && !is.na(n_shards)
if (sharded) {
  cat(sprintf("Shard %d of %d\n", shard_idx, n_shards))
}

# --- load simulation data ---------------------------------------------------

cat("Loading simulate.rda...\n")
load(sim_rda)
# sim.counts.mat: transcripts x 2 conditions (expected read counts)
# iso.dtu:        named logical; TRUE for manipulated DTU transcripts (2 per gene)
# iso.dte:        named logical; TRUE for manipulated DTE transcripts (1 per gene)
# txdf:           GENEID, TXNAME, ntx
# dge.genes / dte.genes / dtu.genes

get_sim_type <- function(gene) {
  if (gene %in% dte.genes) return("DTE")
  if (gene %in% dtu.genes) return("DTU")
  if (gene %in% dge.genes) return("DGE")
  return("Background")
}

counts1 <- sim.counts.mat[, 1]
counts2 <- sim.counts.mat[, 2]

cat(sprintf("  Simulation universe: %d transcripts\n", length(counts1)))
cat(sprintf("  Manipulated DTU transcripts: %d\n", sum(iso.dtu)))
cat(sprintf("  Manipulated DTE transcripts: %d\n", sum(iso.dte)))

# --- structural GT helpers --------------------------------------------------

# manipulated transcripts of a gene, from a named logical iso vector
get_manip_tx <- function(gene, iso_vec) {
  all_tx <- txdf$TXNAME[txdf$GENEID == gene]
  hit <- intersect(all_tx, names(iso_vec))
  if (length(hit) == 0) return(character(0))
  names(which(iso_vec[hit]))
}

# structural GT for one event; returns gt flag + manip-side diagnostics.
# called after inc_tx / skip_tx are determined for the event.
#
# Unifying principle for both rules: GT-positive iff the event REFLECTS THE
# DESIGNED MANIPULATION (changes the inc:skip ratio by design, not incidentally).
# The two rules differ in form only because the manipulations differ:
#
#   DTE: one transcript's abundance changes -> reflected iff that transcript is
#        on a side of the event => "manip participates (inc or skip)".
#
#   DTU: a pair of transcripts swaps in opposite directions -> the swap changes
#        the ratio by design only if the event DISTINGUISHES the two isoforms
#        (puts them in different positions) => "manip span > 1 of inc/skip/neither".
#        If both manip sit on the SAME side, the event cannot tell them apart, so
#        the swap is resolved at some other boundary, not here -> negative.
#
# Note the deliberate asymmetry: "span > 1" is STRICTER than "any PSI change".
# Under an ideal conserving swap a same-side event has exactly zero dPSI; with
# the simulation's reciprocal-but-different-baseline folds it has a small,
# incidental dPSI, which we treat as an artifact (negative), not signal. The
# looser parallel-to-DTE alternative would be "at least one manip participates"
# (NOT both-neither), which would flip both-inc / both-skip events to positive.
# We use separation (the swap-faithful definition) on purpose.
gt_from_structure <- function(gene, inc_tx, skip_tx) {
  stype <- get_sim_type(gene)
  if (stype == "DTU") {
    manip <- get_manip_tx(gene, iso.dtu)
  } else if (stype == "DTE") {
    manip <- get_manip_tx(gene, iso.dte)
  } else {
    manip <- character(0)
  }

  n_manip_inc  <- sum(manip %in% inc_tx)
  n_manip_skip <- sum(manip %in% skip_tx)

  if (stype == "DTU") {
    # GT-positive iff manip transcripts span more than one side.
    # neither = not in inc_tx nor skip_tx; counts as its own side.
    sides <- ifelse(manip %in% inc_tx, "inc",
             ifelse(manip %in% skip_tx, "skip", "neither"))
    gt <- length(unique(sides)) > 1
  } else if (stype == "DTE") {
    gt <- any(manip %in% inc_tx) || any(manip %in% skip_tx)
  } else {
    gt <- FALSE
  }

  # NOTE: gt_positive is the pure structural truth (did a manipulated
  # transcript participate in this event). Whether the event can actually be
  # SCORED is a separate question -- handled downstream by excluding events
  # with an empty inc or skip side (n_inc_tx == 0 | n_skip_tx == 0), since PSI
  # is constant/undefined there. Those are untestable (often our exact-exon
  # matching just missed an isoform rMATS sees in reads), not true negatives.

  list(gt = gt, n_manip_inc = n_manip_inc, n_manip_skip = n_manip_skip)
}

# --- load GTF exon ranges ---------------------------------------------------

cat("Loading GTF exon features (this may take a minute)...\n")
exons_gr <- import(gtf_file, feature.type = "exon",
                   colnames = c("transcript_id", "gene_id"))
cat(sprintf("  Loaded %d exon features\n", length(exons_gr)))

# build per-transcript exon list (sorted by start)
tx_exons <- split(exons_gr, exons_gr$transcript_id)
cat(sprintf("  %d transcripts with exon annotations\n", length(tx_exons)))

# helper: for a transcript, does it contain an exon matching [s, e] exactly?
has_exon <- function(tx_gr, s, e) {
  any(start(tx_gr) == s & end(tx_gr) == e)
}

# compute PSI for inclusion vs skip transcript sets (diagnostic only)
psi <- function(inc_tx, skip_tx, strand, counts1, counts2) {
  inc_tx <- intersect(inc_tx, names(counts1))
  skip_tx <- intersect(skip_tx, names(counts1))
  if (length(inc_tx) == 0 && length(skip_tx) == 0)
    return(c(NA_real_, NA_real_))
  c1_inc  <- sum(counts1[inc_tx],  na.rm = TRUE)
  c2_inc  <- sum(counts2[inc_tx],  na.rm = TRUE)
  c1_skip <- sum(counts1[skip_tx], na.rm = TRUE)
  c2_skip <- sum(counts2[skip_tx], na.rm = TRUE)
  psi1 <- if ((c1_inc + c1_skip) > 0) c1_inc / (c1_inc + c1_skip) else NA_real_
  psi2 <- if ((c2_inc + c2_skip) > 0) c2_inc / (c2_inc + c2_skip) else NA_real_
  c(psi1, psi2)
}

# --- process SE events ------------------------------------------------------

process_SE <- function(events_df, tx_exons, counts1, counts2) {
  rows <- lapply(seq_len(nrow(events_df)), function(i) {
    ev  <- events_df[i, ]
    # coordinates (0-based start -> 1-based by adding 1)
    uES <- ev$upstreamES + 1;   uEE <- ev$upstreamEE
    sES <- ev$exonStart_0base + 1; sEE <- ev$exonEnd
    dES <- ev$downstreamES + 1; dEE <- ev$downstreamEE
    gene <- gsub('"', '', ev$GeneID)

    gene_tx <- txdf$TXNAME[txdf$GENEID == gene]
    gene_tx_gr <- tx_exons[intersect(gene_tx, names(tx_exons))]

    if (length(gene_tx_gr) == 0) {
      inc_tx <- skip_tx <- character(0)
    } else {
      has_up   <- vapply(gene_tx_gr, has_exon, logical(1), uES, uEE)
      has_skip <- vapply(gene_tx_gr, has_exon, logical(1), sES, sEE)
      has_dn   <- vapply(gene_tx_gr, has_exon, logical(1), dES, dEE)
      inc_tx  <- names(gene_tx_gr)[has_up & has_skip & has_dn]
      skip_tx <- names(gene_tx_gr)[has_up & !has_skip & has_dn]
    }

    ps <- psi(inc_tx, skip_tx, ev$strand, counts1, counts2)
    g  <- gt_from_structure(gene, inc_tx, skip_tx)
    data.frame(
      event_type = "SE", ID = ev$ID, GeneID = gene,
      chr = ev$chr, strand = ev$strand,
      coord = sprintf("%d-%d|%d-%d|%d-%d", uES, uEE, sES, sEE, dES, dEE),
      n_inc_tx = length(inc_tx), n_skip_tx = length(skip_tx),
      n_manip_inc = g$n_manip_inc, n_manip_skip = g$n_manip_skip,
      PSI1 = round(ps[1], 4), PSI2 = round(ps[2], 4),
      dPSI = round(ps[2] - ps[1], 4),
      gt_positive = g$gt,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

# --- process A3SS / A5SS events ---------------------------------------------

process_AltSS <- function(events_df, etype, tx_exons, counts1, counts2) {
  # long form uses longExon; short form uses shortE (shortES:shortEE)
  rows <- lapply(seq_len(nrow(events_df)), function(i) {
    ev <- events_df[i, ]
    lES <- ev$longExonStart_0base + 1; lEE <- ev$longExonEnd
    sES <- ev$shortES + 1;             sEE <- ev$shortEE
    fES <- ev$flankingES + 1;          fEE <- ev$flankingEE
    gene <- gsub('"', '', ev$GeneID)

    gene_tx <- txdf$TXNAME[txdf$GENEID == gene]
    gene_tx_gr <- tx_exons[intersect(gene_tx, names(tx_exons))]

    if (length(gene_tx_gr) == 0) {
      long_tx <- short_tx <- character(0)
    } else {
      has_long  <- vapply(gene_tx_gr, has_exon, logical(1), lES, lEE)
      has_short <- vapply(gene_tx_gr, has_exon, logical(1), sES, sEE)
      has_flank <- vapply(gene_tx_gr, has_exon, logical(1), fES, fEE)
      long_tx  <- names(gene_tx_gr)[has_long  & has_flank]
      short_tx <- names(gene_tx_gr)[has_short & has_flank & !has_long]
    }

    ps <- psi(long_tx, short_tx, ev$strand, counts1, counts2)
    g  <- gt_from_structure(gene, long_tx, short_tx)
    data.frame(
      event_type = etype, ID = ev$ID, GeneID = gene,
      chr = ev$chr, strand = ev$strand,
      coord = sprintf("%d-%d|%d-%d|%d-%d", lES, lEE, sES, sEE, fES, fEE),
      n_inc_tx = length(long_tx), n_skip_tx = length(short_tx),
      n_manip_inc = g$n_manip_inc, n_manip_skip = g$n_manip_skip,
      PSI1 = round(ps[1], 4), PSI2 = round(ps[2], 4),
      dPSI = round(ps[2] - ps[1], 4),
      gt_positive = g$gt,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

# --- process RI events ------------------------------------------------------

process_RI <- function(events_df, tx_exons, counts1, counts2) {
  rows <- lapply(seq_len(nrow(events_df)), function(i) {
    ev <- events_df[i, ]
    # riES/riEE = the full retained intron exon (spanning the intron)
    rES <- ev$riExonStart_0base + 1; rEE <- ev$riExonEnd
    uES <- ev$upstreamES + 1;        uEE <- ev$upstreamEE
    dES <- ev$downstreamES + 1;      dEE <- ev$downstreamEE
    gene <- gsub('"', '', ev$GeneID)

    gene_tx <- txdf$TXNAME[txdf$GENEID == gene]
    gene_tx_gr <- tx_exons[intersect(gene_tx, names(tx_exons))]

    if (length(gene_tx_gr) == 0) {
      inc_tx <- skip_tx <- character(0)
    } else {
      has_retained <- vapply(gene_tx_gr, has_exon, logical(1), rES, rEE)
      has_up       <- vapply(gene_tx_gr, has_exon, logical(1), uES, uEE)
      has_dn       <- vapply(gene_tx_gr, has_exon, logical(1), dES, dEE)
      inc_tx  <- names(gene_tx_gr)[has_retained]
      skip_tx <- names(gene_tx_gr)[has_up & has_dn & !has_retained]
    }

    ps <- psi(inc_tx, skip_tx, ev$strand, counts1, counts2)
    g  <- gt_from_structure(gene, inc_tx, skip_tx)
    data.frame(
      event_type = "RI", ID = ev$ID, GeneID = gene,
      chr = ev$chr, strand = ev$strand,
      coord = sprintf("%d-%d|%d-%d|%d-%d", rES, rEE, uES, uEE, dES, dEE),
      n_inc_tx = length(inc_tx), n_skip_tx = length(skip_tx),
      n_manip_inc = g$n_manip_inc, n_manip_skip = g$n_manip_skip,
      PSI1 = round(ps[1], 4), PSI2 = round(ps[2], 4),
      dPSI = round(ps[2] - ps[1], 4),
      gt_positive = g$gt,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

# --- load events from fromGTF (full universe) and JCEC (tested subset) ------

read_rmats <- function(f) {
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", comment.char = "")
}

# collect tested IDs per event type from JCEC files
jcec_ids <- list()
for (etype in c("SE", "A3SS", "A5SS", "RI")) {
  f <- file.path(rmats_dir, paste0(etype, ".MATS.JCEC.txt"))
  if (file.exists(f)) {
    df <- read_rmats(f)
    df <- df[, !duplicated(names(df))]
    jcec_ids[[etype]] <- as.character(df$ID)
    cat(sprintf("  JCEC %s: %d tested events\n", etype, length(jcec_ids[[etype]])))
  }
}

# load each fromGTF universe, strip quotes from GeneID, and split by gene so a
# gene's events can be processed together and written to its own file.
load_split_by_gene <- function(fname) {
  f <- file.path(rmats_dir, fname)
  if (!file.exists(f)) return(list())
  df <- read_rmats(f)
  df$GeneID <- gsub('"', '', df$GeneID)
  split(df, df$GeneID)
}

cat("\nLoading fromGTF event universes...\n")
se_by_gene <- load_split_by_gene("fromGTF.SE.txt")
a3_by_gene <- load_split_by_gene("fromGTF.A3SS.txt")
a5_by_gene <- load_split_by_gene("fromGTF.A5SS.txt")
ri_by_gene <- load_split_by_gene("fromGTF.RI.txt")

all_genes <- sort(unique(c(names(se_by_gene), names(a3_by_gene),
                           names(a5_by_gene), names(ri_by_gene))))
cat(sprintf("  %d genes with at least one annotated event\n", length(all_genes)))

# --- per-gene processing ----------------------------------------------------

COLS <- c("event_type", "ID", "GeneID", "chr", "strand", "coord",
          "n_inc_tx", "n_skip_tx", "n_manip_inc", "n_manip_skip",
          "PSI1", "PSI2", "dPSI", "sim_type",
          "gt_positive", "rmats_tested")

# build the full GT table for a single gene (all its event types)
process_gene <- function(gene) {
  parts <- list()
  if (!is.null(se_by_gene[[gene]]))
    parts$SE   <- process_SE(se_by_gene[[gene]], tx_exons, counts1, counts2)
  if (!is.null(a3_by_gene[[gene]]))
    parts$A3SS <- process_AltSS(a3_by_gene[[gene]], "A3SS", tx_exons, counts1, counts2)
  if (!is.null(a5_by_gene[[gene]]))
    parts$A5SS <- process_AltSS(a5_by_gene[[gene]], "A5SS", tx_exons, counts1, counts2)
  if (!is.null(ri_by_gene[[gene]]))
    parts$RI   <- process_RI(ri_by_gene[[gene]], tx_exons, counts1, counts2)

  df <- bind_rows(parts)
  if (nrow(df) == 0) return(df)

  df$sim_type    <- get_sim_type(gene)
  df$rmats_tested <- mapply(function(etype, id) {
    ids <- jcec_ids[[etype]]
    !is.null(ids) && as.character(id) %in% ids
  }, df$event_type, df$ID)

  df[, COLS]
}

# atomic per-gene write: write to .tmp then rename, so a present .txt means
# the gene is fully done (safe for resume + parallel shards).
write_gene <- function(gene, df) {
  final <- file.path(out_dir, paste0(gene, ".txt"))
  tmp   <- paste0(final, ".tmp")
  write.table(df, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  file.rename(tmp, final)
}

# select this shard's genes (parallel runs), else all genes
genes <- all_genes
if (sharded) {
  genes <- all_genes[(seq_along(all_genes) - 1) %% n_shards == (shard_idx - 1)]
  cat(sprintf("  Shard handles %d genes\n", length(genes)))
}

n_done <- 0L; n_skip <- 0L
t0 <- Sys.time()
for (k in seq_along(genes)) {
  gene  <- genes[k]
  final <- file.path(out_dir, paste0(gene, ".txt"))
  if (file.exists(final)) { n_skip <- n_skip + 1L; next }  # resume: already done

  df <- process_gene(gene)
  if (nrow(df) > 0) write_gene(gene, df)
  n_done <- n_done + 1L

  if (k %% 200 == 0 || k == length(genes)) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  [%d/%d] done=%d skipped=%d (%.1fs)\n",
                k, length(genes), n_done, n_skip, el))
  }
}

cat(sprintf("\nProcessed %d genes (%d newly written, %d skipped/already done)\n",
            length(genes), n_done, n_skip))
cat(sprintf("Per-gene files in: %s\n", out_dir))
cat("Combine into a single table with shell (head + tail), not R.\n")
