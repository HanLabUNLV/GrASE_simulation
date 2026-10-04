#!/usr/bin/env Rscript
#
# scripts/infer_bipartition_gt.R
#
# GT_rule -- the DEFAULT ground truth for evaluation, on GrASE's unit.
# Handles BOTH the exonic bipartition unit (default) and the merged exon+SJ unit
# (--merged). These were two scripts until 2026-09-14; 44% of their code was
# identical, including the rule itself, so a correction applied to one and not
# the other would have produced two disagreeing ground truths.
#
# THE RULE (identical in both modes)
#   y_c(b)    = len(b) * sum of TPM_c over transcripts covering b
#   y_m_c(b)  = len(b) * sum of TPM_c over (transcripts covering b) INTERSECT manipulated
#   y_c(X)    = sum of y_c(b)   over b in X          (extend additively to a bin set)
#   y_m_c(X)  = sum of y_m_c(b) over b in X
#   pi_c      = y_c(D) / ( y_c(D) + y_c(S) )
#   positive  = y_m_1(D) != y_m_2(D)  AND  |pi_2 - pi_1| > 0
#
# Both conjuncts are needed, for different reasons:
#   mass changed  excludes manipulations OUTSIDE the tested D that move pi through
#                 the denominator (a DTE transcript in the shared reference, or on
#                 the complementary path). 16.4% of labelled tests are this case.
#   pi moved      excludes manipulations inside D that contribute proportionally to
#                 D and S, leaving pi unmoved -- a true null.
#
# WHAT --merged CHANGES: only how y_c(D) and y_m_c(D) are built, and only for a
# side whose exonic distinct set is empty. Such a side has no D bins, so in
# exonic mode it is dropped (`if (!length(Dp)) return(NULL)`); merge_exon_sj_counts.R
# gives it an exclusive junction instead, and here the junction stands in for the
# bins:
#   exon-sourced side  y_c(D) = sum over b in D of len(b) * sum TPM_c(tx covering b)
#   sj-sourced side    y_c(D) = sum over j in J of  a     * sum TPM_c(tx using j)
# The denominator S is UNCHANGED in both, which is what puts pi on a mixed scale
# for substituted sides (see --min_dpi_sj in exontest.R).
#
# `a` -- the per-junction analogue of bin length -- is the number of fragment
# positions that yield a split read across a junction, L - 2*overhang + 1.
# Default 91 (100 nt mates, STAR alignSJoverhangMin 5); --a overrides it.
# The BINARY label does not depend on `a`: mass-changed is |a*(X2-X1)| > 0 and a>0
# cancels; pi_2 == pi_1 <=> yD1*yS2 == yD2*yS1, where a cancels. Only the reported
# magnitude true_dpi is on the `a` scale, and that is a stratification covariate.
# len(b) does NOT cancel -- it determines which side carries mass -- so uniform
# coverage along a transcript is the load-bearing assumption here.
#
# Membership caution: the structural condition uses the SAME per-bin transcript
# membership that builds y_c(D). transcripts1/transcripts2 list only transcripts
# traversing that source->sink route; a manipulated transcript can cover a D bin
# and move y_c(D) without appearing in either. Using the path sets here mislabels
# thousands of genuine positives as drift. They are kept only as diagnostics
# (manip_path_mass1/2).
#
# Binary labels; nothing excluded. NO detectability floor: magnitude is a
# stratification covariate only. Zero-mass tests are labelled, not dropped.
# A condition with no mass on D + S is given pi_c = 0, so: empty in both
# conditions -> true_dpi = 0, negative; empty in one -> the bubble appears or
# disappears, true_dpi = pi_other, positive iff D was manipulated.
#
# Computed entirely from design parameters: tpms from simulate.rda (NOT
# sim.counts.mat, which is length-weighted and gives the wrong scale), bin length
# and per-bin transcript lists from dexseq.gff. No reads, quantification or
# fitted values.
#
# Usage:
#   Rscript scripts/infer_bipartition_gt.R                 # exonic unit
#   Rscript scripts/infer_bipartition_gt.R --merged        # merged exon+SJ unit
#   GT_RULE_STRANDED=1 Rscript scripts/infer_bipartition_gt.R --merged
# Out:
#   results/gt/bipartition_gt[.stranded].txt          (default)
#   results/gt/bipartition_merged_gt[.stranded].txt   (--merged)
#   column GT_rule_bipartition; gt_positive kept as an alias
suppressPackageStartupMessages({library(dplyr); library(parallel); library(optparse)})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--merged", action = "store_true", default = FALSE,
              help = "score the MERGED exon+SJ unit: junction-sourced sides get a junction numerator instead of being dropped [default: exonic unit]"),
  make_option("--a",     type = "double",  default = 91,
              help = "fragment positions per junction, the analogue of bin length; --merged only [default 91]"),
  make_option("--cores", type = "integer", default = 16))))
MERGED <- isTRUE(opt$merged)
A_JUNC <- opt$a
## igraph/grase are needed only to map junctions to transcripts, which only the
## merged mode does. Keeping them out of the exonic path means the default run
## has no graph dependency, as it did before the two scripts were folded.
if (MERGED) suppressPackageStartupMessages({library(igraph); library(grase)})
cat(sprintf("mode: %s unit%s\n", if (MERGED) "MERGED exon+SJ" else "exonic",
            if (MERGED) sprintf("  (a = %g)", A_JUNC) else ""))

BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
GFF  <- file.path(BASE, "dexseq.gff")
GML  <- file.path(BASE, "graphml.v28")
OUT  <- file.path(BASE, "results/gt")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
load(file.path(BASE, "swimdown/simulate/data/simulate.rda"))
p1 <- tpms[,1]; p2 <- tpms[,2]

get_st <- function(g) ifelse(g %in% dte.genes, "DTE", ifelse(g %in% dtu.genes, "DTU",
                      ifelse(g %in% dge.genes, "DGE", "Background")))
parse_parts <- function(x) { v <- strsplit(as.character(x), ",")[[1]]
  trimws(v[nchar(trimws(v)) > 0 & trimws(v) != "NA"]) }
parse_tx <- function(x) { v <- strsplit(as.character(x), "[,+]")[[1]]
  trimws(v[nchar(trimws(v)) > 0 & trimws(v) != "NA"]) }

parse_gff <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ep <- ln[grepl("\texonic_part\t", ln, fixed = TRUE)]
  if (!length(ep)) return(NULL)
  f  <- strsplit(ep, "\t", fixed = TRUE)
  st <- as.numeric(vapply(f, `[`, character(1), 4))
  en <- as.numeric(vapply(f, `[`, character(1), 5))
  at <- sub("^([^\t]*\t){8}", "", ep)
  bn <- paste0("E", sub('.*exonic_part_number "([^"]+)".*', "\\1", at))
  tx <- strsplit(sub('.*transcripts "([^"]+)".*', "\\1", at), "+", fixed = TRUE)
  list(len = setNames(en - st + 1, bn), tx = setNames(tx, bn))
}

## junction id -> transcripts, keyed exactly as label_bipartition_introns builds
## them: chr:min(from_pos,to_pos):max(from_pos,to_pos) over the "in" edges.
junction_tx <- function(gene, chr) {
  p <- file.path(GML, paste0(gene, ".graphml"))
  if (!file.exists(p)) return(list())
  ge <- precompute_gene_graph(read_graph(p, format = "graphml"))
  ii <- which(ge$ex_or_in == "in" & !is.na(ge$from_pos) & !is.na(ge$to_pos))
  if (!length(ii)) return(list())
  key <- paste0(chr, ":", pmin(ge$from_pos[ii], ge$to_pos[ii]),
                     ":", pmax(ge$from_pos[ii], ge$to_pos[ii]))
  out <- lapply(seq_along(ii), function(k)
    ge$tx_cols[ge$tx_mat[ii[k], ] > 0L])
  names(out) <- key
  ## the same junction can appear on more than one edge row; union them
  if (anyDuplicated(key)) out <- tapply(out, key, function(z) unique(unlist(z)))
  out
}

## --- merged tests -----------------------------------------------------------
## GT_RULE_STRANDED=1 scores the strand-reconstructed merged tests instead.
## Provenance (which side used a junction) is unchanged -- the merge rule and the
## sjcnt inputs are identical; only the exonic counts differ.
STRANDED <- nzchar(Sys.getenv("GT_RULE_STRANDED"))
files <- if (MERGED) {
  if (STRANDED)
    c(internal = "bipartition.merged.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt",
      TSSTTS   = "bipartition.merged.TSSTTS.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt")
  else
    c(internal = "bipartition.merged.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt",
      TSSTTS   = "bipartition.merged.TSSTTS.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt")
} else {
  if (STRANDED)
    c(internal = "bipartition.internal.stranded.test.EBapprox/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
      TSSTTS   = "bipartition.TSSTTS.stranded.test.EBapprox/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
  else
    c(internal = "bipartition.test.fulldesign/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
      TSSTTS   = "bipartition.test.fulldesign/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
}
GT_OUT <- sprintf("bipartition%s_gt%s.txt",
                  if (MERGED) "_merged" else "",
                  if (STRANDED) ".stranded" else "")
files <- files[file.exists(file.path(BASE, files))]
if (!length(files)) stop("no annotated test files found for this mode")
cat("reading:", paste(names(files), collapse = ", "), "\n")
d <- bind_rows(lapply(names(files), function(s) {
  x <- read.table(file.path(BASE, files[[s]]), header = TRUE, sep = "\t", quote = "",
                  comment.char = "", stringsAsFactors = FALSE)
  x$src <- s
  x[, c("src","gene","event","comparison","ref_ex_part","setdiff1","setdiff2",
        "transcripts1","transcripts2","padj","lfc_diff_net","delta_pi")]
}))

## --- provenance (which side used a junction, and which junction) -------------
## Merged mode only. The exonic tests have no junction-sourced sides by
## construction, so every side is exon-sourced and the row branch below takes
## the bin path unconditionally -- identical to the pre-fold exonic script.
if (MERGED) {
  ## Provenance comes from the SIDECAR, by design, not from the annotated output.
  ## exontest.R carries diff1_source/diff2_source per row (4 characters, needed
  ## on every row by add_significant), but deliberately NOT intron_distinct1/2:
  ## those average 123 characters and the annotated output repeats each event
  ## ~14x on a 10-contrast run, so carrying them would add 1.47 GB to DICE's
  ## 4.30 GB TSS/TTS file to store 422,081 facts fourteen times. This script is
  ## the only consumer, and it needs them once per event -- which is exactly the
  ## shape scripts/extract_merge_provenance.py writes.
  prov <- bind_rows(lapply(names(files), function(s_) {
    p <- file.path(OUT, sprintf("merge_provenance.%s.txt", s_))
    if (!file.exists(p)) stop("missing provenance file: ", p,
                              "\n  run scripts/extract_merge_provenance.py first")
    read.table(p, header = TRUE, sep = "\t", quote = "", comment.char = "",
               stringsAsFactors = FALSE, colClasses = "character")
  }))
  prov$event <- as.character(prov$event)
  d$event    <- as.character(d$event)
  d <- left_join(d, prov, by = c("src","gene","event"))
  cat(sprintf("merged tests: %d  (provenance matched %d)\n",
              nrow(d), sum(!is.na(d$diff1_source))))
} else {
  d$event <- as.character(d$event)
  d$diff1_source <- "exon"; d$diff2_source <- "exon"
  d$intron_distinct1 <- NA_character_; d$intron_distinct2 <- NA_character_
  cat(sprintf("exonic tests: %d\n", nrow(d)))
}

d$sim_type <- get_st(d$gene)
d <- d[d$sim_type %in% c("DTE","DTU"), ]
cat(sprintf("tests in signal genes: %d  (%d genes)\n", nrow(d), length(unique(d$gene))))

## chromosome per gene, for junction keys (merged mode only)
chr_map <- if (MERGED)
  parse_gff_chr_map(file.path(BASE, "ref/gencode.v28.dexseq.bygene.gff")) else list()

gv <- function(v, u) sum(v[intersect(u, names(v))])
genes <- unique(d$gene)

res <- mclapply(genes, function(g) {
  gf <- file.path(GFF, paste0(g, ".dexseq.gff"))
  if (!file.exists(gf)) return(NULL)
  ann <- parse_gff(gf); if (is.null(ann)) return(NULL)
  gtx <- txdf$TXNAME[txdf$GENEID == g]
  tp1 <- p1[intersect(gtx, names(p1))]
  tp2 <- p2[intersect(gtx, names(p2))]
  manip_vec <- iso.dtu | iso.dte
  M <- names(which(manip_vec[intersect(gtx, names(manip_vec))]))
  bins <- names(ann$len)
  m1  <- vapply(bins, function(b) ann$len[[b]] * gv(tp1, ann$tx[[b]]), numeric(1))
  m2  <- vapply(bins, function(b) ann$len[[b]] * gv(tp2, ann$tx[[b]]), numeric(1))
  mm1 <- vapply(bins, function(b) ann$len[[b]] * gv(tp1, intersect(ann$tx[[b]], M)), numeric(1))
  mm2 <- vapply(bins, function(b) ann$len[[b]] * gv(tp2, intersect(ann$tx[[b]], M)), numeric(1))

  sub <- d[d$gene == g, ]
  jtx <- NULL   # built lazily; most genes have no sj-sourced side, and in
                # exonic mode no gene does, so this never fires there
  if (MERGED && any(sub$diff1_source == "sj" | sub$diff2_source == "sj", na.rm = TRUE)) {
    chr <- chr_map[[g]]
    jtx <- if (is.null(chr) || is.na(chr)) list() else junction_tx(g, chr)
  }

  rows <- lapply(seq_len(nrow(sub)), function(k) {
    r    <- sub[k, ]
    one  <- grepl("diff1", r$comparison)
    srcD <- if (one) r$diff1_source else r$diff2_source
    if (is.na(srcD) || srcD == "NA") return(NULL)
    Sp <- intersect(parse_parts(r$ref_ex_part), bins)
    if (!length(Sp)) return(NULL)
    yS1 <- sum(m1[Sp]); yS2 <- sum(m2[Sp])

    if (srcD == "sj") {
      J <- parse_parts(if (one) r$intron_distinct1 else r$intron_distinct2)
      J <- J[J %in% names(jtx)]
      if (!length(J)) return(NULL)
      txJ <- lapply(J, function(j) jtx[[j]])
      yD1 <- A_JUNC * sum(vapply(txJ, function(u) gv(tp1, u), numeric(1)))
      yD2 <- A_JUNC * sum(vapply(txJ, function(u) gv(tp2, u), numeric(1)))
      mD1 <- A_JUNC * sum(vapply(txJ, function(u) gv(tp1, intersect(u, M)), numeric(1)))
      mD2 <- A_JUNC * sum(vapply(txJ, function(u) gv(tp2, intersect(u, M)), numeric(1)))
      TD  <- unique(unlist(txJ))
      nD  <- length(J)
    } else {
      Dp <- intersect(parse_parts(if (one) r$setdiff1 else r$setdiff2), bins)
      if (!length(Dp)) return(NULL)
      yD1 <- sum(m1[Dp]);  yD2 <- sum(m2[Dp])
      mD1 <- sum(mm1[Dp]); mD2 <- sum(mm2[Dp])
      TD  <- unique(unlist(ann$tx[Dp]))
      nD  <- length(Dp)
    }

    z1 <- (yD1 + yS1) <= 0; z2 <- (yD2 + yS2) <= 0
    ## A condition with no mass on D + S has no proportion; define pi_c = 0 there.
    ## This one convention covers both boundary cases: empty in both conditions
    ## gives true_dpi = 0 (negative); empty in one gives true_dpi = pi_other, which
    ## is > 0 exactly when D carries mass in the other condition, and otherwise
    ## manip_mass_changed is already FALSE because manipulated mass <= total mass.
    pi1 <- if (z1) 0 else yD1/(yD1+yS1)
    pi2 <- if (z2) 0 else yD2/(yD2+yS2)
    TS  <- unique(unlist(ann$tx[Sp]))
    Ml  <- intersect(M, union(TD, TS))
    inD <- intersect(Ml, TD); inS <- intersect(Ml, TS)
    tpath <- if (one) parse_tx(r$transcripts1) else parse_tx(r$transcripts2)
    opath <- if (one) parse_tx(r$transcripts2) else parse_tx(r$transcripts1)
    Mp <- intersect(M, tpath)
    data.frame(gene = g, sim_type = r$sim_type, src = r$src, event = r$event,
               comparison = r$comparison, D_source = srcD, n_D = nD,
               true_dpi = abs(pi2 - pi1), obs_dpi = abs(r$delta_pi),
               mass_zero_c1 = z1, mass_zero_c2 = z2,
               padj = r$padj, lfc_diff_net = r$lfc_diff_net,
               n_manip_gene = length(M), n_manip_scope = length(Ml),
               n_manip_tested_path = length(Mp),
               n_manip_other_path = length(intersect(M, opath)),
               manip_massD1 = mD1, manip_massD2 = mD2,
               manip_mass_changed = abs(mD2 - mD1) > 1e-10,
               manip_path_mass1 = gv(tp1, Mp), manip_path_mass2 = gv(tp2, Mp),
               manip_in_D = length(inD), manip_in_S = length(inS),
               manip_D_only = length(setdiff(inD, inS)),
               manip_S_only = length(setdiff(inS, inD)),
               manip_both   = length(intersect(inD, inS)),
               stringsAsFactors = FALSE)
  })
  bind_rows(rows)
}, mc.cores = opt$cores)

r <- bind_rows(res)
cat(sprintf("\nscored %d tests\n", nrow(r)))

## --- the rule, verbatim from infer_bipartition_gt.R -------------------------
EPS <- 1e-10
## true_dpi is already complete (pi_c = 0 for an empty condition, set above).
## true_dpi_eff is kept as an alias because metric_level_comparison_gtrule.R
## bins on it. It used to force 1 for a bubble empty in one condition, which put
## every appearing bubble in the top stratum; it is now the actual |pi_2 - pi_1|.
## The binary label is unchanged either way.
r$true_dpi_eff <- r$true_dpi
r$gt_positive  <- r$manip_mass_changed & !is.na(r$true_dpi_eff) & r$true_dpi_eff > EPS
r$GT_rule_bipartition <- r$gt_positive

write.table(r, file.path(OUT, GT_OUT), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d tests -> %s\n", nrow(r), file.path(OUT, GT_OUT)))

cat("\n=== GT totals by D provenance ===\n")
print(table(D_source = r$D_source, gt_positive = r$gt_positive))
cat("\n=== GT totals by bipartition type ===\n")
print(table(src = r$src, D_source = r$D_source, gt_positive = r$gt_positive))
cat("\n=== by sim_type ===\n")
print(table(sim_type = r$sim_type, gt_positive = r$gt_positive))
cat("\n=== positives stratified by realized effect (reporting only) ===\n")
bs <- cut(r$true_dpi_eff, c(-Inf,0.05,0.1,0.2,Inf), c("<0.05","0.05-0.1","0.1-0.2",">=0.2"))
print(table(D_source = r$D_source[r$gt_positive], stratum = bs[r$gt_positive]))
