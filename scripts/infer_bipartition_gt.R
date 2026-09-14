#!/usr/bin/env Rscript
#
# scripts/infer_bipartition_gt.R
#
# GT_rule -- the DEFAULT ground truth for evaluation, on GrASE's unit.
#
# Stage 2 of ground-truth generation. Stage 1 (infer_diff_exons_gt.R) labels
# exonic PARTS at three stringency levels (level i/ii/iii). This labels
# bipartition TESTS and needs
# GrASE's test output (setdiff1/setdiff2/ref_ex_part), so it must run after the
# tests. GT_rule and the part-level levels live on different units -- do not
# join them as if they shared one.
#
# Rule:
#   mm_c[b] = len(b) * sum of TPM_c over (transcripts covering b) INTERSECT manipulated
#   mD(c)   = sum of mm_c[b] over b in D            manipulated contribution to y_D
#   y_X(c)  = sum of len(b) * sum of TPM_c over transcripts covering b, over b in X
#   pi(c)   = y_D(c) / (y_D(c) + y_S(c))
#   GT_rule_bipartition     = mD(1) != mD(2)  AND  |pi(2) - pi(1)| > 0
#
# Both conjuncts are needed, for different reasons:
#   mD changed   excludes manipulations OUTSIDE the tested D that move pi through
#                the denominator (a DTE transcript in the shared reference, or on
#                the complementary path).
#   pi moved     excludes manipulations inside D that contribute proportionally to
#                D and S, leaving pi unmoved -- a true null.
#
# Binary labels; nothing excluded. NO detectability floor: magnitude is a
# stratification covariate only (metric_level_comparison.R:37,
# gene_level_metrics.R:39, transcript_level_metrics.R:48). Zero-mass tests are
# labelled, not dropped -- zero in both conditions is negative; zero in one only
# means the bubble appears/disappears, positive iff D was manipulated.
#
# Membership caution: the structural condition uses the SAME per-bin transcript
# membership that builds y_D. transcripts1/transcripts2 list only transcripts
# traversing that source->sink route; a manipulated transcript can cover a D part
# and move y_D without appearing in either. Using the path sets here mislabels
# thousands of genuine positives as drift. They are kept only as diagnostics
# (manip_path_mass1/2).
#
# Computed entirely from design parameters: tpms from simulate.rda (NOT
# sim.counts.mat, which is length-weighted and gives the wrong scale), bin length
# and per-bin transcript lists from dexseq.gff. No reads, quantification or
# fitted values. One assumption: uniform coverage along a transcript, which makes
# TPM the correct weight and effective length cancel.
#
# Usage: Rscript scripts/infer_bipartition_gt.R
# Out:   results/gt/bipartition_gt.txt   (column GT_rule_bipartition; gt_positive kept as alias)

suppressPackageStartupMessages(library(dplyr))
BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
GFF  <- file.path(BASE, "dexseq.gff")
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

# parser from infer_diff_exons_gt.R, extended to keep coordinates
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

## Source tests are overridable so the same rule can be applied to the
## strand-reconstructed run. GT_RULE_STRANDED=1 points at the stranded tests and
## writes bipartition_gt.stranded.txt; unset reproduces the original exactly.
STRANDED <- nzchar(Sys.getenv("GT_RULE_STRANDED"))
files <- if (STRANDED) {
  c("bipartition.internal.stranded.test.EBapprox/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
    "bipartition.TSSTTS.stranded.test.EBapprox/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
} else {
  c("bipartition.test.fulldesign/test_bipartition.internal_betabinom_EBapprox.annotated.txt",
    "bipartition.test.fulldesign/test_bipartition.TSSTTS_betabinom_EBapprox.annotated.txt")
}
GT_OUT <- if (STRANDED) "bipartition_gt.stranded.txt" else "bipartition_gt.txt"
d <- bind_rows(lapply(files, function(f) {
  x <- read.table(file.path(BASE, f), header = TRUE, sep = "\t", quote = "",
                  comment.char = "", stringsAsFactors = FALSE)
  x$src <- ifelse(grepl("TSSTTS", f), "TSSTTS", "internal")
  x[, c("src","gene","event","comparison","ref_ex_part","setdiff1","setdiff2",
        "transcripts1","transcripts2","padj","lfc_diff_net","delta_pi")]
}))
d$sim_type <- get_st(d$gene)
d <- d[d$sim_type %in% c("DTE","DTU"), ]
cat(sprintf("tests in signal genes: %d  (%d genes)\n", nrow(d), length(unique(d$gene))))

gv <- function(v, u) sum(v[intersect(u, names(v))])
res <- vector("list", length(unique(d$gene))); i <- 0L
for (g in unique(d$gene)) {
  i <- i + 1L
  gf <- file.path(GFF, paste0(g, ".dexseq.gff"))
  if (!file.exists(gf)) next
  ann <- parse_gff(gf); if (is.null(ann)) next
  gtx <- txdf$TXNAME[txdf$GENEID == g]
  tp1 <- p1[intersect(gtx, names(p1))]          # restrict to this gene once
  tp2 <- p2[intersect(gtx, names(p2))]
  # iso.dtu and iso.dte are aligned over the same transcript universe, so OR
  # them element-wise. c() would concatenate and give duplicate names, where
  # name-indexing silently returns only the iso.dtu half.
  manip_vec <- iso.dtu | iso.dte
  M <- names(which(manip_vec[intersect(gtx, names(manip_vec))]))
  bins <- names(ann$len)
  # bin mass does not depend on the test: compute once per gene
  m1 <- vapply(bins, function(b) ann$len[[b]] * gv(tp1, ann$tx[[b]]), numeric(1))
  m2 <- vapply(bins, function(b) ann$len[[b]] * gv(tp2, ann$tx[[b]]), numeric(1))
  # manipulated CONTRIBUTION to each bin, same membership as m1/m2 above.
  # (transcripts1/transcripts2 list only the transcripts traversing that
  #  source->sink route; a manipulated transcript can cover a D part and move
  #  y_D without appearing in either path set -- using the path sets here was
  #  inconsistent with how y_D is built.)
  mm1 <- vapply(bins, function(b) ann$len[[b]] * gv(tp1, intersect(ann$tx[[b]], M)), numeric(1))
  mm2 <- vapply(bins, function(b) ann$len[[b]] * gv(tp2, intersect(ann$tx[[b]], M)), numeric(1))
  sub <- d[d$gene == g, ]
  rows <- lapply(seq_len(nrow(sub)), function(k) {
    r  <- sub[k, ]
    Dp <- if (grepl("diff1", r$comparison)) parse_parts(r$setdiff1) else parse_parts(r$setdiff2)
    Sp <- parse_parts(r$ref_ex_part)
    Dp <- intersect(Dp, bins); Sp <- intersect(Sp, bins)
    if (!length(Dp) || !length(Sp)) return(NULL)
    yD1 <- sum(m1[Dp]); yS1 <- sum(m1[Sp])
    yD2 <- sum(m2[Dp]); yS2 <- sum(m2[Sp])
    z1 <- (yD1 + yS1) <= 0; z2 <- (yD2 + yS2) <= 0
    pi1 <- if (z1) NA_real_ else yD1/(yD1+yS1)
    pi2 <- if (z2) NA_real_ else yD2/(yD2+yS2)
    TD <- unique(unlist(ann$tx[Dp])); TS <- unique(unlist(ann$tx[Sp]))
    Ml <- intersect(M, union(TD, TS))
    inD <- intersect(Ml, TD); inS <- intersect(Ml, TS)
    # PATH-based placement: T(S) contains alternative-path transcripts too, so
    # bin-derived sets cannot distinguish reference from the other path.
    tpath <- if (grepl("diff1", r$comparison)) parse_tx(r$transcripts1)
             else parse_tx(r$transcripts2)
    opath <- if (grepl("diff1", r$comparison)) parse_tx(r$transcripts2)
             else parse_tx(r$transcripts1)
    Mp  <- intersect(M, tpath)          # path-set membership, kept as diagnostic
    mD1 <- sum(mm1[Dp]); mD2 <- sum(mm2[Dp])   # manipulated contribution to y_D
    data.frame(gene = g, sim_type = r$sim_type, src = r$src, event = r$event,
               comparison = r$comparison,
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
  res[[i]] <- bind_rows(rows)
  if (i %% 400 == 0) cat(".", i, "")
}
cat("\n")
r <- bind_rows(res)
cat(sprintf("\nscored %d tests\n", nrow(r)))

## --- the claim under test: manipulated tx outside the bipartition -> invisible
cat("\n=== tests where NO manipulated transcript is in D or S ===\n")
oo <- r[r$n_manip_scope == 0, ]
cat(sprintf("  n = %d   with true |delta_pi| >= 0.1 : %d   >= 0.01 : %d   max %.4f\n",
    nrow(oo), sum(oo$true_dpi >= 0.1), sum(oo$true_dpi >= 0.01),
    ifelse(nrow(oo), max(oo$true_dpi), NA)))

## --- AGREED RULE (final) ----------------------------------------------------
## gt_positive = manipulated mass on the TESTED PATH changed AND the tested
## ratio actually moved. No detectability floor: magnitude is a STRATIFICATION
## covariate only, per the project convention in metric_level_comparison.R:37,
## gene_level_metrics.R:39 and transcript_level_metrics.R:48.
## Binary labels only -- nothing is excluded from scoring.
EPS <- 1e-10
both_zero <- r$mass_zero_c1 & r$mass_zero_c2
one_zero  <- xor(r$mass_zero_c1, r$mass_zero_c2)
# bubble carries mass in one condition only -> the ratio is undefined but the
# unit did change; positive iff the tested path was the manipulated one.
r$true_dpi_eff <- ifelse(one_zero, 1, ifelse(both_zero, 0, r$true_dpi))
r$gt_positive  <- r$manip_mass_changed & !is.na(r$true_dpi_eff) & r$true_dpi_eff > EPS
r$GT_rule_bipartition <- r$gt_positive   # canonical name; gt_positive kept as an alias

cat(sprintf("\nzero-mass tests: both conditions %d, one condition %d\n",
            sum(both_zero), sum(one_zero)))
write.table(r, file.path(OUT, GT_OUT), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d tests -> %s (GT_rule_bipartition)\n", nrow(r), file.path(OUT, GT_OUT)))

cat("\n=== GT totals (binary, no floor, nothing excluded) ===\n")
print(table(sim_type = r$sim_type, gt_positive = r$gt_positive))
cat("\n=== positives STRATIFIED by realized effect (reporting only) ===\n")
BIN_BREAKS <- c(-Inf, 0.05, 0.1, 0.2, Inf); BIN_LABELS <- c("<0.05","0.05-0.1","0.1-0.2",">=0.2")
bs <- cut(r$true_dpi_eff, BIN_BREAKS, BIN_LABELS)
print(table(sim_type = r$sim_type[r$gt_positive], stratum = bs[r$gt_positive]))
cat("\n=== why each negative is negative ===\n")
why <- ifelse(r$gt_positive, "positive",
       ifelse(!r$manip_mass_changed & r$true_dpi_eff > EPS,
              "moved, but tested path not manipulated",
       ifelse(!r$manip_mass_changed, "tested path not manipulated, no move",
              "tested path manipulated, ratio unmoved")))
print(table(sim_type = r$sim_type, reason = why))
cat("\n=== sanity ===\n")
cat(sprintf("  tests moving with no manip tx in D or S: %d\n",
    sum(r$true_dpi > 0 & r$n_manip_scope == 0, na.rm = TRUE)))
cat(sprintf("  DTU pairs both on tested path called positive: %d (must be 0)\n",
    sum(r$gt_positive & r$sim_type == "DTU" & r$n_manip_tested_path == 2 &
        r$n_manip_other_path == 0)))
