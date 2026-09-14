#!/usr/bin/env Rscript
#
# scripts/infer_lsv_junction_gt.R
#
# GT_rule for MAJIQ, at the JUNCTION level -- the same rule as
# scripts/infer_bipartition_gt.R, instantiated on MAJIQ's unit.
#
# The rule is a template with two slots: every tool here tests a proportion, so
# each has a numerator set and a reference set.
#
#   GT_rule(unit) = [ manipulated mass in the unit's NUMERATOR changed ]
#               AND
#               [ the unit's tested RATIO moved ]
#
#   tool    unit          numerator                 reference
#   GrASE   bipartition   D                         S
#   MAJIQ   junction      that junction's users     the LSV's other junctions
#
# So here:
#   PSI_j(c) = TPM_c(users of j) / TPM_c(users of any junction in the LSV)
#   GT_rule      = manipulated TPM mass among j's OWN users changed
#              AND |PSI_j(2) - PSI_j(1)| > 0
#
# The structural conjunct is required because PSI is compositional: when the
# manipulated transcript uses junction j, EVERY other junction in that LSV has
# its PSI moved by the shared denominator. Those siblings are the analogue of a
# DTE transcript sitting in a GrASE bipartition's shared reference, which is
# negative under GT_rule. Without it, an LSV with five junctions and one manipulated
# transcript scores two positives instead of one.
#
# Binary labels; nothing excluded; NO detectability floor (magnitude is a
# stratification covariate only -- metric_level_comparison.R:37,
# gene_level_metrics.R:39, transcript_level_metrics.R:48).
#
# Computed from tpms in simulate.rda, NOT sim.counts.mat: junction-spanning reads
# go as counts(t)/efflen(t), i.e. as TPM, so summing whole-transcript counts
# weights each transcript by its length and is the wrong scale.
#
# NOTE this is MAJIQ LSV junctions. scripts/infer_rmats_junctions_gt.R is the separate
# structural junction GT built for rMATS; scripts/infer_lsv_gt.old.R is the LSV-level
# structural GT. Different units -- do not join them.
#
# In:    results/gt/lsv_junction_keys.txt  (from make_lsv_junction_keys.R)
# Usage: Rscript scripts/infer_lsv_junction_gt.R
# Out:   results/gt/lsv_junction_gt.txt  (column GT_rule_junction, one row per junction)
#        results/gt/lsv_gt.txt           (column GT_rule_lsv, one row per LSV)

suppressPackageStartupMessages({ library(dplyr); library(rtracklayer) })
BASE <- "/mnt/data1/home/mirahan/GrASE_simulation"
EPS  <- 1e-10
load(file.path(BASE, "swimdown/simulate/data/simulate.rda"))
p1 <- tpms[,1]; p2 <- tpms[,2]
manip_vec <- iso.dtu | iso.dte          # aligned; c() would duplicate names
M <- names(which(manip_vec))
cat(sprintf("manipulated transcripts: %d\n", length(M)))

STRANDED <- nzchar(Sys.getenv("STRANDED"))
KEYS_IN <- if (STRANDED) "results/gt/lsv_junction_keys.stranded.txt" else "results/gt/lsv_junction_keys.txt"
j <- read.table(file.path(BASE, KEYS_IN), header = TRUE,
                sep = "\t", quote = "", stringsAsFactors = FALSE)
cat(sprintf("junctions: %d  LSVs: %d\n", nrow(j), length(unique(j$lsv_id))))

cat("Building tx_junc from GTF...\n")
ex <- import(file.path(BASE, "ref/gencode.v28.annotation.gtf"),
             feature.type = "exon", colnames = "transcript_id")
o <- order(ex$transcript_id, start(ex)); tid <- ex$transcript_id[o]
st <- start(ex)[o]; en <- end(ex)[o]; n <- length(tid); sm <- tid[-n] == tid[-1]
tx_junc <- split(paste(en[-n], st[-1], sep = "-")[sm], tid[-n][sm])

gv <- function(v, u) sum(v[intersect(u, names(v))])
j$true_dpsi_tpm <- NA_real_; j$manip_mass1 <- NA_real_; j$manip_mass2 <- NA_real_
j$n_manip_own <- NA_integer_; j$lsv_zero_c1 <- FALSE; j$lsv_zero_c2 <- FALSE
sig_idx <- which(j$sim_type %in% c("DTE","DTU"))
by_lsv  <- split(sig_idx, j$lsv_id[sig_idx])
cat(sprintf("scoring %d signal-gene LSVs...\n", length(by_lsv)))
k <- 0L
for (lsv in names(by_lsv)) {
  k <- k + 1L
  idx <- by_lsv[[lsv]]; g <- j$gene[idx[1]]; lj <- j$junction[idx]
  gtx <- intersect(txdf$TXNAME[txdf$GENEID == g], names(tx_junc))
  if (!length(gtx)) next
  tp1 <- p1[intersect(gtx, names(p1))]; tp2 <- p2[intersect(gtx, names(p2))]
  users <- lapply(lj, function(x) gtx[vapply(tx_junc[gtx], function(v) x %in% v, logical(1))])
  allu <- unique(unlist(users))
  t1 <- gv(tp1, allu); t2 <- gv(tp2, allu)
  z1 <- t1 <= 0; z2 <- t2 <= 0
  j$lsv_zero_c1[idx] <- z1; j$lsv_zero_c2[idx] <- z2
  if (!z1 && !z2)
    j$true_dpsi_tpm[idx] <- vapply(users, function(u)
        abs(gv(tp2,u)/t2 - gv(tp1,u)/t1), numeric(1))
  # manipulated mass among each junction's OWN users
  own <- lapply(users, function(u) intersect(u, M))
  j$n_manip_own[idx] <- vapply(own, length, integer(1))
  j$manip_mass1[idx] <- vapply(own, function(u) gv(tp1,u), numeric(1))
  j$manip_mass2[idx] <- vapply(own, function(u) gv(tp2,u), numeric(1))
  if (k %% 3000 == 0) cat(".", k, "")
}
cat("\n")
null_idx <- j$sim_type %in% c("Background","DGE")
j$true_dpsi_tpm[null_idx] <- 0
j$manip_mass1[null_idx] <- 0; j$manip_mass2[null_idx] <- 0; j$n_manip_own[null_idx] <- 0L

j$manip_mass_changed <- !is.na(j$manip_mass1) & !is.na(j$manip_mass2) &
                        abs(j$manip_mass2 - j$manip_mass1) > EPS
one_zero <- xor(j$lsv_zero_c1, j$lsv_zero_c2)
j$true_dpsi_eff <- ifelse(one_zero, 1, ifelse(is.na(j$true_dpsi_tpm), 0, j$true_dpsi_tpm))
j$gt_positive_sym <- j$manip_mass_changed & j$true_dpsi_eff > EPS
j$GT_rule_junction <- j$gt_positive_sym    # canonical name; gt_positive_sym kept as an alias

OUT <- file.path(BASE, if (STRANDED) "results/gt/lsv_junction_gt.stranded.txt"
                       else "results/gt/lsv_junction_gt.txt")
write.table(j, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s\n", OUT))

## --- LSV-level rollup ------------------------------------------------------
## An LSV is GT_rule-positive iff at least one of its junctions is. This is the unit
## MAJIQ is conventionally evaluated on (an LSV is called if any junction passes
## threshold), so the rollup keeps the junction rule and the LSV rule consistent
## instead of letting them drift apart.
lsv <- j %>%
  group_by(lsv_id, gene, sim_type, lsv_type, K) %>%
  summarise(n_junc_pos     = sum(GT_rule_junction, na.rm = TRUE),
            GT_rule_lsv    = any(GT_rule_junction, na.rm = TRUE),
            n_junc_called  = sum(sig_p95, na.rm = TRUE),
            called_p95     = any(sig_p95, na.rm = TRUE),
            max_true_dpsi  = suppressWarnings(max(true_dpsi_eff, na.rm = TRUE)),
            n_manip_users  = sum(n_manip_own, na.rm = TRUE),
            .groups = "drop")
lsv$max_true_dpsi[!is.finite(lsv$max_true_dpsi)] <- NA_real_
LOUT <- file.path(BASE, if (STRANDED) "results/gt/lsv_gt.stranded.txt"
                        else "results/gt/lsv_gt.txt")
write.table(lsv, LOUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nwrote %d LSVs -> %s (GT_rule, rolled up from junctions)\n", nrow(lsv), LOUT))
cat("\n=== LSV-level GT_rule ===\n")
print(table(sim_type = lsv$sim_type, GT_rule_lsv = lsv$GT_rule_lsv))
cat("\n  positive junctions per positive LSV:\n")
print(table(lsv$n_junc_pos[lsv$GT_rule_lsv]))

## comparison with the previous structural LSV GT
sf <- file.path(BASE, "results/sim_lsv_gt.txt")
if (file.exists(sf)) {
  old <- read.table(sf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cm <- merge(lsv[, c("lsv_id","sim_type","GT_rule_lsv","max_true_dpsi")],
              old[, c("lsv_id","gt_positive")], by = "lsv_id")
  cm <- cm[cm$sim_type %in% c("DTE","DTU"), ]
  cat("\n=== LSV GT_rule vs the previous structural GT (signal genes) ===\n")
  print(table(structural = cm$gt_positive, GT_rule_lsv = cm$GT_rule_lsv))
  cat(sprintf("  structural-positive but GT_rule-negative: %d\n",
              sum(cm$gt_positive & !cm$GT_rule_lsv)))
  cat(sprintf("  GT_rule-positive but structural-negative: %d\n",
              sum(!cm$gt_positive & cm$GT_rule_lsv)))
}

cat("\n=== junction GT ===\n")
cat(sprintf("  TPM-scale, floor 0.01               : %6d\n",
    sum(!is.na(j$true_dpsi_tpm) & j$true_dpsi_tpm >= 0.01)))
cat(sprintf("  SYMMETRIC (own users + moved, no floor): %4d\n", sum(j$gt_positive_sym)))
cat("\nby sim_type:\n"); print(table(sim_type = j$sim_type, gt_sym = j$gt_positive_sym))

cat("\n=== does the compositional doubling go away? ===\n")
s <- j[j$sim_type %in% c("DTE","DTU"), ]
new <- s %>% group_by(lsv_id) %>% summarise(n = sum(gt_positive_sym), .groups="drop")
cat("  GT-positive junctions per LSV, SYMMETRIC:\n")
print(table(new$n[new$n > 0]))

cat("\n=== sibling junctions: PSI moves but own users not manipulated ===\n")
sib <- s$true_dpsi_eff > EPS & !s$manip_mass_changed
cat(sprintf("  junctions moving with no manipulated own-user: %d  (now NEGATIVE)\n", sum(sib, na.rm=TRUE)))

cat("\n=== decision-space FPR under the symmetric GT ===\n")
tab <- j %>% group_by(sim_type) %>%
  summarise(decisions = n(), neg = sum(!gt_positive_sym),
            FP = sum(sig_p95 & !gt_positive_sym),
            FPR = round(sum(sig_p95 & !gt_positive_sym)/max(sum(!gt_positive_sym),1), 5),
            .groups = "drop")
print(as.data.frame(tab), row.names = FALSE)
