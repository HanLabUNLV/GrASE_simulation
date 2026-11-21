#!/usr/bin/env Rscript
#
# Per-side fp_breakdown table for one DTE gene, the post-hoc lfc_diff_net example,
# from the STRANDED MERGED EBapprox results -- the lane every other figure and
# table uses -- under the CURRENT call rule (pr_curves_three_levels_gtrule.R,
# GrASE_merged_default):
#   padj   plain BH over the raw p-values of every test, internal and TSS/TTS
#          pooled (p.adjust "BH"), not the nested BH padj stored in the files
#   called grase::is_significant() at padj < 0.01 with exontest.R's default gates
#          (|effect_size| >= 0.4, no delta_pi gate, d_support >= 10) and the
#          lfc_diff_net gate OFF,
#          so the table holds the calls the filter acts on
#   status retained (lfc_diff_net > 0) | removed_by_lfc_filter
# The two classes number their events independently, so a side is keyed by
# (src, event, comparison).
#
# fp_type: TP | paired_FP | wrong_side_FP | unrelated_FP, on exonic parts: TP when
# the side's distinct set holds a part of the changed transcript.
#
# Usage: Rscript make_fp_breakdown_example.R [gene, default ENSG00000114491.13]
# Output: example/<gene>_fp_breakdown.txt

suppressPackageStartupMessages({ library(dplyr); library(grase) })

args <- commandArgs(trailingOnly = TRUE)
GENE <- if (length(args) >= 1) args[1] else "ENSG00000114491.13"
PADJ_THR <- 0.01
SIM_DIR <- path.expand("~/GrASE_simulation")
files <- c(internal = "bipartition.merged.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt",
           TSSTTS   = "bipartition.merged.TSSTTS.stranded.test.EBapprox/test_bipartition.merged_betabinom_EBapprox.annotated.txt")

parse_exparts <- function(x) {
  if (is.na(x) || trimws(x) == "" || x == "NA") return(character(0))
  trimws(unlist(strsplit(x, ",")))
}
read_ann <- function(f) {
  read.table(file.path(SIM_DIR, f), header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", comment.char = "", na.strings = c("NA", ""),
             colClasses = c(source = "character", sink = "character", comparison = "character",
                            setdiff1 = "character", setdiff2 = "character"))
}

cat("loading test results...\n")
all <- bind_rows(lapply(names(files), function(s) mutate(read_ann(files[[s]]), src = s)))
## plain BH over every test of both classes, as the benchmark computes it
all$padj_bh <- p.adjust(all$p.value, method = "BH")
gate <- function(d, delta) is_significant(d, PADJ_THR, delta, 0, 0, 10, padj = d$padj_bh,
                                          use_perbase = "both", min_effect = 0.4)
tests <- all %>% filter(gene == GENE)
tests$pre   <- gate(tests, -Inf)
tests$final <- gate(tests, 0)
tests <- tests %>% mutate(setdiff_active = ifelse(comparison == "diff1_vs_ref", setdiff1, setdiff2),
                          event_base = paste(src, event))
cat(sprintf("  %d rows for %s\n", nrow(tests), GENE))

gt <- read.table(file.path(SIM_DIR, "results", "sim_exon_info", paste0(GENE, ".exonic_parts_fc.txt")),
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
gt_pos <- gt$exonic_part[!is.na(gt$fold_change)]
cat(sprintf("  GT-positive exons (%d): %s\n", length(gt_pos), paste(gt_pos, collapse = ", ")))

sig <- tests %>% filter(pre) %>%
  mutate(setdiff_exons = lapply(setdiff_active, parse_exparts),
         any_tp = sapply(setdiff_exons, function(ex) any(ex %in% gt_pos)))
cat(sprintf("  %d sides called before the lfc_diff_net filter\n", nrow(sig)))

tp_bases <- unique(sig$event_base[sig$any_tp])
bubble_pos <- tests %>%
  mutate(any_tp = sapply(lapply(setdiff_active, parse_exparts), function(ex) any(ex %in% gt_pos))) %>%
  group_by(event_base) %>% summarise(bubble_has_gt_pos = any(any_tp), .groups = "drop")
sig <- sig %>% left_join(bubble_pos, by = "event_base") %>%
  mutate(fp_type = case_when(any_tp ~ "TP",
                             event_base %in% tp_bases ~ "paired_FP",
                             bubble_has_gt_pos ~ "wrong_side_FP",
                             TRUE ~ "unrelated_FP"))

out <- sig %>%
  mutate(n_exons = sapply(setdiff_exons, length),
         status  = ifelse(final, "retained", "removed_by_lfc_filter"),
         setdiff = sapply(setdiff_exons, function(x) paste(x, collapse = ",")),
         lfc_diff_net = round(lfc_diff_net, 3), padj = signif(padj_bh, 3),
         delta_pi = round(delta_pi, 3), effect_size = round(effect_size, 2)) %>%
  select(gene, src, event, comparison, ref_ex_part, setdiff, lfc_diff_net, padj,
         delta_pi, effect_size, fp_type, n_exons, status)
out <- bind_rows(filter(out, fp_type != "TP"), filter(out, fp_type == "TP"))

cat("\nSummary:\n"); print(table(out$fp_type, out$status))
out_file <- file.path(SIM_DIR, "example", paste0(GENE, "_fp_breakdown.txt"))
write.table(out, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nWrote: %s\n", out_file))
