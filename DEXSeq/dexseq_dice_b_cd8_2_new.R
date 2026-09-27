library(DEXSeq)
library(tidyverse)
library(BiocParallel)

# Set working directory
setwd("/data2/han_lab/carriehe/DEXSeq/DICE/b_vs_cd8")

#load .rds object back into R session
dxd <- readRDS("dexseq_part_1_new.rds")

BPPARAM = MulticoreParam(workers=15)

#Test for differential exon usage
dxd = testForDEU(dxd, BPPARAM = BPPARAM)

print("done testing for differential exon use")

#Estimate exon fold changes
dxd = estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = BPPARAM)

print ("done estimating exon fold change")

#Extract results
dexseq_results = DEXSeqResults(dxd)

print("done extracting results")

# Convert the results object to a standard data frame
results_df <- as.data.frame(dexseq_results)

# Find any columns that are lists and convert them to comma-separated strings
for (col_name in names(results_df)) {
  if (is.list(results_df[[col_name]])) {
    results_df[[col_name]] <- sapply(results_df[[col_name]], paste, collapse = ",")
  }
}

#Save results to text file
write.table(results_df,
            file = "DEXSeq_DICE_B_vs_CD8_results_new.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

#exons with significantly different usage between B and CD8 cells
sig_results = subset(results_df, padj < 0.05)
write.table(sig_results, file = "DEXSeq_DICE_B_vs_CD8_significant_new.txt", sep = "\t", quote = FALSE)





#nohup Rscript dexseq_dice_b_cd8_2_new.R > script_2_new.log 2>&1 &

