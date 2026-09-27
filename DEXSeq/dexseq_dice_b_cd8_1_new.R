library(DEXSeq)
library(tidyverse)
library(BiocParallel)

# Set working directory
setwd("/data2/han_lab/carriehe/DEXSeq/DICE/b_vs_cd8")

# Flattened annotation file
flattenedFile <- "/data2/han_lab/carriehe/gencode.v34.dexseq.bygene.gff"

# Define sample names
sampleNames <- c("B_1", "B_2", "B_3", "B_4", "CD8_1", "CD8_2", "CD8_3", "CD8_4")


# Define conditions (B cells vs CD8)
conditions <- c("B", "B", "B", "B", "CD8", "CD8", "CD8", "CD8")


# Create sample table
sampleData <- data.frame(
  row.names = sampleNames,
  condition = factor(conditions),
  sample = factor(sampleNames)  # Add sample column for design
)


# Get count files from both B and CD8 folders
countFiles_B <- list.files("/data2/han_lab/carriehe/DEXSeq/DICE/b_vs_cd8/count_files/B_new",
                           pattern = "_counts.txt$",
                           full.names = TRUE)

countFiles_CD8 <- list.files("/data2/han_lab/carriehe/DEXSeq/DICE/b_vs_cd8/count_files/CD8_new",
                             pattern = "_counts.txt$",
                             full.names = TRUE)

# Combine count file paths in the correct sample order
countFiles <- c(countFiles_B, countFiles_CD8)


# Create DEXSeq dataset
dxd <- DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData = sampleData,
  design = ~ sample + exon + condition:exon,
  flattenedfile = flattenedFile
)

dxd = estimateSizeFactors(dxd)
BPPARAM = MulticoreParam(workers=15)

#Estimate dispersions
dxd = estimateDispersions(dxd, BPPARAM=BPPARAM)

print("Done running estimate dispersion")

saveRDS(dxd, file = "dexseq_part_1_new.rds")

#nohup Rscript dexseq_dice_b_cd8_1_new.R > script_1_new.log 2>&1 &

