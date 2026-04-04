library(DEXSeq)
library(tidyverse)
library(BiocParallel)
library(optparse)


GFFfile = "../ref/gencode.v28.dexseq.bygene.gff"
countdir = "../DEXSeq/count_files"
outdir = "dexseq_group1_group2"
celltype1 = "group1"
celltype2 = "group2"


option_list = list(
  make_option(c("--gff"), type="character",  
              help="dexseq gff file", metavar="character"),
  make_option(c("--cntdir"), type="character",  
              help="count_files directory", metavar="character"),
  make_option(c("--outdir"), type="character",  
              help="output directory", metavar="character"),
  make_option(c("--cell1"), type="character",  
              help="group 1", metavar="character"),
  make_option(c("--cell2"), type="character", 
              help="group 2", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

## Check if all required options are provided
#required_opts <- c("gff", "cntdir", "outdir", "cell1", "cell2")
#missing_opts <- required_opts[sapply(required_opts, function(x) is.null(opt[[x]]))]
#if (length(missing_opts) > 0) {
#  cat("Error: Missing required options:", paste(missing_opts, collapse=", "), "\n\n")
#  print_help(opt_parser)
#  stop("Please provide all required options")
#}
#
# Assign required options
#GFFfile = opt$gff
#celltype1 = opt$cell1
#celltype2 = opt$cell2
#outdir = opt$outdir
#countdir = opt$cntdir
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Get count files from both celltype1 and celltype2 folders
countFiles1 <- list.files(paste0(countdir,"/",celltype1),
                           pattern = "_counts.txt$",
                           full.names = TRUE)

countFiles2 <- list.files(paste0(countdir,"/",celltype2),
                             pattern = "_counts.txt$",
                             full.names = TRUE)

# samplesIDs
samples1 <- sub("_counts\\.txt$", "", basename(countFiles1))
samples2 <- sub("_counts\\.txt$", "", basename(countFiles2))
sampleNames = c(paste0(celltype1, '_', samples1), paste0(celltype2, '_', samples2))
conditions = c(rep(celltype1,length(samples1)), rep(celltype2, length(samples2)))
 
# Create sample table
sampleData <- data.frame(
  row.names = sampleNames,
  condition = factor(conditions),
  sample_id = factor(sampleNames)  # Add sample column for design
)


# Combine count file paths in the correct sample order
countFiles <- c(countFiles1, countFiles2)
names(countFiles) = sampleNames

dxdfile = paste0(outdir,"/dxd.",celltype1,"_",celltype2,".rds")
if (file.exists(dxdfile)) {
  dxd <- readRDS(dxdfile)
} else {
  # Create DEXSeq dataset
  dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData = sampleData,
    design = ~ sample_id + exon + condition:exon,
    flattenedfile = GFFfile
  )
  saveRDS(dxd, dxdfile)
}
print(dim(dxd))

# remove exons with zero expression
#keep <- rowSums(assay(dxd, "counts")) != 0

# at least 10 reads in at least 3 samples
keep <- rowSums(assay(dxd, "counts") >= 10) >= 3
dxd <- dxd[keep, ]
print(dim(dxd))

# countsum > 0 in both groups
counts_mat <- assay(dxd, "counts")
is_celltype1 <- colData(dxd)$condition == celltype1 
is_celltype2 <- colData(dxd)$condition == celltype2
safe_rows <- rowSums(counts_mat[, is_celltype1]) > 0 & rowSums(counts_mat[, is_celltype2]) > 0
dxd <- dxd[safe_rows, ]
dxd_filteredbyCount <- dxd
print(dim(dxd))

## proportion > 0.05 in at least 5 samples 
#groupIDs <- rowData(dxd)$groupID
#counts_mat <- assay(dxd, "counts")
#gene_sums <- rowsum(counts_mat, groupIDs)
#exon_gene_totals <- gene_sums[match(groupIDs, rownames(gene_sums)), ]
#usage_prop <- counts_mat / exon_gene_totals
#keep_prop <- rowSums(usage_prop >= 0.05, na.rm = TRUE) >= 5
#keep_prop[is.na(keep_prop)] <- FALSE
#dxd <- dxd[keep_prop, ]
#dxd3 <- dxd
#print(dim(dxd))

# filter out genes with only one exon
keep_groups <- names(which(table(rowData(dxd)$groupID) > 1))
dxd <- dxd[rowData(dxd)$groupID %in% keep_groups, ]
dxd_filteredbyCountMultiExon <- dxd
print(dim(dxd))

# Identify rows where variance is practically zero
#counts_mat <- assay(dxd, "counts")
#row_vars <- apply(counts_mat, 1, var)
#problem_rows <- which(row_vars < 0.01 | is.na(row_vars))
#if(length(problem_rows) > 0) {
#    message("Removing ", length(problem_rows), " numerically unstable rows.")
#    dxd <- dxd[-problem_rows, ]
#}
#print(dim(dxd))

Explist = list(dxd_filteredbyCount=dxd_filteredbyCount, dxd_filteredbyCountMultiExon = dxd_filteredbyCountMultiExon)
for (expnum in 1:(length(Explist)) ) {
  dxdfile = paste0(outdir,"/dxd.",celltype1,"_",celltype2,".",names(Explist)[expnum],".rds")
  dxd = Explist[[expnum]]

  # logit(p)=βcondI(cond1)+βcondI(cond2) no intercept.
  # One coefficient per conditiona
  # Each coefficient is a group-specific log-odds

  if (file.exists(dxdfile)) {
    dxd <- readRDS(dxdfile)
  } else {

    dxd = estimateSizeFactors(dxd)
    BPPARAM = MulticoreParam(workers=15)

    #Estimate dispersions
    dxd = estimateDispersions(dxd, BPPARAM=BPPARAM)
    saveRDS(dxd, dxdfile)
  }

  print(paste0(names(Explist)[expnum], ":"))


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

  totalfile = paste0(outdir,"/all.",celltype1,"_",celltype2,".",names(Explist)[expnum],".txt")
  #Save results to text file
  write.table(results_df,
              file = totalfile,
              sep = "\t",
              quote = FALSE,
              row.names = TRUE)

  #exons with significantly different usage 
  DEUfile = paste0(outdir,"/DEU.",celltype1,"_",celltype2,".",names(Explist)[expnum],".txt")
  sig_results = subset(results_df, padj < 0.05)
  write.table(sig_results, file = DEUfile, sep = "\t", quote = FALSE)


}
