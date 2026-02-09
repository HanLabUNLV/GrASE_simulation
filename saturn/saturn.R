library(DEXSeq)
library(tidyverse)
library(BiocParallel)
library(optparse)
library(satuRn)


GFFfile = "../ref/gencode.v28.dexseq.bygene.gff"
countdir = "../DEXSeq/count_files"
outdir = "saturn_group1_group2"
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

# create SummarizedExperiment
exonInfo <- rowData(dxd)
colnames(exonInfo)[1:2] <- c("isoform_id", "gene_id")
exonInfo$isoform_id <- rownames(exonInfo)
sumExp <- SummarizedExperiment(assays = list(counts=featureCounts(dxd)),
                               colData = sampleAnnotation(dxd),
                               rowData = exonInfo)

sumExp0 <- sumExp
print(dim(sumExp))
# remove exons with zero expression
#keep <- rowSums(assay(sumExp, "counts")) != 0

# at least 10 reads in at least 3 samples
keep <- rowSums(assay(sumExp, "counts") >= 10) >= 3
sumExp <- sumExp[keep, ]
print(dim(sumExp))

# countsum > 0 in both groups
counts_mat <- assay(sumExp, "counts")
is_celltype1 <- colData(sumExp)$condition == celltype1 
is_celltype2 <- colData(sumExp)$condition == celltype2
safe_rows <- rowSums(counts_mat[, is_celltype1]) > 0 & rowSums(counts_mat[, is_celltype2]) > 0
sumExp <- sumExp[safe_rows, ]
sumExp_filteredbyCount <- sumExp
print(dim(sumExp))

## proportion > 0.05 in at least 5 samples 
#gene_ids <- rowData(sumExp)$gene_id
#counts_mat <- assay(sumExp, "counts")
#gene_sums <- rowsum(counts_mat, gene_ids)
#exon_gene_totals <- gene_sums[match(gene_ids, rownames(gene_sums)), ]
#usage_prop <- counts_mat / exon_gene_totals
#keep_prop <- rowSums(usage_prop >= 0.05, na.rm = TRUE) >= 5
#keep_prop[is.na(keep_prop)] <- FALSE
#sumExp <- sumExp[keep_prop, ]
#sumExp3 <- sumExp
#print(dim(sumExp))

# filter out genes with only one exon
keep_groups <- names(which(table(rowData(sumExp)$gene_id) > 1))
sumExp <- sumExp[rowData(sumExp)$gene_id %in% keep_groups, ]
sumExp_filteredbyCountMultiExon <- sumExp
print(dim(sumExp))

# Identify rows where variance is practically zero
#counts_mat <- assay(sumExp, "counts")
#row_vars <- apply(counts_mat, 1, var)
#problem_rows <- which(row_vars < 0.01 | is.na(row_vars))
#if(length(problem_rows) > 0) {
#    message("Removing ", length(problem_rows), " numerically unstable rows.")
#    sumExp <- sumExp[-problem_rows, ]
#}
#print(dim(sumExp))

Explist = list(sumExp_filteredbyCount=sumExp_filteredbyCount, sumExp_filteredbyCountMultiExon = sumExp_filteredbyCountMultiExon)
for (expnum in 1:(length(Explist)) ) {
sumexpfile = paste0(outdir,"/sumexp.",celltype1,"_",celltype2,".",names(Explist)[expnum],".rds")

# logit(p)=βcondI(cond1)+βcondI(cond2) no intercept.
# One coefficient per conditiona
# Each coefficient is a group-specific log-odds

if (file.exists(sumexpfile)) {
  sumExp <- readRDS(sumexpfile)
} else {
  sumExp <- satuRn::fitDTU(object = Explist[[expnum]],
                         formula = ~ 0 + condition,
                         parallel = FALSE,
                         BPPARAM = MulticoreParam(workers=15),
                         verbose = TRUE)


  saveRDS(sumExp, sumexpfile)
}

print(paste0(names(Explist)[expnum], ":"))

# This finds the indices where the model is NOT an error
valid_models <- sapply(rowData(sumExp)$fitDTUModels, function(x) {
    # Check if the object is a valid StatModel and has coefficients
    is(x, "StatModel") && !is.null(x@params$coefficients)
})

# 3. Subset the object to ONLY include the successful fits
message("Successful fits: ", sum(valid_models))
message("Failed fits (including those 66 mu errors): ", sum(!valid_models))

# check the fit
table(sapply(rowData(sumExp)$fitDTUModels, function(x) x@type))

#  coefficient names 
#print(head(rowData(sumExp)$fitDTUModels)) 

# Re-create the contrast matrix to match the formula ~ 0 + condition
# If the levels are "CD4" and "CD4_N_STIM", fitDTU likely named them 
# "conditionCD4" and "conditionCD4_N_STIM"

# Try defining the contrast names based on the design matrix used in fitDTU
group_names <- levels(colData(sumExp)$condition)
# Assuming the fit created names like "conditionCD4"
coef_names <- paste0("condition", group_names) 

L <- matrix(0, ncol = 1, nrow = length(coef_names))
rownames(L) <- coef_names
colnames(L) <- "C1"
L

# Adjust these indices based on which celltype is 1 and 2
L[paste0("condition", celltype1), 1] <- 1
L[paste0("condition", celltype2), 1] <- -1
L

sumExp <- satuRn::testDTU(object = sumExp, 
                          contrasts = L, 
                          diagplot1 = TRUE, 
                          diagplot2 = TRUE, 
                          sort = FALSE, 
                          forceEmpirical = TRUE)



# differential exons
count <- counts(dxd)
exondata <- rowData(sumExp)
ids = rownames(exondata)
fitcount <- counts(dxd)[ids,] 
group1_names = rownames(sampleData)[sampleData$condition==celltype1]
group2_names = rownames(sampleData)[sampleData$condition==celltype2] 
group1_cnts = fitcount[,as.character(colData(dxd)$sample)  %in% group1_names]
group2_cnts = fitcount[,as.character(colData(dxd)$sample)  %in% group2_names]
group1_exon_mean = rowMeans(group1_cnts[,1:length(group1_names)])
group1_total_mean = rowMeans(group1_cnts[,(1+length(group1_names)):(length(group1_names)+length(group1_names))])
group2_exon_mean = rowMeans(group2_cnts[,1:length(group2_names)])
group2_total_mean = rowMeans(group2_cnts[,(1+length(group2_names)):(length(group2_names)+length(group2_names))])
stat <- exondata[,c('exonBaseMean', 'exonBaseVar')]
models <- rowData(sumExp)[,"fitDTUModels"]
coef_list <- lapply(ids, function(id) {
  m <- models[[id]]
  if (is.null(m)) return(NULL)
  m@params$coefficients
})
names(coef_list) <-ids 
coef_list <- Filter(Negate(is.null), coef_list)
# union of all coefficient names
all_coef_names <- sort(unique(unlist(lapply(coef_list, names))))
# build matrix (genes x coefficients) Log-odds change in isoform usage relative to the reference
coef_mat <- do.call(rbind, lapply(coef_list, function(co) {
  out <- setNames(rep(NA_real_, length(all_coef_names)), all_coef_names)
  out[names(co)] <- unname(co)
  out
}))

fitresult <- exondata[,'fitDTUResult_C1']   # fit results estimate βCD-βCD4_N_STI
OR = exp(fitresult$estimates)
results <- cbind.data.frame(ids, group1_exon_mean, group1_total_mean, group2_exon_mean, group2_total_mean, stat,  coef_mat, fitresult, OR)
rownames(results) = ids
totalfile = paste0(outdir,"/all.",celltype1,"_",celltype2,".",names(Explist)[expnum],".txt")
write.table(results, totalfile, sep="\t", quote=FALSE, row.names=FALSE)
# filter DEU at FDR < 0.05
DEU <- rownames(exondata[["fitDTUResult_C1"]][which(exondata[["fitDTUResult_C1"]]$empirical_FDR < 0.05),])
DEUfile = paste0(outdir,"/DEU.",celltype1,"_",celltype2,".",names(Explist)[expnum],".txt")
write.table(results[DEU,], DEUfile, sep="\t", quote=FALSE, row.names=FALSE)

}
