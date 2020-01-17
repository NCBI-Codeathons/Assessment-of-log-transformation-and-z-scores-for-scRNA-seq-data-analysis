library(Matrix)
library(SingleCellExperiment)
library(scran)
library(scater)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_prefix <- args[2]

# Read the data
raw_data <- readMM(input_file)

#Give data generic row and column names.
rownames(raw_data) <- paste0("Feature",1:nrow(raw_data))
colnames(raw_data) <- paste0("Cell",1:ncol(raw_data))

#Create single-cell experiment.
sce <- SingleCellExperiment(assays = list(counts = raw_data))

# cluster
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- logNormCounts(sce)

# normalize data
normalized_data <- assays(sce)[["logcounts"]]

print("Data has been normalized, now saving as a feather file")

#Save as mtx file
#write_feather(x=as.data.frame(normalized_data),path=output_prefix)
writeMM(normalized_data,file=paste0(output_prefix, "_scran_normalized.mtx"))
