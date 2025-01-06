## Raw scRNA-seq reference data are downloaded from Bhaduri, A. et al. doi:10.1038/s41586-021-03910-8 (2021).

library(reticulate)
# use_python("/usr/local/bin/python3.10")
library(Seurat)
library(Matrix)
library(data.table)
library(anndata)
library(scuttle)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(celldex)
library(scRNAseq)
library(splitstackshape)
library(SingleR)
library(purrr)
library(dplyr)

gz_file <- gzcon(file("data/kriegstein/matrix.mtx.gz", "rb"))
matrix <- readMM(gz_file)
close(gz_file)

feature <- fread("data/kriegstein/features.tsv.gz", sep = "\t", header = F)
# cell <- fread("data/kriegstein/barcodes.tsv.gz", sep = "\t", header = F)
meta <- read.csv("data/kriegstein/meta.tsv", sep = "\t")
colnames(meta)[colnames(meta) == "CombinedCluster...Final"] <- "cell_type"
matrix <- matrix[, !meta$cell_type %in% c(0, "Other_combo2_14")]
meta <- meta[!meta$cell_type %in% c(0, "Other_combo2_14"), ]
meta[meta$ConsensusCellType...Final == "OPC", ]$cell_type <- "OPC"
meta[meta$ConsensusCellType...Final == "Microglia", ]$cell_type <- "Microglia"
meta[meta$ConsensusCellType...Final == "Vascular", ]$cell_type <- "Vascular"
meta[meta$ConsensusCellType...Final == "CR", ]$cell_type <- "CR"
meta[meta$ConsensusCellType...Final == "Other" & meta$CellType == "IPC", ]$cell_type <- "IPC_other"
meta[meta$ConsensusCellType...Final == "Other" & meta$CellType == "Dividing", ]$cell_type <- "Dividing_other"
meta[meta$ConsensusCellType...Final == "Other" & meta$CellType == "RG", ]$cell_type <- "RG_other"
meta[meta$ConsensusCellType...Final == "Other" & meta$CellType == "OPC", ]$cell_type <- "OPC"
meta[meta$cell_type == "Other_combo2_17", ]$cell_type <- "Microglia"

rownames(matrix) <- feature$V1
colnames(matrix) <- meta$Cell.Name
print(dim(matrix))

ref_data <- SingleCellExperiment(assays = list(counts = matrix), rowData = feature$V2, 
                                colData = DataFrame(Type=meta$cell_type, row.names = meta$Cell.Name))  
ref_data <- logNormCounts(ref_data)

adata_tot <- read_h5ad("../source_data/merscope_integrated_855.h5ad")
adata_tot <- adata_tot[adata_tot$obs$H1_annotation %in% c("EN-IT", "EN-ET", "EN-Mig")]
adata_tot <- AnnData(adata_tot$raw$X, obs = adata_tot$obs, var = adata_tot$var, obsm = list(spatial = adata_tot$obsm$spatial))
set.seed(1234)
ind_sample <- stratified(adata_tot$obs, group = "H2_annotation", size = 0.005, keep.rownames = T)
gw_rn <- ind_sample$rn
adata <- adata_tot[rownames(adata_tot$obs) %in% gw_rn]
print(dim(adata))
merfish_data <- SingleCellExperiment(assays = list(counts = t(adata$X)), colData = rownames(adata$obs), 
                                     rowData = rownames(adata$var))
merfish_data <- logNormCounts(merfish_data)
pred_clust <- SingleR(test = merfish_data, ref = ref_data, labels = meta$cell_type)
pred_clust <- data.frame(pred_clust)
pred_clust <- pred_clust[, (ncol(pred_clust)-1):ncol(pred_clust) ]
pred_clust$H2 <- adata$obs$H2_annotation
pred_clust$H3 <- adata$obs$H3_annotation

write.csv(pred_clust, "result/kriegstein_merscope_neuron.csv")



