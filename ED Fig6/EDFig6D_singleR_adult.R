## Raw human adult data are downloaded from Hodge, R. D. et al.  doi:10.1038/s41586-019-1506-7 (2019).

library(reticulate)
# use_python("/usr/local/bin/python3.10")
library(Seurat)
# library(edgeR)
library(anndata)
library(scuttle)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(celldex)
library(scRNAseq)
library(splitstackshape)
library(SingleR)

pd <- import("pandas")
sc_expression <- pd$read_pickle("data_allenadult/matrix.pickle")
sc_meta <- read.csv(file = 'data_allenadult/metadata.csv', header = TRUE)
rownames(sc_meta) <- sc_meta$sample_name
sc_meta <- sc_meta[, 2:ncol(sc_meta)]
sc_expression <- sc_expression[, sc_meta$cluster_label != ""]
sc_meta <- sc_meta[sc_meta$cluster_label != "", ]


ref_data <- SingleCellExperiment(assays = list(counts = as.matrix(sc_expression)), rowData = rownames(sc_expression), 
                                 colData = DataFrame(Type=sc_meta$cluster_label, row.names = colnames(sc_expression)))  
ref_data <- logNormCounts(ref_data)

adata_tot <- read_h5ad("../source_data/merscope_integrated_855.h5ad")
adata_tot <- AnnData(adata_tot$raw$X, obs = adata_tot$obs, var = adata_tot$var, obsm = list(spatial = adata_tot$obsm$spatial))
set.seed(1234) 
ind_sample <- stratified(adata_tot$obs, group = "H2_annotation", size = 0.1, keep.rownames = T)
gw_rn <- ind_sample$rn
adata <- adata_tot[rownames(adata_tot$obs) %in% gw_rn]


merfish_data <- SingleCellExperiment(assays = list(counts = t(adata$X)), colData = rownames(adata$obs), 
                                     rowData = rownames(adata$var))
merfish_data <- logNormCounts(merfish_data)


pred_clust <- SingleR(test = merfish_data, ref = ref_data, labels = sc_meta$cluster_label)
write.csv(pred_clust, "allenadult_singler.csv")

