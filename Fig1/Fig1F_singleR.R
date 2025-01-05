library(reticulate)
# use_python("/usr/local/bin/python3.10")
library(Seurat)
library(anndata)
library(scuttle)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(celldex)
library(scRNAseq)
library(splitstackshape)
library(SingleR)

# SMARTseq data are downloaded from Nowakowski, T. J. et al. doi:10.1126/science.aap8809 (2017).
sc_expression <- read.table(file = 'SMARTseq_exprMatrix.tsv', sep = '\t', header = TRUE)
rownames(sc_expression) <- sc_expression$gene
sc_expression <- sc_expression[, 2:ncol(sc_expression)]
rownames(sc_expression) <- substring(rownames(sc_expression), 1, unlist(gregexpr('[|]', rownames(sc_expression)))-1 )

sc_meta <- read.table(file = 'SMARTseq_meta.tsv', sep = '\t', header = TRUE)
sc_expression <- sc_expression[, is.na(sc_meta$WGCNAcluster) == FALSE]
sc_meta <- sc_meta[is.na(sc_meta$WGCNAcluster) == FALSE, ]
sc_expression <- sc_expression[, sc_meta$WGCNAcluster != ""]
sc_meta <- sc_meta[sc_meta$WGCNAcluster != "", ]

cluster_exclude <- c("U1", "U2", "U3", "U4", "MGE-div","MGE-IPC1","MGE-IPC2","MGE-IPC3","MGE-RG1","MGE-RG2")
sc_expression <- sc_expression[, !sc_meta$WGCNAcluster %in% cluster_exclude]
sc_meta <- sc_meta[!sc_meta$WGCNAcluster %in% cluster_exclude, ]
cluster_in <- c("nIN1","nIN2","nIN3","nIN4", "nIN5")
sc_meta$WGCNAcluster[sc_meta$WGCNAcluster %in% cluster_in] <- "nIN"
ref_data <- SingleCellExperiment(assays = list(counts = as.matrix(sc_expression)), rowData = rownames(sc_expression), 
                                colData = DataFrame(Type=sc_meta$WGCNAcluster, row.names = colnames(sc_expression)))  
ref_data <- logNormCounts(ref_data)


adata_tot <- read_h5ad("merscope_integrated_855.h5ad")
adata_tot <- AnnData(adata_tot$raw$X, obs = adata_tot$obs, var = adata_tot$var, obsm = list(spatial = adata_tot$obsm$spatial))
set.seed(1234)
ind_sample <- stratified(adata_tot$obs, group = "H2_annotation", size = 0.1, keep.rownames = T)
gw_rn <- ind_sample$rn
adata <- adata_tot[rownames(adata_tot$obs) %in% gw_rn]
merfish_data <- SingleCellExperiment(assays = list(counts = t(adata$X)), colData = rownames(adata$obs), 
                                     rowData = rownames(adata$var))
merfish_data <- logNormCounts(merfish_data)
pred_clust <- SingleR(test = merfish_data, ref = ref_data, labels = sc_meta$WGCNAcluster)
write.csv(pred_clust, "smartseq_singler_h2.csv")




