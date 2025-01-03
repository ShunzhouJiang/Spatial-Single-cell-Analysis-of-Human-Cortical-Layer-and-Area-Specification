library(Seurat)
library(scuttle)
library(scRNAseq)
library(reticulate)
library(anndata)
library(SingleR)
library(splitstackshape)
use_python("/usr/local/bin/python3.10")

sc <- readRDS("SCRNASeq_All_Integrated_annotated.rds")
h2_cluster <- levels(sc)[grepl("EN", levels(sc))]
sc <- sc[ , Idents(sc) %in% h2_cluster]
set.seed(1234)
ind_sample <- stratified(sc@meta.data, group = "integrated_snn_res.1", size = 0.5, keep.rownames = T)
gw_rn <- ind_sample$rn
sc_sub <- sc[, rownames(sc@meta.data) %in% gw_rn]
ref_data <- SingleCellExperiment(assays = list(counts = GetAssayData(sc_sub)), rowData = rownames(sc_sub), 
                                 colData = DataFrame(Type=sc_sub@meta.data$integrated_snn_res.1, row.names = colnames(sc_sub)))  
ref_data <- logNormCounts(ref_data)


adata_tot <- read_h5ad("gw20.h5ad")
adata_tot <- adata_tot[adata_tot$obs$H1_annotation %in% c("EN-ET", "EN-IT", "EN-Mig")]
adata_tot$obs$source <- paste(adata_tot$obs$sample, adata_tot$obs$region, sep = "-")
adata_tot <- adata_tot[adata_tot$obs$source %in% c("FB080-O1c", "FB121-F1") ]
adata_tot$obs$H2_annotation <- as.character(adata_tot$obs$H2_annotation)
adata_tot$obs$H3_annotation <- as.character(adata_tot$obs$H3_annotation)
adata_tot$obs$H3_annotation <- ifelse(adata_tot$obs$H3_annotation == "NA", 
                                      ifelse(adata_tot$obs$H1_annotation %in% c("EN-ET", "EN-IT"), 
                                             paste0(adata_tot$obs$H2_annotation, "-c0"), 
                                             adata_tot$obs$H2_annotation), adata_tot$obs$H3_annotation)
cluster_keep <- names( which(table(adata_tot$obs$H3_annotation) > 1000) )
adata_tot <- adata_tot[adata_tot$obs$H3_annotation %in% cluster_keep]
adata <- adata_tot
merfish_data <- SingleCellExperiment(assays = list(counts = t(adata$raw$X)), colData = adata$obs, 
                                     rowData = rownames(adata$var))
rownames(merfish_data) <- rownames(adata$var)
merfish_data <- logNormCounts(merfish_data)

pred_clust <- SingleR(test = merfish_data, ref = ref_data, labels = Idents(sc_sub) )
pred_clust <- data.frame(pred_clust)
pred_clust <- pred_clust[, (ncol(pred_clust)-1):ncol(pred_clust) ]
pred_clust$raw <- adata$obs$H3_annotation
write.csv(pred_clust, "result/sc_merscope.csv")
