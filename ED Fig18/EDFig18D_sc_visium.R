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



vis <- readRDS("Visium_A1brain_no14_annotated.rds")
sel_cluster <- c("IZ", "SP-V1", "SP-V2", "L5/6-V1", "L5/6-V2", "L4-V1", "L4-V2", "L3-V1", "L2/3-V2", "L2-V1")
vis <- vis[ , Idents(vis) %in% sel_cluster]
merfish_data <- SingleCellExperiment(assays = list(counts = GetAssayData(vis@assays$Spatial)), 
                                     colData = DataFrame(Type=Idents(vis), row.names = colnames(vis@assays$Spatial)), 
                                     rowData = rownames(vis@assays$Spatial))
merfish_data <- logNormCounts(merfish_data)
pred_clust <- SingleR(test = merfish_data, ref = ref_data, labels = Idents(sc_sub) )
pred_clust <- data.frame(pred_clust)
pred_clust <- pred_clust[, (ncol(pred_clust)-1):ncol(pred_clust) ]
pred_clust$raw <- Idents(vis)
write.csv(pred_clust, "result/sc_visium.csv")


