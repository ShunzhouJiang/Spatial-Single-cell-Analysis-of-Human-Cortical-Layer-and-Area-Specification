# Libraries
library(Seurat)
library(CellChat) # CellChat V1
library(patchwork)
library(ComplexHeatmap)
library(reticulate)

##### 0. Universal variables
# Species
species <- "human"

# Cell types of interest
celltypes <- c("EN-IT-L4-V1", "EN-IT-UL-1", "EN-IT-UL-2", "EN-IT-UL-3", "EN-IT-UL-4", "EN-IT-UL-5")

# Pathway of interest
pathway <- "NRXN"

##### 1. Preprocess the data set
# Load the data set
seuratobj <- readRDS("Final_Integration_MDM_100323.rds")

# Create cell type meta data
seuratobj$celltypes <- Idents(seuratobj)
seuratobj <- seuratobj[, seuratobj$orig.ident %in% c("FB080_I", "FB080_J")]

# Merge the clusters into these categories: 
celltypes <- as.character(seuratobj$celltypes)
# 1. EN-ET (-1 to -7)
celltypes[celltypes %in% paste0("EN-ET-", seq(1,7))] <- "EN-ET"
# 2. EN-IT-DL (-1, -2, -3)
celltypes[celltypes %in% paste0("EN-IT-DL-", seq(1,3))] <- "EN-IT-DL"
# 3. EN-Mig (-1, -2, -3)
celltypes[celltypes %in% paste0("EN-Mig-", seq(1,3))] <- "EN-Mig"
# 4. IN-CGE (IN-CGE-1, IN-CGE-2, INP-CGE)
celltypes[celltypes %in% c("IN-CGE-1", "IN-CGE-2", "INP-CGE")] <- "IN-CGE"
# 5. IN-MGE (-1, -2, -3)
celltypes[celltypes %in% paste0("IN-MGE-", seq(1,3))] <- "IN-MGE"

# Update and sort the Idents
seuratobj$celltypes <- celltypes
Idents(seuratobj) <- seuratobj$celltypes
levels(seuratobj) <- c(sort(setdiff(levels(seuratobj), "unknown")), "unknown")

##### 2. Do cellchat analysis on entire object 
ccMetaData <- data.frame(label = Idents(seuratobj), samples = as.factor(seuratobj$orig.ident))
ccX <- createCellChat(seuratobj@assays$RNA@data, meta = ccMetaData, group.by = 'label')
ccDB <- CellChatDB.human
ccX@DB <- ccDB
ccX <- subsetData(ccX)
ccX <- identifyOverExpressedGenes(ccX)
ccX <- identifyOverExpressedInteractions(ccX)
ccX <- computeCommunProb(ccX)
ccX <- filterCommunication(ccX, min.cells = 10)
ccX <- computeCommunProbPathway(ccX)
ccX <- aggregateNet(ccX)
ccX <- netAnalysis_computeCentrality(ccX)
cellchat_object <- ccX

##### 3. CellChat result analysis and visualization
celltypes_all <- levels(seuratobj)

# Signals contributing most to outgoing or incoming signaling of cell types
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_object, pattern = "outgoing", height = 9.5, width = 7, font.size=5, font.size.title = 12)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_object, pattern = "incoming", height = 9.5, width = 7, font.size=5, font.size.title = 12)
png("Figures/ALL_outgoing_incoming_signal.png", height = 2000, width = 2500, res = 300)
print(ht1 + ht2)
dev.off()
pdf("Figures/ALL_outgoing_incoming_signal.pdf", height = 6, width = 8)
print(ht1 + ht2)
dev.off()

# Heatmap
p <- netVisual_heatmap(cellchat_object, signaling = pathway, color.heatmap = "Reds", font.size.title = 16)
png(paste0("Figures/", pathway, "_net_heatmap.png"), res = 300, height = 1750, width = 1750)
print(p)
dev.off()
pdf(paste0("Figures/", pathway, "_net_heatmap.pdf"), height = 6, width = 6)
print(p)
dev.off()

# Bubble plot of LR for the given pathway
idx_cc <- which(celltypes_all %in% celltypes)
p <- netVisual_bubble(cellchat_object, signaling = pathway, sources.use = seq(1, length(celltypes_all)), targets.use = idx_cc, remove.isolate = TRUE, font.size = 7)
png(paste0("Figures/Bubbleplot_", pathway, ".png"), res = 300, height = 1150, width = 3250)
print(p)
dev.off()
pdf(paste0("Figures/Bubbleplot_", pathway, ".pdf"), height = 3.75, width = 10)
print(p)
dev.off()




