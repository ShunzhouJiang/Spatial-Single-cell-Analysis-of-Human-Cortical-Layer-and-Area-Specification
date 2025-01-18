#Author : Monica D Manam
#         Walsh Lab-BCH
#         2024
#Package Credits : Seurat - Satija Lab
#Visium tutorial : https://satijalab.org/seurat/articles/spatial_vignette ####

#Import Libs
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

#### Fig_Ext_16 - A1 Analysis###

#Read from Zenodo upload 
Visium_A1_brain_011124 <- readRDS("Visium_A1_brain_011124.rds")

#Fig_Ext_16-b
levels(x = Visium_A1_brain_011124) <- c("VZ","iSVZ","iSVZ-L","oSVZ-1","oSVZ-L1",
                                        "oSVZ-L2","oSVZ-3","oSVZ-2","IZ","SP-V1"
                                        ,"SP-V2","L5/6-V1","L5/6-V2","L4-V1",
                                        "L4-V2","L3-V1","L2/3-V2","L2-V1")

#Filtering for markers grouped by clusters with an average log2 fold change greater than 1.
Visium_A1_brain_011124.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(Visium_A1_brain_011124, features = top10$gene) + NoLegend()

table(Visium_A1_brain_011124@active.ident, Visium_A1_brain_011124@meta.data)
# Store cluster identities in object@meta.data$my.clusters
Visium_A1_brain_011124[["my.clusters"]] <- Idents(Visium_A1_brain_011124)


#Fig_Ext_16-a
features_list <- c("NPY", "IGFBPL1", "TSHZ3","GABRA5","NEFM",
                   "ETV1","UNC5C","PDE1A")
Fig_Ext_16-a <- SpatialFeaturePlot(Visium_A1_brain_011124, features = features_list ,  
                            ncol = 4, slot = "data", combine = TRUE,
                            image.alpha = 0, stroke = 0.25, 
                            interactive = FALSE, information = NULL)

#Fig_Ext_16_c
topDEGs <- read.csv("topDEGs.csv")
topDEGs_noSP_5genes <- read.csv("topDEGsnoSP_5genes.csv")
topGenes <- topDEGs$Genes
topDEGs_5genes <-topDEGs_noSP_5genes$Genes
subset_cells <- c('L3-V1',	'L4-V1',	'L5/6-V1',	'SP-V1',	'L2/3-V2',	'L4-V2',	'L5/6-V2',	'SP-V2')
subset_cells_5 <- c('L3-V1',	'L4-V1',	'L5/6-V1',	'L2/3-V2',	'L4-V2',	'L5/6-V2')


Visium_A1_brain_011124_subset <- subset(Visium_A1_brain_011124, idents= c('L3-V1',	'L4-V1',	'L5/6-V1',	'SP-V1',	'L2/3-V2',	'L4-V2',	'L5/6-V2',	'SP-V2'))
Visium_A1_brain_011124_subset2 <- subset(Visium_A1_brain_011124, idents= c('L3-V1',	'L4-V1',	'L5/6-V1','L2/3-V2',	'L4-V2',	'L5/6-V2'))

levels(x = Visium_A1_brain_011124_subset2) <- c('L3-V1',	'L4-V1',	'L5/6-V1','L2/3-V2',	'L4-V2'	,'L5/6-V2')
DoHeatmap(Visium_A1_brain_011124_subset2, features = topDEGs_5genes)

