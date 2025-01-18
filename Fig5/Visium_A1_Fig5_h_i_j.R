#Fig5_i-j.R
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

#sctransform normalizes the data, detects high-variance features, 
#and stores the data in the SCT assay.
A1_brain <- readRDS('/Users/monicam/RWorkspace/Seurat_analysis/VISIUM/A1_brain.rds')
A1_brain <- SCTransform(A1_brain, assay = "Spatial", verbose = FALSE)
A1_brain <- RunPCA(A1_brain, assay = "SCT", verbose = FALSE)
A1_brain <- FindNeighbors(A1_brain, reduction = "pca", dims = 1:30)
#Increased to 1.5 to get more clusters 
A1_brain <- FindClusters(A1_brain, resolution = 1.5)
A1_brain <- RunUMAP(A1_brain, reduction = "pca", dims = 1:30)
plot <- VlnPlot(A1_brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

# QC 1: Remove cluster 14 looking at downstream violin counts - has 0 cells 
A1brain_no14 <- subset(x = A1_brain, idents = 14, invert = TRUE)

A1_SpatialHeatMap <- SpatialFeaturePlot(A1brain_no14, 
                                        max.cutoff = "1e+05",
                                        features = "nCount_Spatial",
                                        image.alpha = 0,
                                        pt.size.factor = 2.5) + theme(legend.position = "right")


#2.Spatial_slide
#Fig5_j
p1 <- DimPlot(Visium_A1_brain_011124, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(Visium_A1_brain_011124, label = TRUE, label.size = 3)
A1_slide <- p1 + p2

#3.
#SpatilDimPlot
A1_spatialdimplot <- SpatialDimPlot(
  A1brain_no14, 
  label = TRUE, 
  label.size = 6, 
  image.alpha = 0,
  pt.size.factor = 2.0,
  stroke = NA
)
#multiple
SpatialDimPlot(A1_brain, cells.highlight = CellsByIdentities(A1_brain),
               facet.highlight = TRUE,
               ncol = 5)

#Annotation : 
A1brain_no14_annotated <- RenameIdents(A1brain_no14, `0` = "iSVZ", `1` = "IZ", 
                                       `2` = "oSVZ-1", `3` = "VZ", `4` = "oSVZ-2",
                                       `5` = "SP-V1", `6` = "oSVZ-L1", 
                                       `7` = "oSVZ-L2", `8` = "L4-V2", 
                                       `9` = "L5/6-V1", `10` = "L4-V1", 
                                       `11` = "oSVZ-3", `12` = "L2/3-V2", 
                                       `13` = "SP-V2", `15` = "L5/6-V2", 
                                       `16` = "L2-V1", `17` = "L3-V1", 
                                       `18` = "iSVZ-L")

fig_5_i <- DimPlot(A1brain_no14_annotated, label = TRUE)

#MARKER GENE ANALYSIS 
table(A1brain_no14@active.ident)
A1brain_no14_annotated.markers <- FindAllMarkers(A1brain_no14_annotated, only.pos = TRUE)

#Filtering for markers grouped by clusters with an average log2 fold change greater than 1.
A1brain_no14_annotated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(A1brain_no14_annotated, features = top10$gene) + NoLegend()

table(A1brain_no14_annotated@active.ident, A1brain_no14_annotated@meta.data)
# Store cluster identities in object@meta.data$my.clusters
A1brain_no14_annotated[["my.clusters"]] <- Idents(A1brain_no14_annotated)

#CALL the function to flip the image
bb5=rotateSeuratImage(A1brain_no14_annotated,rotation = "180")

#Plot Spatial 
fig_5_h <- SpatialFeaturePlot(object = bb5, 
                              features = c("AB13BP", "PDZRN4", "TAFA2", 
                                           "CCBE1", "THSD7B", "IL1RAP",
                                           "FLRT2"),
                              ncol = 3,
                              stroke = NA,
                              max.cutoff = "1e+05",
                              image.alpha = 0,
                              pt.size.factor = 2.5)


#Method to flip Image
# Code to Flipping Images, if the slide is flipped.
# see https://github.com/satijalab/seurat/issues/2702 for more.

# flip_angle %in% c(180, "R90", "L90", "Hf", "Vf")

rotimat=function(foo,rotation){
  if(!is.matrix(foo)){
    cat("Input is not a matrix")
    return(foo)
  }
  if(!(rotation %in% c("180","Hf","Vf", "R90", "L90"))){
    cat("Rotation should be either L90, R90, 180, Hf or Vf\n")
    return(foo)
  }
  if(rotation == "180"){
    foo <- foo %>% 
      .[, dim(.)[2]:1] %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "Hf"){
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  
  if(rotation == "Vf"){
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "L90"){
    foo = t(foo)
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "R90"){
    foo = t(foo)
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  return(foo)
}

rotateSeuratImage = function(seuratVisumObject, slide = "slice1", rotation="Vf"){
  if(!(rotation %in% c("180","Hf","Vf", "L90", "R90"))){
    cat("Rotation should be either 180, L90, R90, Hf or Vf\n")
    return(NULL)
  }else{
    seurat.visium = seuratVisumObject
    ori.array = (seurat.visium@images)[[slide]]@image
    img.dim = dim(ori.array)[1:2]/(seurat.visium@images)[[slide]]@scale.factors$lowres
    new.mx <- c()  
    # transform the image array
    for (rgb_idx in 1:3){
      each.mx <- ori.array[,,rgb_idx]
      each.mx.trans <- rotimat(each.mx, rotation)
      new.mx <- c(new.mx, list(each.mx.trans))
    }
    
    # construct new rgb image array
    new.X.dim <- dim(each.mx.trans)[1]
    new.Y.dim <- dim(each.mx.trans)[2]
    new.array <- array(c(new.mx[[1]],
                         new.mx[[2]],
                         new.mx[[3]]), 
                       dim = c(new.X.dim, new.Y.dim, 3))
    
    #swap old image with new image
    seurat.visium@images[[slide]]@image <- new.array
    
    ## step4: change the tissue pixel-spot index
    img.index <- (seurat.visium@images)[[slide]]@coordinates
    
    #swap index
    if(rotation == "Hf"){
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "Vf"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
    }
    
    if(rotation == "180"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "L90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[2]-img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.index$imagerow
    }
    
    if(rotation == "R90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[1]-img.index$imagerow
    }
    
    return(seurat.visium)
  }  
}





