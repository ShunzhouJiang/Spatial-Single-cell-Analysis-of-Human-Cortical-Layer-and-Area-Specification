#ED Fig5 g-j
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

#Clustering was performed on BCH Cluster : E2
#RUN IMPORT LIBS 
#Get a compute node
#source /programs/biogrids.shrc
#export R_X=4.1
#R

#Read in Zenodo Upload 
D1_brain_no15 <- readRDS("D1_brain_no15.rds")

#Side by Side plot 
p1 <- DimPlot(D1_brain_no15, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(D1_brain_no15, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2


SpatialFeaturePlot(D1_brain_no15, 
                   max.cutoff = "1e+05",
                   features = "nCount_Spatial",
                   image.alpha = 0,
                   pt.size.factor = 1) + theme(legend.position = "right")

#HeatMap
SpatialFeaturePlot(D1_brain_no15, 
                   max.cutoff = "1e+05",
                   features = "nCount_Spatial",
                   image.alpha = 0,
                   pt.size.factor = 2.5) + theme(legend.position = "right")


#D1
q2 <- SpatialDimPlot(
  D1_brain_no15, 
  label = TRUE, 
  label.size = 5, 
  image.alpha = 0,
  pt.size.factor = 2.0,
  stroke = 0.001
)
q2
ggsave("D1_spatial_032624.pdf", width = 12, height = 6)

#Genes of Interest with 
seurat_plot <- SpatialFeaturePlot(object = D1_brain_no15, 
                                  features = c("CUX2","CCBE1","MDGA1","CDH13",
                                               "TRPC6","CHRM2",
                                               "RORB","WHRN","NELL1","ETV6",
                                               "IL1RAPL2","ETV1","TOX",
                                               "DEPTOR","TSHZ2","QRFPR","FN1",
                                               "THSD7B","CDH18",
                                               "GALNTL6","BRINP3","HTR4",
                                               "PPARG","PLD5","CCN2"),
                                  ncol = 5,
                                  stroke = NA,
                                  max.cutoff = "1e+05", #min and max for each feature
                                  image.alpha = 0,
                                  pt.size.factor = 2.5)

seurat_plot
ggsave("spatial_feat_D1_feat_part1.pdf", width = 15, height = 15)

#ED_Fig5_J
SpatialFeaturePlot(object = D1_brain_no15,
                                  #Xuyu March 23 2024
                                  features = c("PDZRN3", "DLX6", "DLX6-AS1",
                                               "ADARB2", "ERBB4", "NRXN3", "DLX2", 
                                               "ZNF536", "PRKCA", "THRB", "TSHZ1", 
                                               "PBX3", "MEIS2", "CALB2", "CDCA7L", 
                                               "SYNPR", "SP8", "CASZ1", "FOXP4", "SP8"),
                                  ncol = 5,
                                  stroke = NA,
                                  max.cutoff = "1e+05",
                                  image.alpha = 0,
                                  pt.size.factor = 2.5)




seurat_plot1 <- SpatialFeaturePlot(object = bb5, 
                                   features = c("ETV1", "IGFBPL1", "UNC5C", 
                                                "TSHZ3", "PDE1A", "GABRA5"),
                                   ncol = 3,
                                   stroke = NA,
                                   max.cutoff = "1e+05",
                                   image.alpha = 0,
                                   pt.size.factor = 2.5)

seurat_plot1
#ETV1, IGFBPL1, UNC5C, TSHZ3, PDE1A, GABRA5 - JAN 26'24

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

#CALL the function 
bb5=rotateSeuratImage(A1brain_no14_annotated,rotation = "180")