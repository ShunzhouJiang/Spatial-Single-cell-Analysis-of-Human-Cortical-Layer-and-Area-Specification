#Import Libs
library(Seurat)
require(data.table)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
set.seed(1234)

#### Fig_Ext_data18_a.2 - VISIUM A1 Analysis###

visium_rds <- readRDS('Visium_A1_brain_011124.rds')
DimPlot(visium_rds, label = TRUE)

#Example feauture list
features_list <- c("ABI3BP", "PDZRN4", "TAFA2")
plot3 <- SpatialFeaturePlot(visium_rds, features = features_list ,  ncol = 3, slot = "data", combine = TRUE,image.alpha = 0, stroke = 0.25, interactive = FALSE, information = NULL)
plot3


##### Fig_Ext_data18_c_e_f - scRNA_Seq Analysis #######

#Reading the processed RDS file
Final_Integration_scRNASeq <- readRDS("Final_Integration_scRNASeq_MDM_100323.rds")

#Extended_Fig_18_c
DimPlot(Final_Integration_scRNASeq, label = TRUE)

# Store cluster identities in object@meta.data$my.clusters
Final_Integration_scRNASeq[["my.clusters"]] <- Idents(Final_Integration_scRNASeq)

# Get number of cells per cluster and per sample of origin
cluster_table <- table(Final_Integration_scRNASeq@meta.data$my.clusters, Final_Integration_scRNASeq@meta.data$orig.ident)

# Convert the table to a data frame
cluster_df <- as.data.frame.matrix(cluster_table)
head(cluster_df)

# Convert row names to a new column 'Cell_Type'
cluster_df$Cell_Type <- rownames(cluster_df)

# Reshape the data to long format
library(tidyr)
cluster_df_long <- gather(cluster_df, key = "Sample", value = "Contribution", -Cell_Type)

# Calculate total contribution per sample
totals_per_sample <- cluster_df_long %>%
  group_by(Sample) %>%
  summarise(Total_Contribution = sum(Contribution))

# Merge total contributions with the original data
cluster_df_long <- merge(cluster_df_long, totals_per_sample, by = "Sample")

# Reorder samples based on total contribution (descending order)
ordered_samples <- c("FB121_3A","FB080_I","FB080_J")
custom_colors <- c("#e31a1c", "#a6cee3", "#1f78b4")
# Convert Sample to factor with ordered levels
cluster_df_long$Sample <- factor(cluster_df_long$Sample, 
                                 levels = ordered_samples)

# Create a ggplot stacked bar plot with proportions and without labels
# Extended_Fig_18_e
ggplot(cluster_df_long, aes(x = Cell_Type, y = Contribution, 
                            fill = Sample)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Proportion of Cell Types in Each Sample",
       x = "Cell Type", y = "Proportion") +
  scale_fill_manual(name = "Sample", 
                    values = setNames(custom_colors, 
                                      ordered_samples))

#Extended_Fig_18_f
FeaturePlot(Final_Integration_scRNASeq, features = c(
                                          "NPY","ABI3BP","PDZRN4" ),
            min.cutoff = "q9",
            order = TRUE,pt.size=0.01)
