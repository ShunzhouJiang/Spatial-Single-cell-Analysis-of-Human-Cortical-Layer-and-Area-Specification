#Import Libraries 
library(Seurat)
require(data.table)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)

#Import Macaque_sub_cluster_Ext15_e.rds
Macaque_sub_cluster <- readRDS("Macaque_sub_cluster_Ext15_e.rds")

#ExtendedFig15_e
DimPlot(Macaque_sub_cluster, label = TRUE) # 3 clusters

#ExtendedFig15_f
FeaturePlot(Macaque_sub_cluster, features = c("NPY", "PDZRN4", "ABI3BP"),
            raster=FALSE, min.cutoff = "q9",
            order = TRUE,pt.size=0.01)

#ExtendedFig15_G
# Get number of cells per refined_region 
refined_region_table <- table(Macaque_sub_cluster@meta.data$my.clusters, 
                              Macaque_sub_cluster@meta.data$refined_region)

# Order age by refined region as well 
# Create a contingency table with three variables
# Convert the multi-dimensional table to a data frame for easier manipulation (optional)
age_group_df <- as.data.frame(age_group_table)

# Convert the tables to data frames
refined_region_df <- as.data.frame.matrix(refined_region_table)

# Convert row names to a new column 'Cell_Type'
refined_region_df$Cell_Type <- rownames(refined_region_df)
age_group_df$Cell_Type <- rownames(age_group_df)


# Reshape the data to long format
refined_region_long <- gather(refined_region_df, key = "Refined_Region",
                              value = "Contribution",
                              -Cell_Type)

#convert contribution to numeric 
refined_region_long$Contribution <- as.numeric(refined_region_long$Contribution)

#remove NAs
refined_region_long <- refined_region_long %>% filter(!is.na(Contribution))

# Calculate total contribution per refined_region 
totals_per_refined_region <- refined_region_long %>%
  group_by(Refined_Region) %>%
  summarise(Total_Contribution = sum(Contribution))

# Merge total contributions with the original data
refined_region_long <- merge(refined_region_long, totals_per_refined_region,
                             by = "Refined_Region")

# Reorder refined_region on total contribution (descending order)
# Convert Refined_Region to factor with ordered levels
refined_region_long$Refined_Region <- factor(refined_region_long$Refined_Region)

# Bringing #35 to left
refined_region_long$Cell_Type <- factor(refined_region_long$Cell_Type)

# Move cluster 35 to the leftmost position
refined_region_long$Cell_Type <- factor(refined_region_long$Cell_Type, 
                                        levels = c("35", setdiff(levels(refined_region_long$Cell_Type), "35")))

##### COLORS #####
# Define custom colors for the Refined_Region
custom_colors <- c(
  "A1C" = "blue",
  "Anterior" = "orange",
  "DFC" = "red1",
  "Dorsolateral" = "red3",
  "Frontal" = "brown",
  "Insula" = "pink",
  "IPC" = "navy",
  "ITC" = "gold",
  "M1C" = "lightgreen",
  "MFC" = "coral",
  "Motor-somatosensory" = "lavender",
  "Occipital" = "maroon1",
  "OFC" = "lightblue",
  "PCC" = "turquoise",
  "Posterior" = "khaki",
  "S1C" = "cyan",
  "STC" = "darkolivegreen",
  "Temporal" = "tomato",
  "V1C" = "magenta",
  "VFC" = "green"
)

#EXTENDED_FIG_15G
# Create a ggplot stacked bar plot for refined_region with cluster 35 on the left
plot_refined_region <- ggplot(refined_region_long, aes(x = Cell_Type, y = Contribution, 
                                                       fill = Refined_Region)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(title = "Proportion of Cell Types in Each Refined Region",
       x = "Cell Type", y = "Proportion") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_refined_region
