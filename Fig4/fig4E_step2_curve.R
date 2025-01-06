# Results based on "fig4E_EDfig11c_gene_updown.ipynb"

library(dplyr)
library(ggplot2)

# Define gene and cluster lists
gene_lst <- c('GLRA2', 'TSHZ3', 'FAM107A', 'VAT1L', 'GABRA5', 'ATP2B4', 'NEFM', 'KLHL1', 'MLIP', 
              'NPR3', 'CYP26A1', 'FOXP1', 'SEMA3E', 'CUX2')
cluster_lst <- c("EN-ET-SP|L6b", "EN-ET-L5|6", "EN-IT-L5|6", "EN-IT-L4", "EN-IT-L3", "EN-IT-L2")

zs_dict <- list()


for (cluster in cluster_lst) {
  zs <- read.csv(paste0("result/", cluster, "_zs.csv"), row.names = 1)
  zs <- zs[, colnames(zs) %in% gene_lst]
  zs_dict[[cluster]] <- zs
}

gene_dict <- list()

# Create gene data frames
for (gene in gene_lst) {
  gene_val <- c()
  gw_val <- c()
  cluster_val <- c()
  
  for (cluster in cluster_lst) {
    gene_val <- c(gene_val, unlist(zs_dict[[cluster]][gene]))
    gw_val <- c(gw_val, rownames(zs_dict[[cluster]]))
    cluster_val <- c(cluster_val, rep(cluster, nrow(zs_dict[[cluster]])))
  }
  
  gene_df <- data.frame(GW = gw_val, Cluster = cluster_val, value = gene_val)
  gene_dict[[gene]] <- gene_df
}

# Define color palette
col_map <- c('EN-ET-SP|L6b' = '#9500ff',
             'EN-ET-L5|6' = '#4a00ff',
             'EN-IT-L5|6' = '#0a98ef',
             'EN-IT-L4' = '#0ff008',
             'EN-IT-L3' = '#ffe000',
             'EN-IT-L2' = '#ec2814')

# Plot and save figures
for (gene in gene_lst) {
  p <- ggplot(data = gene_dict[[gene]], aes(x = GW, y = value, group = Cluster, color = Cluster)) +
    geom_line(size = 0.6) +
    scale_color_manual(values = col_map) +
    labs(y = "Scaled Expression", title = gene, color = "Group") +
    theme_linedraw() + 
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 7, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = "white", linetype="solid", colour ="black", linewidth = 0.2))
  
  ggsave(filename = paste0("fig5D_ver2/", gene, ".pdf"), plot = p, width = 5, height = 4, dpi = 400)
}
