## Results are based on files generated by "fig4E_EDfig11c_gene_updown.ipynb"
library(anndata)
library(stringr)
library(Matrix)
library(MASS)
library(reticulate)
library(Seurat)
library(dplyr)
library(splitstackshape)
library(ggplot2)
library(reshape2)
library(cowplot)


align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  # group_name <- group
  # pdf(paste0(paste0(path, group_name, sep = "/"), ".pdf",sep = ""))
  g
}

group_lst <- c("EN-ET-L5|6", "EN-ET-SP|L6b", "EN-IT-L2", "EN-IT-L3", "EN-IT-L4", "EN-IT-L5|6")
topgene_lst <- list()
bottomgene_lst <- list()

topgene_lst[["EN-ET-L5|6"]] <- c("NR4A2", "NPR3", "LPL", "RGS6", "KCTD12", "ZFPM2", "ZEB2", "ST18", "TSHZ3", "KLHL1")
topgene_lst[["EN-ET-SP|L6b"]] <- c("NPR3", "CYP26A1", "ST18", "HCRTR2", "SORCS1", "CDH13", "BACH2", "TSHZ3", "PDE1A", "MLIP")
topgene_lst[["EN-IT-L2"]] <- c("NEUROD6", "NEFM", "NR4A2", "RORB", "FEZF2", "FOXP1", "TAFA1", "SLITRK5", "NPY", "GLRA2")
topgene_lst[["EN-IT-L3"]] <- c("FEZF2", "NR4A2", "NEUROD6", "ZEB2", "SATB2", "BACH2", "NIPBL.DT", "TSHZ3", "PEX5L", "NEFM")
topgene_lst[["EN-IT-L4"]] <- c("NR4A2", "NPY", "GLRA2", "NIPBL.DT", "NEFM", "TAFA1", "LHX2", "B3GALT2", "KIF26B", "SATB2")
topgene_lst[["EN-IT-L5|6"]] <- c("RSPO3", "ACTN2", "KITLG", "GLRA2", "SEMA3E", "NIPBL.DT", "MYT1L", "STK32B", "BACH2", "GPR85")

bottomgene_lst[["EN-ET-L5|6"]] <- c("VAT1L", "HSPA12A", "NEFM", "GABRB2", "ATP2B4", "MLIP", "PCDH17", "CYP26A1", "SUSD4", "SEMA3E")
bottomgene_lst[["EN-ET-SP|L6b"]] <- c("CBLN2", "ATP2B4", "VAT1L", "TRPM3", "GABRB2", "RAP1GAP", "CCN2", "NEFM", "HSPA12A", "FOXP2")
bottomgene_lst[["EN-IT-L2"]] <- c("PENK", "PRSS12", "NEUROD1", "SCG2", "CPNE8", "CUX2", "CDH13", "KCNJ6", "FAM107A", "NPR3")
topgene_lst[["EN-IT-L3"]] <- c("CUX2", "SCG2", "GABRA5", "CAMK2B", "KCNJ6", "FAM107A", "TMOD1", "SLC17A6", "PENK", "CUX1")
topgene_lst[["EN-IT-L4"]] <- c("GABRA5", "SPOCK1", "FOXP1", "NPNT", "RASGRF2", "GPRC5B", "KLHL1", "FBXO32", "FAM107A", "FSTL5")
topgene_lst[["EN-IT-L5|6"]] <- c("GABRA5", "FSTL5", "KLHL1", "VAT1L", "SPOCK1", "MLIP", "AHI1", "ATP2B4", "FBXO32", "RASGRF2")

for (i in 1:length(group_lst)) {
  setwd("result")
  group <- group_lst[i]
  gene_rank <- read.csv(paste0(group, "_rank.csv"), row.names = 1)
  gene_zs <- read.csv(paste0(group, "_zs.csv"), row.names = 1)
  gene_nz <- read.csv(paste0(group, "_nz.csv"), row.names = 1)
  
  setwd("../")
  gene_incre <- topgene_lst[[group]]
  gene_decre <- bottomgene_lst[[group]]

  zs_top <- gene_zs[, gene_incre]
  nz_top <- gene_nz[, gene_incre]
  zs_bottom <- gene_zs[, gene_decre]
  nz_bottom <- gene_nz[, gene_decre]
  
  order_top <- order(as.vector(unlist(nz_top[1,])), decreasing = TRUE)
  order_bottom <- order(as.vector(unlist(nz_bottom[4,])), decreasing = FALSE)
  zs_top <- zs_top[, order_top]
  zs_bottom <- zs_bottom[, order_bottom]
  zs_select <- cbind(zs_top, zs_bottom)
  
  topgene_lst[[group]] <- colnames(zs_top)
  bottomgene_lst[[group]] <- colnames(zs_bottom)
  
  nz_top <- nz_top[, order_top]
  nz_bottom <- nz_bottom[, order_bottom]
  count_select <- cbind(nz_top, nz_bottom)
  count_select <- count_select * 100
  
  count_select$layer <- rownames(count_select)
  zs_select$layer <- rownames(zs_select)
  
  nonzero_sub_melt <- melt(count_select, id = c("layer"))
  expr_sub_melt <- melt(zs_select, id = c("layer"))
  color_map <- expr_sub_melt$value
  mid <- mean(color_map)
  
  col_min <- floor(min(zs_select[,1:(ncol(zs_select)-1)]))
  col_max <- ceiling(max(zs_select[,1:(ncol(zs_select)-1)]))
  
  x = ggplot(nonzero_sub_melt, aes(x = layer, y = variable, size = value, color = color_map)) + 
    # geom_point(aes(size = value, fill = layer), alpha = 1, shape = 21) + 
    geom_point(aes(size = value)) + 
    scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(20, 40, 60, 80)) +
    labs( x= "", y = "", size = "Percentage of expressed cells (%)", fill = "", color = "Mean expression")  +
    theme(legend.title.align = 0.5,
          legend.text.align = 0.5,
          legend.key=element_blank(),
          axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45), 
          axis.text.y = element_text(colour = "black", face = "bold", size = 10), 
          legend.text = element_text(size = 8, face ="bold", colour ="black"), 
          legend.title = element_text(size = 8, face = "bold"), 
          panel.background = element_blank(),  #legend.spacing.x = unit(0.5, 'cm'),
          legend.position = "right",
          legend.box.just = "top",
          legend.direction = "vertical",
          legend.box="vertical",
          legend.justification="center"
    ) +  
    # theme_minimal() +
    # theme(legend.title.align = 0)+
    
    # scale_fill_manual(values = color_map, guide = FALSE) + 
    scale_y_discrete(limits = rev(levels(nonzero_sub_melt$variable))) +
    scale_color_gradient(low="yellow", high="blue", space ="Lab", limits = c(col_min, col_max), position = "bottom") +
    # scale_fill_viridis_c(guide = FALSE, limits = c(-0.5,1.5)) + 
    guides(colour = guide_colourbar(direction = "vertical", barheight = unit(4, "cm"), title.vjust = 4))
  
  
  
  setwd("fig_gw")
  # pdf(paste0(paste(group, "gw", sep = "_"), ".pdf",sep = ""))
  ggdraw(align_legend(x))
  ggsave(paste0(paste(group, "gw", sep = "_"), ".pdf",sep = ""),
         width = 7, height = 8, limitsize = TRUE)
  # dev.off()
  setwd("../")
}

saveRDS(topgene_lst, "topgene_lst_final.rds")
saveRDS(bottomgene_lst, "bottomgene_lst_final.rds")


