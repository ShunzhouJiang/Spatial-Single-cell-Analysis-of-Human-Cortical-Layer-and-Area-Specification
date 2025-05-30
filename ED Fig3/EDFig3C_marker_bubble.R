## Results are based on files generated by "EDFig3C_marker_detect.py"

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
library(readxl)

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

zs_norm <- function(zs_select) {
  for (i in 1:ncol(zs_select)) {
    zs_select[, i] <- (zs_select[, i] - mean(zs_select[, i])) / sd(zs_select[, i])
  }
  return(zs_select)
}


# marker_total <- readRDS("h2_marker.rds")
marker_total <- read.csv("h3_marker.csv", row.names = 1)
marker_num <- 1
marker_total <- marker_total[1:marker_num, ]

# h2_order <- read_excel("smartseq_singler_h3.xlsx")
# h2_order <- h2_order$H3_annotation
h2_order <- c("EC", "Astro-1", "Astro-2", "Astro-early", "Astro-late1", "Astro-late2", 
                "Astro-late3", "OPC", "IN-CGE1", "IN-CGE2", "IN-CGE3", "IN-CGE4", 
                "IN-CGE5", "IN-GE1", "IN-GE2", "IN-VZ/LGE1", "IN-VZ/LGE2", "INP-LGE1", 
                "INP-LGE2", "INP-VZ", "INP-VZ/GE1", "INP-VZ/GE2", "IN-VIP-late", 
                "IN-MGE4", "IN-MGE5", "IN-SST1", "IN-SST2", "IN-SST3", "IN-SST4", 
                "IN-SST-NPY", "IN-MGE1", "IN-MGE2", "IN-MGE3", "tRG-1", "tRG-2", 
                "tRG-early", "vRG-early", "vRG-early1", "vRG-early2", "oRG-1", 
                "oRG-2", "oRG-3", "oRG-4", "vRG-late1", "vRG-late2", "vRG-late3", 
                "oRG-5", "oRG-6", "oRG-7", "oRG-8", "oRG-9", "IPC-SVZ-1", 
                "IPC-VZ/SVZ", "IPC-oSVZ", "IPC-iSVZ", "IPC-SVZ-2", "EN-oSVZ-1", 
                "En-oSVZ-2", "EN-IZ-1", "EN-IZ-2", "EN-IZ-3", "EN-L2-1", 
                "EN-L2-2", "EN-L2-3", "EN-L2-4", "EN-ET-L6-early1", "EN-ET-L6-early2", 
                "EN-ET-L6-early3", "EN-ET-L6-early4", "EN-ET-L6-early5", "EN-ET-SP-1", 
                "EN-ET-SP-2", "EN-ET-SP-3", "EN-ET-SP-4", "EN-ET-SP-5", "EN-ET-SP-A", 
                "EN-ET-L5-1", "EN-ET-L5/6", "EN-ET-L6-A", "EN-ET-L6-P", "EN-ET-L6-V1", 
                "EN-ET-SP-P1", "EN-ET-SP-P2", "EN-ET-SP-V1T1", "EN-ET-SP-V1T2", 
                "EN-ET-SP-early1", "EN-ET-SP-early2", "EN-ET-SP-early3", 
                "EN-ET-SP-early4", "EN-ET-SP-early5", "EN-IT-Hip", "EN-IT-L4/5-1", 
                "EN-IT-L4/5-early", "EN-IT-L4-1", "EN-IT-L2/3-A2", "EN-IT-L2/3-A1", 
                "EN-IT-L3-A", "EN-IT-L3-P", "EN-IT-L3-V1", "EN-IT-L3-late", 
                "EN-IT-L3/4-1", "EN-IT-L3/4-P", "EN-IT-L3/4-P2", "EN-IT-L3/4-T", 
                "EN-IT-L3/4-early", "EN-IT-L4-A", "EN-IT-L4-V1", "EN-IT-L4-late", 
                "EN-IT-L4/5-late", "EN-IT-L5-1", "EN-IT-L5/6-P", "EN-IT-L6-1", 
                "EN-IT-L6-2", "EN-IT-L6-late")

gene_tot <- c()
for (i in 1:length(h2_order)) {
  h2_i <- gsub("-", ".", h2_order[i])
  h2_i <- gsub("/", ".", h2_i)
  marker_i <- marker_total[, h2_i]
  gene_tot <- c(gene_tot, marker_i)
}
gene_tot <- unique(gene_tot)
gene_tot <- gsub("-", ".", gene_tot)


count_avg <- read.csv("result/h3_nz.csv", row.names = 1)
zs_avg <- read.csv("result/h3_zs.csv", row.names = 1)



count_select <- count_avg[, gene_tot]
count_select <- count_select * 100
zs_select <- zs_avg[, gene_tot]
count_select <- count_select[h2_order, ]
zs_select <- zs_select[h2_order, ]
zs_select <- zs_norm(zs_select)

count_select$layer <- rownames(count_select)
zs_select$layer <- rownames(zs_select)

nonzero_sub_melt <- melt(count_select, id = c("layer"))
expr_sub_melt <- melt(zs_select, id = c("layer"))
nonzero_sub_melt$layer <- factor(nonzero_sub_melt$layer, levels = h2_order)
expr_sub_melt$layer <- factor(expr_sub_melt$layer, levels = h2_order)
color_map <- expr_sub_melt$value
mid <- mean(color_map)

# col_min <- floor(min(zs_select[,1:(ncol(zs_select)-1)]))
col_min <- min(zs_select[,1:(ncol(zs_select)-1)]) - 0.5
# col_max <- ceiling(max(zs_select[,1:(ncol(zs_select)-1)]))
col_max <- max(zs_select[,1:(ncol(zs_select)-1)]) + 0.5

x = ggplot(nonzero_sub_melt, aes(x = layer, y = variable, size = value, color = color_map)) + 
  # geom_point(aes(size = value, fill = layer), alpha = 1, shape = 21) + 
  geom_point(aes(size = value)) + 
  scale_size_continuous(limits = c(-0.000001, 100), range = c(1,10), breaks = c(20, 40, 60, 80)) +
  labs( x= "", y = "", size = "Percentage of expressed cells (%)", fill = "", color = "Mean expression")  +
  theme(legend.title.align = 0.5,
        legend.text.align = 0.5,
        legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 90), 
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

if (marker_num == 1) {
  height <- 26
} else if (marker_num == 3) {
  height <- 36
} else if (marker_num == 5) {
  height <- 42
}

ggdraw(align_legend(x))
ggsave(paste0(paste("h3_marker", marker_num, sep = "_"), ".pdf"),
       width = 38, height = height, limitsize = TRUE)
