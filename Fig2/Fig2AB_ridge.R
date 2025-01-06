library(ggplot2)
library(reticulate)
installed <- py_module_available("anndata")

# If not installed, install it
if (!installed) {
  cat("anndata is not installed. Installing...\n")
  py_install("anndata", pip = TRUE)
} else {
  cat("anndata is already installed.\n")
}

anndata = import("anndata")
adata = anndata$read_h5ad("merscope_integrated_855.h5ad")
obs = adata$obs
obs = obs[obs$sample == "UMB5900" & obs$region == "BA9", ]
obs = obs[!is.na(obs$cortical_depth), ]
types = c('EN-L2-4', 'EN-IT-L3-late', 'EN-IT-L4-late', 'EN-ET-L5-1', 'EN-IT-L6-late')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c("#EC2814", "#F0A608", "#0FF008", '#0832F0', '#F008EC')
names(colors) = types
pdf(file = 'ridge_a.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        geom_vline(xintercept=c(0.9416051031562357, 0.7650373084425757, 0.6371757146516088, 0.5307999237922936), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(),legend.text = element_text(size = 5), panel.background = element_blank()))
dev.off()

types = c('EN-ET-L5-1', 'EN-ET-L6-V1', 'EN-ET-L5/6', 'EN-ET-L6-A', 'EN-ET-SP-early3', 'EN-ET-L6-early4', 'EN-ET-L6-P', 'EN-ET-SP-early4', 'EN-ET-SP-early2', 'EN-ET-SP-3', 'EN-ET-SP-A', 'EN-ET-SP-early1', 'EN-ET-L6-early1', 'EN-ET-L6-early3', 'EN-ET-SP-early5', 'EN-ET-SP-2', 'EN-ET-L6-early2', 'EN-ET-SP-1', 'EN-ET-SP-5', 'EN-ET-SP-4', 'EN-ET-SP-P2', 'EN-ET-SP-P1', 'EN-ET-SP-V1T2', 'EN-ET-L6-early5', 'EN-ET-SP-V1T1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939')
names(colors) = types
pdf(file = 'ridge_b1.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        xlim(0,1) +
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.9416051031562357, 0.7650373084425757, 0.6371757146516088, 0.5307999237922936), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L4-1', 'EN-IT-Hip', 'EN-IT-L6-2', 'EN-IT-L5-1', 'EN-IT-L4/5-early', 'EN-IT-L4/5-late', 'EN-IT-L4/5-1', 'EN-IT-L6-1', 'EN-IT-L6-late', 'EN-IT-L5/6-P')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
names(colors) = types
pdf(file = 'ridge_b2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.9416051031562357, 0.7650373084425757, 0.6371757146516088, 0.5307999237922936), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L3-V1', 'EN-IT-L2/3-A2', 'EN-IT-L3-A', 'EN-IT-L3/4-1', 'EN-IT-L2/3-A1', 'EN-IT-L3/4-early', 'EN-IT-L4-A', 'EN-IT-L4-V1', 'EN-IT-L3-late', 'EN-IT-L3/4-P', 'EN-IT-L3/4-P2', 'EN-IT-L3-P', 'EN-IT-L4-late', 'EN-IT-L3/4-T')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2')
names(colors) = types
pdf(file = 'ridge_b3.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) +
        geom_vline(xintercept=c(0.9416051031562357, 0.7650373084425757, 0.6371757146516088, 0.5307999237922936), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

