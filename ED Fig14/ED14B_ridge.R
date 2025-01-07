library(ggplot2)
library(reticulate)
anndata = import("anndata")

adata = anndata$read_h5ad("gw22_umb1932.h5ad")
obs = adata$obs
obs = obs[!is.na(obs$v1_v2_dist), ]

types = c('EN-IT-L2|3-c1', 'EN-IT-L4-c2', 'EN-IT-L4|5-c2', 'EN-IT-L5|6-2-c1', 'EN-ET-L5|6-c1', 'EN-IT-L2|3-c2', 'EN-IT-L4-c1', 'EN-IT-L4|5-c0', 'EN-IT-L5|6-2-c2', 'EN-ET-L5|6-c0')
colors = c('#FF97EE', '#E971DD', '#D24CCC', '#BC26BB', '#A800A8', '#ACFC64', '#81DC4B', '#56BD32', '#2B9D19', '#007E00')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
names(colors) = types
pdf(file = 'umb1932_v1_v2_ridge.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = v1_v2_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('V1-V2 Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        theme(legend.title = element_blank()))
dev.off()

