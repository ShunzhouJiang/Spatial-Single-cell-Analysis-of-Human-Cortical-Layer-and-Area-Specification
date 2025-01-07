library(ggplot2)
library(reticulate)
anndata = import("anndata")

adata = anndata$read_h5ad("adult_umb5958.h5ad")
obs = adata$obs
obs = obs[!is.na(obs$cortical_depth), ]

types = c('EN-L2|3-2-c0', 'EN-L3-2-c2', 'EN-L4-c0', 'EN-L4-c1', 'EN-L4-c2', 'EN-L4-c3', 'EN-L4c-c1', 'EN-L2|3-1-c0', 'EN-L3-2-c1', 'EN-L4c-c0')
colors = c('#ff00ff', '#F07DC8', '#E068D6', '#D14BB2', '#C23CA0', '#B32A92', '#A800A8', '#2CFF00', '#8FE753', '#007E00')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3 = factor(obs1$H3_annotation, levels = types)
names(colors) = types
pdf(file = 'umb5958_ridge.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3, fill =H3 )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        theme(legend.title = element_blank()))
dev.off()

