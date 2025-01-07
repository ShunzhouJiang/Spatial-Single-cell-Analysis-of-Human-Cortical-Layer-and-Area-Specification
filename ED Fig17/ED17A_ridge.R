library(ggplot2)
library(reticulate)
anndata = import("anndata")

adata = anndata$read_h5ad("gw34_umb5900_ba17.h5ad")
obs = adata$obs
obs = obs[!is.na(obs$cortical_depth), ]

types = c('EN-IT-L3|4-c2', 'EN-IT-L3|4-c3', 'EN-IT-L3|4-c4', 'EN-IT-L3|4-c5', 'EN-IT-L3|4-c6', 'EN-IT-L3|4-c7', 'EN-IT-L4-2-c2', 'EN-IT-L4-2-c3', 'EN-ET-L6-c1', 'EN-ET-L6-V1-c3', 'EN-IT-L3|4-c1', 'EN-IT-L3-c0', 'EN-IT-L4-2-c1', 'EN-ET-L5-c1')
colors = c('#FF97EE', '#F586E6', '#EB75DE', '#E264D6', '#D853CE', '#CE43C7', '#C532BF', '#BB21B7', '#B110AF', '#A800A8', '#ACFC64', '#99EE58', '#3BA721', '#007E00')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3 = factor(obs1$H3_annotation, levels = types)
names(colors) = types
pdf(file = 'umb5900_ba17_ridge.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3, fill =H3 )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        theme(legend.title = element_blank()))
dev.off()

