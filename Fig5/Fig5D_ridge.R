library(ggplot2)
library(reticulate)
anndata = import("anndata")

adata = anndata$read_h5ad("gw20.h5ad")
obs = adata$obs
obs = obs[obs$sample == "FB080" & obs$region == "O1c", ]
obs = obs[!is.na(obs$v1_v2_dist), ]
types = c('EN-IT-L3-c0', 'EN-IT-L3/4-c3', 'EN-ET-L5/6-c4', 'EN-ET-SP-2-c1', 'EN-IT-L3-c3', 'EN-IT-L3/4-c0', 'EN-ET-L5/6-c3', 'EN-ET-SP-2-c3')
colors = c('#FF97EE', '#FF00FF', '#BC00BC', '#A800A8', '#ACFC64', '#00FF00', '#00C800', '#007E00')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3 = factor(obs1$H3_annotation, levels = types)
names(colors) = types
pdf(file = 'ridge.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = v1_v2_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('V1-V2 Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        theme(legend.title = element_blank()))
dev.off()

