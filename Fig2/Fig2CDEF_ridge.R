library(ggplot2)
library(reticulate)
anndata = import("anndata")

adata = anndata$read_h5ad("merscope_integrated_855.h5ad")
<<<<<<< HEAD
obs = adata$obs
obs = obs[obs$sample == "FB123" & obs$region == "F1", ]
=======
obs_all = adata$obs
obs = obs_all[obs_all$sample == "FB123" & obs_all$region == "F1", ]
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
obs = obs[!is.na(obs$cortical_depth), ]

types = c('EN-L2-1', 'EN-IT-L3-A', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1')
obs1 = obs[obs$H3_annotation %in% c('EN-L2-1', 'EN-IT-L3-A', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1'),]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = c('EN-L2-1', 'EN-IT-L3-A', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1'))
colors = c("#EC2814", "#F0A608", "#0FF008", '#0832F0', '#F008EC')
names(colors) = types
<<<<<<< HEAD
pdf(file = 'ridge_c1.pdf', width=5.58, height=3.5)
=======
pdf(file = 'ridge_c_gw22_pfc.pdf', width=5.58, height=3.5)
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        geom_vline(xintercept=c(0.7090927649311176, 0.5404574294859591, 0.4336352966282663, 0.3039781775549624), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(),legend.text = element_text(size = 5), panel.background = element_blank()))
dev.off()

types = c('EN-ET-L5-1', 'EN-ET-L6-V1', 'EN-ET-L5/6', 'EN-ET-L6-A', 'EN-ET-SP-early3', 'EN-ET-L6-early4', 'EN-ET-L6-P', 'EN-ET-SP-early4', 'EN-ET-SP-early2', 'EN-ET-SP-3', 'EN-ET-SP-A', 'EN-ET-SP-early1', 'EN-ET-L6-early1', 'EN-ET-L6-early3', 'EN-ET-SP-early5', 'EN-ET-SP-2', 'EN-ET-L6-early2', 'EN-ET-SP-1', 'EN-ET-SP-5', 'EN-ET-SP-4', 'EN-ET-SP-P2', 'EN-ET-SP-P1', 'EN-ET-SP-V1T2', 'EN-ET-L6-early5', 'EN-ET-SP-V1T1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939')
names(colors) = types
<<<<<<< HEAD
pdf(file = 'ridge_g3.pdf', width=5.58, height=3.5)
=======
pdf(file = 'ridge_d_gw22_pfc.pdf', width=5.58, height=3.5)
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        xlim(0,1) +
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7090927649311176, 0.5404574294859591, 0.4336352966282663, 0.3039781775549624), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5), panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L4-1', 'EN-IT-Hip', 'EN-IT-L6-2', 'EN-IT-L5-1', 'EN-IT-L4/5-early', 'EN-IT-L4/5-late', 'EN-IT-L4/5-1', 'EN-IT-L6-1', 'EN-IT-L6-late', 'EN-IT-L5/6-P')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
names(colors) = types
<<<<<<< HEAD
pdf(file = 'ridge_h3.pdf', width=5.58, height=3.5)
=======
pdf(file = 'ridge_e_gw22_pfc.pdf', width=5.58, height=3.5)
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7090927649311176, 0.5404574294859591, 0.4336352966282663, 0.3039781775549624), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5), panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L3-V1', 'EN-IT-L2/3-A2', 'EN-IT-L3-A', 'EN-IT-L3/4-1', 'EN-IT-L2/3-A1', 'EN-IT-L3/4-early', 'EN-IT-L4-A', 'EN-IT-L4-V1', 'EN-IT-L3-late', 'EN-IT-L3/4-P', 'EN-IT-L3/4-P2', 'EN-IT-L3-P', 'EN-IT-L4-late', 'EN-IT-L3/4-T')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2')
names(colors) = types
<<<<<<< HEAD
pdf(file = 'ridge_i3.pdf', width=5.58, height=3.5)
=======
pdf(file = 'ridge_f_gw22_pfc.pdf', width=5.58, height=3.5)
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) +
        geom_vline(xintercept=c(0.7090927649311176, 0.5404574294859591, 0.4336352966282663, 0.3039781775549624), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5), panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()


<<<<<<< HEAD
setwd(paste(dir, '/', 'FB123_O2', sep = ''))
obs = read.csv('/Users/kylecoleman/data/walsh/all/clustering2/annotations_completed_cortical_dist/FB123_R6-2/A-Occi_obs_cp.csv')
obs$cp_dist = obs$cp_dist**0.5
types = c('EN-L2-1', 'EN-IT-L3-A', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c("#EC2814", "#F0A608", "#0FF008", '#0832F0', '#F008EC')
names(colors) = types
pdf(file = 'ridge_c2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        geom_vline(xintercept=c(0.7868768296116182, 0.6575388963412806, 0.5134491993982202, 0.35806069002273433), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(),legend.text = element_text(size = 5), panel.background = element_blank()))
dev.off()

=======
obs = obs_all[obs_all$sample == "FB123" & obs_all$region == "O2" & obs_all$area=='A-Occi', ]
obs = obs[!is.na(obs$cortical_depth), ]
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a

types = c('EN-ET-L5-1', 'EN-ET-L6-V1', 'EN-ET-L5/6', 'EN-ET-L6-A', 'EN-ET-SP-early3', 'EN-ET-L6-early4', 'EN-ET-L6-P', 'EN-ET-SP-early4', 'EN-ET-SP-early2', 'EN-ET-SP-3', 'EN-ET-SP-A', 'EN-ET-SP-early1', 'EN-ET-L6-early1', 'EN-ET-L6-early3', 'EN-ET-SP-early5', 'EN-ET-SP-2', 'EN-ET-L6-early2', 'EN-ET-SP-1', 'EN-ET-SP-5', 'EN-ET-SP-4', 'EN-ET-SP-P2', 'EN-ET-SP-P1', 'EN-ET-SP-V1T2', 'EN-ET-L6-early5', 'EN-ET-SP-V1T1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939')
names(colors) = types
<<<<<<< HEAD
pdf(file = 'ridge_g4.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
=======
pdf(file = 'ridge_d_gw22_v2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        xlim(0,1) +
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7868768296116182, 0.6575388963412806, 0.5134491993982202, 0.35806069002273433), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L4-1', 'EN-IT-Hip', 'EN-IT-L6-2', 'EN-IT-L5-1', 'EN-IT-L4/5-early', 'EN-IT-L4/5-late', 'EN-IT-L4/5-1', 'EN-IT-L6-1', 'EN-IT-L6-late', 'EN-IT-L5/6-P')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
names(colors) = types
<<<<<<< HEAD
pdf(file = 'ridge_h4.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
=======
pdf(file = 'ridge_e_gw22_v2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7868768296116182, 0.6575388963412806, 0.5134491993982202, 0.35806069002273433), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L3-V1', 'EN-IT-L2/3-A2', 'EN-IT-L3-A', 'EN-IT-L3/4-1', 'EN-IT-L2/3-A1', 'EN-IT-L3/4-early', 'EN-IT-L4-A', 'EN-IT-L4-V1', 'EN-IT-L3-late', 'EN-IT-L3/4-P', 'EN-IT-L3/4-P2', 'EN-IT-L3-P', 'EN-IT-L4-late', 'EN-IT-L3/4-T')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2')
names(colors) = types
<<<<<<< HEAD
pdf(file = 'ridge_i4.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
=======
pdf(file = 'ridge_f_gw22_v2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cortical_depth, stat(count), color=H3_annotation, fill =H3_annotation )) + 
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) +
        geom_vline(xintercept=c(0.7868768296116182, 0.6575388963412806, 0.5134491993982202, 0.35806069002273433), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

<<<<<<< HEAD

setwd(paste(dir, '/', 'FB121_F1/Cing', sep = ''))
obs = read.csv('/Users/kylecoleman/data/walsh/all/clustering2/annotations_completed_cortical_dist/FB121_GW20-3A/B-Cing_obs_cp.csv')
obs$cp_dist = obs$cp_dist**0.5
types = c('EN-IT-L2/3-A1', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c("#F0A608", "#0FF008", '#0832F0', '#F008EC')
names(colors) = types
pdf(file = 'ridge_d.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        geom_vline(xintercept=c(0.7190983855177338, 0.5496677731062829, 0.3135250085751593), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(),legend.text = element_text(size = 5), panel.background = element_blank()))
dev.off()


types = c('EN-ET-L5-1', 'EN-ET-L6-V1', 'EN-ET-L5/6', 'EN-ET-L6-A', 'EN-ET-SP-early3', 'EN-ET-L6-early4', 'EN-ET-L6-P', 'EN-ET-SP-early4', 'EN-ET-SP-early2', 'EN-ET-SP-3', 'EN-ET-SP-A', 'EN-ET-SP-early1', 'EN-ET-L6-early1', 'EN-ET-L6-early3', 'EN-ET-SP-early5', 'EN-ET-SP-2', 'EN-ET-L6-early2', 'EN-ET-SP-1', 'EN-ET-SP-5', 'EN-ET-SP-4', 'EN-ET-SP-P2', 'EN-ET-SP-P1', 'EN-ET-SP-V1T2', 'EN-ET-L6-early5', 'EN-ET-SP-V1T1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939')
names(colors) = types
pdf(file = 'ridge_g2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        xlim(0,1) +
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7190983855177338, 0.5496677731062829, 0.3135250085751593), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L4-1', 'EN-IT-Hip', 'EN-IT-L6-2', 'EN-IT-L5-1', 'EN-IT-L4/5-early', 'EN-IT-L4/5-late', 'EN-IT-L4/5-1', 'EN-IT-L6-1', 'EN-IT-L6-late', 'EN-IT-L5/6-P')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
names(colors) = types
pdf(file = 'ridge_h2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7190983855177338, 0.5496677731062829, 0.3135250085751593), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L3-V1', 'EN-IT-L2/3-A2', 'EN-IT-L3-A', 'EN-IT-L3/4-1', 'EN-IT-L2/3-A1', 'EN-IT-L3/4-early', 'EN-IT-L4-A', 'EN-IT-L4-V1', 'EN-IT-L3-late', 'EN-IT-L3/4-P', 'EN-IT-L3/4-P2', 'EN-IT-L3-P', 'EN-IT-L4-late', 'EN-IT-L3/4-T')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2')
names(colors) = types
pdf(file = 'ridge_i2.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) +
        geom_vline(xintercept=c(0.7190983855177338, 0.5496677731062829, 0.3135250085751593), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()


setwd(paste(dir, '/', 'UMB1117_F1a', sep = ''))
obs = read.csv('/Users/kylecoleman/data/walsh/all/clustering2/annotations_completed_cortical_dist/UMB1117_FP/A-PFC_obs_cp.csv')
obs$cp_dist = obs$cp_dist**0.5
types = c('EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c("#0FF008", '#0832F0', '#F008EC')
names(colors) = types
pdf(file = 'ridge_e.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)+
        geom_vline(xintercept=c(0.7094192203468167, 0.5015364548312168), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(),legend.text = element_text(size = 5), panel.background = element_blank()))
dev.off()


types = c('EN-ET-L5-1', 'EN-ET-L6-V1', 'EN-ET-L5/6', 'EN-ET-L6-A', 'EN-ET-SP-early3', 'EN-ET-L6-early4', 'EN-ET-L6-P', 'EN-ET-SP-early4', 'EN-ET-SP-early2', 'EN-ET-SP-3', 'EN-ET-SP-A', 'EN-ET-SP-early1', 'EN-ET-L6-early1', 'EN-ET-L6-early3', 'EN-ET-SP-early5', 'EN-ET-SP-2', 'EN-ET-L6-early2', 'EN-ET-SP-1', 'EN-ET-SP-5', 'EN-ET-SP-4', 'EN-ET-SP-P2', 'EN-ET-SP-P1', 'EN-ET-SP-V1T2', 'EN-ET-L6-early5', 'EN-ET-SP-V1T1')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede', '#637939')
names(colors) = types
pdf(file = 'ridge_g1.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        xlim(0,1) +
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7094192203468167, 0.5015364548312168), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L4-1', 'EN-IT-Hip', 'EN-IT-L6-2', 'EN-IT-L5-1', 'EN-IT-L4/5-early', 'EN-IT-L4/5-late', 'EN-IT-L4/5-1', 'EN-IT-L6-1', 'EN-IT-L6-late', 'EN-IT-L5/6-P')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
names(colors) = types
pdf(file = 'ridge_h1.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) + 
        geom_vline(xintercept=c(0.7094192203468167, 0.5015364548312168), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

types = c('EN-IT-L3-V1', 'EN-IT-L2/3-A2', 'EN-IT-L3-A', 'EN-IT-L3/4-1', 'EN-IT-L2/3-A1', 'EN-IT-L3/4-early', 'EN-IT-L4-A', 'EN-IT-L4-V1', 'EN-IT-L3-late', 'EN-IT-L3/4-P', 'EN-IT-L3/4-P2', 'EN-IT-L3-P', 'EN-IT-L4-late', 'EN-IT-L3/4-T')
obs1 = obs[obs$H3_annotation %in% types,]
obs1$H3_annotation = factor(obs1$H3_annotation, levels = types)
colors = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2')
names(colors) = types
pdf(file = 'ridge_i1.pdf', width=5.58, height=3.5)
print(ggplot(obs1, aes(x = cp_dist, stat(count), color=H3_annotation, fill =H3_annotation )) + 
        geom_density(alpha=0.3)+  xlab('Cortical Depth') + ylab('Count') + 
        coord_flip() + scale_y_reverse() + 
        scale_fill_manual(values=colors) + 
        scale_color_manual(values=colors) +
        geom_vline(xintercept=c(0.7094192203468167, 0.5015364548312168), linetype="dashed", color = "#1f77b4")+
        theme(legend.title = element_blank(), legend.text = element_text(size = 5),panel.background = element_blank())+ 
        guides(color = guide_legend(override.aes = list(size = 0.3))))
dev.off()

=======
>>>>>>> 89902e0c8a64e6f87e8e0335b6700c05721fa69a
