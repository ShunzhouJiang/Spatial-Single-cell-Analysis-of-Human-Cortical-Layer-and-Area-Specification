library(anndata)
library(stringr)
library(Matrix)
library(MASS)
library(reticulate)
library(Seurat)
library(readxl)
library(dplyr)
library(splitstackshape)

use_python("/usr/local/bin/python3.10")
gw15 <- read_h5ad("gw15.h5ad")
gw20 <- read_h5ad("gw20.h5ad")
gw22 <- read_h5ad("gw22.h5ad")
gw34 <- read_h5ad("gw34.h5ad")


en_group <- read_xlsx("../gene_expression/EN groups based on sankey.xlsx", col_names = TRUE)

# Helper function to replace "cluster" with "Cluster"
rep_cluster <- function(name) {
  gsub("cluster", "Cluster", name)
}

# Define key list and data list
key_lst <- c("clusters from GW15", "clusters from GW20", "clusters from GW22", "clusters from GW34")
data_lst <- c("gw15", "gw20", "gw22", "gw34")

for (i in seq_along(key_lst)) {
  key <- key_lst[i]
  
  # Create group dictionary
  group_dict <- en_group %>%
    filter(!is.na(!!sym(key))) %>%
    mutate(key_col = rep_cluster(!!sym(key))) %>%
    select(key_col, Group.name) %>%
    deframe()
  
  # Map and update the corresponding data
  data <- get(data_lst[i])
  data$obs$group <- group_dict[match(data$obs$H3_annotation, names(group_dict))]
  data <- data[!is.na(data$obs$group), ]
  data$obs <- data$obs %>%
    select(gw, source, H1_annotation, H2_annotation, H3_annotation, group)
  
  assign(data_lst[i], data)
}

gw <- concat(list(gw15, gw20, gw22, gw34))

group_lst <- sort(unique(gw$obs$group))
# group_lst <- c("'EN-ET-SP/L6b")

zs_func <- function(x) {
  x[x > quantile(x, 0.99)] <- quantile(x, 0.99)
  return(mean(x))
}


top_lst <- list()
bottom_lst <- list()
coef_lst <- list()

for (num in 1:length(group_lst)) {
  group <- group_lst[num]
  
  if (group == "'EN-IT-L5/6") {
    gw15_sub <- gw[gw$obs$H3_annotation %in% c("'EN-IT-L5/6 Cluster 1", "'EN-IT-L5/6 Cluster 2",
                                       "'EN-IT-L5/6 Cluster 3", "'EN-IT-L5 Cluster 1",
                                       "'EN-IT-L5 Cluster 2", "'EN-IT-L5 Cluster 3")
                   & gw$obs$gw == 15]
    gw_other <- gw[gw$obs$group == group & gw$obs$H1_annotation != 15]
    gw_group <- concat(list(gw15_sub, gw_other))
  } else if (group == "'EN-IT-L4") {
    gw15_sub <- gw[gw$obs$H3_annotation %in% c("'EN-IT-L3/4-1 Cluster 1")
                   & gw$obs$gw == 15]
    gw_other <- gw[gw$obs$group == group & gw$obs$H1_annotation != 15]
    gw_group <- concat(list(gw15_sub, gw_other))
  } else if (group == "'EN-IT-L3") {
    gw15_sub <- gw[gw$obs$H3_annotation %in% c("'EN-IT-L3/4-1 Cluster 1", "'EN-IT-L5/6 Cluster 1")
                   & gw$obs$gw == 15]
    gw_other <- gw[gw$obs$group == group & gw$obs$H1_annotation != 15]
    gw_group <- concat(list(gw15_sub, gw_other))
  } else if (group == "'EN-IT-L2") {
    gw15_sub <- gw[gw$obs$H3_annotation %in% c("'EN-IT-L3/4-1 Cluster 1")
                   & gw$obs$gw == 15]
    gw_other <- gw[gw$obs$group == group & gw$obs$H1_annotation != 15]
    gw_group <- concat(list(gw15_sub, gw_other))
  } else {
    gw_group <- gw[gw$obs$group == group]
  }
  # rownames(gw_group$obs) <- 1:dim(gw_gene)[1]
  # gw_gene <- data.frame(t(scale(t(gw_group$X))))
  gw_gene <- data.frame(gw_group$X)
  rownames(gw_gene) <- 1:dim(gw_gene)[1]
  colnames(gw_gene) <- rownames(gw_group$var)
  gw_gene$week <- gw_group$obs$gw
  
  gene_avg <- gw_gene %>% group_by(week) %>% summarise(across(everything(), list(zs_func )))
  colnames(gene_avg)[2:ncol(gene_avg)] <- rownames(gw_group$var)
  reg_coef <- c()
  week_num <- 1:4
  for (i in 2:(dim(gene_avg)[2])) {
    avg_i <- as.vector(unlist(gene_avg[,i]))
    fit <- lm(avg_i ~ week_num)
    reg_coef <- c(reg_coef, fit$coefficients[2])
  }
  gene_coef_df <- data.frame(reg_coef)
  gene_coef_df$gene <- rownames(gw_group$var)
  top_gene <- gene_coef_df[order(gene_coef_df$reg_coef ),]$gene[1:10]
  bottom_gene <- gene_coef_df[order(gene_coef_df$reg_coef, decreasing = TRUE),]$gene[1:10]
  
  top_lst[[num]] <- top_gene
  bottom_lst[[num]] <- bottom_gene
  coef_lst[[num]] <- gene_coef_df
}

names(top_lst) <- group_lst
names(bottom_lst) <- group_lst
names(coef_lst) <- group_lst

saveRDS(top_lst, "topgene_lst_raw.rds")
saveRDS(bottom_lst, "bottomgene_lst_raw.rds")
