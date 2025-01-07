import os
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
# from multiprocessing import Pool
from itertools import repeat
import warnings
import seaborn as sns
import sys
from anndata import read_h5ad, concat
warnings.filterwarnings('ignore')


def plot_data():
    adata_15 = read_h5ad("data_region/umb1367_o1_clusters.h5ad")
    adata_15.obs['area'] = adata_15.obs['area'].replace({'C-V2': 'B-V2', 'B-V2': 'B-V2_prev'})
    adata_20 = read_h5ad("data_region/fb080_o1c_clusters.h5ad")
    adata_34 = read_h5ad("data_region/umb5900_ba17_clusters.h5ad")
    adata_adult = read_h5ad("data_region/umb5958_clusters.h5ad")
    return adata_15, adata_20, adata_34, adata_adult
    

def preprocess(adata_15, adata_20, adata_34, adata_adult, itet = "EN-IT"):
    adata_15 = adata_15[adata_15.obs['H1_annotation'] == itet].copy()
    sc.pp.normalize_total(adata_15, target_sum=1e4)
    sc.pp.log1p(adata_15)

    adata_20 = adata_20[adata_20.obs['H1_annotation'] == itet].copy()
    sc.pp.normalize_total(adata_20, target_sum=1e4)
    sc.pp.log1p(adata_20)

    adata_34 = adata_34[adata_34.obs['H1_annotation'] == itet].copy()
    sc.pp.normalize_total(adata_34, target_sum=1e4)
    sc.pp.log1p(adata_34)

    adata_adult = adata_adult[adata_adult.obs['H1_annotation'] == itet].copy()
    sc.pp.normalize_total(adata_adult, target_sum=1e4)
    sc.pp.log1p(adata_adult)
    return adata_15, adata_20, adata_34, adata_adult


def plot_gene(adata_15, adata_20, adata_34, adata_adult, itet = "EN-IT", 
              v1v2_type = "V1", top_num = False, plot_inv_gene = False, deg_path = "DEG_time_result", deg_plot_path = "DEG_time_plot"):
    if itet == "EN-IT":
        name = "IT"
    elif itet == "EN-ET":
        name = "ET"
    gw_lst = ["15", "20", "34", "Adult"]
    adata_tot = [adata_15, adata_20, adata_34, adata_adult]
    lfc_tot = []
    os.makedirs(f"{deg_plot_path}/DEG_list", exist_ok=True)
    
    for gw in gw_lst:
        de_df = pd.read_csv(f"{deg_path}/GW{gw}_{name}_{v1v2_type}_DEG.csv", index_col = 0)
        de_df = de_df[(de_df['log_fold_change'] > 0.5) & (de_df['pvals_adj'] <= 0.05)]
        if top_num:
            de_df = de_df.iloc[:top_num, :]     ## ranked by adjusted p-value
        else:
            de_df[['gene']].to_csv(f"{deg_plot_path}/DEG_list/GW{gw}_{name}_{v1v2_type}_DEG_list.csv")
        gene_lst = de_df['gene'].values
        gene_lst = np.intersect1d(gene_lst, adata_20.var.index)
        print(len(gene_lst))

        v1_dict = {}
        v2_dict = {}
        lfc_dict = {}

        for i in range(len(gene_lst)):
            v1_expr = []
            v2_expr = []
            lfc_gene = []
            for j in range(len(adata_tot)):
                adata = adata_tot[j].copy()
                adata_v1 = adata[(adata.obs['H1_annotation'] == itet) & (adata.obs['area'] == 'A-V1')].copy()
                adata_v2 = adata[(adata.obs['H1_annotation'] == itet) & (adata.obs['area'] == 'B-V2')].copy()
                
                v1_expr.append(np.mean(adata_v1[:, gene_lst[i]].X))
                v2_expr.append(np.mean(adata_v2[:, gene_lst[i]].X))
                
                v1_val = np.max([np.mean(adata_v1[:, gene_lst[i]].X), 1e-4])
                v2_val = np.max([np.mean(adata_v2[:, gene_lst[i]].X), 1e-4])
                if v1v2_type == "V1":
                    lfc_gene.append(np.log2( v1_val / v2_val ))
                elif v1v2_type == "V2":
                    lfc_gene.append(np.log2( v2_val / v1_val ))
                
            v1_dict[gene_lst[i]] = v1_expr
            v2_dict[gene_lst[i]] = v2_expr
            lfc_dict[gene_lst[i]] = lfc_gene
        v2_tot = np.mean([values for values in list(lfc_dict.values())], axis=0)
        lfc_tot.append(v2_tot)

        path = f"{deg_plot_path}/{name}_{v1v2_type}/GW{gw}"
        os.makedirs(path, exist_ok=True)

        if plot_inv_gene:
            for gene in gene_lst:
                values1 = v1_dict[gene]
                values2 = v2_dict[gene]
                indices = np.arange(len(values1))  # Create x positions

                width = 0.35
                fig = plt.figure(figsize=(5,5))
                # Create the bar plot
                plt.bar(indices - width / 2, values1, width, label="V1", alpha=0.8, color="magenta")
                plt.bar(indices + width / 2, values2, width, label="V2", alpha=0.8, color="green")

                plt.title(f"{gene}")
                plt.xlabel("Samples")
                plt.ylabel("Mean expression")
                plt.xticks(indices, ["GW15", "GW20", "GW34", "Adult"])
                plt.tight_layout()
                fig.savefig(f"{path}/{gene}.pdf", format="pdf", dpi = 300)
                plt.close()
    return lfc_tot


def plot_tot(itet = "EN-IT", v1v2_type = "V1", top_num = False, plot_inv_gene = False, deg_path = "DEG_time_result", deg_plot_path = "DEG_time_plot"):
    adata_15, adata_20, adata_34, adata_adult = plot_data()
    adata_15, adata_20, adata_34, adata_adult = preprocess(adata_15, adata_20, adata_34, adata_adult, itet=itet)
    lfc_tot = plot_gene(adata_15, adata_20, adata_34, adata_adult, itet=itet, v1v2_type = v1v2_type, top_num = top_num, 
                        plot_inv_gene = plot_inv_gene, deg_path = deg_path, deg_plot_path = deg_plot_path)
    return lfc_tot
