import os
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
# from multiprocessing import Pool
from itertools import repeat
import warnings
import sys
from anndata import read_h5ad, concat
warnings.filterwarnings('ignore')

def get_data_it(gw_num):
    try:
        gw_num = int(gw_num)
    except:
        gw_num = gw_num
    if gw_num == 20:
        adata_20 = read_h5ad("data_region/fb080_o1c_clusters.h5ad")
        adata_20 = adata_20[adata_20.obs['H1_annotation'] == 'EN-IT'].copy()
        adata_20_v1 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-IT-L3-c0', 'EN-IT-L3/4-c3'])].copy()
        adata_20_v2 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-IT-L3-c3', 'EN-IT-L3/4-c0'])].copy()

        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'
        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)
    
    elif gw_num == 15:
        adata_15 = read_h5ad("data_region/umb1367_o1_clusters.h5ad")
        adata_20_v1 = adata_15[(adata_15.obs['area'].isin(['A-V1'])) & (adata_15.obs['H1_annotation'].isin(['EN-IT']))].copy()
        adata_20_v2 = adata_15[(adata_15.obs['area'].isin(['C-V2'])) & (adata_15.obs['H1_annotation'].isin(['EN-IT']))].copy()
        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'
        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)

    elif gw_num == 34:
        adata_20 = read_h5ad("data_region/umb5900_ba17_clusters.h5ad")
        adata_20_v1 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-IT-L3|4-c2', 'EN-IT-L3|4-c3', 'EN-IT-L3|4-c4', 'EN-IT-L3|4-c5', 'EN-IT-L3|4-c6', 
                                                                   'EN-IT-L3|4-c7', 'EN-IT-L4-2-c2', 'EN-IT-L4-2-c3' ])].copy()
        adata_20_v2 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-IT-L3|4-c1', 'EN-IT-L3-c0', 'EN-IT-L4-2-c1'])].copy()
        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'

        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)
        
    elif gw_num == "Adult":
        adata_20 = read_h5ad("data_region/umb5958_clusters.h5ad")
        adata_20_v1 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-L4c-c1', 'EN-L4-c0', 'EN-L4-c1', 'EN-L4-c2', 'EN-L4-c3', 'EN-L3-2-c2' ])].copy()
        adata_20_v2 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-L4c-c0', 'EN-L3-2-c1'])].copy()
        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'
        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)
        
    sc.pp.normalize_total(adata_compare, target_sum=1e4)
    sc.pp.log1p(adata_compare)
    return adata_compare


def get_data_et(gw_num):
    try:
        gw_num = int(gw_num)
    except:
        gw_num = gw_num
    if gw_num == 20:
        adata_20 = read_h5ad("data_region/fb080_o1c_clusters.h5ad")
        adata_20_v1 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-ET-L5/6-c4', 'EN-ET-SP-2-c1'])].copy()
        adata_20_v2 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-ET-L5/6-c3', 'EN-ET-SP-2-c3'])].copy()
        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'
        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)
    
    elif gw_num == 15:
        adata_15 = read_h5ad("data_region/umb1367_o1_clusters.h5ad")
        adata_20_v1 = adata_15[(adata_15.obs['area'].isin(['A-V1'])) & (adata_15.obs['H1_annotation'].isin(['EN-ET']))].copy()
        adata_20_v2 = adata_15[(adata_15.obs['area'].isin(['C-V2'])) & (adata_15.obs['H1_annotation'].isin(['EN-ET']))].copy()
        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'
        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)
        
    elif gw_num == 34:
        adata_20 = read_h5ad("data_region/umb5900_ba17_clusters.h5ad")
        adata_20_v1 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-ET-L6-c1', 'EN-ET-L6-V1-c3'])].copy()
        adata_20_v2 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-ET-L5-c1'])].copy()
        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'
        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)
        
    elif gw_num == "Adult":
        adata_20 = read_h5ad("data_region/umb5958_clusters.h5ad")
        adata_20_v1 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-L2|3-2-c0' ])].copy()
        adata_20_v2 = adata_20[adata_20.obs['H3_annotation'].isin(['EN-L2|3-1-c0'])].copy()
        adata_20_v1.obs['v1v2'] = 'v1'
        adata_20_v2.obs['v1v2'] = 'v2'
        adata_compare = concat([adata_20_v1, adata_20_v2], axis = 0)        
        
    sc.pp.normalize_total(adata_compare, target_sum=1e4)
    sc.pp.log1p(adata_compare)
    return adata_compare


def find_deg(gw_num, v1v2_type, itet = "IT", path = "DEG_time_result"):
    if itet == "IT":
        adata_compare = get_data_it(gw_num)
    elif itet == "ET":
        adata_compare = get_data_et(gw_num)
    sc.tl.rank_genes_groups(adata_compare, groupby="v1v2", method="t-test" #, groups=["v1"] , reference="v2"
    )
    result = adata_compare.uns['rank_genes_groups']
    group = v1v2_type
    de_df = pd.DataFrame({
        'gene': result['names'][group], 
        'log_fold_change': result['logfoldchanges'][group],
        'pvals': result['pvals'][group],
        'pvals_adj': result['pvals_adj'][group]
    })
    os.makedirs(path, exist_ok=True)
    de_df.to_csv(f"{path}/GW{gw_num}_{itet}_{v1v2_type.upper()}_DEG.csv")
    return de_df
