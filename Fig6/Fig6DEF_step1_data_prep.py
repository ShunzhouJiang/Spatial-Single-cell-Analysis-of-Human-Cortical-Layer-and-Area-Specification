import os
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
# from multiprocessing import Pool
from itertools import repeat
import warnings
import sys
from anndata import read_h5ad


adata_855 = read_h5ad("../source_data/merscope_integrated_855.h5ad")
adata_855 = adata_855[~pd.isna(adata_855.obs['area'])].copy()

def fill_na(adata):
    adata.obs['H3_annotation'] = adata.obs['H3_annotation'].astype(object)
    adata.obs['H3_annotation'] = adata.obs.apply(
        lambda row: f"{row['H2_annotation']}-c0" if pd.isna(row['H3_annotation']) and row['H1_annotation'] in ["EN-ET", "EN-IT"]
        else row['H2_annotation'] if pd.isna(row['H3_annotation'])
        else row['H3_annotation'],
        axis=1)
    return adata


os.makedirs("data_region", exist_ok=True)
## GW15 data preparation
adata_15 = read_h5ad("../source_data/gw15.h5ad")
adata_15 = adata_15.raw.to_adata()
adata_15 = adata_15[(adata_15.obs['sample'] == 'UMB1367') & (adata_15.obs['region'] == 'O1')].copy()
adata_15 = adata_15[adata_15.obs.index.isin(adata_855.obs.index)].copy()
adata_15.obs['area'] = adata_855.obs.loc[adata_15.obs.index, 'area']
adata_15 = fill_na(adata_15)
adata_15.write_h5ad("data_region/umb1367_o1_clusters.h5ad")

## GW20 data preparation
adata_20 = read_h5ad("../source_data/gw20.h5ad")
adata_20 = adata_20.raw.to_adata()
adata_20 = adata_20[(adata_20.obs['sample'] == 'FB080') & (adata_20.obs['region'] == 'O1c')].copy()
adata_20.obs.index = adata_20.obs.index.str.split('-').str[0]
adata_20 = adata_20[adata_20.obs.index.isin(adata_855.obs.index)].copy()
adata_20.obs['area'] = adata_855.obs.loc[adata_20.obs.index, 'area']
adata_20 = fill_na(adata_20)
adata_20.write_h5ad("data_region/fb080_o1c_clusters.h5ad")

## GW34 data preparation
adata_34 = read_h5ad("../source_data/gw34_umb5900_ba17.h5ad")
adata_34 = adata_34.raw.to_adata()
adata_34.obs['H1_annotation'] = adata_34.obs['H1_annotation'].replace({'EN-IT-DL': 'EN-IT', 'EN-IT-UL': 'EN-IT'})
adata_34 = fill_na(adata_34)
adata_34.write_h5ad("data_region/umb5900_ba17_clusters.h5ad")

## Adult data preparation
adata_adult = read_h5ad("../source_data/adult_umb5958.h5ad")
adata_adult = adata_adult.raw.to_adata()
adata_adult.obs['H1_annotation'] = adata_adult.obs['H1_annotation'].replace({'EN-UL': 'EN-IT', 'EN-L2': 'EN-IT', 'EN-DL': 'EN-ET'})
adata_adult = fill_na(adata_adult)
adata_adult.write_h5ad("data_region/umb5958_clusters.h5ad")
