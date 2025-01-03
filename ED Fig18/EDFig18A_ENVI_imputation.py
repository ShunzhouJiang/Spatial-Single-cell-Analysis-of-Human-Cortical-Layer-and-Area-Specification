import os
import pandas as pd                                                    
import numpy as np                                                     
import scanpy as sc                                                                                 
from time import time                                                       
import sys

import matplotlib.pyplot as plt
from anndata import AnnData, read_h5ad, concat
from tqdm import tqdm
import scipy
import scipy.stats as ss
import seaborn as sns
from scENVI import ENVI


adata_sc = read_h5ad("sc.h5ad")
adata_sc.X = adata_sc.X.A
adata_tot = read_h5ad("gw20.h5ad")
adata_tot.obs['source'] = adata_tot.obs['sample'].astype(str) + '-' + adata_tot.obs['region'].astype(str)
adata_tot = adata_tot[adata_tot.obs['source'].isin(["FB080-O1c", "FB121-F1"]) ]
adata_tot.obsm['spatial'] = np.array(adata_tot.obs[['center_x', 'center_y']])
print(adata_sc.shape)
print(adata_tot.shape)

k_nn = 100
covet_gene = 50
# sc_gene_lst = ['ABI3BP', 'PDZRN4', 'FLRT2', 'TAFA2', 'NR1D1', 'IL1RAP', 'PDZRN3', 'DLX6', 'DLX6-AS1', 'ADARB2', 'ERBB4', 'NRXN3', 'DLX2', 'ZNF536', 'PRKCA', 
#              'THRB', 'TSHZ1', 'PBX3', 'MEIS2', 'CALB2', 'CDCA7L', 'SYNPR', 'SP8', 'CASZ1', 'FOXP4']
sc_gene_lst = ['ABI3BP', 'PDZRN4', 'FLRT2', 'TAFA2', 'NR1D1', 'IL1RAP', 'CCBE1', 'THSD7B', 'CDH13', 'TRPC6', 'CHRM2', 
               'LUZP2', 'LRP1B', 'LRRTM4', 'HDAC9', 'FBXL7', 'DTNA', 'SYNDIG1', 'SDK1', 'LMO3', 'TRIQK', 'UNC13C', 'CNTNAP2', 'KCNIP4']

sc_gene_121 = ['PDZRN3', 'DLX6', 'DLX6-AS1', 'ADARB2', 'ERBB4', 'NRXN3', 'DLX2', 'ZNF536', 'PRKCA', 'THRB', 'TSHZ1', 'PBX3', 
                          'MEIS2', 'CALB2', 'CDCA7L', 'SYNPR', 'SP8', 'CASZ1', 'FOXP4']
sc_gene_lst.extend(sc_gene_121)


ENVI_Model = ENVI.ENVI(spatial_data = adata_tot, sc_data = adata_sc, num_HVG = 1000, sc_genes = sc_gene_lst, k_nearest = k_nn, num_cov_genes = covet_gene)

ENVI_Model.Train()
ENVI_Model.impute()
ENVI_Model.infer_COVET()

adata_tot.obsm['imputation'] = ENVI_Model.spatial_data.obsm['imputation']
print(ENVI_Model.spatial_data.obsm['imputation'].shape)

adata_tot.write_h5ad("mers_imputed_1000.h5ad")
