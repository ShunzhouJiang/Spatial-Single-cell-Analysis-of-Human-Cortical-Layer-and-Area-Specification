import pandas as pd                                                    
import numpy as np                                                     
import scanpy as sc 
import os
from anndata import read_h5ad

adata = read_h5ad("../source_data/merscope_integrated_855.h5ad")
adata = adata.raw.to_adata()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, 'H3_annotation', method='t-test')
result = adata.uns['rank_genes_groups']
de_gene = pd.DataFrame(np.array(result['names'][:10]))
de_gene.to_csv("h3_marker.csv")

sc.pp.scale(adata, max_value=10)
adata_df = pd.DataFrame(adata.X)
adata_df.columns = adata.var.index
adata_df['cluster'] = adata.obs['H3_annotation'].values
adata_zs = adata_df.groupby(by='cluster').agg('mean')
adata_nz = adata_df.groupby(by='cluster').agg(lambda x: np.mean(x != np.min(x)))
os.makedirs("result", exist_ok=True)
adata_zs.to_csv("result/h3_zs.csv")
adata_nz.to_csv("result/h3_nz.csv")
