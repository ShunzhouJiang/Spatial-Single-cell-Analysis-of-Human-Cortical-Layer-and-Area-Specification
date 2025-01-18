import scanpy as sc
import pandas as pd
import umap
import seaborn as sns
import matplotlib.pyplot as plt

adata = sc.read('merscope_integrated_855.h5ad')
adata1 = adata[(adata.obs.H1_annotation.isin(['EN-ET', 'EN-IT'])) & ~(pd.isnull(adata.obs.area))].copy()
adata1.obs['section'] = adata1.obs['sample'].astype('str') + '-' + adata1.obs.region.astype('str')

grouped = adata1.obs.groupby(["section", "H3_annotation"]).size().reset_index(name="count")
grouped["proportion"] = grouped.groupby("section")["count"].transform(lambda x: x / x.sum())
result = grouped.pivot(index="section", columns="H3_annotation", values="proportion").fillna(0)

umap_model = umap.UMAP()
umap_coords = umap_model.fit_transform(result)
result['UMAP1'] = umap_coords[:, 0]
result['UMAP2'] = umap_coords[:, 1]
result['area'] = result.index
sns.scatterplot(data=result,x='UMAP1',y='UMAP2',hue='area',palette='rainbow',alpha=0.8); plt.savefig('umap_area_en.png'); plt.clf()


