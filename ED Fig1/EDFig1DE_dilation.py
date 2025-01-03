import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
from tqdm import tqdm
from anndata import AnnData, concat, read_h5ad
import scanpy as sc
from sklearn.cluster import KMeans
from scipy import sparse
from multiprocessing import Pool

def findBlank(name):
    return "Blank" in name

def read_data(path):
    df = pd.read_csv(f"../segment_data/{path}/cell_by_gene.csv", index_col = 0)
    df.index = df.index.astype(str)
    meta = pd.read_csv(f"../segment_data/{path}/cell_metadata.csv", index_col = 0)
    adata = AnnData(df, obs=meta)
    ind_blank = np.array(list(map(findBlank, df.columns)))
    adata = adata[:, ~ind_blank]
    return adata

dila = read_data("202303231348_FB080-BA17-3_VMSC02901")
dila_no = read_data("202303231348_FB080-BA17-3_VMSC02901_nodilation")
dila.obs['source'] = "Dilation"
dila_no.obs['source'] = "No dilation"
adata = concat([dila, dila_no], axis=0)
adata.X = sparse.csr_matrix(adata.X)
adata.write_h5ad("dila_data.h5ad")

# adata = read_h5ad("dila_data.h5ad")
adata.raw = adata.copy()

thre = np.sum(adata.X.A, axis = 1)
print(np.quantile(thre, 0.1))
sc.pp.filter_cells(adata, min_counts=np.quantile(thre, 0.1))

sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=6)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata=adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)


## KM clustering
adata.obs['Cluster'] = KMeans(30, random_state=100).fit_predict(adata.X)

sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=800, figsize=(6,6), format='png')
adata.obs['Cluster'] = adata.obs['Cluster'].astype('category')
sc.pl.umap(adata, color=['source',  "Cluster"], save="dila_comparison_cluster.png")
adata.write_h5ad("dila_cluster.h5ad")

## Compute cluster proportion bar plot (ED Fig 1d)
clusters = np.sort(np.unique(adata.obs['Cluster']))
tab = pd.crosstab(adata.obs['Cluster'], adata.obs['source'])
tab = tab/np.sum(tab,axis=1)[:, None]
tab.to_csv("cluster_prop.csv")
source_A = tab.iloc[:, 0]
source_B = tab.iloc[:, 1]

x = np.arange(len(clusters))
fig, ax = plt.subplots(figsize=(7,4))
ax.bar(x, source_A, label=tab.columns[0], color='skyblue')
ax.bar(x, source_B, bottom=source_A, label=tab.columns[1], color='salmon')
ax.set_xlabel('Clusters')
ax.set_ylabel('Proportion')
ax.set_title('Cluster proportion')
ax.set_xticks(x)
ax.set_xticklabels(clusters, rotation=90, fontsize=8)
ax.set_xticklabels(clusters)
ax.axhline(y=0.5, color='black', linestyle='--', linewidth = 0.8)
ax.legend(title="Source", bbox_to_anchor=(1.05, 0.8))
plt.tight_layout()
fig.savefig("figures/dila_prop.pdf", dpi=400, format='pdf')
plt.show()


## Compute transcript count distribution (ED Fig 1e)

def process_row(i):
    return {
        'Cluster': adata.obs['Cluster'][i],
        'Source': adata.obs['source'][i],
        'value': np.sum(adata.X[i, ])
    }


if __name__ == '__main__':
    adata = read_h5ad("dila_cluster.h5ad")
    adata.X = adata.raw.to_adata().X
    cluster_uniq = np.sort(np.unique(adata.obs['Cluster']))
    source_uniq = np.sort(np.unique(adata.obs['source']))
    with Pool(processes=4) as pool:
        data = list(tqdm(pool.imap(process_row, range(adata.shape[0])), total=adata.shape[0]))

    data = pd.DataFrame(data)
    palette = {source_uniq[0]: 'skyblue', source_uniq[1]: 'salmon'}

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.boxplot(x='Cluster', y='value', hue='Source', data=data, palette=palette, fliersize=0.3)
    plt.title('Distribution of transcriptomics counts')
    plt.xlabel('Clusters')
    plt.ylabel('Transcriptomics counts')

    ax.legend(title="Source", bbox_to_anchor=(1.05, 0.8))
    plt.tight_layout()
    fig.savefig("dila_dist.pdf", dpi=400, format='pdf')
    plt.close()