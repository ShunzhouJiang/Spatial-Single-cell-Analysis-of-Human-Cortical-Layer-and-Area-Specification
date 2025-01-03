## Raw data are downloaded from Trevino, A. E. et al. doi:10.1016/j.cell.2021.07.039 (2021).
## This file is to transform the raw data into h5ad format

from anndata import AnnData
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd

rna_meta = pd.read_csv("data_greenleaf/GSE162170_rna_cell_metadata.txt", sep="\t", index_col = 0)
ct_dict = {"c1": "CGEIN", "c8": "CycProg", "c6": "earlyRG", "c20": "EC", "c2": "GluN1", "c5": "GluN2",
           "c9": "GluN3", "c4": "GluN4", "c0": "GluN5", "c12": "GluN6", "c7": "GluN7", "c15": "GluN8",
           "c10": "lateRG", "c3": "MGEIN", "c11": "mGPC", "c16": "MG", "c14": "nIPC", "c17": "OPC_Oligo",
           "c19": "Peric", "c21": "RBC", "c13": "SP", "c18": "tRG", "c22": "VLMC"}
rna_meta['cell_type'] = rna_meta['seurat_clusters'].map(ct_dict)

rna = pd.read_csv("data_greenleaf/GSE162170_rna_counts.tsv.gz", sep='\t', compression='gzip')
gene_name = pd.read_csv("data_greenleaf/gene_name.csv", index_col=0)
gene_name = gene_name[~pd.isna(gene_name['hgnc_symbol'])]
gene_count = gene_name['ensembl_gene_id'].value_counts()
gene_count = gene_count[gene_count > 1]
gene_name = gene_name[~gene_name['ensembl_gene_id'].isin(gene_count.index)]

ind = np.isin(rna.index, gene_name['ensembl_gene_id'])
rna_filter = rna[ind]
gene_name.index = gene_name['ensembl_gene_id']
index_name = gene_name.loc[rna_filter.index, 'hgnc_symbol']
rna_filter.index = index_name


adata = AnnData(rna_filter.T)
adata.X = csr_matrix(adata.X)

adata.obs = rna_meta.loc[adata.obs.index, ["cell_type", "Sample.ID", "Age",	"Tissue.ID", "Sample.Type", "seurat_clusters"]]
gene_filter = pd.Series(adata.var.index).value_counts()[pd.Series(adata.var.index).value_counts()>1].index
adata = adata[:, ~adata.var.index.isin(gene_filter)].copy()
adata.write_h5ad("data_greenleaf/adata_greenleaf.h5ad")
