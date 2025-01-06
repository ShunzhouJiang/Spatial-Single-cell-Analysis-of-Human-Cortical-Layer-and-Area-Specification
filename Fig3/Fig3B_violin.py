import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool

adata = sc.read('merscope_integrated_855.h5ad')
adata.X = np.exp(adata.X.A)-1
sc.pp.normalize_total(adata)

def find_genes(sample, region, area):
  adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region) & (adata.obs.area==area) & (~adata.obs['cortical_depth'].isna()) & (adata.obs.H1_annotation.isin(['EN-ET', 'EN-IT', 'EN-Mig']))].copy() 
  ge = adata1.X
  ge1 = ge.round()
  ge1 = ge1.astype('int')
  dist_all = [np.concatenate([adata1.obs.cortical_depth[i]*np.ones(ge1[:,j][i]) for i in range(len(adata1.obs.cortical_depth)) if ge1[:,j][i]>0]) for j in range(ge1.shape[1])]
  dist40 = [dist_all[i] for i in np.array([np.quantile(j,0.75)-np.quantile(j,0.25) for j in dist_all]).argsort()[:40]]
  genes = adata1.var.index[np.array([np.quantile(i,0.75)-np.quantile(i,0.25) for i in dist_all]).argsort()[:40]]
  dict1 = dict(zip(genes, dist40))
  genes1 = [[i]*len(dict1[i]) for i in dict1.keys()]
  genes1 = [x for xs in genes1 for x in xs]
  df1 = pd.DataFrame(genes1)
  df1.columns = ['gene']
  df1['cortical_depth'] = np.hstack(dict1.values())
  genes2 = list(df1.groupby('gene').aggregate('median').sort_values(by='cortical_depth').index)
  return genes2

def make_violin(sample, region, area, genes):
  adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region) & (adata.obs.area==area) & (~adata.obs['cortical_depth'].isna()) & (adata.obs.H1_annotation.isin(['EN-ET', 'EN-IT', 'EN-Mig']))].copy() 
  adata1 = adata1[:,genes].copy()
  ge = adata1.X
  ge1 = ge.round()
  ge1 = ge1.astype('int')
  dist40 = [np.concatenate([adata1.obs.cortical_depth[i]*np.ones(ge1[:,j][i]) for i in range(len(adata1.obs.cortical_depth)) if ge1[:,j][i]>0]) for j in range(ge1.shape[1])]
  #dist40 = [dist_all[i] for i in
  #dist40 = [dist_all[i] for i in np.array([np.quantile(j,0.75)-np.quantile(j,0.25) for j in dist_all]).argsort()[:40]]  
  #genes = adata1.var.index[np.array([np.quantile(i,0.75)-np.quantile(i,0.25) for i in dist_all]).argsort()[:40]]
  dict1 = dict(zip(genes, dist40))
  genes1 = [[i]*len(dict1[i]) for i in dict1.keys()]
  genes1 = [x for xs in genes1 for x in xs]
  df1 = pd.DataFrame(genes1)
  df1.columns = ['gene']
  df1['cortical_depth'] = np.hstack(dict1.values())
  layer_types = ['EN-L2-1', 'EN-IT-L3-A', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1']
  l2_3 = (np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[0]].cortical_depth, 0.25)+np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[1]].cortical_depth, 0.75))/2
  l3_4 = (np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[1]].cortical_depth, 0.25)+np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[2]].cortical_depth, 0.75))/2
  l4_5 = (np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[2]].cortical_depth, 0.25)+np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[3]].cortical_depth, 0.75))/2
  l5_6 = (np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[3]].cortical_depth, 0.25)+np.quantile(adata1.obs[adata1.obs.H3_annotation==layer_types[4]].cortical_depth, 0.75))/2
  l = [l2_3,l3_4,l4_5,l5_6]
  order = genes
  plt.figure(figsize=(25,5));
  plot = sns.violinplot(x='gene', y='cortical_depth', hue='gene', data=df1, order=order, density_norm='width', inner = None, dodge=False, cut=0); plot.legend().remove();
  [plot.axhline(i, linestyle = '--') for i in l];
  plt.xticks(rotation=90, fontsize=9); plt.yticks(fontsize=9); plot.set(xlabel=None); plot.set_ylabel('Cortical Depth', fontsize=20); plt.ylim(0,1); plt.tight_layout();
  #plt.savefig(sample + '_' + region + '_' + area + '_gene_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0)
  plt.show()

def main():
  genes = find_genes('FB123', 'F1', 'A-PFC')
  make_violin('FB123', 'F1', 'A-PFC', genes)
  make_violin('FB123', 'O2', 'A-Occi', genes)

if __name__=="__main__":
    main()

  
