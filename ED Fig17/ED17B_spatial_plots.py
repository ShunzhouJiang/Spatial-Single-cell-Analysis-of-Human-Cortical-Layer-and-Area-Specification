import numpy as np
import cv2
from PIL import Image
import scanpy as sc
import sys
import cv2
from scipy.ndimage.morphology import binary_fill_holes
import warnings
import matplotlib.pyplot as plt
import os
from itertools import compress
import glob
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import pandas as pd
import seaborn as sns
import matplotlib.colors
from itertools import repeat


adata = sc.read('/Users/kylecoleman/data/walsh/all/tif/umb5958/umb5958_clusters.h5ad')
#adata.obs.H3_cluster[~np.isnan(adata.obs.H3_cluster)] = adata.obs.H3_cluster[~np.isnan(adata.obs.H3_cluster)].astype('int').astype('str')
#adata.obs['H3'] = adata.obs[['H2_annotation', 'H3_cluster']].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
adata.obs.H3_annotation = adata.obs.H3_annotation.astype('str')
ncol=2

def make_plot(gw, i, s):
    if os.path.exists('/Users/kylecoleman/data/walsh/all/tif/'+ gw + '/increase_dimensions.txt'):
      d = {}
      with open('/Users/kylecoleman/data/walsh/all/tif/'+ gw + '/increase_dimensions.txt') as f:
        for line in f:
          (key, val) = line.split(':')
          d[key] = int(val)
      if 'height' in d:
        height = d['height']
      else:
        height=0
      if 'width' in d:
        width=d['width']
      else: 
        width=0
    else:
      height=0
      width=0

    types = ['EN-L2|3-2-c0', 'EN-L3-2-c2', 'EN-L4-c0', 'EN-L4-c1', 'EN-L4-c2', 'EN-L4-c3', 'EN-L4c-c1', 'EN-L2|3-1-c0', 'EN-L3-2-c1', 'EN-L4c-c0']
    colors = ['#FF97EE', '#F07DC8', '#E068D6', '#D14BB2', '#C23CA0', '#B32A92', '#A800A8', '#ACFC64', '#8FE753', '#007E00']
    if i<2:
      types = [types[i], types[i+7]]
      colors = [colors[i], colors[i+7]]
    elif i==2:
        types = types[i:(i+4)]
        colors = colors[i:(i+4)]

    elif i>2:
      types = [types[i+3], types[i+6]]
      colors = [colors[i+3], colors[i+6]]         

    color_dict = dict(zip(types, colors))
    color_dict.update(dict(zip(np.setdiff1d(adata.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata.obs.H3_annotation.unique(), types)))))
    plot = sc.pl.embedding(adata, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=s, palette = color_dict, alpha=1); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types];    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = ncol, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    #plot.legend(loc = 'center', fontsize=3, ncol = 1, bbox_to_anchor=(0.95,1), markerscale=0.5);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(gw + '_v1_v2_clusters' + str(i) + '_s' + str(s) + '.png', dpi=500); plt.clf()


def main():
    
      make_plot('umb5958', s)

def main():
    for s in [3,4,5]:
      with Pool(4) as pool:
        pool.starmap(make_plot, zip(repeat('umb5958'), range(4), repeat(s)))

if __name__=="__main__":
    main()

