import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import seaborn as sns
from itertools import repeat

adata = sc.read('gw18_umb1759.h5ad')
adata.obs.H3_annotation = adata.obs.H3_annotation.astype('str')
adata.obs.H3_annotation[adata.obs.H3_annotation=='nan'] = adata.obs.H2_annotation[adata.obs.H3_annotation=='nan'] 
adata.obs.H3_annotation = adata.obs.H3_annotation.astype('category')

h1_types = ['EN-Mig', 'EN-ET', 'EN-IT']

def make_plot(h1_type):
    types = adata.obs[adata.obs.H1_annotation==h1_type].H3_annotation.unique()
    colors = sns.color_palette('tab20', n_colors=len(types))
    color_dict = dict(zip(types, colors))
    color_dict.update(dict(zip(np.setdiff1d(adata.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata.obs.H3_annotation.unique(), types)))))
    plot = sc.pl.embedding(adata, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=2, alpha=1, palette = color_dict); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    plt.legend(handles, labels, loc = 'center', fontsize=2, ncol = 2, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig('h3_annotation_' + h1_type + '.png', dpi=800); plt.savefig('h3_annotation_' + h1_type + '.pdf');
    plt.clf()


def main():
  with Pool(3) as pool:  
    pool.map(make_plot, h1_types)

if __name__=="__main__":
    main()

