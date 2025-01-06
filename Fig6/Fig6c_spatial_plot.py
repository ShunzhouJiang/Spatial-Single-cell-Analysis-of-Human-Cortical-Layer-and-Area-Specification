import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from itertools import compress
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar


adata = sc.read('adult_umb5958.h5ad')

def make_plot(gw):
    types = ['EN-L2|3-2-c0', 'EN-L3-2-c2', 'EN-L4-c0', 'EN-L4-c1', 'EN-L4-c2', 'EN-L4-c3', 'EN-L4c-c1', 'EN-L2|3-1-c0', 'EN-L3-2-c1', 'EN-L4c-c0']
    colors = ['#ff00ff', '#F07DC8', '#E068D6', '#D14BB2', '#C23CA0', '#B32A92', '#A800A8', '#2CFF00', '#8FE753', '#007E00']
    color_dict = dict(zip(types, colors))
    unique_annotations = adata.obs['H3_annotation'].dropna().unique()
    remaining_types = np.setdiff1d(unique_annotations, types)
    color_dict.update(dict(zip(remaining_types, ['grey'] * len(remaining_types))))
    plot = sc.pl.embedding(adata, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=3, palette = color_dict, alpha=1); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types];    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = 2, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(gw + '_v1_v2_clusters.png', dpi=500); plt.clf()



def main():
      make_plot('umb5958')

if __name__=="__main__":
    main()

