import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from itertools import compress
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar

adata = sc.read('gw34_umb5900_ba17.h5ad')

def make_plot(gw):
    types = ['EN-IT-L3|4-c2', 'EN-IT-L3|4-c3', 'EN-IT-L3|4-c4', 'EN-IT-L3|4-c5', 'EN-IT-L3|4-c6', 'EN-IT-L3|4-c7', 'EN-IT-L4-2-c2', 'EN-IT-L4-2-c3', 'EN-ET-L6-c1', 'EN-ET-L6-V1-c3', 'EN-IT-L3|4-c1', 'EN-IT-L3-c0', 'EN-IT-L4-2-c1', 'EN-ET-L5-c1']
    colors = ['#FF97EE', '#F586E6', '#EB75DE', '#E264D6', '#D853CE', '#CE43C7', '#C532BF', '#BB21B7', '#B110AF', '#A800A8', '#ACFC64', '#99EE58', '#3BA721', '#007E00']
    color_dict = dict(zip(types, colors))
    unique_annotations = adata.obs['H3_annotation'].dropna().unique()
    remaining_types = np.setdiff1d(unique_annotations, types)
    color_dict.update(dict(zip(remaining_types, ['grey'] * len(remaining_types))))
    plot = sc.pl.embedding(adata, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=2, palette = color_dict, alpha=1); plt.axis('off');
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
    make_plot('umb5900_ba17')

if __name__=="__main__":
    main()

