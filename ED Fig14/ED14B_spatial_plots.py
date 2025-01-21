import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
from itertools import repeat

adata = sc.read('gw21_umb1932.h5ad')

def make_plot_combined():
    types = ['EN-IT-L2|3-c1', 'EN-IT-L4-c2', 'EN-IT-L4|5-c2', 'EN-IT-L5|6-2-c1', 'EN-ET-L5|6-c1', 'EN-IT-L2|3-c2', 'EN-IT-L4-c1', 'EN-IT-L4|5-c0', 'EN-IT-L5|6-2-c2', 'EN-ET-L5|6-c0']
    colors = ['#FF97EE', '#E971DD', '#D24CCC', '#BC26BB', '#A800A8', '#ACFC64', '#81DC4B', '#56BD32', '#2B9D19', '#007E00']
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
    plt.savefig('umb1932_v1_v2_clusters.png', dpi=500); plt.clf()


def make_plot_individual(i):
    types = ['EN-IT-L2|3-c1', 'EN-IT-L4-c2', 'EN-IT-L4|5-c2', 'EN-IT-L5|6-2-c1', 'EN-ET-L5|6-c1', 'EN-IT-L2|3-c2', 'EN-IT-L4-c1', 'EN-IT-L4|5-c0', 'EN-IT-L5|6-2-c2', 'EN-ET-L5|6-c0']
    colors = ['#FF97EE', '#E971DD', '#D24CCC', '#BC26BB', '#A800A8', '#ACFC64', '#81DC4B', '#56BD32', '#2B9D19', '#007E00']
    if i<2:
      types = [types[i], types[i+5]]
      colors = [colors[i], colors[i+5]]         
    elif i==2:
        types = [types[i], types[i+1], types[i+5], types[i+6]]
        colors = [colors[i], colors[i+1], colors[i+5], colors[i+6]]
    elif i>2:
      types = [types[i+1], types[i+6]]
      colors = [colors[i+1], colors[i+6]]         
    color_dict = dict(zip(types, colors))
    unique_annotations = adata.obs['H3_annotation'].dropna().unique()
    remaining_types = np.setdiff1d(unique_annotations, types)
    color_dict.update(dict(zip(remaining_types, ['grey'] * len(remaining_types))))
    plot = sc.pl.embedding(adata, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=2, palette = color_dict, alpha=1); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types];    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = 2, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    #plot.legend(loc = 'center', fontsize=3, ncol = 1, bbox_to_anchor=(0.95,1), markerscale=0.5);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig('umb1932_v1_v2_clusters' + str(i) + '.png', dpi=500); plt.clf()




def main():
    make_plot_combined()
    for i in range(4):
        make_plot_individual(i)


if __name__=="__main__":
    main()

