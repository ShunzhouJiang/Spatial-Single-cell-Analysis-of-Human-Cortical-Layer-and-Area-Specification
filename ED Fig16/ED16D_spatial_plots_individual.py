import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib
import matplotlib.colors as mcolors
import seaborn as sns
from itertools import repeat

adata = sc.read('merscope_integrated_855.h5ad')


def make_plot(sample, region, h1_type):
    adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region)].copy()
    types = adata1.obs[adata1.obs.H1_annotation==h1_type].H3_annotation.unique()
    colors = sns.color_palette('tab20', n_colors=len(types))
    color_dict = dict(zip(types, colors))
    color_dict.update(dict(zip(np.setdiff1d(adata1.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata1.obs.H3_annotation.unique(), types)))))
    plot = sc.pl.embedding(adata1, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=2, alpha=1, palette = color_dict); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types];    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = 2, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(sample + '_' + region + '_' + h1_type + '.png', dpi=800); plt.savefig(sample + '_' + region + '_' + h1_type + '.pdf'); plt.clf()


sample = 'FB080-O1c'
h1_types = ['EN-Mig', 'IPC', 'RG']

def main():
  with Pool(3) as pool:
    pool.starmap(make_plot, zip(repeat(sample.split('-')[0]), repeat(sample.split('-')[1]), h1_types))

if __name__=="__main__":
    main()

