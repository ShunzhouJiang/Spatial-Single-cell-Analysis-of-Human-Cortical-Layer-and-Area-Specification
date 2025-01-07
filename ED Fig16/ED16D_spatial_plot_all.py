import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib
import matplotlib.colors as mcolors
from itertools import repeat

adata = sc.read('merscope_integrated_855.h5ad')

def make_plot(sample, region):
    adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region)].copy()
    types = ['EN-IT-L3/4-P', 'EN-IT-L4-V1']
    types1 = list(set(adata1.obs.H3_annotation.unique()) - set(types))
    types_all = [types[0], types[1]] + types1
    types_all1 = types1 + [types[0], types[1]]
    colormap = plt.cm.rainbow
    num_colors = len(types1)
    colors = [mcolors.to_hex(colormap(i)) for i in np.linspace(0, 1, num_colors)]
    cmap1 = dict(zip(types1, colors))
    cmap1[types[0]] = '#00FF00'
    cmap1[types[1]] = '#FF00FF'
    fig, ax1 = plt.subplots(figsize=(8, 6))
    plot = sc.pl.embedding(adata1, basis="spatial", color = 'H3_annotation', groups = types1, na_color='white', palette = cmap1, show = False, s=0.5, alpha=0, ax=ax1); plt.axis('off');
    mask1 = adata1.obs['H3_annotation'].isin(types1);
    sc.pl.embedding(adata1[mask1], basis="spatial", color = 'H3_annotation', palette = cmap1, show = False, s=0.5, alpha=1, ax=ax1);
    mask2 = adata1.obs['H3_annotation'] == types[0];
    sc.pl.embedding(adata1[mask2], basis="spatial", color = 'H3_annotation', palette = cmap1, show = False, s=0.5, alpha=1, ax=ax1);
    mask3 = adata1.obs['H3_annotation'] == types[1];
    sc.pl.embedding(adata1[mask3], basis="spatial", color = 'H3_annotation', palette = cmap1, show = False, s=0.5, alpha=1, ax=ax1);  
    plot.set_aspect('equal');
    plot.legend().remove();
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types_all];
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.savefig(sample+'_'+region+'.png', dpi=800); plt.savefig(sample+'_'+region+'.pdf');
    plt.clf()
    legend = plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=3, ncol = 4, bbox_to_anchor=(0.5, 0.5), markerscale=0.5);
    plt.axis('off');
    legend.figure.savefig('legend.png', dpi = 400, bbox_inches='tight')



sample = 'FB080'
region = 'O1c'

def main():
  make_plot(sample, region)

if __name__=="__main__":
    main()

