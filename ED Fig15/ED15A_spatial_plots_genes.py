import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib
from itertools import repeat

#gw22_umb1932.h5ad has a typo, experiment is from GW21
adata = sc.read('gw22_umb1932.h5ad')

def make_plot(gene):
    sc.pp.scale(adata, zero_center=True, max_value=6)
    vmin = adata.X.min()
    vmax = adata.X.max()
    cmap = 'YlGnBu'
    plot = sc.pl.embedding(adata, basis="spatial", color = gene, show = False, s=3, color_map = cmap, alpha=1, vmin=vmin, vmax=vmax, colorbar_loc = None);
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); 
    sm.set_array([]);
    plt.colorbar(sm, location = 'top', orientation = 'horizontal', label = gene, shrink = 0.3);
    plt.axis('off');
    plot.set_aspect('equal');
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig('umb1932_' + gene+ '.png', dpi=500); plt.clf()



genes = ['NPY', 'IGFBPL1', 'TSHZ3', 'GABRA5', 'NEFM', 'ETV1', 'UNC5C', 'PDE1A', 'ABI3BP', 'PDZRN4', 'TAFA2', 'THSD7B', 'IL1RAP', 'FLRT2', 'CDH13', 'TRPC6', 'CHRM2']


def main():
    with Pool(3) as pool:
      pool.map(make_plot, genes)

if __name__=="__main__":
    main()

