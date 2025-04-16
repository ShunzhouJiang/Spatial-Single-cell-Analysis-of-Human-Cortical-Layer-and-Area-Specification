import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.colors
from itertools import repeat
import matplotlib


adata = sc.read('gw20.h5ad')

def make_plot(sample, region, gene, cutoff, cmap):
    adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region)].copy()
    sc.pp.scale(adata1, zero_center=True, max_value=cutoff)
    vmin = adata1.X.min()
    vmax = adata1.X.max()
    plot = sc.pl.embedding(adata1, basis="spatial", use_raw = False, color = gene, show = False, s=2, color_map = cmap, alpha=1, vmin=vmin, vmax=vmax, colorbar_loc = None);
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
    plt.savefig(sample + '_' + region + '_'+gene+ '.png', dpi=500); plt.clf()


#sections = ['FB080-O1c', 'FB080-O1b', 'FB080-O1d', 'FB121-O1']

genes = ['NEFM', 'ETV1', 'UNC5C', 'PDE1A', 'NPY', 'IGFBPL1', 'TSHZ3', 'GABRA5']
cutoffs = [6,6,8,6,8,6,6,6]

cmap = ['YlGnBu']*8

def main():
    with Pool(8) as pool:
      pool.starmap(make_plot, zip(repeat('FB080'), repeat('O1c'), genes, cutoffs, cmap))

if __name__=="__main__":
    main()

