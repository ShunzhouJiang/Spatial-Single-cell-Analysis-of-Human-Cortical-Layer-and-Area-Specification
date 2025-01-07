import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.colors
from itertools import repeat
import matplotlib


adata = sc.read('merscope_integrated_855.h5ad')


def make_plot(sample, region, gene):
    adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region)].copy()
    sc.pp.scale(adata1, zero_center = True, max_value = 6)
    plot = sc.pl.embedding(adata1, basis="spatial", use_raw = False, color = gene, show = False, s=2, color_map = cmap, alpha=1, colorbar_loc = None);
    norm = matplotlib.colors.Normalize(vmin=adata1[:,gene].X.min(), vmax=adata1[:,gene].X.max());
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); 
    sm.set_array([]);
    plt.colorbar(sm, location = 'top', orientation = 'horizontal', label = gene, shrink = 0.3);
    plt.axis('off');
    plot.set_aspect('equal');
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(sample+'_'+region + '_'+gene+ '.png', dpi=500); plt.clf()


samples = ['FB121-F1', 'FB080-F2a']

genes = ['GAD2', 'SOX2-OT', 'SP9', 'KITLG', 'BEST3', 'NKX2-1']

cmap = 'YlGnBu'

def main():
  for sample in samples:
    with Pool(8) as pool:
      pool.starmap(make_plot, zip(repeat(sample.split('-')[0]), repeat(sample.split('-')[1]), genes))

if __name__=="__main__":
    main()

