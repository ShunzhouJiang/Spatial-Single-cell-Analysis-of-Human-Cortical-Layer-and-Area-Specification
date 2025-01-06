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
    #vmin = adata1.X.min()
    #vmax = adata1.X.max()
    plot = sc.pl.embedding(adata1, basis="spatial", use_raw = False, color = gene, show = False, s=2, color_map = cmap, alpha=1, colorbar_loc = None);
    #plot = sc.pl.embedding(adata1, basis="spatial", color = gene, show = False, s=2, color_map = cmap, alpha=1, vmin=vmin, vmax=vmax, colorbar_loc = None);
    #norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax);
    norm = matplotlib.colors.Normalize(vmin=adata1[:,gene].X.min(), vmax=adata1[:,gene].X.max());
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); 
    sm.set_array([]);
    plt.colorbar(sm, location = 'top', orientation = 'horizontal', label = gene, shrink = 0.3);
    plt.axis('off');
    plot.set_aspect('equal');
    #handles, labels = plt.gca().get_legend_handles_labels();
    #order = [labels.index(i) for i in types];    
    #plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = ncol, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    #plot.legend(loc = 'center', fontsize=3, ncol = 1, bbox_to_anchor=(0.95,1), markerscale=0.5);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(sample+'_'+region + '_'+gene+ '.png', dpi=500); plt.clf()


samples = ['FB123-F1', 'FB123-F2', 'FB123-P1', 'FB123-O2', 'FB121-F1', 'UMB1117-F1a', 'UMB5900-BA9', 'FB080-O1c']

genes = ['CBLN2', 'SRM', 'RASGRF2', 'B3GALT2', 'SYT6', 'ETV1', 'PENK', 'CUX2', 'GLRA3', 'CUX1', 'RORB', 'PTK2B', 'FOXP1', 'TOX', 'FOXP2', 'TLE4','CALN1', 'SYBU', 'CPNE8', 'CYP26A1',
         'STK32B', 'VSTM2L', 'CRYM', 'GNAL', 'PCDH17', 'FSTL5', 'NEUROG2', 'SLN', 'SOX5', 'TAFA1', 'ARFGEF3', 'OPCML', 'NEFM', 'NFIB', 'PPP3CA','B3GNT2', 'SORCS1', 'TRPM3', 'LPL']


cmap = 'YlGnBu'
def main():
  for sample in samples:
    with Pool(8) as pool:
      pool.starmap(make_plot, zip(repeat(sample.split('-')[0]), repeat(sample.split('-')[1]), genes))

if __name__=="__main__":
    main()

