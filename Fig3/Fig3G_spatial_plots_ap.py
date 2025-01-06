import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
from itertools import repeat


adata = sc.read('merscope_integrated_855.h5ad')

def make_plot(sample, region, types, i):
    adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region)].copy()
    colors = ['red', 'blue']
    color_dict = dict(zip(types, colors))
    color_dict.update(dict(zip(np.setdiff1d(adata.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata.obs.H3_annotation.unique(), types)))))
    plot = sc.pl.embedding(adata1, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=2, palette = color_dict, alpha=1); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types];    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = 2, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(sample + '_' + region + '_' + str(i) + '.png', dpi=500); plt.clf()


samples = ['FB123-F1', 'FB123-F2', 'FB123-P1', 'FB123-O2']
types_all = [('EN-IT-L2/3-A2', 'EN-IT-L3-P'), ('EN-IT-L4-A', 'EN-IT-L4-late'), ('EN-IT-L4/5-1', 'EN-IT-L5/6-P'), ('EN-ET-L6-A', 'EN-ET-L6-P'), ('EN-ET-SP-A', 'EN-ET-SP-P1')]



def main():
  with Pool(5) as pool:
    for sample in samples:
      pool.starmap(make_plot, zip(repeat(sample.split('-')[0]), repeat(sample.split('-')[1]), types_all, range(5)))

if __name__=="__main__":
    main()

