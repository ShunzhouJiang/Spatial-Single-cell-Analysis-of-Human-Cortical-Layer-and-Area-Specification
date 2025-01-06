import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar


adata = sc.read('gw20.h5ad')

def make_plot(sample, region):
    adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region)].copy()
    types = ['EN-IT-L3-c0', 'EN-IT-L3/4-c3', 'EN-ET-L5/6-c4', 'EN-ET-SP-2-c1', 'EN-IT-L3-c3', 'EN-IT-L3/4-c0', 'EN-ET-L5/6-c3', 'EN-ET-SP-2-c3']
    colors = ['#FF97EE', '#FF00FF', '#BC00BC', '#A800A8', '#ACFC64', '#00FF00', '#00C800', '#007E00']
    color_dict = dict(zip(types, colors))
    unique_annotations = adata.obs['H3_annotation'].dropna().unique()
    remaining_types = np.setdiff1d(unique_annotations, types)
    color_dict.update(dict(zip(remaining_types, ['grey'] * len(remaining_types))))
    plot = sc.pl.embedding(adata1, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=2, palette = color_dict, alpha=1); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types];    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = 2, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(sample + '_' + region + '.png', dpi=500); plt.clf()


sections = ['FB080-O1c', 'FB080-O1b', 'FB080-O1d', 'FB121-O1']


def main():
  with Pool(4) as pool:
    pool.starmap(make_plot, zip([i.split('-')[0] for i in sections], [i.split('-')[1] for i in sections]))

if __name__=="__main__":
    main()

