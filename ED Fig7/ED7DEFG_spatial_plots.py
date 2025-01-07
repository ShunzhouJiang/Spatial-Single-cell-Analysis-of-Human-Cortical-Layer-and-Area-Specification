import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
from itertools import repeat

adata = sc.read('merscope_integrated_855.h5ad')

gw15 = ['UMB1367', 'UMB1117']
gw20 = ['FB080', 'FB121']
gw22 = ['FB123']
gw34 = ['UMB5900']


def make_plot(sample, region):
    adata1 = adata[(adata.obs['sample']==sample) & (adata.obs.region==region)].copy()
    ncol=1
    if sample in gw15:
        types = ['EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1']
        colors = ['#0FF008', '#0832F0', '#F008EC']
        color_dict = dict(zip(types, colors))
        color_dict.update(dict(zip(np.setdiff1d(adata.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata.obs.H3_annotation.unique(), types)))))
    elif sample in gw20:
        types = ['EN-IT-L2/3-A1', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1']
        colors = ["#F0A608", "#0FF008", '#0832F0', '#F008EC']
        color_dict = dict(zip(types, colors))
        color_dict.update(dict(zip(np.setdiff1d(adata.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata.obs.H3_annotation.unique(), types)))))
    elif sample in gw22:
        types = ['EN-L2-1', 'EN-IT-L3-A', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1']
        colors = ["#EC2814", "#F0A608", "#0FF008", '#0832F0', '#F008EC']
        color_dict = dict(zip(types, colors))
        color_dict.update(dict(zip(np.setdiff1d(adata.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata.obs.H3_annotation.unique(), types)))))
    elif sample in gw34:
        types = ['EN-L2-4', 'EN-IT-L3-late', 'EN-IT-L4-late', 'EN-ET-L5-1', 'EN-IT-L6-late']
        colors = ["#EC2814", "#F0A608", "#0FF008", '#0832F0', '#F008EC']
        color_dict = dict(zip(types, colors))
        color_dict.update(dict(zip(np.setdiff1d(adata.obs.H3_annotation.unique(), types), ['grey']*len(np.setdiff1d(adata.obs.H3_annotation.unique(), types)))))
    plot = sc.pl.embedding(adata1, basis="spatial", color = 'H3_annotation', groups = types, show = False, s=2, palette = color_dict, alpha=1); plt.axis('off');
    plot.set_aspect('equal');
    handles, labels = plt.gca().get_legend_handles_labels();
    order = [labels.index(i) for i in types];    
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=2, ncol = ncol, bbox_to_anchor=(1.0,1.0), markerscale=0.25);
    plot.get_figure().gca().set_title('');
    scalebar = ScaleBar(1, "um", fixed_value=500, location = 'lower right');
    plot.add_artist(scalebar);
    plt.tight_layout();
    plt.savefig(sample+'_'+region + '/' + sample+'_'+region+'_'+fig+'.png', dpi=500); plt.clf()

fig = 'layer_def'


sections = ['UMB1367-O1', 'FB080-O1c', 'FB123-O2', 'UMB5900-BA18']

def main():
  for i in sections:
        make_plot(i.split('-')[0], i.split('-')[1])

if __name__=="__main__":
    main()

