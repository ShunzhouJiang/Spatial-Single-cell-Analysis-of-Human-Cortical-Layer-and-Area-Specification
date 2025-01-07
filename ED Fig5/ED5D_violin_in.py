import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np


adata = sc.read('merscope_integrated_855.h5ad')
obs_all = adata.obs
gw15 = ['UMB1367', 'UMB1117']
gw20 = ['FB080', 'FB121']
gw22 = ['FB123']
gw34 = ['UMB5900']

order = np.array(['IN-VZ/LGE1', 'IN-CGE5', 'IN-GE2', 'IN-VZ/LGE2', 'INP-VZ',
       'INP-VZ/GE1', 'IN-MGE3', 'IN-CGE4', 'IN-GE1', 'IN-CGE2',
       'INP-LGE1', 'INP-LGE2', 'IN-CGE3', 'IN-MGE2', 'IN-MGE4', 'IN-MGE1',
       'IN-SST4', 'IN-SST2', 'INP-VZ/GE2', 'IN-CGE1', 'IN-SST3',
       'IN-SST1', 'IN-SST-NPY', 'IN-MGE5', 'IN-VIP-late'])


def make_plot(sample, region, area):
      obs = obs_all[(obs_all['sample']==sample) & (obs_all.region==region) & (obs_all.area==area)]
      obs = obs[obs['relative_height'].notna()]
      obs = obs[obs.H1_annotation=='IN']
      obs.H3_annotation = obs.H3_annotation.astype('category').cat.set_categories(order)
      type_rm = list(obs.value_counts('H3_annotation').index[(obs.value_counts('H3_annotation')<10)])
      if len(type_rm)>0:
        obs1 = obs[~obs.H3_annotation.isin(type_rm)]
        obs2 = obs[obs.H3_annotation.isin(type_rm)]
        plt.figure(figsize=(20,5));
        plot = sns.violinplot(x='H3_annotation', y='relative_height', hue='H3_annotation', data=obs1, order=order, palette='rainbow', density_norm='width', inner = None, dodge=False, cut=0);
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=12); plt.yticks(fontsize=12); plot.set(xlabel=None); plot.set_ylabel('Laminar Depth', fontsize=20); plt.ylim(0,1);
        sns.stripplot(x='H3_annotation', y='relative_height', hue='H3_annotation', data=obs2, order=order, palette='rainbow');
        plt.savefig(sample + '-' + region + '-' + area + 'in_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0)l; plt.clf()     
      else:
        plt.figure(figsize=(20,5));
        plot = sns.violinplot(x='H3_annotation', y='relative_height', hue='H3_annotation', data=obs, order=order, palette='rainbow', density_norm='width', inner = None, dodge=False, cut=0);
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=12); plt.yticks(fontsize=12); plot.set(xlabel=None); plot.set_ylabel('Laminar Depth', fontsize=20); plt.ylim(0,1);
        plt.savefig(sample + '-' + region + '-' + area + 'in_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0); plt.clf()




samples_regions_areas = ['FB121_F1_A-PFC', 'FB080_O1c_B-V2', 'UMB1117_F2a_A-PFC']        

samples = [i.split('_')[0] for i in samples_regions_areas]
regions = [i.split('_')[1] for i in samples_regions_areas]
areas = [i.split('_')[2] for i in samples_regions_areas]

          
def main():
  with Pool(12) as pool:
    pool.starmap(make_plot, zip(samples, regions, areas))

if __name__=="__main__":
    main()
