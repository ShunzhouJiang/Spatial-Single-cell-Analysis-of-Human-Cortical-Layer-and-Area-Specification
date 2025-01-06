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

order = np.array(['IN-CGE4', 'IN-CGE1', 'INP-LGE1', 'IN-CGE5', 'IN-CGE2', 'IN-VZ/LGE2', 'IN-CGE3', 'IN-MGE4', 'IN-MGE2', 'INP-VZ/GE1', 'IN-SST-NPY', 'IN-VIP-late', 'INP-VZ', 'IN-GE1', 'IN-MGE5', 'IN-GE2', 'INP-LGE2', 'IN-SST3', 'IN-VZ/LGE1', 'IN-MGE1', 'IN-SST1', 'INP-VZ/GE2', 'IN-SST4', 'IN-MGE3'])

def make_plot(sample, region, area):
  obs = obs_all[(obs_all['sample']==sample) & (obs_all.region==region) & (obs_all.area==area)]
  obs = obs[obs['cortical_depth'].notna()]
  obs_en = obs[obs.H1_annotation.isin(['EN-IT', 'EN-Mig', 'EN-ET'])]
  obs_in = obs[obs.H1_annotation=='IN']
  type_rm = list(obs.value_counts('H3_annotation').index[(obs.value_counts('H3_annotation')<10)])
  if sample in gw15:
    layer_types = ['EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1']
    l4_5 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[0]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.75))/2
    l5_6 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[2]].cortical_depth, 0.75))/2
    l = [l4_5,l5_6]
  elif sample in gw20:
    layer_types = ['EN-IT-L2/3-A1', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1']
    l3_4 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[0]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.75))/2
    l4_5 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[2]].cortical_depth, 0.75))/2
    l5_6 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[2]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[3]].cortical_depth, 0.75))/2
    l = [l3_4,l4_5,l5_6]
  elif sample in gw22:
    layer_types = ['EN-L2-1', 'EN-IT-L3-A', 'EN-IT-L4-1', 'EN-ET-L5-1', 'EN-IT-L6-1']
    l2_3 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[0]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.75))/2
    l3_4 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[2]].cortical_depth, 0.75))/2
    l4_5 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[2]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[3]].cortical_depth, 0.75))/2
    l5_6 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[3]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[4]].cortical_depth, 0.75))/2
    l = [l2_3,l3_4,l4_5,l5_6]
  elif sample in gw34:
    layer_types = ['EN-L2-4', 'EN-IT-L3-late', 'EN-IT-L4-late', 'EN-ET-L5-1', 'EN-IT-L6-late']
    l2_3 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[0]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.75))/2
    l3_4 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[1]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[2]].cortical_depth, 0.75))/2
    l4_5 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[2]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[3]].cortical_depth, 0.75))/2
    l5_6 = (np.quantile(obs_en[obs_en.H3_annotation==layer_types[3]].cortical_depth, 0.25)+np.quantile(obs_en[obs_en.H3_annotation==layer_types[4]].cortical_depth, 0.75))/2
    l = [l2_3,l3_4,l4_5,l5_6]
  obs_in.H3_annotation = obs_in.H3_annotation.astype('category').cat.set_categories(order)
  if ((len(type_rm)>0) & (len(type_rm)<len(obs_in.H3_annotation.unique()))):
        obs_in_notrm = obs_in[~obs_in.H3_annotation.isin(type_rm)]
        obs_in_rm = obs_in[obs_in.H3_annotation.isin(type_rm)]
        #order = list(obs.groupby('H3_annotation').aggregate('median').sort_values(by='cp_dist').index)
        plt.figure(figsize=(35,7));
        plot = sns.violinplot(x='H3_annotation', y='cortical_depth', hue='H3_annotation', data=obs_in_notrm, order=order, palette='rainbow', density_norm='width', inner = None, dodge=False, cut=0);
        sns.stripplot(x='H3_annotation', y='cortical_depth', hue='H3_annotation', data=obs_in_rm, order=order, palette='rainbow'); 
        #plot.axhline(y=l2_3); plot.axhline(y=l3_4); plot.axhline(y=l4_5); plot.axhline(y=l5_6);
        [plot.axhline(i, linestyle = '--') for i in l];
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=9); plt.yticks(fontsize=9); plot.set(xlabel=None); plot.set_ylabel('Cortical Depth', fontsize=20); plt.ylim(0,1);
        plt.tight_layout();
        plt.savefig(sample + '-' + region + '-' + area + 'in_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0)

    elif ((len(type_rm)>0) & (len(type_rm)==len(obs_in.H3_annotation.unique()))):
        #order = list(obs.groupby('H3_annotation').aggregate('median').sort_values(by='cp_dist').index)
        plt.figure(figsize=(35,7));
        #plot = sns.violinplot(x='H3_annotation', y='cortical_depth', hue='H3_annotation', data=obs2, order=order, palette='rainbow', density_norm='width', inner = None, dodge=False, cut=0);
        plot = sns.stripplot(x='H3_annotation', y='cortical_depth', hue='H3_annotation', data=obs_in, order=order, palette='rainbow'); 
        #plot.axhline(y=l2_3); plot.axhline(y=l3_4); plot.axhline(y=l4_5); plot.axhline(y=l5_6);
        [plot.axhline(i, linestyle = '--') for i in l];
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=9); plt.yticks(fontsize=9); plot.set(xlabel=None); plot.set_ylabel('Cortical Depth', fontsize=20); plt.ylim(0,1);
        plt.tight_layout();
        plt.savefig(sample + '-' + region + '-' + area + 'in_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0)
        
    else:
        #order = list(obs.groupby('H3_annotation').aggregate('median').sort_values(by='cp_dist').index)
        plt.figure(figsize=(35,7));
        plot = sns.violinplot(x='H3_annotation', y='cortical_depth', hue='H3_annotation', data=obs_in, order=order, palette='rainbow', density_norm='width', inner = None, dodge=False, cut=0);
        #plot.axhline(y=l2_3); plot.axhline(y=l3_4); plot.axhline(y=l4_5); plot.axhline(y=l5_6);
        [plot.axhline(i, linestyle = '--') for i in l];
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=9); plt.yticks(fontsize=9); plot.set(xlabel=None); plot.set_ylabel('Cortical Depth', fontsize=20); plt.ylim(0,1);
        plt.tight_layout();
        plt.savefig(sample + '-' + region + '-' + area + 'in_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0)

samples_regions_areas = (adata.obs['sample'].astype('str') + '_' + adata.obs.region.astype('str')+ '_' + adata.obs.area.astype('str')).unique()        

samples = [i.split('_')[0] for i in samples_regions_areas]
regions = [i.split('_')[1] for i in samples_regions_areas]
areas = [i.split('_')[2] for i in samples_regions_areas]

          
def main():
  with Pool(12) as pool:
    pool.starmap(make_plot, zip(samples, regions, areas))

if __name__=="__main__":
    main()
