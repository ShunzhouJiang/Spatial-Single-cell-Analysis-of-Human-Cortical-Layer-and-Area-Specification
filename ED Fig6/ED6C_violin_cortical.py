import scanpy as sc
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

h2_colors = ['#92D050','#99D35B','#A0D666','#ADDB7C','#C8E6A7','#00B0F0','#10B5F1','#D8ADFF','#40C3F4','#7FD6F7','#00B050','#10B55B','#20BA66','#40C37C','#7FD6A7','#FFC000','#FFC820','#FFCF40','#FFD760','#FFDE7F','#209AD2','#FEAAA8','#B870FF','#C8E61B','#FF0000','#CC99FF','#7030A0','#8C57B2','#305496','#4A69A3','#647EB0','#7E93BD','#97A8CA']
h2_types = ['RG1', 'oRG1','Astro-late1','tRG','vRG-late',"EN-ET-SP-early",'EN-ET-SP-P','EN-ET-L5/6','EN-ET-L6-early','EN-ET-SP','IPC-oSVZ','IPC-SVZ-1','IPC-iSVZ', 'IPC-SVZ-2', 'IPC-VZ/SVZ','EN-IZ-1','EN-L2','EN-IZ-2','EN-oSVZ-1','En-oSVZ-2','EN-IT-L6','EN-IT-L4','EN-IT-L4/5','EN-IT-L3/4','EN-IT-L2/3','EC','Astro-1', 'OPC','IN-SST', 'IN-CGE', 'INP-VZ/GE', 'IN-VZ/GE', 'IN-MGE']
h2_dict = dict(zip(h2_types, h2_colors))


order = np.array(['EN-ET-L6-early5', 'EN-ET-SP-P2', 'EN-ET-SP-4', 'EN-ET-L6-early3',
       'EN-ET-SP-A', 'EN-ET-SP-5', 'EN-ET-SP-2', 'EN-ET-L6-early2',
       'EN-ET-SP-3', 'EN-ET-L6-early1', 'EN-ET-SP-early4', 'EN-ET-SP-P1',
       'EN-ET-SP-early3', 'EN-ET-SP-early5', 'EN-ET-SP-early2',
       'EN-ET-SP-early1', 'EN-ET-L6-early4', 'EN-ET-SP-1',
       'EN-ET-SP-V1T2', 'EN-ET-SP-V1T1', 'EN-ET-L6-A', 'EN-IT-L6-late',
       'EN-IT-L6-1', 'EN-ET-L6-V1', 'EN-IT-L6-2', 'EN-IT-L5/6-P',
       'EN-ET-L5/6', 'EN-ET-L6-P', 'EN-IT-L5-1', 'EN-IT-L4/5-early',
       'EN-IT-L4/5-1', 'EN-IT-Hip', 'EN-ET-L5-1', 'EN-IT-L4/5-late',
       'EN-IT-L3/4-early', 'EN-IT-L4-1', 'EN-IT-L3-P', 'EN-IT-L3/4-P',
       'EN-IT-L3/4-P2', 'EN-IT-L3/4-T', 'EN-IT-L4-A', 'EN-IT-L4-late',
       'EN-oSVZ-1', 'EN-IT-L4-V1', 'EN-IT-L3-A', 'EN-IT-L3/4-1',
       'EN-IT-L2/3-A1', 'EN-IT-L3-late', 'EN-IT-L3-V1', 'EN-IZ-2',
       'EN-IT-L2/3-A2', 'EN-L2-3', 'En-oSVZ-2', 'EN-IZ-3', 'EN-L2-1',
       'EN-IZ-1', 'EN-L2-2', 'EN-L2-4'])

def make_plot(sample, region, area):
  obs = obs_all[(obs_all['sample']==sample) & (obs_all.region==region) & (obs_all.area==area)]
  obs = obs[obs['cortical_depth'].notna()]
  obs_en = obs[obs.H1_annotation.isin(['EN-IT', 'EN-Mig', 'EN-ET'])]
  type_rm = list(obs_en.value_counts('H3_annotation').index[(obs_en.value_counts('H3_annotation')<50)])
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
  obs_en.H3_annotation = obs_en.H3_annotation.astype('category').cat.set_categories(order)
  if len(type_rm)>0:
        obs_en_notrm = obs_en[~obs_en.H3_annotation.isin(type_rm)]
        obs_en_rm = obs_en[obs_en.H3_annotation.isin(type_rm)]
        plt.figure(figsize=(35,7));
        plot = sns.violinplot(x='H3_annotation', y='cortical_depth', hue='H2_annotation', data=obs_en_notrm, order=order, palette=h2_dict, density_norm='width', inner = None, dodge=False, cut=0);
        sns.stripplot(x='H3_annotation', y='cortical_depth', hue='H2_annotation', data=obs_en_rm, order=order, palette=h2_dict); 
        [plot.axhline(i, linestyle = '--') for i in l];
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=9); plt.yticks(fontsize=9); plot.set(xlabel=None); plot.set_ylabel('Cortical Depth', fontsize=20); plt.ylim(0,1);
        plt.tight_layout();
        plt.savefig(sample + '-' + region + '-' + area + '_in_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0)
  else:
        plt.figure(figsize=(35,7));
        plot = sns.violinplot(x='H3_annotation', y='cortical_depth', hue='H2_annotation', data=obs_en, order=order, palette=h2_dict, density_norm='width', inner = None, dodge=False, cut=0);
        [plot.axhline(i, linestyle = '--') for i in l];
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=9); plt.yticks(fontsize=9); plot.set(xlabel=None); plot.set_ylabel('Cortical Depth', fontsize=20); plt.ylim(0,1);
        plt.tight_layout();
        plt.savefig(sample + '-' + region + '-' + area + '_in_violin.png', dpi=200, bbox_to_inches = 'tight', pad_inches=0)

samples_regions_areas = (adata.obs['sample'].astype('str') + '_' + adata.obs.region.astype('str')+ '_' + adata.obs.area.astype('str')).unique()        

samples = [i.split('_')[0] for i in samples_regions_areas]
regions = [i.split('_')[1] for i in samples_regions_areas]
areas = [i.split('_')[2] for i in samples_regions_areas]

          
def main():
  with Pool(12) as pool:
    pool.starmap(make_plot, zip(samples, regions, areas))

if __name__=="__main__":
    main()
