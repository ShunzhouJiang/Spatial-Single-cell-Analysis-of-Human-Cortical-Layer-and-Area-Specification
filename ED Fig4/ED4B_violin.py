import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np


adata = sc.read('merscope_integrated_855.h5ad')
h2_colors = ['#92D050','#99D35B','#A0D666','#ADDB7C','#C8E6A7','#00B0F0','#10B5F1','#D8ADFF','#40C3F4','#7FD6F7','#00B050','#10B55B','#20BA66','#40C37C','#7FD6A7','#FFC000','#FFC820','#FFCF40','#FFD760','#FFDE7F','#209AD2','#FEAAA8','#B870FF','#C8E61B','#FF0000','#CC99FF','#7030A0','#8C57B2','#305496','#4A69A3','#647EB0','#7E93BD','#97A8CA']
h2_types = ['RG1', 'oRG1','Astro-late1','tRG','vRG-late',"EN-ET-SP-early",'EN-ET-SP-P','EN-ET-L5/6','EN-ET-L6-early','EN-ET-SP','IPC-oSVZ','IPC-SVZ-1','IPC-iSVZ', 'IPC-SVZ-2', 'IPC-VZ/SVZ','EN-IZ-1','EN-L2','EN-IZ-2','EN-oSVZ-1','En-oSVZ-2','EN-IT-L6','EN-IT-L4','EN-IT-L4/5','EN-IT-L3/4','EN-IT-L2/3','EC','Astro-1', 'OPC','IN-SST', 'IN-CGE', 'INP-VZ/GE', 'IN-VZ/GE', 'IN-MGE']
h2_dict = dict(zip(h2_types, h2_colors))

obs_all = adata.obs

order = np.array(['tRG', 'vRG-late', 'IN-VZ/GE', 'INP-VZ/GE', 'IN-MGE', 'IN-CGE', 'IPC-SVZ-2', 'IPC-VZ/SVZ', 'IPC-iSVZ', 'RG1', 'IPC-SVZ-1', 'Astro-1', 'oRG1', 'IPC-oSVZ', 'IN-SST', 'En-oSVZ-2', 'EN-oSVZ-1', 'EN-IZ-1', 'EN-IZ-2', 'EC', 'Astro-late1', 'OPC', 'EN-L2', 'EN-ET-L6-early', 'EN-ET-SP', 'EN-ET-SP-P', 'EN-ET-SP-early', 'EN-ET-L5/6', 'EN-IT-L6', 'EN-IT-L4/5', 'EN-IT-L4', 'EN-IT-L3/4', 'EN-IT-L2/3'])

def make_plot(sample, region, area):
  obs = obs_all[(obs_all['sample']==sample) & (obs_all.region==region) & (obs_all.area==area)]
  obs = obs[obs['relative_height'].notna()]
  type_rm = list(obs.value_counts('H2_annotation').index[(obs.value_counts('H2_annotation')<50)])
  if len(type_rm)>0:
        obs1 = obs[~obs.H2_annotation.isin(type_rm)]
        obs2 = obs[obs.H2_annotation.isin(type_rm)]
        plt.figure(figsize=(30,15));
        plot = sns.violinplot(x='H2_annotation', y='relative_height', hue='H2_annotation', data=obs1, order=order, palette=h2_dict, density_norm='width', inner = None, dodge=False, cut=0);
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=12); plt.yticks(fontsize=12); plot.set(xlabel=None); plot.set_ylabel('Laminar Depth', fontsize=20); plt.ylim(0,1);
        sns.stripplot(x='H2_annotation', y='relative_height', hue='H2_annotation', data=obs2, order=order, palette=h2_dict);
        plt.show()
        #plt.savefig(sample + '_' + region + '_' + area + '_violin.png', dpi=500, bbox_to_inches = 'tight'); plt.clf()
  else:
        plt.figure(figsize=(30,15));
        plot = sns.violinplot(x='H2_annotation', y='relative_height', hue='H2_annotation', data=obs, order=order, palette=h2_dict, density_norm='width', inner = None, dodge=False, cut=0);
        plot.legend().remove(); plt.xticks(rotation=90, fontsize=12); plt.yticks(fontsize=12); plot.set(xlabel=None); plot.set_ylabel('Laminar Depth', fontsize=20); plt.ylim(0,1);
        plt.show()
        #plt.savefig(sample + '_' + region + '_' + area + '_violin.png', dpi=500, bbox_to_inches = 'tight'); plt.clf()


samples_regions_areas = (adata.obs['sample'].astype('str') + '_' + adata.obs.region.astype('str')+ '_' + adata.obs.area.astype('str')).unique()        

samples = [i.split('_')[0] for i in samples_regions_areas]
regions = [i.split('_')[1] for i in samples_regions_areas]
areas = [i.split('_')[2] for i in samples_regions_areas]


          
def main():
  with Pool(12) as pool:
    pool.starmap(make_plot, zip(samples, regions, areas))

if __name__=="__main__":
    main()
