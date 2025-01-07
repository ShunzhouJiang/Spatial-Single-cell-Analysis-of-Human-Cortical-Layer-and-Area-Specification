import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import umap
import seaborn


#adata = sc.read('merscope_integrated_855.h5ad')
#sc.pp.scale(adata, zero_center=True, max_value=6)
#data = adata.X
#np.random.seed(100)
#samples = np.random.choice(data.shape[0],1000000)
#obs = adata.obs.iloc[samples]
#del adata
#data1 = data[samples,:]
#embedder = umap.UMAP(init = 'random', random_state=100).fit(data1)
#merscope_umap = embedder.embedding_
obs = pd.read_csv('umap_obs.csv', index_col = 0)

color_dict = {'RG1':'#FADBD8', 'oRG1':'#F1948A','Astro-late1':'#E74C3C','tRG':'#B03A2E','vRG-late':'#78281F',
              "EN-ET-SP-early":'#E8DAEF','EN-ET-SP-P':'#BB8FCE','EN-ET-L5/6':'#8E44AD','EN-ET-L6-early':'#6C3483','EN-ET-SP':'#4A235A',
              'IPC-oSVZ':'#D6EAF8','IPC-SVZ-1':'#7FB3D5','IPC-iSVZ':'#2980B9', 'IPC-SVZ-2':'#1F618D', 'IPC-VZ/SVZ':'#154360',
              'EN-IZ-1':'#D4EFDF','EN-L2':'#7DCEA0','EN-IZ-2':'#27AE60','EN-oSVZ-1':'#1E8449','En-oSVZ-2':'#145A32',
              'EN-IT-L6':'#FDEBD0','EN-IT-L4':'#F8C471','EN-IT-L4/5':'#F39C12','EN-IT-L3/4':'#B9770E','EN-IT-L2/3':'#7E5109',
              'EC':'#979A9A',
              'Astro-1':'#D6DBDF', 'OPC':'#ABB2B9',
              'IN-SST':'#FCF3CF', 'IN-CGE':'#F7DC6F', 'INP-VZ/GE':'#F1C40F', 'IN-VZ/GE':'#B7950B', 'IN-MGE':'#7D6608'}              
h2_colors = ['#92D050','#99D35B','#A0D666','#ADDB7C','#C8E6A7','#00B0F0','#10B5F1','#D8ADFF','#40C3F4','#7FD6F7','#00B050','#10B55B','#20BA66','#40C37C','#7FD6A7','#FFC000','#FFC820','#FFCF40','#FFD760','#FFDE7F','#209AD2','#FEAAA8','#B870FF','#C8E61B','#FF0000','#CC99FF','#7030A0','#8C57B2','#305496','#4A69A3','#647EB0','#7E93BD','#97A8CA']
h2_color_dict = dict(zip(color_dict.keys(), h2_colors))
#obs['x'] = merscope_umap[:,0]
#obs['y'] = merscope_umap[:,1]

plt.figure(figsize = (10,10)); seaborn.scatterplot(data = obs, x = 'umap_x', y = 'umap_y', hue = 'H2_annotation', palette = h2_color_dict, s=0.5); plt.axis('off'); handles, labels = plt.gca().get_legend_handles_labels();
order = [labels.index(i) for i in color_dict.keys()]; plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=3, ncol = 3, bbox_to_anchor=(0.8,0.8), markerscale=0.5); 
plt.savefig('umap_1000000_h2.png', dpi=500); plt.clf()

h1_colors = ['#92D050','#00B0F0','#00B050','#FFC000','#FF0000','#CC99FF','#7030A0','#305496']
h1_types = ['RG', 'EN-ET', 'IPC', 'EN-Mig', 'EN-IT', 'EC', 'Glia', 'IN']
h1_dict = dict(zip(h1_types, h1_colors))
plt.figure(figsize = (10,10)); seaborn.scatterplot(data = obs, x = 'umap_x', y = 'umap_y', hue = 'H1_annotation', palette = h1_dict, s=0.5); plt.axis('off'); handles, labels = plt.gca().get_legend_handles_labels(); order = [labels.index(i) 
for i in h1_dict.keys()]; plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc = 'center', fontsize=3, ncol = 1, bbox_to_anchor=(0.8,0.8), markerscale=0.5); plt.savefig('umap_1000000_h1.png', dpi=500); plt.clf()

