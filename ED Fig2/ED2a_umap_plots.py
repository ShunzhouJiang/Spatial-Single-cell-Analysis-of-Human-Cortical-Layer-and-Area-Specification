import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


obs = pd.read_csv('umap_obs.csv', index_col = 0)

plt.figure(figsize = (10,10)); sns.scatterplot(data = obs, x = 'umap_x', y = 'umap_y', hue = 'sample', s=0.5); plt.axis('off'); plt.savefig('umap_sample.png', dpi=500); plt.clf()
plt.figure(figsize = (10,10)); sns.scatterplot(data = obs, x = 'umap_x', y = 'umap_y', hue = 'gw', s=0.5); plt.axis('off'); plt.savefig('umap_gw.png', dpi=500); plt.clf()
plt.figure(figsize = (10,10)); sns.scatterplot(data = obs, x = 'umap_x', y = 'umap_y', hue = 'cortical_area', s=0.5); plt.axis('off'); plt.savefig('umap_cortical_area.png', dpi=500); plt.clf()



