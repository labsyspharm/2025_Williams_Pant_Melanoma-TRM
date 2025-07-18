# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 10:10:31 2025

@author: shp918
"""

adata = ad.read_h5ad(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\python\adata_BPO.h5ad")
#subset to mdata (100000)
mdata = sc.pp.subsample(adata, fraction=None, n_obs=300000, random_state=0, copy=True)
mdata = sm.tl.cluster(mdata, method='leiden', resolution=0.3, use_raw=False, log=False)
mdata.obs['leiden'].value_counts()
sm.pl.heatmap(mdata, groupBy='leiden', standardScale='column', figsize=(10,8), showPrevalence=True)

sc.pp.neighbors(mdata, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(mdata) # Build a UMAP to visualize the neighbourhood graph


plt.rcParams['figure.facecolor'] = 'white'         # Background of the figure
plt.rcParams['axes.facecolor'] = 'white'   
sc.pl.umap(mdata, color=['imageid']) # View the clustering



sc.pl.umap(mdata, color=['leiden']) # View the clustering                                                                                                                                                                                                                                                                                                        
sc.pl.umap(mdata, color=['imageid']) # View the clustering
sc.pl.umap(mdata, color=['Time_point']) # View the clustering
sc.pl.umap(mdata, color=['Tissue_ID']) # View the clustering

sc.pl.umap(mdata, color=['Tissue']) # View the clustering
sc.pl.umap(mdata, color=['distance_from_epidermis2500_binned']) # View the clustering



default_color = '#d3d3d3'
customColors = {
    'Trms' : '#8559a5',
    'CD103-': default_color,  # <--- Added quotes here
    'Tregs' : default_color,
    'Endothelial cells' : default_color,
    'T helper': default_color,
    'DCs' : default_color,
    'Macrophage' : default_color,
    'Immune cells': default_color,
    'Keratinocytes': default_color,
    'Unknown': default_color,
    'Melanocytes': default_color,
    'B cells' : default_color,
    'T cells' : default_color,
}

# --- Set Matplotlib parameters *before* the plotting function call ---
plt.rcParams['figure.facecolor'] = 'white' # Background of the figure
plt.rcParams['axes.facecolor'] = 'white'   # Background of the plot area
sc.pl.umap(mdata, color=['phenotype_L2'], palette=customColors) # View the clustering
