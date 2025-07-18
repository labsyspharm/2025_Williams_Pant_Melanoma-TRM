# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 09:15:16 2025

@author: shp918
"""

import anndata as ad
import pandas as pd
import scimap as sm
import scanpy as sc
import numpy as np
import sys
import os
import seaborn as sns; sns.set(color_codes=True)
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

adata = ad.read_h5ad(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\PCA_phaseII\cycif_analysis_materials\data_p135_e9\p135e9_cy13_qcV2.6merged.h5ad")


adata = sm.hl.dropFeatures(adata, drop_markers=['5hmc', 'Langerin', 'LaminB1', 'CDKN1A', 'CD1c', 'pH3','GZMB', 'aSMA', 'PDL1', 'HLA/DRB1', 'aSMA', 'HLA/ABC', 'CD19', 'PRAME', 'PCNA', 'CD163', 'L1CAM', 'VIM', 'HLA/ABC', 'CD66b', 'MITF', 'CD16', 'B2M'])


# Markers in dataset
adata.var.index
adata.obs['phenotype_L2'].value_counts()
adata.obs['imageid'].value_counts()
adata.n_obs
adata.obs['distance_from_epidermis2500_binned'].value_counts()
adata.obs['histopath_roi'].value_counts()
adata.obs
adata.obs.columns
adata.var.head()

phenotype = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\Phenotype\phenotype_human_L1.csv")
adata = sm.tl.phenotype_cells (adata, phenotype, gate=0.5, label='phenotype_L2', imageid='imageid', pheno_threshold_percent=None, pheno_threshold_abs=None)
adata.obs['phenotype_L2'].value_counts()
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotype_L2', dendrogram=False, use_raw=False, cmap="vlag", standard_scale='var')


#Dotplot Works
plt.rcParams['figure.facecolor'] = 'white'         # Background of the figure
plt.rcParams['axes.facecolor'] = 'white'  
save_directory = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II"
os.makedirs(save_directory, exist_ok=True)
sc.settings.figdir = save_directory 
sc.pl.dotplot(adata, var_names= adata.var.index, groupby='phenotype_L2', dendrogram=False, use_raw=False, cmap="vlag", standard_scale='var', save='dotplot_phenotype_L2_human.pdf' )

sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L1', method='percent', plot_tool='matplotlib', figsize=(10, 10))

sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotypeSHP', dendrogram=False, use_raw=False, cmap="vlag", standard_scale='var')


sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L1', subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T hellper', 'DCs', 'T cells', 'B cells', 'Immune cells'], method='percent', plot_tool='matplotlib', figsize=(10, 10))

sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L1', subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T hellper', 'DCs', 'T cells', 'B cells', 'Immune cells'], method='percent', plot_tool='matplotlib', figsize=(10, 10))

sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L1', method='percent', plot_tool='matplotlib', figsize=(10, 10))


sm.pl.stacked_barplot(adata, x_axis='distance_from_epidermis2500_binned', y_axis='phenotype_L2', subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'], method='percent', plot_tool='matplotlib', figsize=(10, 10))

distance_from_epidermis2500_binned

sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L1', method='percent', plot_tool='matplotlib', figsize=(10, 10))


subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],

adata.write(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\PCA_phaseII\cycif_analysis_materials\data_p135_e9\p135e9_cy13_allcells_gated_SMP_V5.h5ad")


fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(21, 7))
metric_col='phenotype_L1'

sns.scatterplot(
    data=adata.obs,
    x="X_centroid",
    y="Y_centroid",
    hue=metric_col,
    palette='colorblind',
    alpha=0.8,
    s=5,
    ax=ax1,
)

plt.show()

Plot the stacked barplot with the color map
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='percent', subset_xaxis = ['T1', 'T2', 'T3'],
                    subset_yaxis = ['Tregs', 'T helper', 'Keratinocytes', 'Melanocytes', 'Macrophage', 'Endothelial cells', 'Immune cells', 'DCs', 'Trms', 'CD103-'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results",
                      fileName='cell_proportions_N.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

adata.obs['histopath_roi'] = adata.obs['histopath_roi'].fillna('Other')



rename= {'VGP': ['VGP'], 
         'Rest': ['normal','precursor','Other']}


adata = sm.hl.rename (adata, rename, from_column='histopath_roi', to_column='histopath_2')
adata.obs['histopath_2'].value_counts()
