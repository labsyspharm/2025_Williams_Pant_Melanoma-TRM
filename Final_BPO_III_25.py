iI# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 11:52:01 2025

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

feature_table_path = [r"N:\lsp-analysis\cycif-production\150-BPO\MousepcaII_24\LSP17991\quantification\LSP17991--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\BPO_PCA_2025\LSP27996\quantification\LSP27996--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\MousepcaII_24\LSP28001\quantification\LSP28001--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\MousepcaII_24\LSP28006\quantification\LSP28006--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\BPO_PCA_2025\LSP28011\quantification\LSP28011--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\MousepcaII_24\LSP28021\quantification\LSP28021--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\MousepcaII_24\LSP28026\quantification\LSP28026--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\MousepcaII_24\LSP28031\quantification\LSP28031--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\BPO_PCA_2025\LSP28036\quantification\LSP28036--unmicst_cellRing.csv",
                      r"N:\lsp-analysis\cycif-production\150-BPO\BPO_PCA_2025\LSP28041\quantification\LSP28041--unmicst_cellRing.csv" ]

adata = sm.pp.mcmicro_to_scimap (feature_table_path, remove_dna=True, log =True,drop_markers = ["DAPI1", "ab1", "ab2", "DAPI2","DAPI3", "DAPI4", "DAPI5", "DAPI6", "PD-L1", "F480", "Ki67"])

manual_gate = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\Gates\gates.csv")
adata = sm.pp.rescale(adata, gate=manual_gate, log=True, imageid='imageid', failed_markers=None, method='all', random_state=0)

# Markers in dataset
adata.var.index
adata.obs['LOOSE'].value_counts()
adata.obs['imageid'].value_counts()
adata.n_obs
adata.obs['Tissue_ID'].value_counts()
adata.obs
adata.obs.columns
adata.var.head()

#Remove_ROIs
roi1 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP17991_remove.csv")
adata=sm.hl.addROI_omero(adata, roi1, label='LOOSE', subset='LSP17991--unmicst_cellRing', naming_column='Text')

roi2 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28001_remove.csv")
adata=sm.hl.addROI_omero(adata, roi2, label='LOOSE', subset='LSP28001--unmicst_cellRing', naming_column='Text', overwrite=False)

roi3 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28011_remove.csv")
adata=sm.hl.addROI_omero(adata, roi3, label='LOOSE', subset='LSP28011--unmicst_cellRing',  naming_column='Text', overwrite=False)

roi4 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28021_remove.csv")
adata=sm.hl.addROI_omero(adata, roi4, label='LOOSE', subset='LSP28021--unmicst_cellRing', naming_column='Text', overwrite=False)

roi5 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28026_remove.csv")
adata=sm.hl.addROI_omero(adata, roi5, label='LOOSE', subset='LSP28026--unmicst_cellRing', naming_column='Text', overwrite=False)

roi6 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28036_remove.csv")
adata=sm.hl.addROI_omero(adata, roi6, label='LOOSE', subset='LSP28036--unmicst_cellRing', naming_column='Text', overwrite=False)

roi7 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28041_remove.csv")
adata=sm.hl.addROI_omero(adata, roi7, label='LOOSE', subset='LSP28041--unmicst_cellRing', naming_column='Text', overwrite=False)

roi01 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP27996_remove.csv")
adata=sm.hl.addROI_omero(adata, roi01, label='LOOSE', subset='LSP27996--unmicst_cellRing', naming_column='Text', overwrite=False)

roi02 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28006_remove.csv")
adata=sm.hl.addROI_omero(adata, roi02, label='LOOSE', subset='LSP27996--unmicst_cellRing', naming_column='Text', overwrite=False)

roi03 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28031_remove.csv")
adata=sm.hl.addROI_omero(adata, roi03, label='LOOSE', subset='LSP28031--unmicst_cellRing', naming_column='Text', overwrite=False)


area = ['Remove']
adata=adata[(~adata.obs['LOOSE'].isin(area))] 

adata.write(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\python\adata_BPO_III_final.h5ad")


#Tissue_ID
roi8 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP17991_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi8, label='Tissue_ID', subset='LSP17991--unmicst_cellRing')

roi9 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28001_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi9, label='Tissue_ID', subset='LSP28001--unmicst_cellRing', overwrite=False)

roi10 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28011_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi10, label='Tissue_ID', subset='LSP28011--unmicst_cellRing', overwrite=False)

roi11 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28021_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi11, label='Tissue_ID', subset='LSP28021--unmicst_cellRing', overwrite=False)

roi12 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28026_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi12, label='Tissue_ID', subset='LSP28026--unmicst_cellRing', overwrite=False)

roi13 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28036_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi13, label='Tissue_ID', subset='LSP28036--unmicst_cellRing', overwrite=False)

roi14 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28041_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi14, label='Tissue_ID', subset='LSP28041--unmicst_cellRing', overwrite=False)

roi04 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP27996_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi04, label='Tissue_ID', subset='LSP27996--unmicst_cellRing', overwrite=False)

roi05 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28006_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi05, label='Tissue_ID', subset='LSP28006--unmicst_cellRing', overwrite=False)

roi06 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28031_sample_id.csv")
adata=sm.hl.addROI_omero(adata, roi06, label='Tissue_ID', subset='LSP28031--unmicst_cellRing', overwrite=False)


#Location
roi15 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP17991_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi15, label='Location', subset='LSP17991--unmicst_cellRing')

roi16 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28001_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi16, label='Location', subset='LSP28001--unmicst_cellRing', overwrite=False)

roi17 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28011_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi17, label='Location', subset='LSP28011--unmicst_cellRing', overwrite=False)

roi18 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28021_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi18, label='Location', subset='LSP28021--unmicst_cellRing', overwrite=False)

roi19 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28026_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi19, label='Location', subset='LSP28026--unmicst_cellRing', overwrite=False)

roi20 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28036_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi20, label='Location', subset='LSP28036--unmicst_cellRing', overwrite=False)

roi21 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28041_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi21, label='Location', subset='LSP28041--unmicst_cellRing', overwrite=False)

roi07 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP27996_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi07, label='Location', subset='LSP27996--unmicst_cellRing', overwrite=False)

roi08 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28006_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi08, label='Location', subset='LSP28006--unmicst_cellRing', overwrite=False)

roi09 = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\ROIs\LSP28031_upper_lower.csv")
adata=sm.hl.addROI_omero(adata, roi09, label='Location', subset='LSP28031--unmicst_cellRing', overwrite=False)


phenotype = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\Phenotype\phenotype.csv")
adata = sm.tl.phenotype_cells (adata, phenotype, gate=0.5, label='phenotype', imageid='imageid', pheno_threshold_percent=None, pheno_threshold_abs=None)
adata.obs['phenotype'].value_counts()
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotype', dendrogram=False, use_raw=False, cmap="vlag", standard_scale='var')

phenotype = pd.read_csv(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\Phenotype\phenotype_L1.csv")
adata = sm.tl.phenotype_cells (adata, phenotype, gate=0.5, label='phenotype_L1', imageid='imageid', pheno_threshold_percent=None, pheno_threshold_abs=None)
adata.obs['phenotype_L1'].value_counts()
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='phenotype_L1', dendrogram=False, use_raw=False, cmap="vlag", standard_scale='var')

adata.obs['Location'].value_counts()
adata.obs['Tissue_ID'].value_counts()
adata.obs['Time_point'].value_counts()

adata.obs['Tissue_ID'] = adata.obs['Tissue_ID'].fillna('Other')

adata.obs['Location'] = adata.obs['Location'].fillna('Other')

sm.pl.spatial_scatterPlot (adata=adata,
                 colorBy = ['phenotype'],
                 subset = 'LSP28041--unmicst_cellRing',
                 figsize=(4,4),
                 s=0.5)


#Dotplot Works
plt.rcParams['figure.facecolor'] = 'white'         # Background of the figure
plt.rcParams['axes.facecolor'] = 'white'  
save_directory = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II"
os.makedirs(save_directory, exist_ok=True)
sc.settings.figdir = save_directory 
sc.pl.dotplot(adata, var_names= adata.var.index, groupby='phenotype_L1', dendrogram=False, use_raw=False, cmap="vlag", standard_scale='var', save='dotplot_phenotype_L1.pdf' )


save_directory = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results"
os.makedirs(save_directory, exist_ok=True)
sc.settings.figdir = save_directory
sc.pl.umap(mdata, color=['phenotype'], save='dotplot_phenotyle_II.pdf')


adata.write(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\python\adata_BPO.h5ad")


adata = ad.read_h5ad(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\python\adata_BPO.h5ad")


sm.pl.stacked_barplot(adata, x_axis='Tissue_ID', y_axis='phenotype_L1', method='percent', plot_tool='matplotlib', figsize=(10, 10))

T1 = ['5619d17PL','5600d17SK','5708d10SK', '5816d17PL', '5706d10SK', '5574d17PL', '5707d10SK', '5705d18SK', '5619d7SK', '5711d7SK', '5711d17SK'] #Early
T2 = ['5600d24PL','5708d21PL','5816d24PL', '5574d24PL', '5707d21PL', '5810d21PL', '5706d21PL'] #Mid
T3 = ['5708d44T1','5706d44T2', '5574d44T1', '5707d44T2', '5705d44SK', '5711d44T2', '5619d44T1', '5810d44T1', '5600d41T2'] #Late

rename= {'T1': ['5619d17PL','5600d17SK','5708d10SK', '5816d17PL', '5706d10SK', '5574d17PL', '5707d10SK', '5705d18SK', '5619d7SK', '5711d7SK', '5711d17SK'], 
         'T2': ['5600d24PL','5708d21PL','5816d24PL', '5574d24PL', '5707d21PL', '5810d21PL', '5706d21PL'],
         'T3':['5708d44T1', '5706d44T2', '5574d44T1', '5707d44T2', '5705d44SK', '5711d44T2', '5619d44T1', '5810d44T1', '5600d41T2']}


adata = sm.hl.rename (adata, rename, from_column='Tissue_ID', to_column='Time_point')

adata.obs['Tissue_ID'] = adata.obs['Tissue_ID'].replace('5574d44T1', '5574d41T2')

#typo correction
#adata.obs['Location'] = adata.obs['Location'].replace('lower_5816_d24PL', 'lower_5816d24PL')

adata.obs['Tissue_ID'].value_counts()
adata.obs['Location'].value_counts()



upper = ['upper_5619d17PL', 'upper_5600d17SK', 'upper_5708d10SK', 'upper_5816d17PL', 'upper_5706d10SK', 
         'upper_5574d17PL', 'upper_5707d10SK', 'upper_5705d18SK', 'upper_5619d7SK', 'upper_5711d7SK', 
         'upper_5711d17SK', 'upper_5600d24PL', 'upper_5708d21PL', 'upper_5816d24PL', 'upper_5574d24PL', 'upper_5707d21PL', 
         'upper_5810d21PL', 'upper_5706d21PL', 'upper_5800d24PL', 'upper_5708d44T1', 'upper_5706d44T2', 'upper_5574d41T1', 'upper_5707d44T2', 
         'upper_5705d44SK', 'upper_5711d44T2', 'upper_5810d44T1', 'upper_5619d44T1', 'upper_5810d44T1', 'upper_5600d41T2'] #upper

lower = ['lower_5619d17PL', 'lower_5600d17SK', 'lower_5708d10SK', 'lower_5816d17PL', 'lower_5706d10SK', 
         'lower_5574d17PL', 'lower_5707d10SK', 'lower_5705d18SK', 'lower_5619d7SK', 'lower_5711d7SK', 'lower_5711d17SK', 
         'lower_5600d24PL', 'lower_5708d21PL', 'lower_5816d24PL', 'lower_5574d24PL', 'lower_5707d21PL', 'lower_5810d21PL', 
         'lower_5706d21PL', 'lower_5800d24PL', 'lower_5708d44T1', 'lower_5706d44T2', 'lower_5574d41T1', 'lower_5707d44T2', 'lower_5705d44SK', 
         'lower_5711d44T2', 'lower_5810d44T1', 'lower_5619d44T1', 'lower_5810d44T1', 'lower_5600d41T2'] #lower


rename = {'upper': ['upper_5574d41T2', 'upper_5711d44T2', 'upper_5707d44T2', 'upper_5708d44T1', 'upper_5706d44T2',         
              'upper_5619d17PL', 'upper_5600d17SK', 'upper_5810d21PL', 'upper_5706d21PL', 'upper_5708d10SK', 
              'upper_5816d17PL', 'upper_5706d10SK', 'upper_5574d17PL', 'upper_5707d10SK', 'upper_5600d24PL',
              'upper_5705d18SK', 'upper_5619d7SK', 'upper_5816d24PL', 'upper_5707d21PL', 'upper_5708d21PL',
              'upper_5711d17SK', 'upper_5705d44SK', 'upper_5711d7SK', 'upper_5574d24PL', 'upper_5810d44T1', 'upper_5800d24PL',
              'upper_5619d44T1', 'upper_5810d44T1', 'upper_5600d41T2'],
          'lower': ['lower_5574d41T2', 'lower_5711d44T2', 'lower_5707d44T2', 'lower_5708d44T1', 'lower_5706d44T2',
                        'lower_5619d17PL', 'lower_5600d17SK', 'lower_5810d21PL', 'lower_5706d21PL', 'lower_5708d10SK',
                        'lower_5816d17PL', 'lower_5706d10SK', 'lower_5574d17PL', 'lower_5707d10SK', 'lower_5600d24PL',
                        'lower_5705d18SK', 'lower_5619d7SK', 'lower_5816d24PL', 'lower_5707d21PL', 'lower_5708d21PL',
                        'lower_5711d17SK', 'lower_5705d44SK', 'lower_5711d7SK', 'lower_5574d24PL', 'lower_5810d44T1', 
                        'lower_5619d44T1', 'lower_5810d44T1', 'lower_5600d41T2', 'lower_5800d24PL']}


adata = sm.hl.rename(adata, rename, from_column='Location', to_column='Location_combined')


sm.pl.stacked_barplot(adata, x_axis='Location_combined', y_axis='phenotype', method='percent', plot_tool='matplotlib', figsize=(10, 10))

sm.pl.stacked_barplot(adata, x_axis='Location_combined', y_axis='phenotype', method='percent', subset_xaxis = ['upper', 'lower'], plot_tool='matplotlib', figsize=(10, 10))
sm.pl.stacked_barplot(adata, x_axis='Location_combined', y_axis='phenotype', method='absolute', subset_xaxis = ['upper', 'lower'], plot_tool='matplotlib', figsize=(10, 10))


sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='percent', plot_tool='matplotlib', figsize=(10, 10))

sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='percent', plot_tool='matplotlib', figsize=(10, 10))


sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='percent', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['Tregs', 'T helper', 'CD103-Gzmb-', 'CD103-Gzmb+', 'CD103+Gzmb-', 'CD103+Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))

sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='absolute', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['Tregs', 'T helper', 'CD103-Gzmb-', 'CD103-Gzmb+', 'CD103+Gzmb-', 'CD103+Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))


sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='percent', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['Tregs', 'T helper', 'Keratinocytes', 'Melanocytes', 'Macrophage', 'Endothelial cells', 'Immune cells', 'DCs', 'Trms', 'CD103-'], plot_tool='matplotlib', figsize=(10, 10))
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='absolute', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['Tregs', 'T helper', 'Keratinocytes', 'Melanocytes', 'Macrophage', 'Endothelial cells', 'Immune cells', 'DCs', 'Trms', 'CD103-'], plot_tool='matplotlib', figsize=(10, 10))



#Stackbar_plot_publication
# Create a dictionary with your colors
cell_colors = {
    'Melanocyte': '#fcba79', 
    'Endothelial cells': '#8d574c', 
    'Macrophage': '#b5bd61', 
    'Keratinocytes': '#2077b5', 
    'Trms': '#d62729', 
    'T helper': '#f69696', 
    'DCs': '#8559a5',
    'Tregs': '#c4afd3', 
    'Immune cells': '#f57e20', 
    'CD103-': '#d87ab1',
}



cell_colors = {
    'Melanocyte': '#d87ab1', 
    'Endothelial cells': '#269d68', 
    'Macrophage': '#8d574c', 
    'Keratinocytes': '#8559a5', 
    'Trms': '#afc7e5', 
    'T helper': '#9dd089', 
    'DCs': '#f57e20',
    'Tregs': '#1ebcce', 
    'Immune cells': '#d62729', 
    'CD103-': '#2077b5',
}



# Create a colormap from your color dictionary
cmap = mcolors.ListedColormap([cell_colors[key] for key in cell_colors.keys()])


# Plot the stacked barplot with the color map
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='percent', subset_xaxis = ['T1', 'T2', 'T3'],
                    subset_yaxis = ['Tregs', 'T helper', 'Keratinocytes', 'Melanocytes', 'Macrophage', 'Endothelial cells', 'Immune cells', 'DCs', 'Trms', 'CD103-'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results",
                      fileName='cell_proportions_N.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})


import matplotlib.colors as mcolors

cell_colors = {
    'Melanocytes': '#d87ab1', 
    'Endothelial cells': '#269d68', 
    'Macrophage': '#8d574c', 
    'Keratinocytes': '#8559a5', 
    'Trms': '#afc7e5', 
    'T helper': '#9dd089', 
    'DCs': '#f57e20',
    'Tregs': '#1ebcce', 
    'Immune cells': '#d62729', 
    'CD103-': '#2077b5',
}

# Define the order in which you want to show categories
subset_y = ['Tregs', 'T helper', 'Keratinocytes', 'Melanocytes', 'Macrophage',
            'Endothelial cells', 'Immune cells', 'DCs', 'Trms', 'CD103-']

# Extract colors in the same order
ordered_colors = [cell_colors[cell] for cell in subset_y]

# Create a colormap
cmap = mcolors.ListedColormap(ordered_colors)

# Plot
sm.pl.stacked_barplot(
    adata,
    x_axis='Time_point',
    y_axis='phenotype',
    method='percent',
    subset_xaxis=['T1', 'T2', 'T3'],
    subset_yaxis=subset_y,
    plot_tool='matplotlib',
    matplotlib_cmap=cmap,  # ✅ Correct argument
    saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results",
    fileName='cell_proportions_N.pdf',
    figsize=(6, 4),
    **{'edgecolor': 'none'}
)






# Plot the stacked barplot with the color map
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='percent', subset_xaxis = ['T1', 'T2', 'T3'],
                    subset_yaxis = ['Tregs', 'T helper', 'Keratinocytes', 'Melanocytes', 'Macrophage', 'Endothelial cells', 'Immune cells', 'DCs', 'Trms', 'CD103-'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results",
                      fileName='cell_proportions.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype', method='percent', subset_xaxis = ['T1', 'T2', 'T3'],
                    subset_yaxis = ['Tregs', 'T helper', 'Keratinocytes', 'Melanocytes', 'Macrophage', 'Endothelial cells', 'Immune cells', 'DCs', 'Trms', 'CD103-'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results",
                      #fileName='cell_proportions.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

subset_xaxis=['5619d17PL','5600d17SK','5708d10SK', '5816d17PL', '5706d10SK', '5574d17PL', '5707d10SK', 
             '5705d18SK', '5619d7SK', '5711d7SK', '5711d17SK', '5600d24PL','5708d21PL','5816d24PL', '5574d24PL', 
             '5707d21PL', '5810d21PL', '5706d21PL', '5708d44T1','5706d44T2', '5574d44T1', '5707d44T2', '5705d44SK', 
             '5711d44T2'], order_xaxis=['5619d17PL','5600d17SK','5708d10SK', '5816d17PL', '5706d10SK', '5574d17PL', '5707d10SK', 
             '5705d18SK', '5619d7SK', '5711d7SK', '5711d17SK', '5600d24PL','5708d21PL','5816d24PL', '5574d24PL', 
             '5707d21PL', '5810d21PL', '5706d21PL', '5708d44T1','5706d44T2', '5574d44T1', '5707d44T2', '5705d44SK', 
             '5711d44T2'],


sm.pl.stacked_barplot(adata, x_axis='Tissue_ID', y_axis='phenotype_L1', method='percent', subset_yaxis = ['Tregs', 'T helper', 'CD103-Gzmb-', 'CD103-Gzmb+', 'CD103+Gzmb-', 'CD103+Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))

sm.pl.stacked_barplot(adata, x_axis='Tissue_ID', y_axis='phenotype_L1', method='absolute',
                      subset_xaxis=['5574d44T1', '5711d44T2', '5707d44T2', '5708d44T1', '5706d44T2',
                                    '5619d17PL', '5600d17SK', '5810d21PL', '5706d21PL', '5708d10SK',
                                    '5816d17PL', '5706d10SK', '5574d17PL', '5707d10SK', '5600d24PL',
                                    '5705d18SK', '5619d7SK', '5816d24PL', '5707d21PL', '5708d21PL',
                                    '5711d17SK', '5705d44SK', '5711d7SK', '5574d24PL'],
                      subset_yaxis = ['Tregs', 'T helper', 'CD103-Gzmb-', 'CD103-Gzmb+', 'CD103+Gzmb-', 'CD103+Gzmb+'], 
                      plot_tool='matplotlib', figsize=(10, 10))



sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='percent', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['CD103-Gzmb-', 'CD103-Gzmb+', 'CD103+Gzmb-', 'CD103+Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='absolute', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['CD103-Gzmb-', 'CD103-Gzmb+', 'CD103+Gzmb-', 'CD103+Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))

#CD8 T cells
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='percent', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['CD103-Gzmb-', 'CD103-Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='absolute', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['CD103-Gzmb-', 'CD103-Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))

sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='percent', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['CD103+Gzmb-', 'CD103+Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='absolute', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['CD103+Gzmb-', 'CD103+Gzmb+'], plot_tool='matplotlib', figsize=(10, 10))
 
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='percent', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['Tregs', 'T helper'], plot_tool='matplotlib', figsize=(10, 10))
sm.pl.stacked_barplot(adata, x_axis='Time_point', y_axis='phenotype_L1', method='absolute', subset_xaxis = ['T1', 'T2', 'T3'], subset_yaxis = ['Tregs', 'T helper'], plot_tool='matplotlib', figsize=(10, 10))


column_to_filter = 'Time_point'
mask = adata.obs[column_to_filter] != 'Other'
bdata = adata[mask, :].copy()
print(bdata.obs['Time_point'].value_counts())


column_to_filter = 'Location_combined'
mask = adata.obs[column_to_filter] != 'Other'
bdata = adata[mask, :].copy()
print(bdata.obs['Location_combined'].value_counts())

adata.write(r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\python\adata_BPO_061025.h5ad")

np.expm1(2)*0.32
np.expm1(4)*0.32
np.expm1(6)*0.32
np.expm1(8)*0.32
np.expm1(10)*0.32

For depth
2500*0.32
5000*0.32