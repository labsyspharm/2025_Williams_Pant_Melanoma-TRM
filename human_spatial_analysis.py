# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 16:38:52 2025

@author: shp918
"""

adata = sm.tl.spatial_distance (adata, phenotype='phenotype_L2')
sm.pl.spatial_distance (adata, phenotype='phenotype_L2',figsize=(5,4))

sm.pl.spatial_distance (adata, method='distribution',distance_from='Tregs',
                        distance_to = 'Trms', imageid='histopath_2', log=True, 
                        height=3,phenotype='phenotype_L2', 
                        saveDir=r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/spatial_analysis", 
                        fileName='Distance_Treg_to_Trm_histopath_2_distribution.pdf', 
                        aspect=9/8)


#Human plotted
distance_to = ['CD103-', 'Trms']
sm.pl.spatial_distance (adata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, phenotype='phenotype_L2',
                        log=True, imageid='histopath_2',
                        saveDir=r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/spatial_analysis", 
                        fileName='Distance_Treg_to_CD8_histopath_2_distance_numeric.pdf', 
                        height=3, aspect=18/16)

######RETURN DATA FOR T TEST######
distance_to = ['CD103-', 'Trms']
returned_data = sm.pl.spatial_distance (adata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, phenotype='phenotype_L2',
                        log=True, imageid='histopath_2',
                        saveDir=r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/spatial_analysis", 
                        fileName='Distance_Treg_to_CD8_histopath_2_distance_numeric.pdf', 
                        height=3, aspect=18/16,  return_data=True)

# Convert the returned data to a DataFrame
df_histopath_2 = pd.DataFrame(returned_data)

file_path = r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/df_histopath_2.csv"

# Save the DataFrame to a CSV file
df_histopath_2.to_csv(file_path, index=False)

#######XXXXX#####################


distance_to = ['CD103-', 'Trms']
sm.pl.spatial_distance (adata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, phenotype='phenotype_L2',
                        log=True, imageid='RCNs',
                        #saveDir=r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/spatial_analysis", 
                        #fileName='Distance_T helper_to_CD8_RCNs_numeric.pdf', 
                        height=3, aspect=18/16)


distance_to = ['CD103-', 'Trms']
sm.pl.spatial_distance (adata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, phenotype='phenotype_L2',
                        log=True, imageid='histopath_roi',
                        #saveDir=r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/spatial_analysis", 
                        #fileName='Distance_Treg_to_CD8_histopath_roi_distance_numeric.pdf', 
                        height=3, aspect=18/16)


sm.pl.spatial_distance(adata, method='heatmap', phenotype='phenotype_L2', imageid='RCNs')


#subset distance
distance_to = ['CD103-', 'Trms']
sm.pl.spatial_distance (adata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, phenotype='phenotype_L2',
                        log=True, imageid='distance_from_epidermis2500_binned',
                        subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)','[10000, 12500)',],
                        order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)','[10000, 12500)', ],
                        #saveDir=r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/spatial_analysis", 
                        #fileName='Distance_Treg_to_CD8_histopath_roi_distance_numeric.pdf', 
                        height=3, aspect=18/16)

subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)','[10000, 12500)',],
order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)','[10000, 12500)', ],


adata = sm.tl.spatial_count(adata, phenotype='phenotype_L2', method='radius', radius=80, label='spatial_count')
adata = sm.tl.spatial_cluster(adata, df_name='spatial_count', method='kmeans', k=5, label='neigh_kmeans')
# adata = sm.tl.spatial_cluster(adata, df_name='spatial_count', method='kmeans', k=30, label='neigh_kmeans')
sm.pl.stacked_barplot (adata, x_axis='neigh_kmeans', y_axis='phenotype_L2')


cell_type_colors = {
    'Melanocytes': '#259DD0',
    'Endothelial cells': '#6BA292',
    'Macrophage': '#8C4843',
    'Keratinocytes': '#2077b5',
    'Trms': '#d62729',
    'T helper': '#f69696',
    'DCs': '#8559a5',
    'Tregs': '#c4afd3',
    'Immune cells': '#f57e20',
    'CD103-': '#d87ab1',
    'B cells': '#1b9e77',
    'T cells': '#7570b3',
    'Unknown': '#999999',
}

# Use the correct dictionary name
cmap = mcolors.ListedColormap([cell_type_colors[key] for key in cell_type_colors])

sm.pl.stacked_barplot(adata, x_axis='neigh_kmeans', y_axis='phenotype_L2', method='percent', 
                      subset_yaxis= [ 'Melanocytes', 'Endothelial cells', 'Keratinocytes', 'B cells', 'Immune cells', 'DCs','Macrophage','Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\spatial_analysis",
                      #fileName='neigh_kmeans_barplot_L2.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='neigh_kmeans', y_axis='phenotype_L2', method='percent', 
                      subset_yaxis= ['B cells', 'Immune cells', 'DCs','Macrophage','Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\spatial_analysis",
                      #fileName='neigh_kmeans_barplot_L2_immune cells.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})


sm.pl.groupCorrelation(adata, groupBy='phenotype_L2', condition='neigh_kmeans', figsize=(9,6))

rename_dict_adata = {'RCN1': ['0', '13', '5', '11'],#Melanocytes
                     'RCN2': ['10', '4'],#Keratinocytes
                     'RCN3': ['2', '6', '9', '14', '18'],#Macrophjage
                     'RCN4': ['1', '12'],#Endothelial cells
                     'RCN5': ['17'],#B cells
                     'RCN6': ['8'],#DCs
                     'RCN7': ['16'],#CD103-
                     'RCN8': ['7'],#Trms
                     'RCN9': ['3', '19'],#CD4
                     'RCN10': ['15'],#Immune cells
                     }

rename_dict_adata = {'RCN1': ['0', '13', '5', '11'],#Melanocytes
                     'RCN2': ['10', '4'],#Keratinocytes
                     'RCN3': ['2', '6', '9', '14', '18'],#Macrophjage
                     'RCN4': ['1', '12'],#Endothelial cells
                     'RCN5': ['17'],#B cells
                     'RCN6': ['8'],#DCs
                     'RCN7': ['16', '7', '3', '19'],#T cells
                     'RCN8': ['15'],#Immune cells
                     }

rename_dict_adata = {'RCN1': ['0',],#Immune
                     'RCN2': ['1'],#Keratinocytes
                     'RCN3': ['2', ],#Melanocytes
                     'RCN4': ['3', ],#CD8 T cells
                     'RCN5': ['4'],#Endothelial cells
                     }

  
adata = sm.hl.rename(adata, rename=rename_dict_adata, from_column='neigh_kmeans', to_column='RCNs')

sm.pl.pie(adata, phenotype='phenotype_L2', group_by='RCNs', ncols=3)

sm.pl.stacked_barplot (adata, x_axis='imageid', y_axis='RCNs',
                       **{'edgecolor': 'none'})

sm.pl.pie(adata, phenotype='RCNs', group_by='imageid', ncols=3)

sm.pl.stacked_barplot (adata, x_axis='distance_from_epidermis2500_binned', y_axis='RCNs',
                       subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                                       '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                       '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                      order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)', '[10000, 12500)', 
                                       '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                      '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                      subset_yaxis=['RCN5', 'RCN6', 'RCN7', 'RCN8', 'RCN9', 'RCN10'],
                      order_yaxis=['RCN5', 'RCN6', 'RCN7', 'RCN8', 'RCN9', 'RCN10'],
                       **{'edgecolor': 'none'})

sm.pl.stacked_barplot (adata, x_axis='distance_from_epidermis2500_binned', y_axis='RCNs',
                       subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)','[10000, 12500)',],
                       order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)','[10000, 12500)', ],
                       #subset_yaxis=['RCN5', 'RCN6', 'RCN7', 'RCN8', 'RCN9', 'RCN10'],
                       #order_yaxis=['RCN5', 'RCN6', 'RCN7', 'RCN8', 'RCN9', 'RCN10'],
                       **{'edgecolor': 'none'})




#######
from scipy import stats
