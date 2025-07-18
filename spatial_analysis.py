# -*- coding: utf-8 -*-
"""
Created on Thu May 29 16:19:55 2025

@author: shp918
"""

adata = sm.tl.spatial_distance (adata, phenotype='phenotype_L2')
sm.pl.spatial_distance (adata, phenotype='phenotype_L2',figsize=(5,4))

sm.pl.spatial_distance (adata, method='distribution',distance_from='Tregs',distance_to = 'CD103-', phenotype='phenotype_L2',imageid='imageid', log=True, height=3, aspect=9/8)
sm.pl.spatial_distance (adata, method='distribution',distance_from='Tregs',distance_to = 'Trms', phenotype='phenotype_L2',imageid='imageid', log=True, height=3, aspect=9/8)

sm.pl.spatial_distance (bdata, method='distribution',distance_from='Tregs',distance_to = 'Trms',
                        imageid='Location_combined', log=True, height=3, 
                        saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results", 
                        fileName='Distance_Treg_to_Trms.pdf', aspect=9/8)


sm.pl.spatial_distance (bdata, method='distribution',distance_from='Tregs',distance_to = 'Trms', 
                        imageid='Time_point', log=True, height=3, 
                        saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results", 
                        fileName='Distance_Treg_to_Trms_Timepoint.pdf', 
                        aspect=9/8)

subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)']

sm.pl.spatial_distance (bdata, method='distribution',distance_from='Tregs',distance_to = 'Trms', imageid='Location_combined', log=False, height=3, aspect=9/8)



#Human plotted
distance_to = ['CD103-', 'Trms']
sm.pl.spatial_distance (adata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, phenotype='phenotype_L2',
                        log=True, imageid='histopath_roi',
                        saveDir=r"C:/Users/shp918/HMS Dropbox/Shishir Pant/Shishir-Jason/BPO_III/results/Results_II/spatial_analysis", 
                        fileName='Distance_Treg_to_CD8_histopath_roi_distance_numeric.pdf', 
                        height=3, aspect=18/16)



distance_to = ['CD103-', 'Trms']
sm.pl.spatial_distance (adata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, phenotype='phenotype_L2',
                        log=True, imageid='histopath_roi',
                        #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results", 
                        #fileName='Distance_Treg_to_Trms_Location_combined_numeric.pdf', 
                        height=3, aspect=9/8)





distance_to = ['CD103-', 'Trms']
sm.pl.spatial_distance (bdata, method='numeric', 
                        distance_from='Tregs', distance_to=distance_to, 
                        log=True, imageid='Time_point', 
                        saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results", 
                        fileName='Distance_Treg_to_Trms_Timepoint_numeric.pdf', 
                        height=3, aspect=9/8)


adata = sm.tl.spatial_count(adata, phenotype='phenotype', method='radius', radius=80, label='spatial_count')
adata = sm.tl.spatial_cluster(adata, df_name='spatial_count', method='kmeans', k=3, label='neigh_kmeans',)
# adata = sm.tl.spatial_cluster(adata, df_name='spatial_count', method='kmeans', k=30, label='neigh_kmeans')
sm.pl.stacked_barplot (adata, x_axis='neigh_kmeans', y_axis='phenotype')

sm.pl.stacked_barplot (adata, x_axis='neigh_kmeans', y_axis='phenotype',  subset_yaxis = ['Tregs', 'T helper', 'Trms', 'CD103-', 'Macrophage', 'DCs', 'Immune cells'],
    plot_tool='matplotlib', matplotlib_cmap=cmap, 
 #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Project_Analysis\iBIP2\plots", 
 **{'edgecolor': 'none'})


sm.pl.groupCorrelation(adata, groupBy='phenotype', condition='neigh_kmeans', figsize=(9,6))

rename_dict_adata = {'RCN1': ['0'], #Macrophjage
                     'RCN2': ['1'],#T cells
                     'RCN3': ['2'],#Tregs
                     }

 adata = sm.hl.rename(adata, rename=rename_dict_adata, from_column='neigh_kmeans', to_column='RCNs')

sm.pl.pie(adata, phenotype='RCNs', group_by='Time_point', ncols=3)

sm.pl.stacked_barplot (adata, x_axis='Time_point', y_axis='RCNs',
                       **{'edgecolor': 'none'})


column_to_filter = 'Time_point'
mask = adata.obs[column_to_filter] != 'Other'
adata = adata[mask, :].copy()
print(bdata.obs['Time_point'].value_counts())


column_to_filter = 'Location_combined'
mask = adata.obs[column_to_filter] != 'Other'
adata = adata[mask, :].copy()
print(bdata.obs['Location_combined'].value_counts())
