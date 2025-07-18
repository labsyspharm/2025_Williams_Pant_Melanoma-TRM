# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 15:44:53 2025

@author: shp918
"""


from matplotlib import colors as mcolors

cell_type_colors = {
    'Melanocytes': '#259DD0',
    'Endothelial cells': '#6BA292',
    'Macrophage': '#8C4843',
    'Keratinocytes': '#2077b5',
    'Trms': '#d62729',
    'T helper': '#f69696',
    'DCs': '#8559a5',
    'Keratinocyte': '#afc7e5',
    'Tregs': '#c4afd3',
    'Immune cells': '#f57e20',
    'CD103-': '#d87ab1',
    'B cells': '#1b9e77',
    'T cells': '#7570b3',
    'Unknown': '#999999',
}

# Use the correct dictionary name
cmap = mcolors.ListedColormap([cell_type_colors[key] for key in cell_type_colors])


# Plot the stacked barplot with the color map
sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L2', method='percent', 
                      subset_yaxis= ['Immune cells', 'Macrophage', 'T cells', 'DCs', 'Tregs',  'CD103-', 'T helper', 'Trms','B cells', ],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='Immune_cell_proportions_abs_human_imageid.pdf',
                      figsize=(12, 6),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='histopath_2', y_axis='phenotype_L2', method='percent', 
                      #subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      #fileName='cell_proportions_abs_human_imageid_2.pdf',
                      figsize=(12, 6),  **{'edgecolor': 'none'})




sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L2', method='percent', 
                      #subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      #fileName='cell_proportions_percent_human_imageid.pdf',
                      figsize=(12, 6,  **{'width': 0.7, 'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L2', method='percent', 
                      subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='T_cell_subset_percent_human_imageid.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L2', method='absolute', 
                      subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='T_cell_subset_absolute_human_imageid.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='histopath_roi', y_axis='phenotype_L2', method='percent', 
                      subset_xaxis= ['normal', 'precursor',  'VGP'],
                      #subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      order_xaxis= ['normal', 'precursor',  'VGP'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='cell_proportions_percent_human_histopathroi.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='histopath_roi', y_axis='phenotype_L2', method='absolute', 
                      subset_xaxis= ['normal', 'precursor',  'VGP'],
                      #subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      order_xaxis= ['normal', 'precursor',  'VGP'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='cell_proportions_absolute_human_histopathroi.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='histopath_roi', y_axis='phenotype_L2', method='absolute', 
                      #subset_xaxis= ['normal', 'precursor',  'VGP'],
                      subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      order_xaxis= ['normal', 'precursor',  'VGP', 'Other'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='T_cell_subset_absolute_histopathroi.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='histopath_roi', y_axis='phenotype_L2', method='percent', 
                      #subset_xaxis= ['normal', 'precursor',  'VGP'],
                      subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                      order_xaxis= ['normal', 'precursor',  'VGP'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='T_cell_subset_percent_histopathroi.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='distance_from_epidermis2500_binned', y_axis='phenotype_L2', method='absolute', 
                      subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                                      '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                      '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                     order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                                     '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                     '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='cell_proportions_abs_human_distance_from_epidermis2500_binned.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='distance_from_epidermis2500_binned', y_axis='phenotype_L2', method='percent', 
                      subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', 
                                      '[7500, 10000)'],
                     order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)',
                                     '[7500, 10000)'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='cell_proportions_abs_human_percent_subset_binned.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})

sm.pl.stacked_barplot(adata, x_axis='distance_from_epidermis2500_binned', y_axis='phenotype_L2', method='percent', 
                      subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                                      '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                      '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                     order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                                     '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                     '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap,
                      saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      fileName='cell_proportions_percent_human_distance_from_epidermis2500_binned.pdf',
                      figsize=(6, 4),  **{'edgecolor': 'none'})
             
sm.pl.stacked_barplot(adata, x_axis='distance_from_epidermis2500_binned', y_axis='phenotype_L2', method='absolute', 
                       subset_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                                       '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                       '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                      order_xaxis = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[10000, 12500)', 
                                      '[7500, 10000)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)',
                                      '[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)'],
                       subset_yaxis= ['Trms', 'Tregs',  'CD103-', 'T helper', 'T cells'],
                       plot_tool='matplotlib', matplotlib_cmap=cmap,
                       saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                       fileName='T_cell_subset_absolute_human_distance_from_epidermis2500_binned.pdf',
                       figsize=(6, 4),  **{'edgecolor': 'none'})              


# Plot the stacked barplot with the color map
sm.pl.stacked_barplot(adata, x_axis='imageid', y_axis='phenotype_L2', method='percent',
                      #subset_yaxis= ['Trms', 'Tregs', 'CD103-', 'T helper', 'T cells'],
                      plot_tool='matplotlib', matplotlib_cmap=cmap, # Use your defined 'cmap'
                      #saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II",
                      #fileName='cell_proportions_abs_human_imageid.pdf',
                      figsize=(12, 46,
                      # --- CRUCIAL CHANGE: Adjust the 'width' parameter ---
                      **{'edgecolor': 'none', 'width': 0.7} # Smaller width (default is usually 0.8 or 1.0)
                     )


