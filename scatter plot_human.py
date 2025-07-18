# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 10:03:56 2025

@author: shp918
"""

default_color = '#d3d3d3'

customColors = {
    'Trms' : '#8559a5',
    'CD103-': '#259DD0',  # <--- Added quotes here
    'Tregs' : '#d62729',
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


# The rest of your plotting code
sm.pl.spatial_scatterPlot (adata,
    colorBy = ['phenotype_L2'],
    topLayer = ['Trms', 'Tregs', 'CD103-'],
    subset = 'LSP15219',
    figsize=(10,10,),
    s=0.5,
    plotLegend=True,
    fontsize=3,
    dpi=300,
    vmin=0,
    vmax=1,
    customColors=customColors,
    fileName='ScatterPlot_LSP15219.png',
    saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\scatterplots"
    )



Melanocytes: #259DD0 - A vibrant sky blue / light blue #
Endothelial cells: #6BA292 - A muted blue-green / seafoam green #
Macrophage: #8C4843 - A deep reddish-brown / brick red #
Keratinocytes: #2077b5 - A medium to dark true blue #
Trms: #d62729 - A strong, bright red #
T helper: #f69696 - A soft, light pink #
DCs: #8559a5 - A medium to dark purple / plum #
Keratinocyte: #afc7e5 - A light, somewhat muted blue #
Tregs: #c4afd3 - A light lavender / pale purple #
Immune cells: #f57e20 - A bright, bold orange #
CD103-: #d87ab1 - A medium pinkish-purple / rosy mauve #
B cells: #1b9e77 - A rich, deep teal / blue-green #
T cells: #7570b3 - A medium blue-purple / muted violet #
Unknown: #999999 - A neutral medium gray #
