# -*- coding: utf-8 -*-
"""
Created on Fri May 30 15:51:25 2025

@author: shp918
"""


default_color = '#d3d3d3'

customColors = {
    'Trms' : '#06d6a0',
    'CD103-': default_color,
    'Tregs': '#ef476f',
    'Endothelial cells' : default_color,
    'T helper': default_color,
     'DCs' : default_color,
    'Macrophage' : default_color,
    'Immune cells': default_color,
    'Keratinocytes': default_color,
    'Unknown': default_color,
    'Melanocytes': default_color,
}


#LSP28001--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP28001--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP28001.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")


#LSP27996--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP27996--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP27996.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")

#LSP17991--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP17991--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP17991.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")


#LSP28036--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP28036--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP28036.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")

#LSP28021--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP28021--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP28021.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")

#LSP28026--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP28026--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP28026.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")

#LSP28011--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP28011--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP28011.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")

#LSP28041--unmicst_cellRing
sm.pl.spatial_scatterPlot (adata, 
                 colorBy = ['phenotype'],
                 topLayer = ['Trms', 'Tregs', 'T helper', 'CD8103-'],
                 subset = 'LSP28041--unmicst_cellRing',
                 figsize=(10,10,),
                 s=0.5,
                 plotLegend=True,
                 fontsize=3,
                 dpi=300,
                 vmin=0,
                 vmax=1,
                 customColors=customColors,
                 fileName='scimapScatterPlot_LSP28041.png',
                 saveDir=r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results")

