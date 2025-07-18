# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 15:16:15 2025

@author: shp918
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
import os

# --- Define the specific imageid and target cell type for THIS plot ---
target_image_id = 'LSP13179'
target_cell_type = 'Tregs' # *** CHANGED FOR Tregs ***

# --- Plotting parameters (ADJUST THESE TO FINE-TUNE) ---
# These are the "knobs" you'd turn to change the plot's appearance.
n_bins = 200
sigma_smoothing_target = 10
num_contour_levels = 15
contour_min_density_threshold = 0.005

scatter_dot_size_all_cells = 0.1
all_cells_scatter_color = 'lightgrey'
all_cells_scatter_alpha = 0.6

contour_line_color = 'black'
contour_line_width = 0.7

scatter_dot_size_target = 0.4
target_cell_scatter_color = 'red' # *** CHANGED COLOR FOR Tregs, you can pick any color! ***

# --- Save directory (make sure this path exists on your system) ---
save_base_dir = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\contour_plots_abundance"


# --- 1. Prepare Data for the target image ---
adata_image_data = adata[adata.obs['imageid'] == target_image_id].copy()

# Get ALL spatial coordinates for the subsetted image to define the canvas
if 'X_centroid' in adata_image_data.obs and 'Y_centroid' in adata_image_data.obs:
    all_x_coords = adata_image_data.obs['X_centroid'].values
    all_y_coords = adata_image_data.obs['Y_centroid'].values
else:
    raise ValueError("Spatial coordinates (adata.obs['X_centroid'] and/or 'Y_centroid') not found in the AnnData object.")

# Determine the overall grid boundaries for the entire image
x_min, x_max = all_x_coords.min(), all_x_coords.max()
y_min, y_max = all_y_coords.min(), all_y_coords.max()

# --- Calculate aspect ratio for dynamic figsize ---
x_range = x_max - x_min
y_range = y_max - y_min

if y_range < 0:
    y_range = abs(y_range)

base_width = 10 # inches
figure_width = base_width
figure_height = base_width * (y_range / x_range)


# --- 2. Create the "intensity" / presence array for target cells ---
is_target_cell = (adata_image_data.obs['phenotype_L2'] == target_cell_type).astype(int)

# Get coordinates specifically for target cells for the scatter overlay
x_coords_target = adata_image_data.obs['X_centroid'][is_target_cell == 1].values
y_coords_target = adata_image_data.obs['Y_centroid'][is_target_cell == 1].values

# --- Check if there are any target cells ---
if len(x_coords_target) == 0:
    print(f"No {target_cell_type} cells found in image {target_image_id}. Skipping plot generation.")
    # Exit or handle gracefully if no target cells are found
    # You might want to remove this `exit()` if you want the script to continue
    # even if one plot can't be generated.
    exit()


# --- 3. Calculate 2D Histogram for Target Cell Abundance (Weighted) ---
density_map_target, x_edges, y_edges = np.histogram2d(all_x_coords, all_y_coords,
                                                      bins=[n_bins, n_bins],
                                                      range=[[x_min, x_max], [y_min, y_max]],
                                                      weights=is_target_cell)

# --- 4. Smooth the Target Cell Density Map ---
density_map_target_smoothed = gaussian_filter(density_map_target.T, sigma=sigma_smoothing_target)

# Create meshgrid for plotting contours
X_grid, Y_grid = np.meshgrid((x_edges[:-1] + x_edges[1:]) / 2,
                             (y_edges[:-1] + y_edges[1:]) / 2)


# --- 5. Plotting the Combined Visual ---

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
fig, ax = plt.subplots(figsize=(figure_width, figure_height))

# --- Plot the background tissue outline using ALL cells as very small scatter points ---
ax.scatter(all_x_coords, all_y_coords, c=all_cells_scatter_color,
           s=scatter_dot_size_all_cells, alpha=all_cells_scatter_alpha,
           rasterized=True)


# --- Plot the Contour Lines on top ---
contour_levels_raw = np.linspace(density_map_target_smoothed.min(),
                                 density_map_target_smoothed.max(),
                                 num_contour_levels)

contour_levels = contour_levels_raw[contour_levels_raw >= contour_min_density_threshold]

if len(contour_levels) > 1:
    ax.contour(X_grid, Y_grid, density_map_target_smoothed,
               levels=contour_levels,
               colors=contour_line_color,
               linewidths=contour_line_width,
               alpha=0.9)
else:
    print(f"Warning: Not enough distinct density values above {contour_min_density_threshold} to plot meaningful contours for {target_cell_type}. Consider adjusting smoothing, binning, or `contour_min_density_threshold`.")


# --- Overlay original Target Cell scatter points (small dots) ---
ax.scatter(x_coords_target, y_coords_target, c=target_cell_scatter_color,
           s=scatter_dot_size_target, alpha=0.7, label=f'{target_cell_type} Cells', rasterized=True)


# --- Plot Customization ---
ax.invert_yaxis()

ax.set_title(f'{target_cell_type} Cell Abundance (Image {target_image_id})', fontsize=16)
ax.set_xlabel('X Coordinate', fontsize=14)
ax.set_ylabel('Y Coordinate', fontsize=14)

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_max, y_min)

ax.set_aspect('equal', adjustable='box')

plt.tight_layout()

ax.legend(loc='upper right', frameon=True)


# --- Saving the plot ---
os.makedirs(save_base_dir, exist_ok=True)
file_name = f'{target_cell_type}_Abundance_ContoursScatter_{target_image_id}_single_plot.png' # Unique filename
file_path = os.path.join(save_base_dir, file_name)

plt.savefig(file_path, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {file_path}")
plt.show()
plt.close(fig) # Close the figure