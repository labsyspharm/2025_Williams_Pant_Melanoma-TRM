# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 13:18:23 2025

@author: shp918
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
import os

# Assuming 'adata' is your loaded AnnData object

# --- Define the specific imageid and target cell type ---
target_image_id = 'LSP13179' # Current image ID
target_cell_type = 'Trms' # The cell type you want to show abundance for

# --- Plotting parameters (ADJUST THESE TO FINE-TUNE) ---
n_bins = 200        # Resolution of the density grid (e.g., 100-400). Higher = finer detail.
sigma_smoothing_trms = 10 # Smoothing for Trms density. Experiment with 5-20.

num_contour_levels = 15 # INCREASED: Number of contour lines (e.g., 10-25).
contour_min_density_threshold = 0.0 # Adjusted: Minimum density value for contours to appear. Adjust as needed.

# Background (All Cells) Scatter Parameters
scatter_dot_size_all_cells = 0.1 # Very very small dots for overall tissue background (e.g., 0.05-0.2)
all_cells_scatter_color = 'lightgrey' # Color for the tissue background
all_cells_scatter_alpha = 0.6 # Transparency for the tissue background scatter

# Trms Heatmap Parameters
heatmap_cmap = 'Reds' # Colormap for the Trms abundance heatmap
heatmap_alpha = 0.8   # Adjusted: Transparency of the Trms abundance heatmap (e.g., 0.7-1.0)
heatmap_vmin_manual = 0.0001 # **CRUCIAL CHANGE:** Start heatmap color from this density. ADJUST THIS VALUE!
                          # Try 0.05, 0.1, 0.2, 0.5 depending on your density scale.
                          # Values below this will be white.

# Contour Line Parameters
contour_line_color = 'black' # Color of the contour lines
contour_line_width = 0.7 # Width of the contour lines

# Trms Cells Overlay Scatter Parameters
scatter_dot_size_trms = 0.4 # Size of dots for target cells (Trms) (e.g., 0.3-0.8)
trms_scatter_color = 'Purple' # Color for the overlaid Trms cells


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
x_coords_trms = adata_image_data.obs['X_centroid'][is_target_cell == 1].values
y_coords_trms = adata_image_data.obs['Y_centroid'][is_target_cell == 1].values


# --- 3. Calculate 2D Histogram for Target Cell Abundance (Weighted) ---
density_map_target, x_edges, y_edges = np.histogram2d(all_x_coords, all_y_coords,
                                                      bins=[n_bins, n_bins],
                                                      range=[[x_min, x_max], [y_min, y_max]],
                                                      weights=is_target_cell)

# --- 4. Smooth the Target Cell Density Map ---
density_map_target_smoothed = gaussian_filter(density_map_target.T, sigma=sigma_smoothing_trms)

# Create meshgrid for plotting contours (based on the density map grid)
X_grid, Y_grid = np.meshgrid((x_edges[:-1] + x_edges[1:]) / 2,
                             (y_edges[:-1] + y_edges[1:]) / 2)


# --- 5. Plotting the Combined Visual ---

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.figure(figsize=(figure_width, figure_height))
ax = plt.gca()

# --- Plot the background tissue outline using ALL cells as very small scatter points ---
ax.scatter(all_x_coords, all_y_coords, c=all_cells_scatter_color,
           s=scatter_dot_size_all_cells, alpha=all_cells_scatter_alpha,
           rasterized=True)

# --- Plot the Target Cell Abundance Heatmap (Colored Overlay) ---
heatmap_plot = ax.imshow(density_map_target_smoothed,
                         extent=[x_min, x_max, y_max, y_min],
                         origin='upper',
                         cmap=heatmap_cmap,
                         alpha=heatmap_alpha,
                         vmin=heatmap_vmin_manual, # This controls the white background!
                         vmax=density_map_target_smoothed.max()
                        )

# Add color bar for the Trms abundance
plt.colorbar(heatmap_plot, ax=ax, label=f'{target_cell_type} Cell Density', shrink=0.7)


# --- Plot the Contour Lines on top ---
# Calculate levels based on the smoothed map's range, but filter by threshold
contour_levels_raw = np.linspace(density_map_target_smoothed.min(),
                                 density_map_target_smoothed.max(),
                                 num_contour_levels)

contour_levels = contour_levels_raw[contour_levels_raw >= contour_min_density_threshold]

if len(contour_levels) > 1:
    contour_lines = ax.contour(X_grid, Y_grid, density_map_target_smoothed,
                               levels=contour_levels,
                               colors=contour_line_color,
                               linewidths=contour_line_width,
                               alpha=0.9)
else:
    print(f"Warning: Not enough distinct density values above {contour_min_density_threshold} to plot meaningful contours. Consider adjusting smoothing, binning, or `contour_min_density_threshold`.")


# --- Overlay original Target Cell scatter points (small dots) ---
ax.scatter(x_coords_trms, y_coords_trms, c=trms_scatter_color,
           s=scatter_dot_size_trms, alpha=0.7, label=f'{target_cell_type} Cells', rasterized=True)


# --- Plot Customization ---
ax.invert_yaxis()

plt.title(f'{target_cell_type} Cell Abundance (Image {target_image_id})', fontsize=16)
plt.xlabel('X Coordinate', fontsize=14)
plt.ylabel('Y Coordinate', fontsize=14)

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_max, y_min)

ax.set_aspect('equal', adjustable='box')

plt.tight_layout()

plt.legend(loc='upper right', frameon=True)
plt.show()

# --- Saving the plot ---
contour_save_dir = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\contour_plots_abundance"
os.makedirs(contour_save_dir, exist_ok=True)
file_name = f'{target_cell_type}_Abundance_Heatmap_Contour_WhiteBG_{target_image_id}.png' # New filename
file_path = os.path.join(contour_save_dir, file_name)

plt.savefig(file_path, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {file_path}")
plt.show()






############
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
import os

# Assuming 'adata' is your loaded AnnData object

# --- Define the specific imageid and target cell type ---
target_image_id = 'LSP13179' # Current image ID
target_cell_type = 'Trms' # The cell type you want to show abundance for

# --- Plotting parameters (ADJUST THESE TO FINE-TUNE) ---
n_bins = 200        # Resolution of the density grid (e.g., 100-400). Higher = finer detail.
sigma_smoothing_trms = 10 # Smoothing for Trms density. Experiment with 5-20.

num_contour_levels = 15 # Number of contour lines (e.g., 10-25).
contour_min_density_threshold = 0.001 # Minimum density value for contours to appear. Adjust as needed.

# Background (All Cells) Scatter Parameters
scatter_dot_size_all_cells = 0.1 # Very very small dots for overall tissue background (e.g., 0.05-0.2)
all_cells_scatter_color = 'lightgrey' # Color for the tissue background
all_cells_scatter_alpha = 0.6 # Transparency for the tissue background scatter

# Contour Line Parameters
contour_line_color = 'black' # Color of the contour lines
contour_line_width = 0.7 # Width of the contour lines

# Trms Cells Overlay Scatter Parameters
scatter_dot_size_trms = 0.4 # Size of dots for target cells (Trms) (e.g., 0.3-0.8)
trms_scatter_color = 'Purple' # Color for the overlaid Trms cells


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
x_coords_trms = adata_image_data.obs['X_centroid'][is_target_cell == 1].values
y_coords_trms = adata_image_data.obs['Y_centroid'][is_target_cell == 1].values


# --- 3. Calculate 2D Histogram for Target Cell Abundance (Weighted) ---
density_map_target, x_edges, y_edges = np.histogram2d(all_x_coords, all_y_coords,
                                                      bins=[n_bins, n_bins],
                                                      range=[[x_min, x_max], [y_min, y_max]],
                                                      weights=is_target_cell)

# --- 4. Smooth the Target Cell Density Map ---
density_map_target_smoothed = gaussian_filter(density_map_target.T, sigma=sigma_smoothing_trms)

# Create meshgrid for plotting contours
X_grid, Y_grid = np.meshgrid((x_edges[:-1] + x_edges[1:]) / 2,
                             (y_edges[:-1] + y_edges[1:]) / 2)


# --- 5. Plotting the Combined Visual ---

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.figure(figsize=(figure_width, figure_height))
ax = plt.gca()

# --- Plot the background tissue outline using ALL cells as very small scatter points ---
ax.scatter(all_x_coords, all_y_coords, c=all_cells_scatter_color,
           s=scatter_dot_size_all_cells, alpha=all_cells_scatter_alpha,
           rasterized=True)

# --- Removed: Target Cell Abundance Heatmap (imshow) ---
# This layer is now excluded as per your request.
# The `ax.imshow(...)` and `plt.colorbar(...)` lines are removed.


# --- Plot the Contour Lines on top ---
contour_levels_raw = np.linspace(density_map_target_smoothed.min(),
                                 density_map_target_smoothed.max(),
                                 num_contour_levels)

contour_levels = contour_levels_raw[contour_levels_raw >= contour_min_density_threshold]

if len(contour_levels) > 1:
    contour_lines = ax.contour(X_grid, Y_grid, density_map_target_smoothed,
                               levels=contour_levels,
                               colors=contour_line_color,
                               linewidths=contour_line_width,
                               alpha=0.9)
else:
    print(f"Warning: Not enough distinct density values above {contour_min_density_threshold} to plot meaningful contours. Consider adjusting smoothing, binning, or `contour_min_density_threshold`.")


# --- Overlay original Target Cell scatter points (small dots) ---
ax.scatter(x_coords_trms, y_coords_trms, c=trms_scatter_color,
           s=scatter_dot_size_trms, alpha=0.7, label=f'{target_cell_type} Cells', rasterized=True)


# --- Plot Customization ---
ax.invert_yaxis()

plt.title(f'{target_cell_type} Cell Abundance (Image {target_image_id})', fontsize=16)
plt.xlabel('X Coordinate', fontsize=14)
plt.ylabel('Y Coordinate', fontsize=14)

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_max, y_min)

ax.set_aspect('equal', adjustable='box')

plt.tight_layout()

plt.legend(loc='upper right', frameon=True)



# --- Saving the plot ---
contour_save_dir = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\contour_plots_abundance"
os.makedirs(contour_save_dir, exist_ok=True)
file_name = f'{target_cell_type}_Abundance_ContoursOnly_{target_image_id}.png' # New filename
file_path = os.path.join(contour_save_dir, file_name)

plt.savefig(file_path, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {file_path}")
plt.show()