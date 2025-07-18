# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 13:00:26 2025

@author: shp918
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# --- Data Preparation Steps (REQUIRED if not already run in current session) ---
# Assuming 'adata' is your loaded AnnData object and it has 'imageid', 'distance_from_epidermis2500_binned', and 'phenotypeSHP' in adata.obs

# 1. Calculate absolute counts and initial proportions per imageid and phenotypeSHP
data_for_counts = adata.obs[['imageid', 'distance_from_epidermis2500_binned', 'phenotype_L2']].copy()
data_for_counts = data_for_counts.dropna(subset=['imageid', 'distance_from_epidermis2500_binned', 'phenotype_L2'])

# Corrected groupby order for unstacking phenotypeSHP into columns
absolute_counts_per_image = data_for_counts.groupby(['imageid', 'distance_from_epidermis2500_binned', 'phenotype_L2']).size().unstack(fill_value=0)
total_cells_per_image = data_for_counts.groupby(['imageid', 'distance_from_epidermis2500_binned']).size()

# Corrected: Removed 'level' parameter to resolve MultiIndex ambiguity
initial_proportions_df = absolute_counts_per_image.div(total_cells_per_image, axis=0)


# --- 2. Normalize these proportions to TUMOR cells ---
normalized_proportions_to_plot = initial_proportions_df.copy()

# Define the new reference cell type
reference_cell_type = 'Melanocytes'

# Ensure the new reference_cell_type column exists in initial_proportions_df for normalization
if reference_cell_type not in normalized_proportions_to_plot.columns:
    raise ValueError(f"The '{reference_cell_type}' cell type is not present in your data's phenotype_L2 for normalization.")

# Get a list of all cell type columns to normalize (excluding the reference_cell_type itself)
cell_type_columns_to_normalize = [col for col in normalized_proportions_to_plot.columns if col != reference_cell_type]

for col_name in cell_type_columns_to_normalize:
    new_column_name = f"{col_name}_norm_to_{reference_cell_type}"
    normalized_proportions_to_plot[new_column_name] = normalized_proportions_to_plot[col_name] / normalized_proportions_to_plot[reference_cell_type]
    # Handle division by zero (where reference_cell_type proportion was 0)
    normalized_proportions_to_plot[new_column_name] = normalized_proportions_to_plot[new_column_name].replace([np.inf, -np.inf], np.nan)

normalized_proportions_to_plot = normalized_proportions_to_plot.drop(columns=initial_proportions_df.columns, errors='ignore')


# --- 3. Prepare the data for plotting (Melt for Seaborn) ---
df_melted_for_plotting = normalized_proportions_to_plot.reset_index()

df_melted_for_plotting = df_melted_for_plotting.melt(
    id_vars=['imageid', 'distance_from_epidermis2500_binned'],
    var_name='Cell_Type',
    value_name=f'Proportion_Relative_to_{reference_cell_type}' # Dynamic value name
)

df_melted_for_plotting['Cell_Type'] = df_melted_for_plotting['Cell_Type'].str.replace(f'_norm_to_{reference_cell_type}', '')
df_melted_for_plotting = df_melted_for_plotting.dropna(subset=[f'Proportion_Relative_to_{reference_cell_type}'])


# --- Define the subset of cell types to plot (as in previous interaction) ---
subset_yaxis_types = ['Tregs', 'Trms', 'CD103-']

# Filter the melted DataFrame for this cell type subset
df_subset_plotting = df_melted_for_plotting[
    df_melted_for_plotting['Cell_Type'].isin(subset_yaxis_types)
].copy()

#['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)','[10000, 12500)', '[12500, 15000)', '[15000, 17500)', '[17500, 20000)','[20000, 22500)', '[22500, 25000)', '[25000, 27500)', '[27500, 30000)',
    

# --- Define and filter the x-axis order for your desired range ---
full_binned_order = ['[0, 2500)', '[2500, 5000)', '[5000, 7500)', '[7500, 10000)']


# Filter df_subset_plotting to include only these bins
df_subset_plotting_filtered_xaxis = df_subset_plotting[
    df_subset_plotting['distance_from_epidermis2500_binned'].isin(full_binned_order)
].copy()

# Ensure the x-axis column is categorical with the desired order
df_subset_plotting_filtered_xaxis['distance_from_epidermis2500_binned'] = pd.Categorical(
    df_subset_plotting_filtered_xaxis['distance_from_epidermis2500_binned'],
    categories=full_binned_order,
    ordered=True
)


# --- Plotting the Box Plot (Filtered x-axis range and improved colors) ---

fig_width = len(full_binned_order) * len(subset_yaxis_types) * 0.2 + 4
plt.figure(figsize=(fig_width, 7))
sns.set_theme(style="whitegrid")

color_palette = 'Dark2' # Or 'tab10', 'Set1', 'Paired'

ax = sns.boxplot(
    data=df_subset_plotting_filtered_xaxis,
    x='distance_from_epidermis2500_binned',
    y=f'Proportion_Relative_to_{reference_cell_type}', # Dynamic y-axis
    hue='Cell_Type',
    order=full_binned_order,
    hue_order=subset_yaxis_types,
    palette=color_palette,
    fliersize=0
)

sns.stripplot(
    data=df_subset_plotting_filtered_xaxis,
    x='distance_from_epidermis2500_binned',
    y=f'Proportion_Relative_to_{reference_cell_type}', # Dynamic y-axis
    hue='Cell_Type',
    order=full_binned_order,
    hue_order=subset_yaxis_types,
    dodge=True,
    size=5,
    color='black',
    alpha=0.6,
    jitter=0.15,
    ax=ax,
    legend=False
)

# --- Plot Customization ---
plt.title(f'Proportion of Selected Cell Types Relative to {reference_cell_type} by Distance from Epidermis', fontsize=16)
plt.xlabel('Distance From Epidermis (microns)', fontsize=14)
plt.ylabel(f'Proportion Relative to {reference_cell_type}', fontsize=14) # Dynamic y-label

plt.xticks(rotation=45, ha='right', fontsize=10)
ax.tick_params(axis='y', labelsize=12)

handles, labels = ax.get_legend_handles_labels()
num_hues = len(subset_yaxis_types)
plt.legend(handles[:num_hues], labels[:num_hues],
           title='Cell Type',
           bbox_to_anchor=(1.03, 1), loc='upper left',
           title_fontsize='13', fontsize='11')

plt.tight_layout(rect=[0, 0, 0.85, 1])


# --- Optional: Saving the plot ---
#import os
save_dir = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\boxplot"
if not os.path.exists(save_dir):
 os.makedirs(save_dir)
file_name = 'boxplot_subset_norm_to_{melanocytes}_up_to_12500}.pdf'
file_path = os.path.join(save_dir, file_name)
plt.savefig(file_path, bbox_inches='tight')
print(f"Plot saved to: {file_path}")
plt.show()









########FOR HISTOROI#####
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os # Make sure to import os for saving plots

# --- Data Preparation Steps for 'histopath_roi' ---
# Assuming 'adata' is your loaded AnnData object and it has 'imageid', 'histopath_roi', and 'phenotype_L2' in adata.obs


# Create a copy for 'histopath_roi' specific processing
data_for_counts_histopath = adata.obs[['imageid', 'histopath_roi', 'phenotype_L2']].copy()
data_for_counts_histopath = data_for_counts_histopath.dropna(subset=['imageid', 'histopath_roi', 'phenotype_L2'])

# Corrected groupby order for unstacking phenotype_L2 into columns
absolute_counts_per_image_histopath = data_for_counts_histopath.groupby(['imageid', 'histopath_roi', 'phenotype_L2']).size().unstack(fill_value=0)
total_cells_per_image_histopath = data_for_counts_histopath.groupby(['imageid', 'histopath_roi']).size()

# Corrected: Removed 'level' parameter to resolve MultiIndex ambiguity
initial_proportions_df_histopath = absolute_counts_per_image_histopath.div(total_cells_per_image_histopath, axis=0)


# --- 2. Normalize these proportions to TUMOR cells for 'histopath_roi' ---
normalized_proportions_to_plot_histopath = initial_proportions_df_histopath.copy()

# Define the reference cell type (can be the same as before if desired)
reference_cell_type_histopath = 'Melanocytes'

# Ensure the new reference_cell_type column exists
if reference_cell_type_histopath not in normalized_proportions_to_plot_histopath.columns:
    raise ValueError(f"The '{reference_cell_type_histopath}' cell type is not present in your data's phenotype_L2 for normalization (for histopath_roi analysis).")

# Get a list of all cell type columns to normalize (excluding the reference_cell_type itself)
cell_type_columns_to_normalize_histopath = [col for col in normalized_proportions_to_plot_histopath.columns if col != reference_cell_type_histopath]

for col_name in cell_type_columns_to_normalize_histopath:
    new_column_name = f"{col_name}_norm_to_{reference_cell_type_histopath}"
    normalized_proportions_to_plot_histopath[new_column_name] = normalized_proportions_to_plot_histopath[col_name] / normalized_proportions_to_plot_histopath[reference_cell_type_histopath]
    # Handle division by zero (where reference_cell_type proportion was 0)
    normalized_proportions_to_plot_histopath[new_column_name] = normalized_proportions_to_plot_histopath[new_column_name].replace([np.inf, -np.inf], np.nan)

normalized_proportions_to_plot_histopath = normalized_proportions_to_plot_histopath.drop(columns=initial_proportions_df_histopath.columns, errors='ignore')


# --- 3. Prepare the data for plotting (Melt for Seaborn) for 'histopath_roi' ---
df_melted_for_plotting_histopath = normalized_proportions_to_plot_histopath.reset_index()

df_melted_for_plotting_histopath = df_melted_for_plotting_histopath.melt(
    id_vars=['imageid', 'histopath_roi'], # Changed 'distance_from_epidermis2500_binned' to 'histopath_roi'
    var_name='Cell_Type',
    value_name=f'Proportion_Relative_to_{reference_cell_type_histopath}' # Dynamic value name
)

df_melted_for_plotting_histopath['Cell_Type'] = df_melted_for_plotting_histopath['Cell_Type'].str.replace(f'_norm_to_{reference_cell_type_histopath}', '')
df_melted_for_plotting_histopath = df_melted_for_plotting_histopath.dropna(subset=[f'Proportion_Relative_to_{reference_cell_type_histopath}'])


# --- Define the subset of cell types to plot (same as previous interaction) ---
subset_yaxis_types_histopath = ['Tregs', 'Trms', 'CD103-']

# Filter the melted DataFrame for this cell type subset
df_subset_plotting_histopath = df_melted_for_plotting_histopath[
    df_melted_for_plotting_histopath['Cell_Type'].isin(subset_yaxis_types_histopath)
].copy()

# --- Define and filter the x-axis order for 'histopath_roi' ---
# Use the unique values from your 'histopath_roi' column, you can define a specific order if desired.
# For example, if you want 'normal' and 'precursor' first, then 'VGP', then 'Other':
histopath_roi_order = ['normal', 'precursor', 'VGP', 'Other'] # Adjust this order as needed
# Or dynamically get unique values:
# histopath_roi_order = df_subset_plotting_histopath['histopath_roi'].unique().tolist()
# Note: If dynamically getting unique values, the order might not be intuitive.

# Filter df_subset_plotting_histopath to include only these categories (if any are missing)
df_subset_plotting_histopath_filtered_xaxis = df_subset_plotting_histopath[
    df_subset_plotting_histopath['histopath_roi'].isin(histopath_roi_order)
].copy()

# Ensure the x-axis column is categorical with the desired order
df_subset_plotting_histopath_filtered_xaxis['histopath_roi'] = pd.Categorical(
    df_subset_plotting_histopath_filtered_xaxis['histopath_roi'],
    categories=histopath_roi_order,
    ordered=True
)


# --- Plotting the Box Plot for 'histopath_roi' ---

fig_width_histopath = len(histopath_roi_order) * len(subset_yaxis_types_histopath) * 0.2 + 4
plt.figure(figsize=(fig_width_histopath, 7))
sns.set_theme(style="whitegrid")

color_palette_histopath = 'Dark2' # Or 'tab10', 'Set1', 'Paired'

ax_histopath = sns.boxplot(
    data=df_subset_plotting_histopath_filtered_xaxis,
    x='histopath_roi', # Changed x-axis to 'histopath_roi'
    y=f'Proportion_Relative_to_{reference_cell_type_histopath}', # Dynamic y-axis
    hue='Cell_Type',
    order=histopath_roi_order, # Used the new order
    hue_order=subset_yaxis_types_histopath,
    palette=color_palette_histopath,
    fliersize=0
)

sns.stripplot(
    data=df_subset_plotting_histopath_filtered_xaxis,
    x='histopath_roi', # Changed x-axis to 'histopath_roi'
    y=f'Proportion_Relative_to_{reference_cell_type_histopath}', # Dynamic y-axis
    hue='Cell_Type',
    order=histopath_roi_order, # Used the new order
    hue_order=subset_yaxis_types_histopath,
    dodge=True,
    size=5,
    color='black',
    alpha=0.6,
    jitter=0.15,
    ax=ax_histopath,
    legend=False
)

# --- Plot Customization ---
plt.title(f'Proportion of Selected Cell Types Relative to {reference_cell_type_histopath} by Histopathology ROI', fontsize=16) # Updated title
plt.xlabel('Histopathology ROI', fontsize=14) # Updated x-label
plt.ylabel(f'Proportion Relative to {reference_cell_type_histopath}', fontsize=14) # Dynamic y-label

plt.xticks(rotation=45, ha='right', fontsize=10)
ax_histopath.tick_params(axis='y', labelsize=12)

handles_histopath, labels_histopath = ax_histopath.get_legend_handles_labels()
num_hues_histopath = len(subset_yaxis_types_histopath)
plt.legend(handles_histopath[:num_hues_histopath], labels_histopath[:num_hues_histopath],
           title='Cell Type',
           bbox_to_anchor=(1.03, 1), loc='upper left',
           title_fontsize='13', fontsize='11')

plt.tight_layout(rect=[0, 0, 0.85, 1])


# --- Optional: Saving the plot ---
save_dir_histopath = r"C:\Users\shp918\HMS Dropbox\Shishir Pant\Shishir-Jason\BPO_III\results\Results_II\boxplot_histopath_roi" # New save directory
if not os.path.exists(save_dir_histopath):
   os.makedirs(save_dir_histopath)
file_name_histopath = f'boxplot_subset_norm_to_{reference_cell_type_histopath}_histopath_roi.pdf' # New file name
file_path_histopath = os.path.join(save_dir_histopath, file_name_histopath)
plt.savefig(file_path_histopath, bbox_inches='tight')
print(f"Plot saved to: {file_path_histopath}")
plt.show()