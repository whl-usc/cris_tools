"""
contact:    wlee9829@gmail.com
date:       2024_08_02
python:     python3.10
script:     dge_barplot.py

This Python script summarizes fold change significance of log2(x+1) RSEM or
DESEQ2 normalized gene expression counts downloaded from the UCSC Xena web 
latform.
"""
# Define version
__version__= "1.1.0"

# Version notes
__update_notes__ = """
1.1.0 
    - Remove unnecessary comments
    - Remove hard-coded paths and replace with argparse logic.

1.0.0
    - Initial commit, set up function logic and styling.
"""

# Import Packages
from datetime import datetime
from matplotlib.colors import Normalize
from scipy.stats import ttest_ind
import arg parse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys
import textwrap

###########################################################################
# 1. Set up functions

def read_input(file_path):
    data = pd.read_csv(file_path)
    
    # Count the n=0 for each tissue type, show warning if n=0 > 50% total count.
    total_counts = len(data)
    zero_counts = data.apply(lambda col: (col == 0).sum())
    zero_proportions = zero_counts / total_counts
    high_zero =zero_proportions[zero_proportions > 0.50]
    if not high_zero.empty:
        zero_counts_df = pd.DataFrame(high_zero, columns=['Zero Proportion'])
        
        print(f"Warning: These tissue types have '0' values exceeding 50% of "
            " total counts:")
        print(zero_counts_df)

    return data

def data_cleanup(dataframe):

    relative_expression = []

    for tissue_type in tissue_types:
        # Extract normal and tumor values
        normal_values = data[f'{tissue_type} Normal']
        tumor_values = data[f'{tissue_type} Tumor']

        # Drop NaN values
        combined_data = pd.DataFrame({'Normal': normal_values, 'Tumor': tumor_values}).dropna()
        normal_cleaned = combined_data['Normal']
        tumor_cleaned = combined_data['Tumor']

        # Convert log2 values back to original scale
        normal_original = 2**normal_cleaned - 1
        tumor_original = 2**tumor_cleaned - 1

        # Calculate means and fold change
        mean_normal = np.mean(normal_original)
        mean_tumor = np.mean(tumor_original)
        fold_change = mean_tumor / mean_normal

        # Perform t-test
        _, p_value = ttest_ind(normal_original, tumor_original, equal_var=False)

        # Calculate -log10(p-value)
        neg_log_p_value = -np.log10(p_value) if p_value <= 0.05 else 0

        # Determine significance stars
        if p_value <= 0.001:
            significance = '***'
        elif p_value <= 0.01:
            significance = '**'
        elif p_value <= 0.05:
            significance = '*'
        else:
            significance = ''

        # Append results
        relative_expression.append([tissue_type, fold_change, neg_log_p_value, significance])

    # Convert to DataFrame for plotting
    expression_data = pd.DataFrame(relative_expression, columns=['Tissue Type', 'Fold Change', 'Significance', 'Stars'])
    expression_data = expression_data.set_index('Tissue Type')

    # Ensure x-axis values are categorical
    expression_data.reset_index(inplace=True)
    expression_data['Tissue Type'] = pd.Categorical(expression_data['Tissue Type'], categories=tissue_types)

    # Define color map and norm for fold change
    cmap = plt.get_cmap('coolwarm')
    norm = Normalize(vmin=expression_data['Fold Change'].min(), vmax=expression_data['Fold Change'].max())

    return norm

# Set up figure and axes
fig, ax = plt.subplots(figsize=(10, 2.5))  # Adjusted figure height

# Define tissue types
tissue_types = [
    'Adrenal Gland', 'Bile Duct', 'Bladder', 'Brain', 'Breast',
    'Cervix', 'Colon', 'Esophagus', 'Head And Neck', 'Kidney',
    'Liver', 'Lung', 'Ovary', 'Pancreas', 'Prostate', 'Rectum',
    'Skin', 'Stomach', 'Testis', 'Thyroid', 'Uterus'
]

# Plot alternating background colors
for i, tissue_type in enumerate(tissue_types):
    if i % 2 == 0:
        ax.axvspan(i - 0.5, i + 0.5, facecolor='lightgray', alpha=0.25)

# Plot small boxes with color representing fold change
sns.scatterplot(
    x='Tissue Type', 
    y=0,  # Constant y-value for small boxes
    hue='Fold Change', 
    palette=cmap, 
    data=expression_data,
    s=600,  
    marker='s', 
    edgecolor=None, 
    legend=None,
    ax=ax
)

# Add color bar for fold change
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.005, pad=0)  # Adjusted color bar size
cbar.set_label('Fold Change')

# Annotate significance stars above each marker
for index, row in expression_data.iterrows():
    ax.text(index, 0.03,  # Adjusted y-value to place stars directly above squares
            row['Stars'], 
            ha='center', va='bottom', fontsize=12, color='black')

# Plot labels and title for figure
ax.set_xlabel('Tissue Type')
ax.set_ylabel('-log$_{10}$(p)')
ax.set_title('Tumor vs. Normal (Fold Change)', fontsize=12, fontweight='bold')

# Set x-tick labels with counts
x_labels = expression_data['Tissue Type']
ax.set_xticks(range(len(x_labels)))
ax.set_xticklabels(x_labels, rotation=45, rotation_mode='anchor', ha='right', fontsize=8)

# Remove y-axis ticks and spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.tight_layout()

# Save and show plot
plt.savefig("deseq2_Fold_Change_with_Color_Bar_and_Significance.png", dpi=400)
plt.show()
