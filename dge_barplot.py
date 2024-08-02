import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.stats import ttest_ind

# Load data
file_path = 'RSEM_GLP1R.csv'
data = pd.read_csv(file_path)

# Define tissue types
tissue_types = [
    'Adrenal Gland', 'Bile Duct', 'Bladder', 'Brain', 'Breast',
    'Cervix', 'Colon', 'Esophagus', 'Head And Neck', 'Kidney',
    'Liver', 'Lung', 'Ovary', 'Pancreas', 'Prostate', 'Rectum',
    'Skin', 'Stomach', 'Testis', 'Thyroid', 'Uterus'
]

# Prepare data for plotting
relative_expression = []

for tissue_type in tissue_types:
    normal_col = f'{tissue_type} Normal'
    tumor_col = f'{tissue_type} Tumor'
    
    if normal_col not in data.columns or tumor_col not in data.columns:
        print(f"Skipping {tissue_type} due to missing columns.")
        continue
    
    normal_values = data[normal_col]
    tumor_values = data[tumor_col]
    
    # Drop NaN values from both normal and tumor values
    combined_data = pd.DataFrame({'Normal': normal_values, 'Tumor': tumor_values})
    combined_data = combined_data.dropna()
    
    normal_cleaned = combined_data['Normal']
    tumor_cleaned = combined_data['Tumor']

    # Calculate fold change for all data
    mean_normal = np.mean(normal_cleaned + 1)
    mean_tumor = np.mean(tumor_cleaned + 1)
    fold_change = np.log2(mean_tumor / mean_normal)

    # Determine fold change calculation based on file type
    # if "RSEM" in file_path:
    #     # Calculate fold change with log2 transformation for RSEM data
    #     mean_normal = np.mean(normal_cleaned + 1)
    #     mean_tumor = np.mean(tumor_cleaned + 1)
    #     fold_change = np.log2(mean_tumor / mean_normal)
    
    # elif "DESEQ2" in file_path:
    #     # Calculate fold change without log2 transformation for DESEQ2 data
    #     mean_normal = np.mean(normal_cleaned)
    #     mean_tumor = np.mean(tumor_cleaned)
    #     fold_change = mean_tumor / mean_normal

    # Perform a t-test to get the p-value
    _, p_value = ttest_ind(normal_cleaned, tumor_cleaned, equal_var=False)
    
    # Calculate -log10(p-value) for significance
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

    # Append fold change, significance, and stars
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

# Set up figure and axes
fig, ax = plt.subplots(figsize=(14, 8))

# Plot small boxes with color representing fold change
sns.scatterplot(
    x='Tissue Type', 
    y=[1] * len(expression_data),  # Constant y-value for small boxes
    hue='Fold Change', 
    palette=cmap, 
    data=expression_data,
    s=400,  
    marker='s',  
    legend=None,
    ax=ax
)

# Add color bar for fold change
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', fraction=0.005, pad=0.04)
cbar.set_label('Fold Change')

# Annotate significance stars above each marker
for index, row in expression_data.iterrows():
    ax.text(index, 1.1, 
            row['Stars'], 
            ha='center', va='bottom', fontsize=10, color='black')

# Plot labels and title for figure
ax.set_xlabel('Tissue Type')
ax.set_ylabel('')
ax.set_title('Relative Gene Expression: Tumor vs Normal', fontsize=12, fontweight='bold')

# Set x-tick labels with counts
x_labels = expression_data['Tissue Type']
ax.set_xticks(range(len(x_labels)))
ax.set_xticklabels(x_labels, rotation=45, rotation_mode='anchor', ha='right', fontsize=8)

# Remove y-axis ticks and spines
ax.yaxis.set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.tight_layout()

# Save and show plot
plt.savefig("rsem_Fold_Change_with_Color_Bar_and_Significance.png")
plt.show()