"""
contact:    wlee9829@gmail.com
date:       2024_08_06
python:     python3.10
script:     dge_2genes.py

This Python script compares the distribution of log2(x+1) RSEM normalized data
between two genes.
"""

# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
1.0.0   
    -   Initial commit, set up outline of logic and functions.
"""

# Import Packages
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import mannwhitneyu
from scipy import stats
import seaborn as sns
import sys
import textwrap

###########################################################################
# 1. Set up functions

def read_input(file_path1, file_path2):
    """
    Read the data for RNA expression counts for two specified genes from
    separate CSV files.

    Args:
        file_path1 (str): CSV file containing expression data for first gene.
        file_path2 (str): CSV file containing expression data for second gene.

    Returns:
        tuple: Two DataFrames, one for each gene's tumor data.
    """
    placeholder = -9999

    data1 = pd.read_csv(file_path1).fillna(placeholder)
    data2 = pd.read_csv(file_path2).fillna(placeholder)

    def filter_tumor_columns(data):
        data.columns = [col.strip() for col in data.columns]
        tumor_columns = [col for col in data.columns if col.endswith('Tumor')]
        filtered_data = data[tumor_columns] if tumor_columns else pd.DataFrame()

        return filtered_data

    tumor_data1 = filter_tumor_columns(data1)
    tumor_data2 = filter_tumor_columns(data2)

    return tumor_data1, tumor_data2

def calc_significance(df, gene1_name, gene2_name, placeholder_value):
    """
    Calculate significance between two genes for each tissue type.
    This function computes p-values using a Mann-Whitney U test, ignoring placeholder values.

    Args:
        df (pd.DataFrame): A DataFrame containing gene expression data for two genes.
        gene1_name (str): The name of the first gene to be compared.
        gene2_name (str): The name of the second gene to be compared.
        placeholder_value (float): The placeholder value used to fill NaNs.

    Returns:
        p_values (dict): Dictionary with tissue types as keys and p-values as values.
        significance_levels (dict): Dictionary with tissue types as keys and significance level stars as values.
    """
    placeholder_value = -9999

    p_values = {}
    significance_levels = {}
    
    tissue_types = df['tissue_type'].unique()
    for tissue in tissue_types:
        gene1_expr = df[(df['tissue_type'] == tissue) & (df['gene'] == gene1_name)]['expression']
        gene2_expr = df[(df['tissue_type'] == tissue) & (df['gene'] == gene2_name)]['expression']
        
        # Filter out placeholder values
        gene1_expr = gene1_expr[gene1_expr != placeholder_value]
        gene2_expr = gene2_expr[gene2_expr != placeholder_value]
        
        if len(gene1_expr) > 1 and len(gene2_expr) > 1:
            _, p_value = stats.mannwhitneyu(gene1_expr, gene2_expr, alternative='two-sided')
            p_values[tissue] = p_value
            if p_value < 0.001:
                significance_levels[tissue] = '***'
            elif p_value < 0.01:
                significance_levels[tissue] = '**'
            elif p_value < 0.05:
                significance_levels[tissue] = '*'
            else:
                significance_levels[tissue] = ''
        else:
            significance_levels[tissue] = ''
    
    return p_values, significance_levels

def plot(dataframe1, dataframe2, gene1_name, gene2_name, output_prefix='', placeholder_value=-9999):
    """
    Generates and saves plots of the log2 transformed normalized_counts,
    separated by phenotype (tissue types). This function creates a boxplot
    and a strip plot to visualize the expression levels of two genes across
    different tissue types, ignoring placeholder values.

    Args:
        dataframe1 (pd.DataFrame): A DataFrame containing gene expression data for gene 1.
        dataframe2 (pd.DataFrame): A DataFrame containing gene expression data for gene 2.
        gene1_name (str): The name of the first gene to be plotted.
        gene2_name (str): The name of the second gene to be plotted.
        output_prefix (str, optional): A prefix to be added to the output file name. Defaults to ''.
        placeholder_value (float, optional): The placeholder value used to fill NaNs. Defaults to -9999.

    Returns:
        None. The plot is saved as PNG and SVG files.
    """
    # Rename columns for clarity
    dataframe1.columns = [f'{gene1_name}_{col}' for col in dataframe1.columns]
    dataframe2.columns = [f'{gene2_name}_{col}' for col in dataframe2.columns]
    
    # Combine the two dataframes
    combined_df = pd.concat([dataframe1, dataframe2], axis=1)
    
    # Reshape the data for plotting
    df = pd.melt(combined_df, var_name='tissue_type', value_name='expression')
    
    # Extract gene names from column headers
    df['gene'] = df['tissue_type'].apply(lambda x: x.split('_')[0])
    df['tissue_type'] = df['tissue_type'].apply(lambda x: x.split('_', 1)[1])

    # Filter out placeholder values
    df = df[df['expression'] != placeholder_value]

    # Create color mapping for gene1 and gene2
    unique_genes = df['gene'].unique()
    gene_palette = {gene: 'blue' if gene == gene1_name else 'red' for gene in unique_genes}
    
    # Set figure parameters
    num_tissues = len(df['tissue_type'].unique())
    fig_width = min(num_tissues, 16)
    plt.figure(figsize=(fig_width, 8))
    
    # Define boxplot variables
    ax = sns.boxplot(
            x='tissue_type',
            y='expression',
            data=df,
            hue='gene',
            palette=gene_palette,
            whis=[0, 100],
            linewidth=2,
            fliersize=0.5,
            showcaps=True,
            boxprops={'facecolor': 'none', 'edgecolor': 'black', 'linewidth': 0.75},
            whiskerprops={'color': 'black', 'linewidth': 0.75},
            medianprops={'color': 'black', 'linewidth': 0.75},
            capprops={'color': 'gray', 'linewidth': 0.5}
        )

    # Define stripplot variables
    sns.stripplot(
            x='tissue_type',
            y='expression',
            data=df,
            hue='gene',
            palette=gene_palette,
            jitter=True,
            edgecolor='black',
            size=5,
            alpha=0.25,
            ax=ax,
            dodge=True  # Align strip plots with the boxplots
        )

    # Plot labels and title for figure
    ax.set_xlabel('', fontsize=8, fontweight='bold')
    ax.set_ylabel(f'{gene1_name} and {gene2_name} Expression (Log2(x+1) RSEM Normalized Count)',
        fontsize=8, fontweight='bold')
    ax.set_title(f'{gene1_name} and {gene2_name} Expression by Tissue Type', fontsize=12,
        fontweight='bold')

    # Get unique tissue types for x-axis labels
    tissue_types = df['tissue_type'].unique()

    # Prepare labels with counts for each gene and tissue type
    x_labels = []
    for tissue_type in tissue_types:
        counts_gene1 = df[(df['tissue_type'] == tissue_type) & (df['gene'] == gene1_name)].shape[0]
        counts_gene2 = df[(df['tissue_type'] == tissue_type) & (df['gene'] == gene2_name)].shape[0]
        x_labels.append(f"{tissue_type} ({gene1_name}: n={counts_gene1}, {gene2_name}: n={counts_gene2})")
    
    ax.set_xticklabels(x_labels, rotation=45, rotation_mode='anchor', 
                        ha='right', fontsize=8)

    ax.set_xlim(-0.75, len(tissue_types) - 0.75)
    
    # Remove spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Print statistics on tissue types
    print(f"Summary of datasets:\n")
    unique_tissue_types = df['tissue_type'].unique()
    print(f"\t\t     Unique sample types : {len(unique_tissue_types)}")

    total_count = df.shape[0]
    print(f"\t\tTotal number of datasets : {total_count}\n")
    for tissue_type in unique_tissue_types:
        count = df[df['tissue_type'] == tissue_type].shape[0]
        print(f"{tissue_type.rjust(40)} : {count}")
    print("-" * 60)

    # Add significance annotations above the boxplot groups
    p_values, significance_levels = calc_significance(df, gene1_name, gene2_name, placeholder_value=-9999)
    significance_annotations = {
        '***': {'color': 'black', 'alpha': 1.0},
        '**': {'color': 'black', 'alpha': 1.0},
        '*': {'color': 'black', 'alpha': 1.0},
        '': {'color': 'black', 'alpha': 0.0}
    }

    # Calculate mid-points for annotation
    for tissue_type, significance in significance_levels.items():
        if significance:
            tissue_index = list(tissue_types).index(tissue_type)
            y_coord = df[df['tissue_type'] == tissue_type]['expression'].max()
            
            ax.annotate(
                significance,
                (tissue_index, y_coord),
                xytext=(0, 10),
                textcoords='offset points',
                ha='center',
                va='center',
                color=significance_annotations[significance]['color'],
                fontsize=8,
                fontweight='bold',
                backgroundcolor='none',
                alpha=significance_annotations[significance]['alpha']
            )

            # Add gray background behind the annotation
            ax.axvspan(
                tissue_index - 0.5,
                tissue_index + 0.5,
                alpha=0.1,
                edgecolor='black',
                linewidth=0.5,
                facecolor='gray',
                zorder=-1)

    # Add custom legend
    legend_labels = [f'{gene1_name}', f'{gene2_name}']
    legend_handles = [plt.Line2D([0], [0], 
        marker='o', color='w', markerfacecolor=color, markersize=5) 
        for color in ['blue', 'red']]
    legend = ax.legend(
        legend_handles, 
        legend_labels, 
        loc='upper right',
        fontsize=8,
        frameon=False,
        title='Genes',
        title_fontsize='9'
    )
    
    # Add gray background for each tissue type
    for i in range(len(tissue_types)):
        ax.axvspan(
            i - 0.5,
            i + 0.5,
            alpha=0.1,
            color='lightgray',
            zorder=-1
        )
        
    # Save the plot as PNG and SVG
    output_file = f"{output_prefix}{gene1_name}_{gene2_name}_plot"
    plt.subplots_adjust(bottom=0.1)
    plt.tight_layout()
    plt.savefig(f"{output_file}.png", dpi=400)
    plt.savefig(f"{output_file}.svg", dpi=400)
    plt.close()

    time = str(datetime.now())[:-7]
    print(f"Plot saved as {output_file} on {time}.")

def parse_args():
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog="%(prog)s",
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent("""\
###########################################################################
NOTE: The input strings are positional.

1. file_path1:      Input file (.csv) structured with tissue types as columns
                    and expression count values as new rows for every gene.

2. file_path2:      See above.

3. -o, --output-prefix  An output name preceding 'gene_name'_plot.png.

###########################################################################
"""),
    usage=
"""
    \npython3 %(prog)s file_path1 file_path2
""")
    parser.add_argument('file_path1', type=str, 
        help='Path to the CSV file for the first gene.')
    parser.add_argument('file_path2', type=str, 
        help='Path to the CSV file for the second gene.')
    parser.add_argument(
        '-o', '--output', type=str, help='User given name for output plot.')

    return parser.parse_args()

def main(args):
    """
    Main function to execute the script.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    gene1_file = args.file_path1; gene1 = gene1_file.replace(".csv", "")
    gene2_file = args.file_path2; gene2 = gene2_file.replace(".csv", "")

    tumor_data1, tumor_data2 = read_input(gene1_file, gene2_file)
    
    plot(tumor_data1, tumor_data2, gene1, gene2, 
        output_prefix=args.output if args.output else '')

if __name__ == "__main__":
    args = parse_args()
    main(args)
    sys.exit()