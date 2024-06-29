"""
contact:    wlee9829@gmail.com
date:       2024_06_27
python:     python3.10
script:     DGE_plot.py

This Python script plots the distribution of log2(x+1) RSEM normalized gene 
expression counts downloaded from the UCSC Xena web platform.
"""
# Define version
__version__ = "1.1.2"

# Version notes
__update_notes__ = """
1.1.2
    -   Add the Wilcoxon rank-sum (Mann-Whitney U) test on columns with paired 
        "Normal" and "Tumor" (calc_significance).
    -   Added 'Other' tissue type exclusion flag (-x, --exclude)
    -   Added option to name output figures (-o, --output-prefix)
    -   Added statistics printout (-s, --stats)

1.1.1
    -   Added boxplot and stripplot function, prints tissue_type count 
        for each of the categories. 

1.1.0
    -   Add function to skip processing input files and use CSV instead
        (-csv, --csv-file)

    -   Added CSV output function to increse processing speed.

1.0.1
    -   Combined input file reading functions to convert sample_ids to tissue
        type by phenotype file.

1.0.0   
    -   Initial commit, set up outline of logic and functions.
"""

# Import Packages
import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import ranksums
import seaborn as sns

###########################################################################
# 1. Set up functions

def read_input(input_file, phenotype_file, gene_name):
    """
    Processes the input gene expression and phenotype data files to generate 
    a combined DataFrame with normalized counts and renamed columns.

    This function reads the input expression file and phenotype description 
    file, processes the data to extract and rename columns based on tissue 
    types and sample categories, and returns a structured DataFrame. It 
    supports both plain text and gzip-compressed input files.

    Args:
        input_file (str): Path to the gene expression count data file.
        phenotype_file (str): Path to the phenotype description file.
        gene_name (str): The name of the gene to extract expression data for.

    Returns:
        pd.DataFrame: A transformed DataFrame with normalized counts and 
        renamed columns based on tissue types and sample categories.
    
    The phenotype file is expected to have the following columns:
        - sample
        - detailed_category
        - primary disease or tissue
        - _primary_site
        - _sample_type
        - _gender
        - _study
    """
    # Read input expression file
    if "gzip" in str(input_file):
        with gzip.open(input_file, 'rt') as file:
            count_df = pd.read_csv(file, sep='\t')
    else:
        with open(input_file, 'rt') as file:
            count_df = pd.read_csv(file, sep='\t')
    
    # Check DataFrame shape, search for specified gene
    if count_df.shape[0] == 1 and count_df.iloc[0, 0] == gene_name:
        print(f"Processing data for {gene_name}.")
        norm_count = count_df

    else:
        print(f"Searching data for {gene_name}.")
        norm_count = count_df[count_df.iloc[:, 0] == gene_name]

    norm_count = norm_count.transpose().reset_index()
    norm_count.columns = ['identifier', 'Expression']
    expression = norm_count[1:]

    # Read phenotype description file, process information
    with gzip.open(phenotype_file, 'rt', errors='replace') as file:
        phenotype_df = pd.read_csv(file, sep='\t', encoding='utf-8-sig')

    phenotype_df['_primary_site'] = phenotype_df['_primary_site'].str.title()

    # Merge and update identifier based on phenotype information
    merged_count = pd.merge(expression, phenotype_df, 
        left_on='identifier', right_on='sample', how='left')

    def get_sample_type(sample_type):
        if pd.isna(sample_type):
            return 'Other'
        elif 'Normal' in sample_type:
            return 'Normal'
        elif 'Tumor' in sample_type:
            return 'Tumor'
        else:
            return 'Other'

    merged_count['identifier'] = merged_count.apply(
        lambda row: 
        f"{row['_primary_site']} "
        f"{get_sample_type(row['_sample_type'])}", axis=1)

    # Select and restructure columns
    merged_count = merged_count[['identifier', 'Expression']]
    transformed_df = (merged_count
        .set_index('identifier')
        .stack()
        .groupby(level=0)
        .apply(list)
        .apply(pd.Series)
        .T)
    transformed_df.columns = [
        'Sympathetic Nervous System Other' if 'sympathetic' in col.lower() \
            and 'other' in col.lower() else 
        'Sympathetic Nervous System Tumor' if 'sympathetic' in col.lower() \
            and 'tumor' in col.lower() else 
        'Soft Tissue or Bone Normal' if 'soft tissue,bone' in col.lower() \
            and 'normal' in col.lower() else
        'Soft Tissue or Bone Tumor' if 'soft tissue,bone' in col.lower() \
            and 'tumor' in col.lower() else
        col
        for col in transformed_df.columns
        ]
    
    transformed_df = transformed_df.drop(
        columns=['nan Normal', 'nan Other'],axis=1)
    transformed_df.to_csv(f'{gene_name}.csv', index=False)
    
    # Print the transformed DataFrame
    # print(transformed_df)
    return transformed_df

def calc_significance(dataframe):
    """
    Perform Wilcoxon rank-sum tests between columns with matching names 
    when stripped of "Normal" or "Tumor". Prints the p-value for each matched 
    pair and summarizes significance based on the p-values.

    Args:
        dataframe (pd.DataFrame): A DataFrame containing gene expression data.

    Returns:
        dict: A dictionary containing p-values.
    """
    p_values = {}; significance_levels = {}

    # Iterate over unique base names (stripped of "Normal" and "Tumor")
    for base_name in dataframe.columns.str.replace('Normal|Tumor', '', 
        regex=True).unique():
        normal_col = f"{base_name}Normal"
        tumor_col = f"{base_name}Tumor"
        
        # Check if both Normal and Tumor columns exist
        if normal_col in dataframe.columns and tumor_col in dataframe.columns:
            normal_values = dataframe[normal_col].dropna()
            tumor_values = dataframe[tumor_col].dropna()
            
            if len(normal_values) > 0 and len(tumor_values) > 0:
                stat, p_value = ranksums(normal_values, tumor_values)
                
                # Define significance level
                if p_value < 0.001:
                    significance = '***'
                elif p_value < 0.01:
                    significance = '**'
                elif p_value < 0.05:
                    significance = '*'
                else:
                    significance = ''

                p_values[base_name] = p_value
                significance_levels[base_name] = significance

            else:
                print(f"Skipping {normal_col} vs {tumor_col} due to "
                    f"insufficient data.")

    # Summarize significance findings if any comparisons were made
    sorted_p_values = sorted(p_values.items(), key=lambda x: x[1])

    if p_values:
        print(f"Summary of Significance Calculations:\n")
        for base_name, p_value in sorted_p_values:
            if p_value < 0.001:
                significance = '***' 
            elif p_value < 0.01:
                significance = '**'
            elif p_value < 0.05:
                significance = '*' 
            else:
                significance = ''

            print(f"{base_name.rjust(40)}: (p-value = {p_value:.4f}) "
                f"{significance} ")
        print("-" * 60)

    return p_values, significance_levels

def plot(dataframe, gene_name, output_prefix='', exclude_other=False,
        stats=False):
    """
    Generates and saves a plots of the log2 transformed normalized_counts,
    separated by phenotype (tissue types). This function creates a boxplot
    and a strip plot to visualize the expression levels of a specified gene 
    across different tissue types. It includes options to exclude columns with
    'Other' in their names. The plot is saved as a PNG file.

    Args:
        dataframe (pd.DataFrame): A DataFrame containing gene expression data.
        gene_name (str): The name of the gene to be plotted.
        output_prefix (str, optional): A prefix to be added to the output file 
            name. Defaults to ''.
        exclude_other (bool, optional): If True, exclude columns containing 
            'Other' from the plot. Defaults to False.
        stats (bool, optional): If True, includes statistics on the columns that
            share a common name by varying condition.

    Returns:
        None. The plot is saved as a PNG file.
    """
    # Set up the data.
    if exclude_other:
        dataframe = dataframe.loc[:, ~dataframe.columns.str.contains('Other')]

    # Lower this number if too strict
    row_count_filter = dataframe.count() >= 5
    dataframe = dataframe.loc[:, row_count_filter]

    df = pd.melt(dataframe, var_name='tissue_type', value_name='expression')
    filtered_df = df[df['expression'].notna()]
    counts_dict = filtered_df.groupby('tissue_type').size().to_dict()

    if stats:
        # Print statistics on tissue types, and number of counts for each.
        print(f"Summary of datasets:\n")
        tissue_types = filtered_df['tissue_type'].nunique()
        print(f"\t\t     Unique sample types : {tissue_types}")

        total_count = filtered_df.shape[0]
        print(f"\t\tTotal number of datasets : {total_count}\n")
        for tissue_type, count in counts_dict.items():
            print(f"{tissue_type.rjust(40)} : {count}")
        print("-" * 60)

    # Color mapping for tissue types.
    def set_color(column_name):
        if 'Normal' in column_name:
            return 'blue'
        elif 'Tumor' in column_name:
            return 'red'
        else:
            return 'gray'

    colors = [set_color(col) for col in dataframe.columns 
        if col != 'tissue_type']

    # Set figure parameters
    plt.figure(figsize=(16, 10))

    # Define boxplot variables 
    ax = sns.boxplot(
                x='tissue_type', 
                y='expression', 
                data=filtered_df, 
                hue='tissue_type', 
                palette=colors, 
                whis=[0, 100],
                linewidth=2, 
                fliersize=0.5, 
                showcaps=True, 
                boxprops={'facecolor':'none', 
                    'edgecolor':'black', 
                    'linewidth': 0.75}, 
                whiskerprops={'color':'black', 'linewidth': 0.75}, 
                medianprops={'color':'black', 'linewidth': 0.75}, 
                capprops={'color':'gray', 'linewidth': 0.5}
            )

    # Define stripplot variables
    sns.stripplot(  
                x='tissue_type', 
                y='expression',
                data=filtered_df, 
                hue='tissue_type', 
                palette=colors, 
                jitter=True, 
                edgecolor='black', 
                size=4, 
                alpha=0.5,
            )

    # Plot labels and title for figure
    plt.xlabel('', fontsize=8, fontweight='bold')
    plt.ylabel(f'{gene_name} Expression (Log2(x+1) RSEM Normalized Count)',
        fontsize=8, fontweight='bold')
    plt.title(f'{gene_name} Expression by Tissue Type', fontsize=12,
        fontweight='bold')

    # Add counts to x-axis labels
    x_labels = [f"{label} (n={counts_dict.get(label, 0)})" 
        for label in filtered_df['tissue_type'].unique()]

    plt.xticks(ticks=range(len(x_labels)), labels=x_labels, 
        rotation=90, ha='center', fontsize=8)

    # Remove spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # Add custom legend
    legend_labels = ['Normal', 'Tumor', 'Other']
    legend_handles = [plt.Line2D([0], [0], 
        marker='o', color='w', markerfacecolor=color, markersize=6) 
        for color in ['blue', 'red', 'gray']]
    legend = plt.legend(
        legend_handles, 
        legend_labels, 
        loc='upper left',
        frameon=True,
        handletextpad=0.4,
        labelspacing=0.3,
        borderpad=0.5
        )
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_facecolor('white')

    if stats:
        # Add significance annotations above the boxplot groups
        p_values, significance_levels = calc_significance(dataframe)
        significance_annotations = {
            '***': {'color': 'gray', 'alpha': 0.5},
            '**': {'color': 'gray', 'alpha': 0.5},
            '*': {'color': 'gray', 'alpha': 0.5}
        }
        x_labels = filtered_df['tissue_type'].unique()
        for base_name, significance in significance_levels.items():
            if base_name in x_labels:
                index = np.where(filtered_df['tissue_type'].unique() \
                    == base_name)[0][0]
                ax.annotate(
                    significance, 
                    (index, filtered_df.loc[filtered_df['tissue_type'] \
                        == base_name, 'expression'].max()), 
                    xytext=(0, 10), 
                    textcoords='offset points', 
                    ha='center', 
                    va='center',
                    color=significance_annotations[significance]['color'], 
                    fontsize=8, 
                    fontweight='bold', 
                    backgroundcolor='white',
                    alpha=significance_annotations[significance]['alpha'])

    # Save the plot as PNG
    output_file = f"{output_prefix}{gene_name}_plot.png"
    plt.subplots_adjust(bottom=0.1)
    plt.tight_layout()
    plt.savefig(output_file, dpi=400)
    plt.close()
    print(f"Plot saved as {output_file}")

def main():
    """
    Main function to set up the argument parser, handle input arguments, 
    and call the relevant functions for processing and plotting gene 
    expression data.

    This function supports the following functionalities:
    - Reading gene expression data from a CSV file or input file.
    - Extracting gene expression data for a specific gene.
    - Excluding columns containing "Other" from plotting.
    - Setting a prefix for the output file names.
    - Printing statistics on tissue types and counts.

    The function parses command-line arguments, reads the input data, 
    calculates statistical significance, and generates plots for gene 
    expression data.

    Args:
        -csv, --csv-file (str, optional): Path to a CSV file to read data from.
        input_file (str, optional): Path to the gene expression count data file.
        gene_name (str): The name of the gene to extract expression data for.
        -o, --output-prefix (str, optional): Prefix for the output file names.
        -x, --exclude (bool, optional): Flag to exclude columns containing 
        "Other" from plotting.
        -s, --stats (bool, optional): Flag to print statistics on tissue types 
        and counts.

    Raises:
        argparse.ArgumentError: If the required arguments are not provided 
        when --csv-file is not used.
    """
    parser = argparse.ArgumentParser(
        description='Process, plot gene expression from input files.')

    parser.add_argument(
        '-csv', '--csv-file', type=str,
        help='Optional. Path to CSV file. Reduces memory use if specified.',
        default=None)

    parser.add_argument(
        'input_file', type=str, nargs='?',
        help=('Required if CSV file (-csv, --csv-file) is not provided. Path to'
            ' gene expression count data file for processing.'))

    parser.add_argument(
        'gene_name', type=str,
        help='The name of the gene to extract expression data for.')

    parser.add_argument(
        '-o', '--output-prefix', type=str, default='',
        help='Optional. Prefix for the output file names. Defaults to None.')

    parser.add_argument(
        '-x', '--exclude', action='store_true',
        help='Optional. Exclude columns containing "Other" from plotting.')

    parser.add_argument(
        '-s', '--stats', action='store_true',
        help='Optional. Prints statistics on the tissue types and counts.')

    args = parser.parse_args()
    csv_file = args.csv_file
    gene_name = args.gene_name
    output_prefix = args.output_prefix
    exclude_other = args.exclude
    stats = args.stats

    if csv_file:
        combined_df = pd.read_csv(csv_file)
        print(f"Searching {csv_file} for {gene_name} data.")
        print("-" * 60)

    else:
        if not args.input_file or not args.gene_name:
            parser.error('You must provide input_file and gene_name' \
                ' if csv-file is not provided.')
        input_file = args.input_file

        if os.path.exists('TcgaTargetGTEX_phenotype.txt.gz'):
            phenotype_file = 'TcgaTargetGTEX_phenotype.txt.gz'
        else:
            phenotype_file = '*_phenotype.txt.gz'
            print(f"Attempting to use file {phenotype_file} as identifier.")
            
        count_df = read_input(input_file, phenotype_file, gene_name)
        combined_df = count_df
        print("-" * 60)
    
    if output_prefix:
        output_prefix = str(output_prefix+"_")

    plot(combined_df, gene_name, output_prefix, exclude_other, stats)

if __name__ == "__main__":
    main()