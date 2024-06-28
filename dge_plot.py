"""
contact:    wlee9829@gmail.com
date:       2024_06_27
python:     python3.10
script:     DGE_plot.py

This Python script plots the distribution of log2(x+1) RSEM normalized gene 
expression counts downloaded from the UCSC Xena web platform.
"""
# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
1.1.1
    -   Added boxplot function, including the counting for datasets. 
    -   Function to skip processing input files and use CSV instead

1.1.0
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
import seaborn as sns
import textwrap

###########################################################################
# 1. Set up functions

def read_input(input_file, phenotype_file, gene_name):
    """
    We highly recommended pre-processing the input files;
        Row 1: sample_ids
        Row 2: normalized read counts

    Check the row count of input_file. If =2, use file as is. If >2 rows,
    parse file, search subsequent rows for matching "gene_name". Coverts
    the input_count dataframe sample_ids to names in the phenotype
    description file.

    Returns combined dataframe with normalized_count and renamed columns.
    """
    # Read input expression file
    if "gzip" in str(input_file):
        with gzip.open(input_file, 'rt') as file:
            count_df = pd.read_csv(file, sep='\t')
    else:
        with open(input_file, 'rt') as file:
            count_df = pd.read_csv(file, sep='\t')
    
    # Check the DataFrame shape and search for appropriate gene 
    if count_df.shape[0] == 1 and count_df.iloc[0, 0] == gene_name:
        print(f"Processing data for {gene_name}.")
        norm_count = count_df

    else:
        print(f"Searching data for {gene_name}.")
        norm_count = count_df[count_df.iloc[:, 0] == gene_name]

    norm_count = norm_count.transpose().reset_index()
    norm_count.columns = ['identifier', 'Expression']
    expression = norm_count[1:]

    # Read phenotype description file
    with gzip.open(phenotype_file, 'rt', errors='replace') as file:
        phenotype_df = pd.read_csv(file, sep='\t', encoding='utf-8-sig')

    """
    Phenotype file is separated into these columns:
        1) sample
        2) detailed_category
        3) primary disease or tissue
        4) _primary_site
        5) _sample_type
        6) _gender
        7) _study
    """
    # Process phenotype information
    phenotype_df['_primary_site'] = phenotype_df['_primary_site'].str.title()

    # Merge and update identifier based on phenotype information
    merged_count = pd.merge(expression, phenotype_df, 
        left_on='identifier', right_on='sample', how='left')

    # Need to fix the name filtering. Some issues here.
    merged_count = merged_count.dropna(subset=['_primary_site'])

    merged_count.loc[merged_count['_primary_site'].str.contains(
        'sympathetic', case=False), 'identifier'] = 'Sympathetic Nervous System'
    
    merged_count.loc[merged_count['_primary_site'].str.contains(
        '"Soft Tissue,Bone Tumor"'), 'identifier'] = 'Soft Tissue or Bone Tumor'

    print(merged_count)

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
    
    transformed_df.to_csv(f'{gene_name}.csv', index=False)

    return transformed_df

def plot(dataframe, gene_name, tissue_types):
    """
    Plot a boxplot of the log2 transformed normalized_counts, 
    separated by phenotype (tissue types) and outputs PNG file.
    """
    header_row = list(dataframe.columns)
    print(header_row)

    # Calculate counts of data points for each tissue type
    counts = filtered['Tissue Type'].value_counts()
    counts_dict = counts.to_dict()
    total_count = counts.sum()
    print(f"Total data points: {total_count}")

    # Set figure parameters
    plt.figure(figsize=(20, 8))
    colors = {tissue:'blue' if tissue in tissue_types else 'red' 
        for tissue in filtered['Tissue Type'].unique()}

    # Boxplot variables
    sns.boxplot(
                x='Tissue Type', 
                y='Expression', 
                data=filtered, 
                hue='Tissue Type', 
                palette=colors, 
                whis=[0, 100],
                linewidth=2, 
                fliersize=0.5, 
                showcaps=True, 
                boxprops={'facecolor':'none', 'edgecolor':'black', 
                    'linewidth': 0.75}, 
                whiskerprops={'color':'black', 'linewidth': 0.75}, 
                medianprops={'color':'black', 'linewidth': 0.75}, 
                capprops={'color':'gray', 'linewidth': 0.5}
            )

    # Stripplot variables
    sns.stripplot(  
                x='Tissue Type', 
                y='Expression',
                data=filtered, 
                hue='Tissue Type', 
                palette=colors, 
                jitter=True, 
                edgecolor='black', 
                size=4, 
                alpha=0.25
            )

    # Plot labels and title
    plt.xlabel('Tissue Type', fontsize=6, fontweight='bold')
    plt.ylabel(f'{gene_name} Expression (Log2(x+1) RSEM Normalized Count)',
        fontsize=8, fontweight='bold')
    plt.title(f'{gene_name} Expression by Tissue Type', fontsize=12,
        fontweight='bold')

    # Add counts to x-axis labels
    x_labels = [f"{label}\n(n={counts_dict.get(label, 0)})" 
        for label in filtered['Tissue Type'].unique()]
    plt.xticks(range(len(counts)), x_labels, rotation=90, 
        ha='center', fontsize=8)
    plt.subplots_adjust(bottom=0.2)

    # Remove spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()

    # Save the plot as PNG
    output_file = f"{gene_name}_boxplot.png"
    plt.savefig(output_file)
    plt.close()
    print(f"Plot saved as {output_file}")

def main():
    """
    Main function to set up argument parser and call relevant functions.
    """
    parser = argparse.ArgumentParser(description='Plot gene expression \
        distribution from input files.')
    parser.add_argument('input_file', type=str, help='Gene expression count \
        data file for processing')
    parser.add_argument('gene_name', type=str, help='Gene name to extract \
        expression data for')
    parser.add_argument('working_directory', type=str, help='PATH to working \
        directory; use "." for current directory')
    parser.add_argument('--csv-file', type=str, help='Optional flag, \
        specifies a CSV file to read data from, reduces memory use', 
        default=None)

    args = parser.parse_args()
    input_file=args.input_file
    gene_name = args.gene_name
    working_directory=args.working_directory
    csv_file = args.csv_file

    # Add or remove tissue types as appropriate.
    tissue_types = [
        "Adipose Tissue",
        "Adrenal Gland",
        "Bile duct",
        "Bladder",
        "Blood",
        "Blood Vessel",
        "Bone Marrow",
        "Brain",
        "Breast",
        "Cervix",
        "Cervix Uteri",
        "Colon",
        "Endometrium",
        "Esophagus",
        "Eye",
        "Fallopian Tube",
        "Head and Neck region",
        "Heart",
        "Kidney",
        "Lining of body cavities",
        "Liver",
        "Lung",
        "Lymphatic tissue",
        "Muscle",
        "Nerve",
        "Ovary",
        "Pancreas",
        "Paraganglia",
        "Pituitary",
        "Prostate",
        "Rectum",
        "Salivary Gland",
        "Skin",
        "Small Intestine",
        "Soft tissue,Bone",
        "Spleen",
        "Stomach",
        "Sympathetic Nervous System",
        "Testis",
        "Thymus",
        "Thyroid",
        "Thyroid Gland",
        "Uterus",
        "Vagina",
        "White blood cell"
    ]

    if csv_file:
        combined_df = pd.read_csv(f'{gene_name}.csv')

    else:
        phenotype_file='''TcgaTargetGTEX_phenotype.txt.gz'''
        count_df = read_input(
            input_file, phenotype_file, gene_name) 
    
    csv_df = pd.read_csv(f'{gene_name}.csv', encoding='utf-8-sig')
    plot(csv_df, gene_name, tissue_types)

if __name__ == "__main__":
    main()