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

def read_input(input_file, phenotype_file, gene_name, tissue_types):
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
    with gzip.open(input_file, 'rt') as file:
        count_df = pd.read_csv(file, sep='\t')
    
    # Check the DataFrame shape and search for appropriate gene 
    first_gene = count_df.iloc[0, 0]

    if count_df.shape[0] == 1 and first_gene == gene_name:
        print(f"Plotting data for {first_gene}.")
        norm_count = count_df

    else:
        print(f"Searching data for {gene_name}.")
        norm_count = count_df[count_df.iloc[:, 0] == gene_name]

    # Read phenotype description file
    with gzip.open(phenotype_file, 'rt') as file:
        phenotype_df = pd.read_csv(file, sep='\t')

    # Function to find matching phenotype category for a given row
    def find_matching_category(row):
        for tissue in tissue_types:
            if tissue in row[2]:
                return tissue
            elif tissue in row[1]:
                return tissue
        return row[0]

    phenotype_dict = {row[0]: find_matching_category(row) for index, 
        row in phenotype_df.iterrows()}
    
    norm_count.columns = [phenotype_dict.get(col, col) for col 
        in norm_count.columns]

    norm_df = norm_count.stack().groupby(level=1, 
        group_keys=True).apply(list).apply(pd.Series).T
    norm_df.drop(columns=['sample'], inplace=True)
    norm_df.to_csv(f'{gene_name}.csv', index=False)

    # print(norm_df)
    return norm_df

def plot(dataframe, gene_name, tissue_types):
    """
    Plot a boxplot of the log2 transformed normalized_counts, 
    separated by phenotype (tissue types) and outputs PNG file.
    """
    # Define the order of data based on tissue/cancer type.
    order = [
        'Adipose Tissue', 
        'Adrenal Gland', 'ACC', # Adrenocortical carcinoma
        'Artery',
        'Bladder', 'BLCA', # Bladder urothelial carcinoma
        'Blood', 'DLBC', # Diffuse large B-cell lymphoma
        'Blood Vessel',
        'Bone Marrow', 'LAML', # Acute myeloid leukemia
        'Brain', 'GBM', # Glioblastoma multiforme
        'Breast', 'BRCA', # Breast invasive carcinoma
        'Cells', 'PCPG', # Pheochromocytoma and paraganglioma
        'Cervix', 'CESC', # Cervical & endocervical cancer
        'Colon', 'COAD', # Colon adenocarcinoma
        'Esophagus', 'ESCA', # Esophageal carcinoma
        'Fallopian Tube',
        'Heart',
        'Kidney', 'KIRC', # Kidney clear cell carcinoma
        'Liver', 'LIHC', # Liver hepatocellular carcinoma
        'Lung', 'LUAD', # Lung adenocarcinoma
        'Muscle',
        'Nerve',
        'Ovary', 'OV', # Ovarian serous cystadenocarcinoma
        'Pancreas', 'PAAD', # Pancreatic adenocarcinoma
        'Pituitary',
        'Prostate', 'PRAD', # Prostate adenocarcinoma
        'Salivary Gland', 
        'SARC', # Sarcoma
        'Skin', 'SKCM', # Skin cutaneous melanoma
        'Small Intestine',
        'Spleen',
        'Stomach', 'STAD', # Stomach adenocarcinoma
        'Testis', 'TGCT', # Testicular germ cell tumor
        'Thyroid', 'THCA', # Thyroid carcinoma
        'Uterus', 'UCEC', # Uterine corpus endometrioid carcinoma
        'Vagina'
        ]

    data = dataframe.reindex(order, axis=1)
    filtered = data.dropna(subset=['Expression'])
    print(filtered)

    filtered['Tissue Type'] = pd.Categorical(filtered.index,
        categories=order, ordered=True)

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
    parser.add_argument('gene_name', type=str, help='Gene name to extract \
        expression data for')
    parser.add_argument('working_directory', type=str, help='PATH to working \
        directory that contains all relevant files')
    parser.add_argument('--csv-file', type=str, help='Optional flag, \
        specifies a CSV file to read data from', default=None)

    args = parser.parse_args()
    gene_name = args.gene_name
    working_directory=args.working_directory
    csv_file = args.csv_file

    # Add or remove tissue types as appropriate.
    tissue_types = [
    "Adipose Tissue",
    "Adrenal gland",
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
        input_file='''TcgaTargetGtex_gene_expected_count.gz'''
        phenotype_file='''TcgaTargetGTEX_phenotype.txt.gz'''
        count_df = read_input(input_file, phenotype_file, 
            gene_name, tissue_types) 
        csv_df = pd.read_csv(f'{gene_name}.csv')

    plot(csv_df, gene_name, tissue_types)

if __name__ == "__main__":
    main()