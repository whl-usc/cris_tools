"""
contact:    wlee9829@gmail.com
date:       2024_06_25
python:     python3.10
script:     DGE_plot.py

This Python script plots the distribution of log2(x+1) normalized gene 
expression counts downloaded from the UCSC Xena web platform.
"""
# Define version
__version__ = "1.2.0"

# Version notes
__update_notes__ = """
    
1.2.0
    -   Changed the plot parameters for readability and consistency.
    -   Added optional argument to directly plot the generated expression 
        csv without having to redo data transformation.

1.1.0
    -   Wrote boxplot function, including count for datasets.
    -   Function to write-out the expression to csv for faster processing

1.0.1
    -   Combined input file reading function to also convert sample_ids
        to the tissue type by phenotype file (GTEX) or the cancer type by
        name (TCGA).

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

def read_gtex(input_file, phenotype_file, gene_name):
    """
    Highly recommended to pre-process the input;
        Row 1: sample_ids
        Row 2: normalized read counts

    Check the row count of input_file. If =2, use file as is. If >2 rows,
    parse file, search subsequent rows for matching "gene_name". Coverts
    gtex_count dataframe sample_ids to names in the phenotype
    description file.

    Returns combined dataframe with gtex_count and renamed columns.
    """
    # Read input expression file
    with gzip.open(input_file, 'rt') as file:
        exp_df = pd.read_csv(file, sep='\t')
    
    # Check the DataFrame shape and search for appropriate gene 
    first_gene = exp_df.iloc[0, 0]

    if exp_df.shape[0] == 1 and first_gene == gene_name:
        print(f"Plotting data for {first_gene}.")
        gtex_count = exp_df

    else:
        print(f"Searching data for {gene_name}.")
        gtex_count = exp_df[exp_df.iloc[:, 0] 
            == gene_name]

    # Read phenotype description file
    with gzip.open(phenotype_file, 'rt') as file:
        phenotype_df = pd.read_csv(file, sep='\t')

    tissue_types = [
        'Adipose Tissue', 'Adrenal Gland', 'Artery', 'Bladder', 'Blood',
        'Blood Vessel', 'Bone Marrow', 'Brain', 'Breast', 'Cells', 'Cervix',
        'Colon', 'Esophagus', 'Fallopian Tube', 'Heart', 'Kidney', 'Liver',
        'Lung', 'Muscle', 'Nerve', 'Ovary', 'Pancreas', 'Pituitary', 
        'Prostate', 'Salivary Gland', 'Skin', 'Small Intestine', 'Spleen', 
        'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina'
    ]

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
    
    gtex_count.columns = [phenotype_dict.get(col, col) for col 
        in gtex_count.columns]

    gtex_df = gtex_count.stack().groupby(level=1, 
        group_keys=True).apply(list).apply(pd.Series).T
    gtex_df.drop(columns=['sample'], inplace=True)
    gtex_df.to_csv(f'{gene_name}.csv', index=False)
    # print(gtex_df)

    return tissue_types

def read_tcga(input_file, gene_name):
    """
    Pre-process the input as with before.

    Returns tcga_count dataframe with renamed columns.
    """
    # Read input expression file
    with gzip.open(input_file, 'rt') as file:
        exp_df = pd.read_csv(file, sep='\t')
        
    # Identify the first gene in the DataFrame
    first_gene = exp_df.iloc[0, 0]

    if exp_df.shape[0] == 1 and first_gene == gene_name:
        print(f"Plotting data for {first_gene}.")
        tcga_count = exp_df
    else:
        print(f"Searching {input_file} data for {gene_name}.")
        tcga_count = exp_df[exp_df.iloc[:, 0] == gene_name]

    cancer_types = [
        'ACC',       # Adrenocortical carcinoma
        'BLCA',      # Bladder urothelial carcinoma
        'BRCA',      # Breast invasive carcinoma
        # 'BRCA-Basal',# Breast invasive carcinoma - Basal
        # 'BRCA-Her2', # Breast invasive carcinoma - Her2
        # 'BRCA-LumA', # Breast invasive carcinoma - LumA
        # 'BRCA-LumB', # Breast invasive carcinoma - LumB
        'CESC',      # Cervical & endocervical cancer
        'CHOL',      # Cholangiocarcinoma
        'COADREAD',  # Colorectal adenocarcinoma
        'COAD',      # Colon adenocarcinoma
        'DLBC',      # Diffuse large B-cell lymphoma
        'ESCA',      # Esophageal carcinoma
        'GBM',       # Glioblastoma multiforme
        'HNSC',      # Head and neck squamous cell carcinoma
        # 'HNSC-HPV+', # Head and neck squamous cell carcinoma - HPV positive
        # 'HNSC-HPV-', # Head and neck squamous cell carcinoma - HPV negative
        'KICH',      # Kidney chromophobe
        'KIRC',      # Kidney clear cell carcinoma
        'KIRP',      # Kidney papillary cell carcinoma
        'LAML',      # Acute myeloid leukemia
        'LGG',       # Brain lower grade glioma
        'LIHC',      # Liver hepatocellular carcinoma
        'LUAD',      # Lung adenocarcinoma
        'LUSC',      # Lung squamous cell carcinoma
        'MESO',      # Mesothelioma
        'OV',        # Ovarian serous cystadenocarcinoma
        'PAAD',      # Pancreatic adenocarcinoma
        'PCPG',      # Pheochromocytoma and paraganglioma
        'PRAD',      # Prostate adenocarcinoma
        'READ',      # Rectum adenocarcinoma
        'SARC',      # Sarcoma
        'SKCM',      # Skin cutaneous melanoma
        # 'SKCM/Metastasis',  # Skin cutaneous melanoma - Metastasis
        'STAD',      # Stomach adenocarcinoma
        'TGCT',      # Testicular germ cell tumor
        'THCA',      # Thyroid carcinoma
        'THYM',      # Thymoma
        'UCEC',      # Uterine corpus endometrioid carcinoma
        'UCS',       # Uterine carcinosarcoma
        'UVM'        # Uveal melanoma
    ]

    cancer_abbr = input_file.split('.')[2]
    if cancer_abbr in cancer_types:
        new_columns = [tcga_count.columns[0]] + [cancer_abbr] * (
            len(tcga_count.columns) - 1)
        tcga_count.columns = new_columns

    tcga_df = tcga_count.stack().groupby(level=1).apply(list) \
        .apply(pd.Series).T
    tcga_df.drop(columns=['sample'], inplace=True)
#    print(tcga_df)
    return tcga_df

def plot(dataframe, gene_name, tissue_types):
    """
    Plot a boxplot of the log2 transformed normalized_counts, 
    separated by phenotype (tissue types) and outputs PNG file.
    """
    # Prepare data for seaborn
    melted_data = dataframe.melt(var_name='Tissue Type', 
        value_name='Expression')
    filtered = melted_data.dropna(subset=['Expression'])
    print(filtered)

    # Calculate counts of data points for each tissue type
    counts = filtered['Tissue Type'].value_counts()
    counts_dict = counts.to_dict()
    total_count = counts.sum()
    print(f"Total data points: {total_count}")

    # Setting the figure size
    plt.figure(figsize=(18, 10))

    # Determine colors for each tissue type
    color_palette = {tissue: 'blue' if tissue in tissue_types else 'red' 
                     for tissue in filtered['Tissue Type'].unique()}

    # Set dot transparency and thickness of lines
    dot_alpha = 0.25  # 50% transparency
    line_width = 0.5  # Adjust as needed

    # Boxplot with custom color for each tissue type
    sns.boxplot(x='Tissue Type', y='Expression', data=filtered, 
                palette=color_palette, whis=[0, 100], linewidth=2, 
                fliersize=0, showcaps=False, boxprops={'facecolor':'none', 
                'edgecolor':'black', 'linewidth': 2}, 
                whiskerprops={'color':'gray', 'linewidth': 2}, 
                medianprops={'color':'black', 'linewidth': 2}, 
                capprops={'color':'gray', 'linewidth': 2})

    # Stripplot with custom color for each tissue type
    sns.stripplot(x='Tissue Type', y='Expression', data=filtered, 
                jitter=True, color='red', edgecolor='black', 
                size=4, alpha=0.5)

    # Customize plot labels and title
    plt.xlabel('Tissue Type', fontsize=6, fontweight='bold')
    plt.ylabel(f'{gene_name} Expression (Log2(x+1) Normalized Count)',
        fontsize=6, fontweight='bold')
    plt.title(f'{gene_name} Expression by Tissue Type', fontsize=8,
        fontweight='bold')

    # Remove spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
 
    # Add counts to x-axis labels
    x_labels = [f"{label}\n(n={counts_dict.get(label, 0)})" 
        for label in filtered['Tissue Type'].unique()]
    plt.xticks(range(len(counts)), x_labels, rotation=90, 
        ha='center', fontsize=6)
    plt.subplots_adjust(bottom=0.1)

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

    args = parser.parse_args()
    gene_name = args.gene_name
    working_directory=args.working_directory

    gtex_input='''gtex_glp1r.gz'''
    gtex_phenotype='''GTEX_phenotype.gz'''
    gtex_df = read_gtex(gtex_input, gtex_phenotype, gene_name)
    
    TCGA_data = [
        'TCGA.ACC.sampleMap%2FHiSeqV2.gz',  # ACC
        'TCGA.BLCA.sampleMap%2FHiSeqV2.gz', # BLCA
        'TCGA.BRCA.sampleMap%2FHiSeqV2.gz', # BRCA
        'TCGA.CESC.sampleMap%2FHiSeqV2.gz', # CESC
        'TCGA.CHOL.sampleMap%2FHiSeqV2.gz', # CHOL
        'TCGA.COADREAD.sampleMap%2FHiSeqV2.gz', # COADREAD
        'TCGA.COAD.sampleMap%2FHiSeqV2.gz', # COAD
        'TCGA.DLBC.sampleMap%2FHiSeqV2.gz', # DLBC
        'TCGA.ESCA.sampleMap%2FHiSeqV2.gz', # ESCA
        'TCGA.GBM.sampleMap%2FHiSeqV2.gz', # GBM
        'TCGA.HNSC.sampleMap%2FHiSeqV2.gz', # HNSC
        'TCGA.KICH.sampleMap%2FHiSeqV2.gz', # KICH
        'TCGA.KIRC.sampleMap%2FHiSeqV2.gz', # KIRC
        'TCGA.KIRP.sampleMap%2FHiSeqV2.gz', # KIRP
        'TCGA.LAML.sampleMap%2FHiSeqV2.gz', # LAML
        'TCGA.LGG.sampleMap%2FHiSeqV2.gz', # LGG
        'TCGA.LIHC.sampleMap%2FHiSeqV2.gz', # LIHC
        'TCGA.LUAD.sampleMap%2FHiSeqV2.gz', # LUAD
        'TCGA.LUSC.sampleMap%2FHiSeqV2.gz', # LUSC
        'TCGA.MESO.sampleMap%2FHiSeqV2.gz', # MESO
        'TCGA.OV.sampleMap%2FHiSeqV2.gz', # OV
        'TCGA.PAAD.sampleMap%2FHiSeqV2.gz', # PAAD
        'TCGA.PCPG.sampleMap%2FHiSeqV2.gz', # PCPG
        'TCGA.PRAD.sampleMap%2FHiSeqV2.gz', # PRAD
        'TCGA.READ.sampleMap%2FHiSeqV2.gz', # READ
        'TCGA.SARC.sampleMap%2FHiSeqV2.gz', # SARC
        'TCGA.SKCM.sampleMap%2FHiSeqV2.gz', # SKCM
        'TCGA.STAD.sampleMap%2FHiSeqV2.gz', # STAD
        'TCGA.TGCT.sampleMap%2FHiSeqV2.gz', # TGCT
        'TCGA.THCA.sampleMap%2FHiSeqV2.gz', # THCA
        'TCGA.THYM.sampleMap%2FHiSeqV2.gz', # THYM
        'TCGA.UCEC.sampleMap%2FHiSeqV2.gz', # UCEC
        'TCGA.UCS.sampleMap%2FHiSeqV2.gz', # UCS
        'TCGA.UVM.sampleMap%2FHiSeqV2.gz'  # UVM
    ]

    gtex_df = pd.read_csv(f'{gene_name}.csv')
    tcga_df_list = []

    # Iterate over TCGA data paths and call read_tcga function
    for tcga_file in TCGA_data:
        tcga_count = read_tcga(f'{working_directory}/{tcga_file}', gene_name)
        tcga_df_list.append(tcga_count)

    combined_df = pd.concat([gtex_df] + tcga_df_list, axis=1)
    combined_df.to_csv(f'{gene_name}.csv', index=False)

    plot(combined_df, gene_name, gtex_df)
    
if __name__ == "__main__":
    main()