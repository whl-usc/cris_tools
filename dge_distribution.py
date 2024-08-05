"""
contact:    wlee9829@gmail.com
date:       2024_08_02
python:     python3.10
script:     dge_distribution.py

This Python script summarizes sample distribution for expression of gene
expression counts downloaded from the UCSC Xena web platform.
"""

# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
1.0.0
    -   Initial commit, set up function logic and styling.
"""

# Import Packages
from datetime import datetime
from scipy.stats import mannwhitneyu
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import textwrap

###########################################################################
# 1. Set up functions

def plot_expression_histograms(file_path1, file_path2, separated, output):
    """
    Plot histograms of RNA expression counts for two specified genes from
    separate CSV files.

    Args:
        file_path1 (str): CSV file containing expression data for first gene.
        file_path2 (str): CSV file containing expression data for second gene.
    """
    data1 = pd.read_csv(file_path1)
    data2 = pd.read_csv(file_path2)

    gene1 = file_path1.replace("RSEM_", "").replace("DSEQ2_", 
        "").replace(".csv", "")
    gene2 = file_path2.replace("RSEM_", "").replace("DSEQ2_", 
        "").replace(".csv", "")

    melted_data1 = data1.melt(var_name='Tissue', 
        value_name='Expression Count').dropna()
    melted_data2 = data2.melt(var_name='Tissue', 
        value_name='Expression Count').dropna()

    mw_stat, mw_p_value = mannwhitneyu(melted_data1['Expression Count'],
        melted_data2['Expression Count'], alternative='two-sided')

    print(f"Mann-Whitney U Test: Statistic={mw_stat:.2e}," 
        f" p-value={mw_p_value:.2e}")

    plt.figure(figsize=(10, 6))

    sns.histplot(melted_data1['Expression Count'], edgecolor=None, bins=30,
        color='grey', label=gene1, alpha=0.85)
    sns.histplot(melted_data2['Expression Count'], edgecolor=None, bins=30,
        color='black', label=gene2, alpha=0.85)
    sns.despine()
    plt.title('Gene Expression Counts')
    plt.xlabel('Gene Expression (log2(x+1) Count)')
    plt.ylabel('Frequency (Number of Samples)')
    plt.legend()

    # Add p-value text to the plot
    plt.text(0.25, 0.95, f'Mann-Whitney p-value={mw_p_value:.2e}', 
             transform=plt.gca().transAxes, fontsize=8, 
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))

    plt.tight_layout()

    if output:
        outname = str(output)
    else: 
        outname = str(f"{gene1}_{gene2}_distribution.png")

    plt.savefig(outname, dpi=400)
    time = str(datetime.now())[:-7]
    print(f"Plot saved as {outname} on {time}.")

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
    parser.add_argument('-s', '--separated', type=str, 
        help='Pool sample types or keep separate. Default is none, '
            'which pools all tissue types ',
        default=False)
    parser.add_argument(
        '-o', '--output', type=str, help='Filename for saving the plot')

    return parser.parse_args()

def main():
    """
    Main function to execute the script.
    """
    args = parse_args()
    plot_expression_histograms(args.file_path1, args.file_path2, 
        args.separated, args.output)

if __name__ == "__main__":
    main()