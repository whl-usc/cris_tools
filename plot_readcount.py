"""
contact:    wlee9829@gmail.com
date:       2024_07_03
python:     python3.10
script:     plot_readcount.py

This Python script plots gene counts in user-defined styles. 
Graph can be log(2) transformed to find normal distribution.

"""
# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
1.0.0   
    -   Initial commit, set up functions.
"""

# Import Packages
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

###########################################################################
def readfile(file_path):
    """
    Opens coverage text file and parses information.
    """
    coverage = pd.read_csv(file_path, delimiter='\t')

    print(coverage)
    return coverage

def plot_gene_counts(dataframe, output_file):
    """
    Plots log2 transformed gene counts from the dataframe.
    """
    df['log_count'] = np.log2(df['count'] + 1)

def parse_args():
    """
    Main function to set up argument prasing, handles input arguments, calls
    relevant functions for processing and plotting read count data.
    
    Args:

    """
    parser = argparse.ArgumentParser(
        prog="plot_readcount.py",        
        description="Plot gene counts from a tab-separated file.")
    
    parser.add_argument("input_file", 
        help="Path to the input file containing gene count data.")
    
    parser.add_argument("output_file", 
        help="Path to save the output plot as a PNG file.")

def main():
    input_file = args.input_file
    output_file = args.output_file

    readfile(input_file)

if __main__=="__main__":
    args = parse_args()
    main()
    sys.exit() 