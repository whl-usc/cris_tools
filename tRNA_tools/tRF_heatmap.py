# Import Packages
from textwrap import dedent
import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import re

###########################################################################
# 1. Set up functions
def read_input(file_path, log2scale=False):
    """
    Read and process the input file to create the matrix of tRNA x tRNA.

    Args:
        file_path (str): Path to the input file.
        log2scale (bool): Whether to apply log2 transformation.

    Returns:
        dict: Dictionary with (start_pos, end_pos) as keys and maxRPM as values.
        int: Length of RNA sequences.
    """
    with open(file_path, 'r') as infile:
        lines = infile.readlines()

    records = [lines[n:n+4] for n in range(0, len(lines), 4)]

    position_rpm = {}
    max_RNA_length = 0

    for record in records:
        sequence = record[0].strip()
        if not sequence:
            continue

        maxRPM = float(record[2].strip())
        start_pos = len(re.match(r"^(\s*)", record[0]).group(0)) + 1
        end_pos = start_pos + len(sequence) - 1

        max_RNA_length = max(max_RNA_length, len(record[0]))

        # Add the RPM value to the appropriate (start, end) pair
        if (start_pos, end_pos) not in position_rpm:
            position_rpm[(start_pos, end_pos)] = 0
        position_rpm[(start_pos, end_pos)] += maxRPM

    # Apply log2 scaling if needed
    if log2scale:
        for key in position_rpm:
            if position_rpm[key] > 0:
                position_rpm[key] = math.log2(position_rpm[key])

    # Print example data for verification
    print(f"Max RNA length: {max_RNA_length}")

    return position_rpm, max_RNA_length

def plot_heatmap(position_rpm, RNAlen, log2scale, output_file, title):
    """
    Plot a heatmap with start positions on the x-axis and end positions on the y-axis.

    Args:
        position_rpm (dict): Dictionary with (start_pos, end_pos) as keys and RPM values as values.
        RNAlen (int): Length of RNA sequences.
        log2scale (bool): Whether the data is log2 transformed.
        output_file (str): Base name for the output file.
        title (str): Title of the heatmap.
    """
    # Initialize an empty matrix of zeros
    matrix = np.zeros((RNAlen, RNAlen))
    
    # Populate the matrix with maxRPM values from the dictionary
    for (start_pos, end_pos), maxRPM in position_rpm.items():
        if 1 <= start_pos <= RNAlen and 1 <= end_pos <= RNAlen:
            matrix[end_pos, start_pos] = maxRPM
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')
    heatmap = ax.imshow(matrix, cmap='Blues', interpolation='nearest', origin='lower')

    # Add a colorbar with custom ticks
    nmax = matrix.max()
    tick1 = int(str(int(nmax / 5))[0] + (len(str(int(nmax / 5))) - 1) * '0') if nmax > 0 else 1
    tickslist = list(range(0, int(nmax) + tick1, tick1))
    colorbar = plt.colorbar(mappable=heatmap, ax=ax, shrink=0.25, aspect=10, ticks=tickslist)
    colorbar_label = 'Log2(Max RPM)' if log2scale else 'Max RPM'
    colorbar.set_label(colorbar_label, fontsize=8)

    # Set axis labels
    ax.set_xlabel('Start Position', fontsize=8)
    ax.set_ylabel('End Position', fontsize=8)

    # Add major tick marks every 2 positions
    major_ticks = np.arange(0, RNAlen + 2, 2)
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)

    # Rotate x-axis labels
    plt.xticks(rotation=90, fontsize=8)

    # Adjust tick parameters
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.tick_params(axis='both', which='minor', length=0)
    
    # Add gridlines to encapsulate each cell
    ax.set_xticks(np.arange(0.5, RNAlen, 1), minor=True)
    ax.set_yticks(np.arange(0.5, RNAlen, 1), minor=True)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=0.3, zorder=1)

    # Adjust the limits to match the heatmap data bounds
    ax.set_xlim(-0.5, RNAlen - 0.5)
    ax.set_ylim(-0.5, RNAlen - 0.5)
    ax.set_aspect('equal')
    
    # Add title to the heatmap
    ax.set_title(title, fontsize=10)

    # Save the heatmap plot
    plt.tight_layout()
    plt.savefig(f"{output_file}_heatmap{'_log2' if log2scale else ''}.png", dpi=300)
    plt.close()
    print(f"Plot saved as {output_file}_heatmap{'_log2' if log2scale else ''}.png")

def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Plot heatmap of tRF positions using MINTbase data."
    )
    parser.add_argument(
        'file_path', type=str,
        help='Path to the input text file containing tRF information.'
    )
    parser.add_argument(
        'log2scale', type=int,
        help='Perform log2 transformation: 1 for yes, 0 for no.'
    )
    parser.add_argument(
        'output_file', type=str, nargs='?',
        default=None,
        help='Base filename for saving the plots. If not provided, a default name will be generated.'
    )
    return parser.parse_args()

def main():
    """
    Main function to execute the script.
    """
    args = parse_args()
    log2scale = bool(args.log2scale)
    
    # Read input and get RPM values for start and end positions
    position_rpm, RNAlen = read_input(args.file_path, log2scale)

    # Generate default output file name if not provided
    if args.output_file is None:
        base_name = os.path.basename(args.file_path)
        base_name = base_name.replace('tRNA-', '').replace('_MINTbase.txt', '')
        output_file = base_name
    else:
        output_file = args.output_file
    
    # Use the base name for the title
    title = output_file
    
    # Plot heatmap with start and end RPM values
    plot_heatmap(position_rpm, RNAlen, log2scale, output_file, title)

if __name__ == "__main__":
    main()