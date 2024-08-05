"""
contact:    wlee9829@gmail.com
date:       2000_00_00
python:     python3.10
script:     SCRIPT_NAME.py

DESCRIPTION HERE
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
from matplotlib.colors import Normalize
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys
import textwrap

###########################################################################
# 1. Set up functions
def function(arg):
    """
    DOCSTRING

    Args:

    Returns:

    """
    return

def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
#        prog="<SCRIPT_NAME>",
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent("""\
###########################################################################
NOTE: NOTES HERE

DESCRIPTION, NUMBERED

###########################################################################
"""),
    usage=
"""
    \npython3 %(prog)s <ARGS>
""")
    parser.add_argument(
        'ARG', type=str,
        help='HELP DESCRIPTION.')

    return parser.parse_args()

def main(args):
    """
    Main function to execute the script.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """

if __name__ == "__main__":
    args = parse_args()
    main(args)