"""
contact:    wlee9829@gmail.com
date:       2024_03_15
python:     python3.10
script:     crssant_benchmark.py

This script is standalone tool for benchmarking CRSSANT analysis based on 
previously generated data. 
"""
# Define version
__version__="1.0"

# Version notes
__update_notes__="""
Version 1.0:
- Benchmarking based on multiple regression.

To Do:
"""
# Import Packages
from datetime import datetime
import argparse
import os
import sys
import textwrap
import time
from sklearn.linear_model import LinearRegression
import numpy as np

###########################################################################
# 1. Perform multiple regression analysis for the datasets.

def timenow():
    """
    Returns the current timestamp as a string.

    Returns:
        str: Current timestamp in format 'YYY-MM-DD HH:MM:SS'.
    """
    time = str(datetime.now())[:-7]

    return time

def skip():	
	"""
	Performs a multiple regression analysis using data from the CRSSANT run times that "skipped" the hs45S chromosome that contains rRNA reads. The data performed here is the mean between multiple datasets.
	"""

	# Set up arrays from multiple runs.

	#######################################################
	# SHARC-exo skip
	reads = [26747, 46079, 17821, 2066, 39674, 22777, 
			6732, 76884, 69212, 23079, 6245]
	
	genes = [18, 22, 13, 0, 2, 12, 
			4, 32, 19, 12, 4]
	
	runtime = [47, 125, 47, 28, 32, 38, 
			36, 389, 47, 43, 28]

	variables = np.array([reads, genes]).T
	outcome = np.array(runtime) 
	model = LinearRegression().fit(variables, outcome)

	print('Coefficients:', model.coef_)
	print('Intercept:', model.intercept_)

	#######################################################
	# PARIS-exo skip
	reads = [55860, 5128, 19209, 40823]
	genes = [24, 2, 3, 12]
	runtime = [72, 66, 67, 61]

	variables = np.array([reads, genes]).T
	outcome = np.array(runtime) 
	model = LinearRegression().fit(variables, outcome)

	print('Coefficients:', model.coef_)
	print('Intercept:', model.intercept_)

	#######################################################
	# PARIS2-Virus skip
	reads = []
	genes = []
	runtime = []

	#return model.coef_, model.intercept_

# 3. Main, define accepted arguments. 
def main():
    parser = argparse.ArgumentParser(
        prog="benchmark.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
###########################################################################
Pass
"""),    
usage='''
\npython3 %(prog)s
''')

    parser.add_argument('-V', '--version', action='version', 
        version=f'%(prog)s {__version__}')

    args = parser.parse_args()

    skip()

if __name__ == "__main__":
    main()
sys.exit()