"""
contact:    wlee9829@gmail.com
date:       2024_03_15
python:     python3.10
script:     crssant_eval.py

This script is multifunctional tool for assessing CRSSANT resource use and
generates a multiple regression analysis based on the generated data. 
"""
# Define version
__version__="1.0.1"

# Version notes
__update_notes__="""

1.0.1
    - Wrote function to read CRSSANT output files, write to Log.out file.

1.0.0
    - Calculations for run speed based on multiple regression.

To Do:

    - Improve error handling and logging.
    - Add CSV output for file types

"""

# Import Packages
from datetime import datetime
import argparse
import os
import pandas as pd
import subprocess
import sys
import textwrap
import time
from sklearn.linear_model import LinearRegression
import numpy as np

###########################################################################
# 1. Define common functions for the data analysis.

def timenow():
    """
    Returns the current timestamp as a string.

    Returns:
        str: Current timestamp in format 'YYY-MM-DD HH:MM:SS'.
    """
    time = str(datetime.now())[:-7]

    return time

def read_outfile(minreads, outfile):
    """
    Extracts data from the CRSSANT read statistics from the log.out files.

    input:  crssant_skip.out, crssant_low.out, crssant_all.out
    """
    out_data = {'Sample_Name':[], 'rRNA_reads':[], 'Reads':[], 'Genes':[]}
    run_times = {'Reading':[], 'Clustering':[], 'Assembly':[], 'Run_Time':[]}
    is_low = False

    with open(outfile, 'r') as file:
        for line in file:
            if '.cliques.t_o0.2.sam' in line:
                file_names = (line.split()[1]
                    .replace('.cliques.t_o0.2.sam', ''))
                out_data['Sample_Name'].append(file_names)
            elif "Subsampling from" in line and "low" in outfile:
                rRNA_reads = line.split()[2]
                is_low = True
                out_data['rRNA_reads'].append(rRNA_reads)
            elif 'Number of reads:' in line:
                num_reads = line.split()[3]
                out_data['Reads'].append(num_reads)
            elif 'Number of genes:' in line:
                num_genes = line.split()[3]
                out_data['Genes'].append(num_genes)
            
            #Read job resources to file.
            elif "SLURM_JOB_ID" in line:
                job_id = line.split()[2]
                subprocess.run(f'jobinfo {job_id} >> resources_{outfile}',
                    shell=True)

            #Read DG run statistics.    
            elif 'Completed reading input in:' in line:
                reading = line.split()[4].replace('s', '')
                run_times['Reading'].append(reading)
            elif 'Completed clustering in:' in line:
                clustering = line.split()[3].replace('s', '')
                run_times['Clustering'].append(clustering)
            elif 'Completed assembly in:' in line:
                assembly = line.split()[3].replace('s', '')
                run_times['Assembly'].append(assembly)
            elif 'real' in line:
                time = line.split()[1].replace('s', '')
                minutes = time.split("m")[0]
                seconds = time.split("m")[1]
                real = int(minutes)*60+float(seconds)
                run_times['Run_Time'].append(real)

        if not is_low:
            for _ in range(len(out_data['Sample_Name'])):
                out_data['rRNA_reads'].append("0")
    
    resources = {'CPUs':[], 'Tasks':[], 'Mem_Req':[], 'Mem_Used':[],
                'Req_Time':[], 'Wall_Time':[], 'CPU_Time':[]}    
    with open(f'resources_{outfile}', 'r') as file:
        for line in file:
            if 'CPUs' in line:
                cpus = line.split()[2]
                resources['CPUs'].append(cpus)
            elif 'Tasks' in line:
                tasks = line.split()[2]
                resources['Tasks'].append(tasks)
            elif 'Reserved memory' in line:
                mem_req = line.split()[3].replace('G', '')
                resources['Mem_Req'].append(mem_req)
            elif 'Max memory used' in line:
                mem_used = line.split()[4].replace('G', '')
                resources['Mem_Used'].append(mem_used)
            elif 'Reserved walltime' in line:
                req_time = line.split()[3]
                resources['Req_Time'].append(req_time)
            elif 'Used walltime' in line:
                walltime = line.split()[3]
                resources['Wall_Time'].append(walltime)
            elif 'Used CPU time' in line:
                cputime = line.split()[4]
                resources['CPU_Time'].append(cputime)

    read_stats = pd.DataFrame(out_data)
    cpu_stats = pd.DataFrame(resources)
    resources = [cpu_stats] * len(read_stats)
    resource_stats = pd.concat(resources, ignore_index=True)
    time_stats = pd.DataFrame(run_times)
    log_out = pd.concat([read_stats, resource_stats, time_stats], 
        axis=1).reset_index(drop=True)

    os.remove(f'resources_{outfile}')
    try:
        if os.path.exists("Log.out"):
            existing_df = pd.read_csv("Log.out", sep='\t')
            combine = pd.concat([existing_df, log_out], ignore_index=True)
            combine.to_csv("Log.out", sep='\t', index=False)
        else:
            log_out.to_csv("Log.out", sep='\t', index=False)
    except Exception as e:
        print(f"An error occured:", {e})

    return log_out

def skip_read(logfile): 
    """
    Performs a multiple regression analysis using data from the CRSSANT run
    times that "skipped" the hs45S chromosome that contains rRNA reads. The
    data analyzed here is the mean between multiple datasets.
    """

    # Set up arrays from multiple runs.
    #######################################################
    df = pd.read_csv(logfile, sep='\t')
    skipped_df = df[~df['Sample_Name'].str.contains('all') & df['Sample_Name'].str.contains('low')]
    
    runtime = df['Assembly'].tolist()
    reads = df['Reads'].tolist()
    genes = df['Genes'].tolist()

    variables = np.array([reads, genes]).T
    outcome = np.array(runtime) 
    model = LinearRegression().fit(variables, outcome)

    print('Coefficients:', model.coef_)
    print('Intercept:', model.intercept_)

    #return model.coef_, model.intercept_

###########################################################################
# 2. Main, define accepted arguments. 
def main():
    parser = argparse.ArgumentParser(
        prog="benchmark.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        formatter_class=argparse.RawTextHelpFormatter,      
        description=textwrap.dedent("""\
###########################################################################
Pass
###########################################################################
"""),    
usage='''

- For reading the output and writing to a log file only: 
\npython3 %(prog)s [-L]

''')

    parser.add_argument('-L', '--log-out', action='store_true',
        help='Parses crssant output files, writes resource use and statistics'
            ' to a log file.')
    parser.add_argument('-min', '--min_coverage', 
        help="Min_coverage is defined as minimum number of reads for a region"
        " before it is considered well covered. Defaults to [5]")
    parser.add_argument('-V', '--version', action='version', 
        version=f'%(prog)s {__version__}\n{__update_notes__}', 
        help='Print version + update notes and exit.')
    ##########################################################################

    args = parser.parse_args()

    if args.log_out:
        if args.min_coverage:
            minreads = args.min_coverage
        else:
            minreads = "5"

    outfiles = ["crssant_skip.out", "crssant_low.out", "crssant_all.out"]
    for outfile in outfiles:
        try:
            print(f'{outfile}')
            read_outfile(minreads, outfile)
        except Exception as e:
            print(f"An error occured:", {e})

    skip_read("Log.out")

if __name__ == "__main__":
    main()
sys.exit()
