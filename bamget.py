"""
contact:    wlee9829@gmail.com
date:       2024_02_09
python:     python3.10
script:     bamget.py

This script is tool for preparing a BED file for CRSSANT analysis based on 
either the minimum coverage per annotated gene regions or nucleotide. Files
should be named according to the shell scripts in the rna2d3d repository. 
See the options associated with the analysis types.

Pedersen, B. S.; Quinlan, A. R.; Mosdepth: quick coverage calculation for
genomes and exomes, Bioinformatics, Volume 34, Issue 5, March 2018, 867–868

Danecek, P.; Bonfield, J. K.; Liddle, J.; Marshall, J.; Ohan, V.; Pollard, 
M. O.; Whitwham, A.; Keane, T.; McCarthy, S. A.; Davies, R. M.; Li, H.
GigaScience, Volume 10, Issue 2, February 2021

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification
of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019).
"""
# Define version
__version__ = "1.0.0"

# Version notes
__update_notes__ = """
1.0.0
    -   Gene-region based analysis function and output to BED6 file is complete.

To Do:
    -   Add position-based analysis
    -   Add support for non homo-sapiens (HS) genome files.
"""

# Import Packages
from datetime import datetime
import argparse
import glob
import gzip
import os
import random
import subprocess
import sys
import textwrap
import time

###########################################################################
# 1. Basic functions that common to either analyses function.

def timenow():
    """
    Returns the current timestamp as a string.

    Returns:
        str: Current timestamp in format 'YYY-MM-DD HH:MM:SS'.
    """
    time = str(datetime.now())[:-7]

    return time

def collapse_gene_regions(annotation_file, skip_chromosome=None):
    """
    Prepares a "collapsed" gene region containing file by parsing through a 
    modified GTF annotation file. This function was written with the hg38 
    annotation in mind, but will be updated for other gene annotation files
    in future iterations. See help section for details about preparing the 
    modified annotation.bed from a Gencode GTF file.

    Note that this function is hardcoded to filter for any genes that fall 
    into the ENCODE blacklist region. In this way, problematic regions of the 
    genome are excluded.

    Args:
        annotation_file: input annotation file with gene information.

    Returns:
        gene_regions (dict): gene names as keys, chromosome, start, stop, 
        and strand information as values.        
    """
    valid_chromosomes = set([f"chr{i}" for i in range(1, 23)] + 
        ["chrX", "chrY", "chrM"] + [])
    # Add extra "chromosomes" as needed into the last list. This function
    # will eventually be updated to be able to handle more unique cases.

    # Blacklist is from ENCODE, regions that have anomalous, unstructured, or 
    # unordinarily high signal in multiple NGS experiments.
    # Defined blacklisted regions as a list tuples (chromosome, start, stop):
    blacklist_regions = [
        ("chr10", 38528030, 38529790), ("chr10", 42070420, 42070660),
        ("chr16", 34571420, 34571640), ("chr16", 34572700, 34572930),
        ("chr16", 34584530, 34584840), ("chr16", 34585000, 34585220),
        ("chr16", 34585700, 34586380), ("chr16", 34586660, 34587100),
        ("chr16", 34587060, 34587660), ("chr16", 34587900, 34588170),
        ("chr16", 34593000, 34593590), ("chr16", 34594490, 34594720),
        ("chr16", 34594900, 34595150), ("chr16", 34595320, 34595570),
        ("chr16", 46380910, 46381140), ("chr16", 46386270, 46386530),
        ("chr16", 46390180, 46390930), ("chr16", 46394370, 46395100),
        ("chr16", 46395670, 46395910), ("chr16", 46398780, 46399020),
        ("chr16", 46400700, 46400970), ("chr1", 124450730, 124450960),
        ("chr20", 28513520, 28513770), ("chr20", 31060210, 31060770),
        ("chr20", 31061050, 31061560), ("chr20", 31063990, 31064490),
        ("chr20", 31067930, 31069060), ("chr20", 31069000, 31069280),
        ("chr21", 8219780, 8220120), ("chr21", 8234330, 8234620),
        ("chr2", 90397520, 90397900), ("chr2", 90398120, 90398760),
        ("chr3", 93470260, 93470870), ("chr4", 49118760, 49119010),
        ("chr4", 49120790, 49121130), ("chr5", 49601430, 49602300),
        ("chr5", 49657080, 49657690), ("chr5", 49661330, 49661570)]

    # Check if the output BED file already exists; check column formatting. 
    if os.path.exists('annotation.bed'):
        try:
            print(f"\nAttempting to read regions from: 'annotation.bed'")
            gene_regions = {}
            with open('annotation.bed', 'r') as f:
                for line in f:
                    parts = line.split("\t")
                    chromosome = parts[0].strip("'\"")
                    if skip_chromosome and chromosome in skip_chromosome:
                        continue
                    start = int(parts[1])
                    stop = int(parts[2])
                    gene_name = parts[3]
                    strand = parts[4]
                    gene_regions[gene_name] = (chromosome, start, 
                        stop, strand)
        except Exception as e:
            error_message = (f"\nERROR: Failed to parse annotation.bed. Check"
                f" to see if the columns are as defined in the help section.")
            print(error_message)
            sys.exit()

    else:
        try:
            if annotation_file is not None:
                print(f"\nSorting gene regions from {annotation_file}" 
                    f" into 'annotations.bed'...") 
            gene_regions = {}
            with open(annotation_file, 'r') as f:
                for line in f:
                    parts = line.split("\t")
                    chromosome = parts[0].strip("'\"")
                    if skip_chromosome and chromosome in skip_chromosome:
                        continue
                    if chromosome not in valid_chromosomes:
                        chromosome = chromosome[3:]
                    start = int(parts[1])
                    stop = int(parts[2])
                    strand = parts[3]
                    gene_name = parts[4]

                    if (chromosome, start, stop) not in blacklist_regions:
                        if gene_name not in gene_regions:
                            gene_regions[gene_name] = (chromosome,
                                start, stop, strand)
                        else:
                            gene_regions[gene_name] = (chromosome,
                                min(start, gene_regions[gene_name][1]), 
                                max(stop, gene_regions[gene_name][2]), strand)

                # Write gene regions to BED file
                with open("annotation.bed", 'w') as file:
                    for gene_name, gene_info in gene_regions.items():
                            chromosome, start, stop, strand = gene_info
                            file.write(
                                f"{chromosome}\t{start}\t{stop}\t{gene_name}"
                                f"\t{strand}\n")
        except Exception as e:
            if annotation_file == None:
                error_message = (f"\nERROR: No annotation file provided.")
            else:
                error_message = (f"\nERROR: Failed to parse annotation file."
                    f" Check if file exists and is formatted correctly.")
            print(error_message)
            sys.exit()

    return gene_regions

def sort_bam(input_file, max_reads=None):
    """
    Checks the input file type, converts to BAM if necessary, 
    then sorts the BAM file for input into the analysis functions.
    Requires samtools packages. 

    If the argument max_reads is passed, a maximum number of randomized reads
    per given gene region will be pulled from the input file. 
    Requires the bedtools package and annotation file.

    Args:
        input_file: either a SAM or BAM alignment mapping file.
        
    (Optional):
        max_reads: integer, maximum number of reads per gene region.  

    Returns:
        file: sorted_bam file using samtools.
    """
    input_prefix = os.path.splitext(input_file)[0]

    try:
        if input_file.endswith(".sam"):
            sam_to_bam = f"{input_prefix}.bam"
            subprocess.run(['samtools', 'view', '-bS', '-o', 
                sam_to_bam, input_file], check=True)

        else:
            sam_to_bam = input_file
            sorted_bam = f"{input_prefix}_sorted.bam"

        if max_reads is not None:
            try:
                header_file = "header.tmp"
                hs45S_reads = "hs45S.tmp"
                chrom = "chrom.tmp"
                with open(header_file, "w") as header:
                    subprocess.run(['samtools', 'view', '-H', sam_to_bam], 
                        stdout=header, check=True)

                is_hs45S = f'awk \'$3 == "hs45S" {{print}}\''
                hs45S_extract = ['samtools', 'view', sam_to_bam, '|', 
                    is_hs45S, '>', hs45S_reads]
                subprocess.run(' '.join(hs45S_extract), 
                    shell=True, check=True)

                not_hs45S = f'awk \'$3 != "hs45S" {{print}}\''
                chrom_extract = ['samtools', 'view', sam_to_bam, '|', 
                    not_hs45S, '>', chrom]
                subprocess.run(' '.join(chrom_extract), 
                    shell=True, check=True)

                with open("hs45S.tmp", 'r') as f:
                    hs45S_reads = f.readlines()
                    hs45S_count = len(hs45S_reads)
                    print(f"\nSubsampling from {hs45S_count} rRNA reads...")
                    if hs45S_count > int(max_reads):
                        hs45S_low = random.sample(hs45S_reads, int(max_reads))
                    else:
                        hs45S_low = hs45S_reads

                with open("hs45S_low.tmp", 'w') as f:
                    f.writelines(hs45S_low)

                sam_to_bam = (f"{input_prefix}_low.bam")
                sorted_bam = f"{input_prefix}_low_sorted.bam"
                
                out_files = ["header.tmp", "hs45S_low.tmp", "chrom.tmp"]
                with open(sam_to_bam, 'w') as output_file:
                    for filename in out_files:
                        with open(filename, 'r') as input_file:
                            for line in input_file:
                                output_file.write(line)

            except Exception as e:
                error_message = (f"\nERROR: Sub-sampling failed."
                    f" Check if file contains reads mapped to 'hs45S'"
                    f" chromosome.")
                print(error_message)

            finally:
                for file in ["header.tmp", "hs45S.tmp", "hs45S_low.tmp",
                    "chrom.tmp"]:
                    if os.path.exists(file):
                        os.remove(file)

        subprocess.run(['samtools', 'sort', '-o', sorted_bam, sam_to_bam])
        subprocess.run(['samtools', 'index', sorted_bam])

    except Exception as e:
        error_message = (f"\nERROR: Failed to parse the provided file input."
            f" Check to see if file format is correct. Details: {e}.")
        print(error_message)
        sys.exit()       
    
    return sorted_bam

###########################################################################
# 2. Gene region based analysis - calculation and writing out. 

def mosdepth_regions(sorted_bam, annotation, min_coverage=None, 
    skip_chromosome=None):
    """
    Function that parses the provided sorted_bam file to determine gene
    regions that meet or exceed minimum coverage as defined by user. Uses
    mosdepth functions to determine reads at each gene region, then collects
    alignments >= min_coverage. Defining num_cores should speed up the depth
    calculations, but benchmarking for mosdepth shows no notable improvement
    past 4 threads.

    Args:
        sorted_bam:         PATH to BAM file.
        annotation:         bed file that defines the gene regions, derived 
                            from collapse_gene_regions() function.
    (Optional)
        min_coverage:       Integer defining minimum number of reads/gene
                            before the region is considered covered.
        skip_chromosome:    comma separated list of chromosomes to skip 
                            from the input file.
    Returns:
        covered_genes (list): genes with coverage that exceed min_coverage
    """
    covered_genes = []
    if os.path.exists('tmp.regions.bed.gz'):
        print(f"\nMean coverage summary exists. Using for calculations...")
        
    else:            
        if os.cpu_count() >= 4:
            num_threads = str(4) 
        else:
            num_threads = os.cpu_count()

        mosdepth_args = ['mosdepth', '-n',
                        '-b', 'annotation.bed', 
                        '-t', str(num_threads), 
                        'tmp', sorted_bam]
        subprocess.run(mosdepth_args, check=True)

    if min_coverage is None:
        min_coverage = 10

    if skip_chromosome:
        skipped = skip_chromosome.split(',')
    else:
        skipped = []

    with gzip.open('tmp.regions.bed.gz', 'rt') as f:
        for line in f:
            chrom, start, end, gene, cov_mean = line.strip().split()
            if chrom in skipped:
                continue
            if float(cov_mean) >= float(min_coverage):
                covered_genes.append(gene)

    return covered_genes

def region_to_bed(covered_genes, gene_regions, output):
    """
    Reads the return from region based analysis to write output to file.

    Args:
        covered_genes (list):   List of genes that have high_coverage_regions.
        gene_regions (dict):    Dictionary of all genomic regions from the
                                annotation file.
        output (str):           PATH to the output BED file. 

    Returns:
        None
    """
    mode = 'a' if os.path.exists(output+'.bed') else 'w'

    existing_genes = set()
    if mode =='a':
        with open(output+'.bed', 'r') as f:
            for line in f:
                gene_name = line.split('\t')[3].strip()
                existing_genes.add(gene_name)

    with open(output+'.bed', mode, newline='') as f:
        for gene_name in covered_genes:
            if gene_name not in existing_genes:
                if "tRNA" not in gene_name: 
                    if gene_name in gene_regions:
                        chrom, start, stop, strand = gene_regions[gene_name]
                        line = (f"{chrom}\t{start}\t{stop}\t "
                            f" {gene_name}\t1000\t{strand}\n")
                        f.write(line)
                        existing_genes.add(gene_name)

###########################################################################
# 2. Position based analysis - calculation and plotting to histograms. 

    # THIS FUNCTION WILL PLOT A HISTOGRAM WITH THE OUTPUT OF MOSDEPTH
    # INSTEAD OF PRODUCING A FILE FOR GENES.BED

    # def mosdepth_positions(sorted_bam, window_size=None):
    #     """
    #     Function that parses the provided sorted_bam file to determine positions
    #     that meet or exceed minimum coverage as defined by user. Uses mosdepth
    #     functions to determine reads at each nucleotide "window", then collects 
    #     alignments >= min_coverage which overlap the gene_regions. Defining
    #     num_cores may speed up depth calculations, with little improvement past 
    #     4 threads.

    #     Args/Returns:
    #         Same as mosdepth_regions() function.
    #     """
    #     if os.path.exists('tmp.regions.bed.gz'):
    #         print(f"WARNING: Does the tmp.regions.bed")
    #         with gzip.open('tmp.regions.bed.gz', 'rt') as f:
    #             for i, line in enumerate(f):
    #                 chrom, start, end, cov_mean = line.strip().split()
    #                 window_size = int(end) - int(start)
    #                 break

    #         print(f"Mean coverage summary for positions already exists with"
    #             f" window size of: {window_size}" 
    #             f" If window is not as expected, remove file and try again.")
    #     else:

    #         if os.cpu_count() >= 4:
    #             num_threads = str(4) 
    #         else:
    #             num_threads = str(os.cpu_count())

    #         # Includes a conditional window_size, defaults to 20 nt.
    #         # This might need to be changed to 15 as a default for smaller genes. 
    #         if window_size is None:
    #             window = 20
    #         else:
    #             window = int(window_size)
      
    #         mosdepth_args = ['mosdepth', '-n', 
    #                         '-b', str(window), 
    #                         '-t', num_threads, 
    #                         'tmp', sorted_bam]
    #         subprocess.run(mosdepth_args, check=True)

    #     covered_positions = []    
    #     with gzip.open('tmp.regions.bed.gz', 'rt') as f:
    #         for line in f:
    #             chrom, start, end, cov_mean = line.strip().split()
    #             if cov_mean >= min_coverage:
    #                 gene = f"{start}_{end}"
    #                 covered_positions = [chrom, start, end, gene]
    #                 covered_positions.append(gene)
    #     return covered_positions

    # RE-WRITE THIS FUNCTION. MAKE IT GENERATE A HISTOGRAM INSTEAD.
    # def position_to_bed(covered_positions, output):
    #     """
    #     Reads the return from  position based analysis to write output to file. 
    #     Args:
    #         covered_positions:  Set of genes that have high coverage positions.
    #         output (str):   PATH to the output BED file. 
    #     Returns:
    #         None
    #     """
    #     mode = 'a' if os.path.exists(output) else 'w'           
    #     existing_genes = set()
    #     if mode =='a':
    #         with open(output+'.bed', 'r') as f:
    #             for line in f:
    #                 existing_genes.add(line.split('\t')[3]).strip()
    #     with open(output+'.bed', mode, newline='') as f:
    #         for gene in covered_positions:
    #             if gene in existing_genes:
    #                 if "tRNA" not in gene_name: 
    #                     line_pos = f"{chrom}\t{start}\t{stop}\t{gene}\t1000\t+\n"
    #                     line_neg = f"{chrom}\t{start}\t{stop}\t{gene}\t1000\t-\n"
    #                     f.write(line_pos); f.write(line_neg)
    #                     existing_genes.add(gene_name)

###########################################################################
# 3. Main, define accepted arguments. 

def main():
    parser = argparse.ArgumentParser(
        prog="bamget.py",
        formatter_class=argparse.RawTextHelpFormatter,
        description=textwrap.dedent("""\
###########################################################################

If used as part of the automated CRSSANT/rna2d3d pipeline, this script must 
be run in the parent directory where *_pri_crssant.bam is located after 
mapping is completed. Use the 'gene' analysis function. It will process the 
BAM file containing gapped_1 and trans reads, determine highly covered 
positions based on coverage per gene region in the annotation file. Outputs 
are written to the BED6 file, which can be directly used for DG assembly.

NOTE: Arguments should only be provided in the following order:

1. analysis_type:       Available options: 'gene', 'position'

                        Gene region-based uses mean coverage of reads for 
                        genomic coordinate regions to determine which genes 
                        are covered above a minimum coverage, writes 'genes' 
                        into a BED6 file for CRSSANT DG assembly. Requires
                        annotation file to be specified (see #4 [-a]).
                        
                        Position-based uses read coverage at nucleotide
                        positions of n=<window_size> to determine regions of 
                        high coverage and plots a histogram for read coverage
                        profile.  

                        See optional arguments for changing default values.

2. input (file path):   Intended as PATH of the *pri_crssant SAM or BAM 
                        file generated after the mapping.sh script from the 
                        rna2d3d pipeline is used. However, any SAM or BAM
                        file can be provided to check for read depth.

3. output (prefix):     Any string to be used for the output file prefix. 
                        Recommended to include min_coverage as part of the 
                        name to specify cutoff values. Supplying the same 
                        output filename will append to existing lines. 

4. -a, --annotation:    Annotation file (path) containing gene regions in
                        modified GTF format:

                        chrom   chrom_start   chrom_end   gene_name   strand

                        The modified GTF can be generated using a gencode 
                        annotation file and the following commands:

                        zcat gencode.v45.basic.annotation.gtf.gz
                        awk -F'\\t '$3 == "gene" && $9 ~ /gene_type 
                        "protein_coding"/ && $9 !~ /gene_name "ENSG0000"/ 
                        {split($9, a, "\\""); print $1 "\\t" $4 "\\t" $5 "\\
                        t" a[6] "\\t" $7}' > annotation.bed

                        If providing a pre-made annotation file, the format 
                        must follow the tab-separated columns described above.

                        NOTE: Specific genomic regions or RNAs are defined as 
                        a separate 'chromosome' in the PARIS/SHARC pipelines. 
                        These must manually be added as separate lines to the 
                        annotation_file. Always check to see if they are in 
                        the output BED6 files generated using this script.

                        List of genes added as a separate 'chromosome':
                        ---------------------------------
                        RN7SK   1       331     RN7SK   +
                        RN7SL   1       288     RN7SL   +
                        RNU7    1       63      RNU7    +
                        RNY     1       112     RNY1    +
                        RNY     213     314     RNY2    +
                        RNY     415     510     RNY3    +
                        RNY     611     694     RNY4    +
                        U13     1       120     U13     +
                        U14AB   1       92      U14A    +
                        U14AB   193     283     U14B    +
                        U17     1       207     U17     +
                        U3      1       217     U3      +
                        U8      1       136     U8      +
                        hs12S   1       954     12S     +
                        hs16S   1       1559    16S     +
                        hs45S   3654    5523    18S     +
                        hs45S   6600    6757    5.8S    +
                        hs45S   7924    12994   28S     +
                        hs5S    1       121     hs5S    +
                        hssnRNA 1       164     U1      +
                        hssnRNA 265     451     U2      +
                        hssnRNA 552     696     U4      +
                        hssnRNA 797     902     U6      +
                        hssnRNA 1003    1118    U5      +
                        hssnRNA 1219    1353    U11     +
                        hssnRNA 1454    1603    U12     +
                        hssnRNA 1704    1833    U4atac  +
                        hssnRNA 1934    2058    U6atac  +
                        ---------------------------------

OPTIONAL ARGUMENTS: SPECIFIC TO GENE-BASED ANALYSIS

-min, --min_coverage:   Defaults to [10] if no value is provided.

                        Integer value to set cutoff of minimum mean number of 
                        reads for a given genomic region. This value should
                        vary based on the dataset. Start with a larger value
                        and decrease to include genes with lower coverage.

OPTIONAL ARGUMENTS: SPECIFIC TO POSITION-BASED ANALYSIS

-w, --window-size:      Defaults to [20] if no value is provided.
                        
                        Window size to check for positions from input file.  
                        Start with a larger value to minimize compute time.

OPTIONAL ARGUMENTS: COMMON TO EITHER ANALYSIS

-max, --max-reads:      Defaults to [10,000] if no value is provided.

                        Parses input BAM file and sets a cutoff for 
                        number of rRNA reads allowed for coverage analysis. 

                        Recommended for analysis that contain mostly rRNA
                        reads, as this will decrease DG assembly time. Note 
                        that the 'annotation' flag [-a] is a paired argument. 

-k, --keep-files:       Retains the mean_coverage analysis files; 
                        
                            *.regions.bed.gz*
                            annotation.bed

                        WARNING: Error may occur if using the same file 
                        between analysis types. By default, -k is FALSE.
                        Provide this argument only when repeating the same
                        analysis type and changing values for min_coverage 
                        or skip_chromosome.

-s, --skip-chromosome:  Provide a list of comma separated chromosomes to 
                        omit when performing the data analysis. Highly 
                        recommended for analyses that do not need rRNA reads,
                        as this will significantly decrease analysis time. 

                            -skip=chr1,chr2,chr3  
                        
-stats, --statistics:   Provides basic statistics of the job as a log file;
                        
                            Initial arguments passed to the script
                            Job timing & memory used
                            Read coverage per gene
 
###########################################################################
"""),
usage='''
\npython3 %(prog)s analysis_type input output [options]

Usage examples:
    
- For gene region-based analysis: 
    %(prog)s gene input output -a=annotation_file \
[-min, -max, -skip, -stats, -k]\n

- For position-based analysis: 
    %(prog)s position input output \
[-w, -max, -skip, -stats, -k]\n
''')

    # Arguments to be provided in the command line.
    parser.add_argument('analysis_type', 
        help="Must be 'gene' or 'position'. See help text for more details.")
    parser.add_argument('input', 
        help="Input file in unsorted BAM or SAM format.")
    parser.add_argument('output', 
        help="File prefix name for analysis outputs.")
    # Flags, specific to the gene-based analysis.
    parser.add_argument('-a', '--annotation', 
        help="Annotation file in the modified GTF format. See help text for"
        " more details.")
    parser.add_argument('-min', '--min_coverage', 
        help="Min_coverage is defined as minimum number of reads for a region"
        " before it is considered well covered. Defaults to [10]")
    # Optional flags, specific to the position-based analysis.
    parser.add_argument('-w', '--window-size',
        help="Optional parameter for position-based analysis."
        " Defaults to [20]")
    # Optional flags, common to either analysis type. Set default values.
    parser.add_argument('-k', '--keep-files', action='store_true', 
        help="Optional parameter to keep the intermediary files.")
    parser.add_argument('-max', '--max-reads', 
        help="Optional parameter to define the maximum number of reads. If"
        " provided, subsamples highly covered gene regions.")
    parser.add_argument('-skip', '--skip-chromosome', nargs='?', 
        help="Optional parameter to skip chromosomes during processing.")
    parser.add_argument('-stats', '--statistics', action='store_true', 
        help="Optional parameter to give statistics for the run. Times the" 
        " processed for benchmarking in CRSSANT analysis.")
    parser.add_argument('-V', '--version', action='version', 
        version=f'%(prog)s {__version__}\n{__update_notes__}', 
        help='Print version + update notes and exit.')

    args = parser.parse_args()

    ##########################################################################

    # Begin timing the job for the --statistics flag. 
    start_time = time.time()

    # Check for optional arguments to modify functions.
    min_coverage = args.min_coverage # Gene-based: Default 10
    window_size = args.window_size # Position-based: Default 20
    keep_files = args.keep_files # Default: False
    max_reads = args.max_reads # Default: None
    skip_chromosome = args.skip_chromosome # Default: None
    statistics = args.statistics # Default: None

    # Check for mandatory positional arguments.
    analysis_type = args.analysis_type # Must be 'gene' or 'position'
    input_file = args.input # Must be in unsorted BAM or SAM format
    output = args.output # Any string, no file extension necessary
    annotation = args.annotation # Annotation file in modified GTF format

    # Set up the I/O and parse required arguments.
    if analysis_type != "gene" and analysis_type != "position": 
        print(f"Provided analysis type not supported."
            f" Use 'gene' or 'position'.")
        sys.exit()

    # Parse input file for analysis, set up the output.
    file_extension = os.path.splitext(input_file)[1].lower()
    if os.path.splitext(input_file)[1] != ".sam" and \
        os.path.splitext(input_file)[1] != ".bam":
        print(f"File type not supported."
            f" Make sure the file is in the SAM or BAM format.")
        sys.exit()
    else:
        sorted_bam = sort_bam(input_file, max_reads)    

    # Clean output prefix, strip file extensions accidentally provided.
    if '.bed' in output or '.pdf' in output:
        output = os.path.splitext(output)[0]

    # Check options for gene-based analysis.
    if analysis_type == "gene" or max_reads == True:
        # Check for gene regions based on the annotation file.
        gene_regions = collapse_gene_regions(annotation)

    if analysis_type == "gene":
        if min_coverage is not None and float(min_coverage) < 10:
            print("WARNING: Using a low min_coverage value increases" 
                " analysis time. CTRL + C to cancel now if desired...")
        try:
            time.sleep(5)
            print("Continuing...")
        except KeyboardInterrupt:
            print("Operation canceled by user.")
            sys.exit(1)
        covered_genes = mosdepth_regions(sorted_bam, 
            annotation, min_coverage, skip_chromosome)
        region_to_bed(covered_genes, gene_regions, output)

    # For coverage by nucleotide position, output to histograms.
    if analysis_type == "position":
        positions = mosdepth_positions(sorted_bam, window_size)
        position_to_pdf(histogram, output)

    # End the job here and provide the statistics.
    end_time = time.time()
    print(f"\nBAMget.py completed on {timenow()}.\n")

    # Set up the statistics 
    if statistics == True:
        elapsed_time = "{:.2f}".format(end_time - start_time)
        if min_coverage and float(min_coverage) < 10:
            elapsed_time = "{:.2f}".format(float(elapsed_time) - 5.00)
        print(f"     Elapsed time: {elapsed_time}s")
        if analysis_type == "gene":
            count = sum("tRNA" not in gene for gene in covered_genes)
            print(f"  Genes collected: {count}")
        elif analysis_type == "position":
            pass
            # Fill in a function here, show histograms or something useful.

    # Remove intermediate files unless -k argument is provided.
    if keep_files == True:
        print(f"WARNING: Keeping intermediate files between analysis types"
            f" may lead to incorrect coverage values. Remove and start over"
            f" if this was unintended.\n")
    else:
        os.remove('tmp.regions.bed.gz')
        os.remove('tmp.regions.bed.gz.csi')

    # Remove unnecessary output files from mosdepth.
    try:
        os.remove('tmp.mosdepth.global.dist.txt')
        os.remove('tmp.mosdepth.region.dist.txt')
        os.remove('tmp.mosdepth.summary.txt')
        os.remove('annotation.bed')
    
    except FileNotFoundError:
        pass #Ignore, if the file does not exist.

if __name__ == "__main__":
    main()
sys.exit()