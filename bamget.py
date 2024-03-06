"""
contact:    wlee9829@gmail.com
date:       2024_02_09
python:     python3.10
script:     bamget-1.0.py

This script is tool for preparing a BED file for CRSSANT analysis based on 
either the minimum coverage per annotated gene regions or nucleotide. Files
should be named according to the shell scripts in the rna2d3d repository. 
See the options associated with the analysis types.

Pedersen, B. S.; Quinlan, A. R.; Mosdepth: quick coverage calculation for
genomes and exomes, Bioinformatics, Volume 34, Issue 5, March 2018, 867–868

Danecek, P.; Bonfield, J. K.; Liddle, J.; Marshall, J.; Ohan, V.; Pollard, 
M. O.; Whitwham, A.; Keane, T.; McCarthy, S. A.; Davies, R. M.; Li, H.
GigaScience, Volume 10, Issue 2, February 2021
"""

# Import Packages
from datetime import datetime
import argparse
import glob
import gzip
import os
import subprocess
import sys
import textwrap
import time

###########################################################################
# 1. Basic functions that can be used for most analyses.

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
    modified GTF annotation file. This function was written with the hg38.bed 
    in mind, but will be updated for other gene annotation files in future
    iterations. See the help section for more details about preparing the 
    modified GTF file from a generic gencode GTF.

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
        ["chrX", "chrY", "chrM"])

    # ### FIX THIS SECTION. SKIPPING CHROMOSOMES IS STILL BROKEN.
    # if skip_chromosome != None:
    #     skip_chromosome = skip_chromosome.split(',')
    #     print(f"Skipping reads that are on chromosomes: {skip_chromosome}")

    # Define blacklisted regions as list of tuples (chromosome, start, stop).
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
    if os.path.exists("annotation.bed"):
        try:
            print(f"\nReading gene regions from 'annotations.bed'.\n")
            gene_regions = {}
            with open("annotation.bed", 'r') as f:
                for line in f:
                    parts = line.split("\t")
                    chromosome = parts[0].strip("'\"")
                    if skip_chromosome and chromosome in skip_chromosome:
                        continue
                    start = int(parts[1])
                    stop = int(parts[2])
                    gene_name = parts[3]
                    strand = parts[4]
                    gene_regions[gene_name] = (chromosome, start, stop, strand)

        except Exception as e:
            error_message = (f"\nError: Failed to parse annotation.bed. Check"
                f" to see if the columns are as defined in the help section.")
            print(error_message)
            sys.exit(error_message)

    else:
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

    return gene_regions

def sort_bam(input_file, annotation=None, max_reads=None, skip_chromosome=None):
    """
    Checks the input file type, converts to BAM if necessary, 
    then sorts the BAM file for input into the analysis functions.
    Requires samtools packages. 

    If the argument max_reads is passed, a maximum number of randomized reads per given gene region will be pulled from the input file. 
    Requires the bedtools package.

    Args:
        input_file: either a SAM or BAM alignment mapping file.

    Returns:
        file: sorted_bam file using samtools.
    """
    input_prefix = os.path.splitext(input_file)[0]

    if max_reads == None: 
        sam_to_bam = f"{input_prefix}.bam"
        sorted_bam = f"{input_prefix}_sorted.bam"
            
    else:

        # # WRITE THIS SO THAT IF MAX-READS IS CALLED, ANNOTATION IS NEEDED.
        # cutoff = int(max_reads)

        # # Check the output of the bedtools intersect and parse for cutoff.
        # subprocess.run([
        #     'bedtools', 'bamtobed', 
        #     '-i', input_file, 'bam_tmp.bed'], shell=True, check=True)

        # subprocess.run([
        #     'bedtools', 'intersect', 
        #     '-a', 'annotation.bed', 
        #     '-b', 'bam_tmp.bed', 
        #     '>', 'overlaps.bed'], check=True)

        # # Re-write the BAM file with randomized number of reads. numpy.random?
        # # samtools seed?
        # subprocess.run([
        #    'awk',
        #    '{count[$4]++}',
        #    'END {for (gene in count) print gene, count[gene]}',
        #    'overlaps.bed'], capture_output=True, text=True)

        # with open ('bedtools_output', 'r') as f:
        #     for line in f:
        #         parts = line.split("\t")
        #         chromosome = parts[0].strip("'\"")

        # sam_to_bam = f"{input_prefix}_max-{cutoff}.bam"
        # sorted_bam = f"{input_prefix}_max-{cutoff}_sorted.bam"

    try:
        if input_file.endswith(".sam"):
            subprocess.run(['samtools', 'view', '-bS', '-o', 
                sam_to_bam, input_file])
            subprocess.run(['samtools', 'sort', '-o', 
                sorted_bam, sam_to_bam])
            subprocess.run(['samtools', 'index', sorted_bam])

        elif input_file.endswith(".bam"):
            subprocess.run(['samtools', 'sort', '-o', 
                sorted_bam, input_file])
            subprocess.run(['samtools', 'index', sorted_bam])

    except Exception as e:
        error_message = f("\nError: Failed to parse the provided file input."
            f" Check to see if file format is correct. Details: {e}.")
        print(error_message)
        sys.exit(error_message)

    return sorted_bam

###########################################################################
# 2. Coverage depth calculations 

def mosdepth_regions(sorted_bam, annotation, min_coverage):
    """
    Function that parses the provided sorted_bam file to determine gene regions
    that meet or exceed minimum coverage as defined by user. Uses mosdepth 
    functions to determine reads at each gene region, then collects alignments 
    >= min_coverage. Defining num_cores should speed up the depth calculations, 
    but benchmarking for mosdepth shows no notable improvement past 4 threads.

    Args:
        sorted_bam:         PATH to BAM file.
        annotation:         bed file that defines the gene regions, derived 
                            from collapse_gene_regions() function.
        min_coverage:       Integer defining minimum number of reads/nt before 
                            the region is considered covered.

    Returns:
        covered_genes (list): genes with coverage that exceed min_coverage
    """
    covered_genes = []
    if os.path.exists('tmp.regions.bed.gz'):
        print(f"Mean coverage summary already exists. Using for calculations.")
        
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

    with gzip.open('tmp.regions.bed.gz', 'rt') as f:
        for line in f:
            chrom, start, end, gene, cov_mean = line.strip().split()
            if cov_mean >= min_coverage:
                covered_genes.append(gene)

    return covered_genes

# RE-WRITE THIS FUNCTION SO THAT IT PLOTS A HISTOGRAM WITH THE OUTPUT
# INSTEAD OF USING THE SAME FILE FOR GENES.BED

# def mosdepth_positions(sorted_bam, annotation, min_coverage, 
#     gene_regions, window_size=None):
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
    mode = 'a' if os.path.exists(output) else 'w'

    existing_genes = set()
    if mode =='a':
        with open(output, 'r') as f:
            for line in f:
                gene_name = line.split('\t')[3].strip()
                existing_genes.add(gene_name)

    with open(output, mode, newline='') as f:
        for gene_name in covered_genes:
            if gene_name not in existing_genes:
                if "tRNA" not in gene_name: 
                    if gene_name in gene_regions:
                        chrom, start, stop, strand = gene_regions[gene_name]
                        line = (f"{chrom}\t{start}\t{stop}\t "
                            f" {gene_name}\t1000\t{strand}")
                        f.write(line)
                        existing_genes.add(gene_name)


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
# 3. Main, defining accepted arguments. 

def main():
    parser = argparse.ArgumentParser(
        prog="bamget-1.0.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
###########################################################################

If used as part of the automated CRSSANT/rna2d3d pipeline, this script would
best be run in the parent directory where *_pri_crssant.bam is located after 
the mapping script completes. It will process the BAM file that contains 
gapped_1 and trans reads, determining highly covered positions based on either 
the read coverage per region defined in a user-provided annotation file, or 
positions that fall within a user-provided window-size. Outputs are written to 
a BED6 file, which can then be directly used for CRSSANT DG assembly. 

NOTE: Arguments should only be provided in the following order:

1. analysis_type:       Options are 'gene' OR 'position'. 

                        Gene region-based analysis uses the mean coverage of 
                        reads for a specific gene region to determine which
                        genes are highly covered and writes to the 'genes.bed'
                        used in CRSSANT DG assembly. 
                        
                        Position-based analysis uses the read coverage at 
                        nucleotides positions of n=<window_size> to determine
                        regions of high coverage and plots a histogram for read 
                        coverage profile. See options [-c, -w, -m]. 

2. input (file path):   Intended to be the PATH of the *pri_crssant SAM or BAM 
                        file generated after the mapping.sh script from the 
                        rna2d3d pipeline is used. However, any SAM or BAM file 
                        can be provided to check for read depth.

3. output (prefix):     Any string to be used for the output file prefix. It is 
                        recommended to include min_coverage as part of the 
                        name to specify the cutoff values. Supplying the 
                        same output filename will append to existing lines. 

4. min_coverage (int):  Integer value defines the minimum number of reads for
                        each analysis_type. Default value of min_coverage is 
                        (10) with default window_size (1), unless option [-w]
                        is specified. 

                        This value is recommended to be changed according to 
                        the dataset. Start with a larger value and decrease to 
                        include genes with lower coverage.

                        For 'gene' based analysis, it defines the minimum mean 
                        number of reads for a genomic region before being
                        considered high coverage.

                        For 'position' based analysis, it sets the cutoff for
                        minimum number of reads at any given nucleotide 
                        position within the window (1 unless [-w]) before 
                        being considered a highly covered position.

4. annotation (file):   Annotation file (path) containing gene regions in a 
                        modified GTF format:

                        chrom   chrom_start   chrom_end   gene_name   strand

                        The modified GTF can be generated using a gencode 
                        annotation file and the following commands:

                        zcat gencode.v45.basic.annotation.gtf.gz
                        awk -F'\\t '$3 == "gene" && $9 ~ /gene_type 
                        "protein_coding"/ && $9 !~ /gene_name "ENSG0000"/ 
                        {split($9, a, "\\""); print $1 "\\t" $4 "\\t" $5 "\\t" 
                        a[6] "\\t" $7}' > annotation.bed

                        If providing a pre-made annotation file, the format 
                        must be correct and the file re-named "annotation.bed".
                        If working with PARIS/SHARC pipeline data, see below.

                        NOTE: Specific genomic regions or RNAs are defined as a
                        separate 'chromosome' in the PARIS/SHARC pipelines. 
                        These must be manually included as separate lines in
                        the annotation_file. Check to see if present in the
                        output BED6 files generated using this script.

                        Example 'genes' as a separate 'chromosome':

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

OPTIONAL ARGUMENTS: SPECIFIC TO POSITION-BASED ANALYSIS

-w, --window-size:      Window size of the positions from the input file.  
                        Start with a larger value to minimize compute time.
                        Defaults to 20 if provided with no value. 

OPTIONAL ARGUMENTS: COMMON

-m, --max-reads:        Parses the provided input BAM file and collects a 
                        maximum number of reads before the coverage analysis. 
                        Recommended for datasets that contain rRNA reads, as
                        this will decrease DG assembly time. Defaults to 1000
                        if provided with no value.

-s, --skip-chromosome:  Provide a list of comma separated chromosomes to 
                        omit when performing the data analysis. Recommended 
                        for analyses that do not need rRNA reads, as this will 
                        significantly decrease DG assembly time. 

                            e.g., chr1,chr2,chr3  
                        
-S, --statistics:       Provides basic statistics of the job as a log file;
                        
                            Initial arguments passed to the script
                            Job timing & memory used
                            Read coverage per gene

-k, --keep-files:       Retains the mean_coverage analysis intermediary files; 
                        
                            *.regions.bed.gz*
                            annotation.bed

                        WARNING: Error may occur if using the same intermediary 
                        file between analysis types. By default, -k is FALSE.
                        Provide this argument only when performing the same
                        analysis type and changing the min_coverage value.
 
###########################################################################
"""),
usage='''
\npython3 %(prog)s analysis_type [options] input_file min_coverage \
    annotation_file output 

Usage examples:
    
- For a gene region-based analysis: 
    %(prog)s gene input_file min_coverage annotation_file output \
[-k, -m, -s, -S]\n

- For a position-based analysis: 
    %(prog)s position input_file min_coverage annotation_file output \
[-w, -k, -m, -s, -S]\n
''')
    # Arguments for the command line.
    parser.add_argument('analysis_type', 
        help="'gene' or 'position' options. See help text for more details.")
    
    parser.add_argument('input_file', 
        help="Input file in either a BAM or SAM format.")
        
    parser.add_argument('output', help="File prefix name for analysis outputs.")

    # Flags, specific to the gene-based analysis.

    parser.add_argument('min_coverage', 
        help="Min_coverage is defined as the minimum number of reads for either"
        " a given region or position before it is considered highly" 
        " covered.")

    parser.add_argument('annotation_file', 
        help="Annotation file in the modified GTF format. See help text for "
        " more details.")


    # Optional flags, specific to the position-based analysis.
    parser.add_argument('-w', '--window-size',
        help="Optional parameter for position-based analysis. Defaults to 20" 
        " nt window size.")

    # Optional flags, common to either analysis type. Set default values.
    parser.add_argument('-m', '--max-reads', 
        help="Optional parameter to define the maximum number of reads. If "
        " provided, subsamples highly covered gene regions.")

    parser.add_argument('-s', '--skip-chromosome', nargs='?', 
        help="Optional parameter to skip specified chromosomes during "
        " processing.")

    parser.add_argument('-S', '--statistics', action='store_true', 
        help="Optional parameter to give statistics for the run. Times the" 
        " processed for benchmarking in CRSSANT analysis.")

    parser.add_argument('-k', '--keep-files', action='store_true', 
        help="Optional parameter to keep the intermediary files.")

    args = parser.parse_args()

    ###########################################################################
    # Parse arguments and perform sub-functions.

    # Time the job for the --statistics flag. 
    start_time = time.time()

    # Check for gene regions based on the annotation file.
    gene_regions = collapse_gene_regions(args.annotation_file)    

    # Check positional arguments, existing files, strip file extensions. 
    analysis_type = args.analysis_type
    if analysis_type != "gene" and analysis_type != "position": 
        print(f"Provided analysis type not supported."
            f" Use 'gene' or 'position'.")
        sys.exit()

    input_file = args.input_file
    min_coverage = args.min_coverage

    if min_coverage <= 10:
        print(f"Warning, using a low min_coverage value takes significantly 
            f" longer for analysis.")

    output = args.output
    if '.bed' in output:
        output = os.path.splitext(output)[0]

    # Check for optional arguments. 
    window_size = args.window_size
    max_reads = args.max_reads
    skip_chromosome = args.skip_chromosome


    file_extension = os.path.splitext(input_file)[1].lower()
    if os.path.splitext(input_file)[1] != ".sam" and \
        os.path.splitext(input_file)[1] != ".bam":
        print(f"File type not supported."
            f" Make sure the file is in the SAM or BAM format.")
        sys.exit()

    # Parse input file for analysis, set up the output.
    sorted_bam = sort_bam(input_file, 'annotation.bed', 
        max_reads, skip_chromosome)    

    # For coverage per gene region, output to bed.
    if args.analysis_type == "gene":
        covered_genes = mosdepth_regions(sorted_bam, 
            'annotation.bed', min_coverage)
        region_to_bed(covered_genes, gene_regions, output)

    # For coverage by nucleotide position, output to histograms.
    elif analysis_type == "position":
        covered_positions = mosdepth_positions(sorted_bam, 'annotation.bed',
            min_coverage, gene_regions, window_size)
        position_to_bed(covered_positions, output)

    end_time = time.time()
    print(f"\nJob completed at {timenow()}.\n")

    if args.statistics:
        elapsed_time = "{:.2f}".format(end_time - start_time)
        print(f"     Elapsed time: {elapsed_time}s")

        if analysis_type = "gene":
            with open(output, 'rt') as f:
                count = len(f.readlines())
            print(f"  Genes collected: {count}")
        elif analysis_type = "position":
            # Fill in a function here, show num_reads or something else useful.

    # Remove the intermediate files unless -k argument is provided.
    if args.keep_files:
        print(f"\nWARNING: Keeping intermediate files between analysis types"
            f" may lead to incorrect coverage values. Remove and start over"
            f" if this was unintended.")
    else:
        os.remove('tmp.regions.bed.gz')
        os.remove('tmp.regions.bed.gz.csi')

    # Remove unnecessary output files from mosdepth.
    if os.path.exists('tmp.mosdepth.global.dist.txt'):
        os.remove('tmp.mosdepth.global.dist.txt')
    if os.path.exists('tmp.mosdepth.region.dist.txt'):
        os.remove('tmp.mosdepth.region.dist.txt')
    if os.path.exists('tmp.mosdepth.summary.txt'):
        os.remove('tmp.mosdepth.summary.txt')

if __name__ == "__main__":
    main()
sys.exit()