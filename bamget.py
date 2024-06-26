"""
contact:    wlee9829@gmail.com
date:       2024_02_09
python:     python3.10
script:     bamget.py

This script is multifunctional tool for preparing file in the CRSSANT 
analysis pipeline. Main functions generate output based on either the minimum 
coverage per annotated gene regions or nucleotide position. Files should be
named according to the shell scripts in the rna2d3d repository. 
See the options associated with the analysis types.

Relevant publications:

Pedersen, B. S.; Quinlan, A. R.; Mosdepth: quick coverage calculation for
genomes and exomes, Bioinformatics, Volume 34, Issue 5, March 2018, 867–868

Danecek, P.; Bonfield, J. K.; Liddle, J.; Marshall, J.; Ohan, V.; Pollard, 
M. O.; Whitwham, A.; Keane, T.; McCarthy, S. A.; Davies, R. M.; Li, H.
GigaScience, Volume 10, Issue 2, February 2021

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification
of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019).
"""
# Define version
__version__ = "2.1.1"

# Version notes
__update_notes__ = """
2.1.1
    -   Added new optional flag to specify which chromosome to plot in the 
        'position' based analysis.
    -   Updated docstrings to reflect new options.

2.1.0
    -   Added support for multi-threaded processing. 
    -   Updated docstrings and formatting.

2.0.0
    -   Added support for position-based analysis function, output modified
        BAM file based on cutoff and window size.

1.1.1
    -   Updated wording and formatting for code readability. 

1.1.0
    -   hg38_blacklist.v2 added to excluded regions.

1.0.0
    -   Gene-region based analysis function and output to BED6 file.

To Do:
    -   Add histogram plotting feature
    -   Add support for non homo-sapiens (HS) genome files.
"""

# Import Packages
from datetime import datetime
import argparse
import glob
import gzip
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import multiprocessing
import os
import pandas as pd
import random
import shutil
import subprocess
import sys
import textwrap
import time

###########################################################################
# 1. Basic functions common to either analyses function.

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
    # unordinarily high signal in multiple NGS experiments. Uses Boyle lab's
    # hg38-blacklist.v2.bed.gz from the GitHub page.

    # Define blacklisted signal artifact regions in the human genome
    # as a list tuples (chromosome, start, stop):
    blacklist_regions = [
    ("chr1", 0, 792500), ("chr1", 91386300, 91388400),
    ("chr1", 103594400, 103760600), ("chr1", 121605200, 124938900),
    ("chr1", 125067600, 125086000), ("chr1", 125130200, 143562200),
    ("chr1", 161423100, 161472400), ("chr1", 168348600, 168349900),
    ("chr1", 224010800, 224017000), ("chr1", 236713000, 236715600),
    ("chr1", 248932700, 248956400), ("chr2", 1221700, 1223900),
    ("chr2", 1594700, 1605200), ("chr2", 3179600, 3182100),
    ("chr2", 4643800, 4648800), ("chr2", 10952800, 10955000),
    ("chr2", 13718700, 13737700), ("chr2", 21903500, 21906400),
    ("chr2", 32865900, 32869900), ("chr2", 32915300, 32918400),
    ("chr2", 33766500, 33768400), ("chr2", 36183000, 36184500),
    ("chr2", 49228700, 49230700), ("chr2", 64359300, 64377000),
    ("chr2", 86655300, 86661100), ("chr2", 86900700, 87078100),
    ("chr2", 87119300, 87189800), ("chr2", 87217000, 87866200),
    ("chr2", 88771000, 88806500), ("chr2", 89235300, 89947100),
    ("chr2", 90246300, 91735500), ("chr2", 91783000, 91924800),
    ("chr2", 91969000, 94569500), ("chr2", 95849400, 96067900),
    ("chr2", 97106300, 97615800), ("chr2", 109198400, 109200700),
    ("chr2", 109744600, 110095200), ("chr2", 110229200, 110633400),
    ("chr2", 111253600, 111500500), ("chr2", 112346200, 112441300),
    ("chr2", 113370100, 113662700), ("chr2", 130496800, 130716400),
    ("chr2", 132201000, 132288900), ("chr2", 132353600, 132364500),
    ("chr2", 148880800, 148882800), ("chr2", 161277700, 161283400),
    ("chr2", 181274800, 181276800), ("chr2", 226108500, 226110400),
    ("chr2", 234889800, 234894400), ("chr2", 239642200, 239645600),
    ("chr2", 240308100, 240310300), ("chr2", 241589300, 241591800),
    ("chr2", 242005900, 242011100), ("chr2", 242110100, 242193500),
    ("chr3", 3895200, 3896700), ("chr3", 4916700, 4922500),
    ("chr3", 14091000, 14092500), ("chr3", 15187200, 15207800),
    ("chr3", 15592100, 15603300), ("chr3", 16176800, 16179200),
    ("chr3", 16679700, 16682500), ("chr3", 19499700, 19504000),
    ("chr3", 19624000, 19627100), ("chr3", 21983200, 21988100),
    ("chr3", 24053500, 24054900), ("chr3", 26384800, 26404100),
    ("chr3", 29993900, 29999900), ("chr3", 36987500, 36995000),
    ("chr3", 38083400, 38085400), ("chr3", 38406100, 38430900),
    ("chr3", 39366700, 39386000), ("chr3", 40219400, 40240500),
    ("chr3", 49671000, 49696700), ("chr3", 51457800, 51462000),
    ("chr3", 57326800, 57328500), ("chr3", 65124100, 65126100),
    ("chr3", 65510000, 65513900), ("chr3", 65697400, 65699300),
    ("chr3", 66273800, 66275200), ("chr3", 68076400, 68077800),
    ("chr3", 69047300, 69053600), ("chr3", 69475300, 69479700),
    ("chr3", 75630100, 75707800), ("chr3", 75736400, 75754600),
    ("chr3", 78948800, 78950500), ("chr3", 80876000, 80894000),
    ("chr3", 89345600, 89370500), ("chr3", 90156400, 90175500),
    ("chr3", 90455400, 91297100), ("chr3", 91516200, 93749200),
    ("chr3", 96616300, 96619300), ("chr3", 97905100, 97923200),
    ("chr3", 101674800, 101698400), ("chr3", 103224300, 103236500),
    ("chr3", 106665700, 106669700), ("chr3", 106975900, 106979600),
    ("chr3", 108751100, 108755100), ("chr3", 111019500, 111024600),
    ("chr3", 121933800, 121936400), ("chr3", 122414300, 122417500),
    ("chr3", 122735500, 122796600), ("chr3", 122837000, 122838700),
    ("chr3", 75630100, 75707800), ("chr3", 75736400, 75754600),
    ("chr3", 78948800, 78950500), ("chr3", 80876000, 80894000),
    ("chr3", 89345600, 89370500), ("chr3", 90156400, 90175500),
    ("chr3", 90455400, 91297100), ("chr3", 91516200, 93749200),
    ("chr3", 96616300, 96619300), ("chr3", 97905100, 97923200),
    ("chr3", 101674800, 101698400), ("chr3", 103224300, 103236500),
    ("chr3", 106665700, 106669700), ("chr3", 106975900, 106979600),
    ("chr3", 108751100, 108755100), ("chr3", 111019500, 111024600),
    ("chr3", 121933800, 121936400), ("chr3", 122414300, 122417500),
    ("chr3", 122735500, 122796600), ("chr3", 122837000, 122838700),
    ("chr3", 133177100, 133179800), ("chr3", 133551500, 133579500),
    ("chr3", 135437200, 135439100), ("chr3", 136954600, 136969200),
    ("chr3", 137168400, 137169900), ("chr3", 138575800, 138595900),
    ("chr3", 139190800, 139194700), ("chr3", 153236200, 153241300),
    ("chr3", 155544100, 155546700), ("chr3", 156279000, 156283500),
    ("chr3", 157080800, 157093400), ("chr3", 158511300, 158513100),
    ("chr3", 160941200, 160948700), ("chr3", 161001900, 161014100),
    ("chr3", 165573100, 165591000), ("chr3", 166228200, 166232400),
    ("chr3", 168012100, 168016800), ("chr3", 170567000, 170569900),
    ("chr3", 170864300, 170881400), ("chr3", 171626600, 171637700),
    ("chr3", 174829200, 174831800), ("chr3", 176828700, 176833000),
    ("chr3", 177660600, 177664000), ("chr3", 178926800, 178941300),
    ("chr3", 183016900, 183019100), ("chr3", 183955400, 183958700),
    ("chr3", 187893900, 187896100), ("chr3", 192739300, 192742700),
    ("chr3", 194323600, 194334900), ("chr3", 195477900, 195507300),
    ("chr3", 195616000, 195750100), ("chr3", 195775500, 195791400),
    ("chr3", 195914100, 196028300), ("chr3", 196249400, 196251900),
    ("chr3", 196897800, 196899800), ("chr3", 197030600, 197035800),
    ("chr3", 197383400, 197428800), ("chr3", 197454700, 197460800),
    ("chr3", 197598400, 197680900), ("chr3", 198099800, 198295500),
    ("chr4", 0, 69200), ("chr4", 554100, 556500),
    ("chr4", 1427000, 1468900), ("chr4", 6002700, 6005700),
    ("chr4", 7863000, 7865000), ("chr4", 9212700, 9369600),
    ("chr4", 40291700, 40318200), ("chr4", 49077200, 51816100),
    ("chr4", 55327200, 55329200), ("chr4", 77994000, 78009600),
    ("chr4", 119274400, 119301700), ("chr4", 146285100, 146305300),
    ("chr4", 162420500, 162422400), ("chr4", 166554300, 166581300),
    ("chr4", 181238800, 181242300), ("chr4", 189232500, 189236300),
    ("chr4", 189834900, 189849700), ("chr4", 189877500, 190023700),
    ("chr4", 190048600, 190214500), ("chr5", 0, 44100),
    ("chr5", 548300, 564100), ("chr5", 647600, 651700),
    ("chr5", 1326100, 1334600), ("chr5", 2144600, 2147800),
    ("chr5", 2489800, 2491700), ("chr5", 3322100, 3325100),
    ("chr5", 6967700, 6971700), ("chr5", 17516800, 17600200),
    ("chr5", 21477600, 21497600), ("chr5", 25381400, 25384300),
    ("chr5", 34177900, 34244800), ("chr5", 45522900, 45525200),
    ("chr5", 45743000, 45744800), ("chr5", 46433900, 46687700),
    ("chr5", 46708100, 50165300), ("chr5", 60759700, 60762500),
    ("chr5", 63320900, 63335500), ("chr5", 69540700, 71359500),
    ("chr5", 71850000, 71852800), ("chr5", 74685400, 74712400),
    ("chr5", 78452400, 78457600), ("chr5", 78848400, 78872800),
    ("chr5", 80649100, 80653100), ("chr5", 85641800, 85662700),
    ("chr5", 93947500, 93948900), ("chr5", 94567100, 94570400),
    ("chr5", 100045500, 100076300), ("chr5", 106425500, 106429500),
    ("chr5", 109259500, 109265400), ("chr5", 111302100, 111308300),
    ("chr5", 114156700, 114158300), ("chr5", 119904000, 119905600),
    ("chr5", 123760300, 123762200), ("chr5", 134922500, 134929400),
    ("chr5", 139005500, 139011600), ("chr5", 146610000, 146615500),
    ("chr5", 153071100, 153077000), ("chr5", 156658300, 156665400),
    ("chr5", 161606000, 161611700), ("chr5", 171083900, 171090200),
    ("chr5", 175904500, 176118000), ("chr5", 176590700, 176593000),
    ("chr5", 177636700, 177684700), ("chr5", 177960500, 177981400),
    ("chr5", 178584600, 178586600), ("chr5", 181172600, 181538200),
    ("chr6", 256500, 382800), ("chr6", 861500, 864200),
    ("chr6", 1052800, 1054800), ("chr6", 26669200, 26832300),
    ("chr6", 33484600, 33486400), ("chr6", 34070600, 34074000),
    ("chr6", 38262000, 38301600), ("chr6", 39455800, 39460500),
    ("chr6", 44043600, 44080000), ("chr6", 44180600, 44182900),
    ("chr6", 51874900, 51901300), ("chr6", 54961900, 54967400),
    ("chr6", 58432200, 60242300), ("chr6", 61321800, 61493000),
    ("chr6", 61573200, 61575100), ("chr6", 61661900, 61673400),
    ("chr6", 103709000, 103715100), ("chr6", 115254900, 115256800),
    ("chr6", 143799900, 143801800), ("chr6", 156035300, 156040100),
    ("chr6", 157309500, 157324800), ("chr6", 160611700, 160647400),
    ("chr6", 170145300, 170147200), ("chr6", 170376900, 170401000),
    ("chr6", 170465400, 170468400), ("chr7", 0, 49600), 
    ("chr7", 224500, 241300), ("chr7", 904700, 907100), 
    ("chr7", 1271400, 1273500), ("chr7", 45251000, 45253000), 
    ("chr7", 56369500, 56375600), ("chr7", 57485300, 57497800), 
    ("chr7", 57611600, 57637700), ("chr7", 58031800, 60997400), 
    ("chr7", 61017800, 61075200), ("chr7", 61102900, 61725200), 
    ("chr7", 62265700, 62409500), ("chr7", 62430000, 62520600), 
    ("chr7", 65488000, 65496500), ("chr7", 100951400, 100968300), 
    ("chr7", 100991500, 101004600), ("chr7", 102474700, 102686400), 
    ("chr7", 142665100, 142668500), ("chr7", 144180800, 144377300), 
    ("chr7", 145996400, 146018600), ("chr7", 152375800, 152435100),
    ("chr7", 158131400, 158156200), ("chr7", 158594300, 158596200), 
    ("chr7", 158893100, 158918100), ("chr7", 159334900, 159345900),
    ("chr8", 7209800, 7914700), ("chr8", 7940500, 8075700),
    ("chr8", 8128200, 8204600), ("chr8", 12136900, 12614300),
    ("chr8", 43236700, 43262600), ("chr8", 43937900, 45969600),
    ("chr8", 46829400, 46832000), ("chr8", 57204900, 57216100),
    ("chr8", 59168700, 59170400), ("chr8", 67584500, 67592700),
    ("chr8", 69688400, 69691100), ("chr8", 71406700, 71412400),
    ("chr8", 75444100, 75448200), ("chr8", 81841500, 81851900),
    ("chr8", 85642100, 85829300), ("chr8", 88685900, 88691700),
    ("chr8", 96171200, 96173100), ("chr8", 99494900, 99496800),
    ("chr8", 105789200, 105793800), ("chr8", 141491400, 141493500),
    ("chr8", 141871100, 141875200), ("chr8", 143641400, 143670500),
    ("chr8", 144124800, 144137600), ("chr9", 319900, 322400), 
    ("chr9", 33656600, 33660000), ("chr9", 35912600, 35915300), 
    ("chr9", 38824200, 39089400), ("chr9", 39846200, 40771100), 
    ("chr9", 40792500, 41323100), ("chr9", 41492300, 41635600), 
    ("chr9", 41661300, 42119600), ("chr9", 42364000, 42410600), 
    ("chr9", 42899400, 42901300), ("chr9", 43263100, 61518900), 
    ("chr9", 61735300, 63548000), ("chr9", 63761400, 64027300), 
    ("chr9", 64135000, 65390600), ("chr9", 65579400, 66874600), 
    ("chr9", 66959000, 68398100), ("chr9", 70037200, 70039600), 
    ("chr9", 76174300, 76176200), ("chr9", 83222900, 83226900), 
    ("chr9", 85071600, 85075100), ("chr9", 85164800, 85166100), 
    ("chr9", 108502000, 108506600), ("chr9", 134164500, 134185500), 
    ("chr9", 137326800, 137330600), ("chr9", 137715200, 137722200),
    ("chr9", 137841200, 137846800), ("chr9", 138222000, 138394700),
    ("chr10", 0, 45700), ("chr10", 38481300, 38596500),
    ("chr10", 38782600, 38967900), ("chr10", 39901300, 41712900),
    ("chr10", 41838900, 42107300), ("chr10", 42279400, 42322500),
    ("chr10", 126946300, 126953400), ("chr10", 133625800, 133797400),
    ("chr11", 0, 194500), ("chr11", 518900, 520700),
    ("chr11", 584400, 586500), ("chr11", 964100, 966000),
    ("chr11", 1015700, 1019300), ("chr11", 1091000, 1098200),
    ("chr11", 3652800, 3655600), ("chr11", 10506900, 10511100),
    ("chr11", 28206300, 28236700), ("chr11", 50813600, 54383000),
    ("chr11", 61084500, 61130400), ("chr11", 70370400, 70372400),
    ("chr11", 73509800, 73511700), ("chr11", 77885600, 77887600),
    ("chr11", 93417500, 93427700), ("chr11", 94232700, 94240400),
    ("chr11", 103408700, 103410600), ("chr11", 121175000, 121187000),
    ("chr11", 131679500, 131681500), ("chr11", 135075600, 135086600),
    ("chr12", 0, 77800), ("chr12", 371800, 422400),
    ("chr12", 2254900, 2257000), ("chr12", 2519800, 2540500),
    ("chr12", 5928900, 5933000), ("chr12", 20550500, 20552400),
    ("chr12", 20768400, 20770300), ("chr12", 29790400, 29834600),
    ("chr12", 34715400, 37269100), ("chr12", 41362700, 41364600),
    ("chr12", 61471100, 61473000), ("chr12", 66473900, 66475800),
    ("chr12", 101147000, 101155000), ("chr12", 113079600, 113081500),
    ("chr12", 124430500, 124440300), ("chr12", 124905900, 124941800),
    ("chr12", 130386400, 130394100), ("chr12", 131475300, 131478600),
    ("chr12", 131576000, 131589700), ("chr12", 132223300, 132243400),
    ("chr12", 132455100, 132465200), ("chr12", 133249000, 133275300),
    ("chr13", 16087600, 16165300), ("chr13", 16226300, 18171400),
    ("chr13", 18211000, 18216100), ("chr13", 57140500, 57172500),
    ("chr13", 109423200, 109425200), ("chr13", 114353300, 114364300),
    ("chr14", 0, 18670900), ("chr14", 18695400, 19724300),
    ("chr14", 23033300, 23098600), ("chr14", 26629300, 26634900),
    ("chr14", 31793800, 31798100), ("chr14", 32483400, 32486000),
    ("chr14", 34537100, 34562600), ("chr14", 35947200, 35950000),
    ("chr14", 37351000, 37356700), ("chr14", 44025100, 44027200),
    ("chr14", 44705100, 44709900), ("chr14", 45477100, 45482500),
    ("chr14", 46865300, 46866500), ("chr14", 54235600, 54240000),
    ("chr14", 57112100, 57118100), ("chr14", 74711700, 74729000),
    ("chr14", 86074000, 86076000), ("chr14", 86593300, 86595200),
    ("chr14", 88443700, 88458100), ("chr14", 100525900, 100527800),
    ("chr14", 101267600, 101272200), ("chr14", 101674400, 101676400),
    ("chr14", 104288100, 104290200), ("chr14", 105215000, 105240900),
    ("chr14", 105568500, 105583900), ("chr14", 105616500, 105618600),
    ("chr14", 106326900, 106367700), ("chr15", 0, 17035000),
    ("chr15", 17058500, 19790100), ("chr15", 20005600, 22606300),
    ("chr15", 23125400, 23357400), ("chr15", 25757700, 25759100),
    ("chr15", 28304900, 28683400), ("chr15", 30066300, 30627500),
    ("chr15", 30844100, 30859900), ("chr15", 32153700, 32626200),
    ("chr15", 54925700, 54932200), ("chr15", 56311200, 56314600),
    ("chr15", 72635200, 72687100), ("chr15", 74068100, 74102000),
    ("chr15", 75254100, 75299800), ("chr15", 77698600, 77700600),
    ("chr15", 82321000, 82374600), ("chr15", 82421200, 82541700),
    ("chr15", 84405300, 84524700), ("chr15", 101752300, 101764800),
    ("chr15", 101892700, 101991100), ("chr16", 29430800, 29566900),
    ("chr16", 34061400, 34121400), ("chr16", 34272000, 34633100),
    ("chr16", 34657200, 34672500), ("chr16", 34694600, 34772000),
    ("chr16", 34832600, 34922100), ("chr16", 34945600, 35072500),
    ("chr16", 36166300, 36202400), ("chr16", 36225200, 46423000),
    ("chr16", 46449700, 46467000), ("chr16", 90100500, 90338300),
    ("chr17", 0, 137600), ("chr17", 294900, 317900),
    ("chr17", 448200, 510900), ("chr17", 1061500, 1066100),
    ("chr17", 1307700, 1312000), ("chr17", 19025700, 19237400),
    ("chr17", 21783300, 22054000), ("chr17", 22520400, 22527300),
    ("chr17", 22745200, 26629800), ("chr17", 26766800, 26987200),
    ("chr17", 43227600, 43324300), ("chr17", 45511500, 45641300),
    ("chr17", 53104900, 53107300), ("chr18", 0, 64600),
    ("chr18", 105200, 113200), ("chr18", 971000, 976500),
    ("chr18", 2841300, 2861500), ("chr18", 15367200, 20940300),
    ("chr18", 46961600, 47031700), ("chr18", 47852300, 47854300),
    ("chr18", 52791800, 52793800), ("chr18", 74615900, 74618100),
    ("chr18", 76966200, 76968500), ("chr18", 78436900, 78438700),
    ("chr18", 79013800, 79040300), ("chr18", 79617800, 79621500),
    ("chr18", 80257400, 80373200), ("chr19", 0, 271200),
    ("chr19", 7019100, 7061300), ("chr19", 7449400, 7452000),
    ("chr19", 8740100, 8800500), ("chr19", 24330100, 27274500),
    ("chr19", 27337600, 27427400), ("chr19", 34386800, 34393500),
    ("chr19", 34860600, 34866200), ("chr19", 36267900, 36313700),
    ("chr19", 37264900, 37304300), ("chr19", 44393300, 44416700),
    ("chr19", 47903000, 47959700), ("chr19", 50090500, 50140400),
    ("chr19", 58538700, 58617600), ("chr20", 0, 67900),
    ("chr20", 26364200, 28916900), ("chr20", 28939400, 29264700),
    ("chr20", 30995400, 31246000), ("chr20", 47893800, 47900200),
    ("chr21", 0, 8679600), ("chr21", 9159900, 9735300),
    ("chr21", 10013900, 10069600), ("chr21", 10094700, 10505100),
    ("chr21", 10650900, 12965800), ("chr21", 43212400, 43280900),
    ("chr21", 46682700, 46709900), ("chr22", 10687700, 11428100),
    ("chr22", 11496900, 11873100), ("chr22", 11976900, 15154400),
    ("chr22", 16258000, 16385800), ("chr22", 18175900, 18947300),
    ("chr22", 20337400, 20343300), ("chr22", 21113500, 21554000),
    ("chr22", 49972700, 49975300), ("chr22", 50642800, 50644900),
    ("chr22", 50786600, 50818400), ("chrX", 0, 329300), 
    ("chrX", 362400, 388500), ("chrX", 456500, 531800), 
    ("chrX", 723800, 739500), ("chrX", 864500, 930400), 
    ("chrX", 1049100, 1054300), ("chrX", 1085100, 1175500), 
    ("chrX", 1200600, 1209400), ("chrX", 1249200, 1269000), 
    ("chrX", 1289500, 1298900), ("chrX", 1365300, 1458700), 
    ("chrX", 1480900, 1492800), ("chrX", 1816200, 1820600), 
    ("chrX", 2223900, 2521900), ("chrX", 2580600, 2751300), 
    ("chrX", 3966700, 3968700), ("chrX", 5481200, 5486100), 
    ("chrX", 6933400, 6938700), ("chrX", 7587600, 7591800), 
    ("chrX", 9403600, 9415100), ("chrX", 10785000, 10809700), 
    ("chrX", 10966600, 10976800), ("chrX", 11218800, 11221100), 
    ("chrX", 11840900, 11848000), ("chrX", 14085100, 14109500), 
    ("chrX", 14286500, 14289300), ("chrX", 16361200, 16366000), 
    ("chrX", 16498100, 16503400), ("chrX", 19940200, 19946300), 
    ("chrX", 21340600, 21345700), ("chrX", 25773300, 25776000), 
    ("chrX", 26176400, 26181400), ("chrX", 30767800, 30772600), 
    ("chrX", 31077600, 31082600), ("chrX", 31511400, 31535800), 
    ("chrX", 34416800, 34425900), ("chrX", 36465200, 36471200), 
    ("chrX", 37628400, 37633500), ("chrX", 42872300, 42910700), 
    ("chrX", 49317500, 49623500), ("chrX", 50019400, 50033700), 
    ("chrX", 50056700, 50066100), ("chrX", 51202300, 51268100), 
    ("chrX", 51427500, 51432400), ("chrX", 52175000, 52228100), 
    ("chrX", 52442800, 52538100), ("chrX", 53761700, 53789500), 
    ("chrX", 55180400, 55184500), ("chrX", 56754900, 56781100), 
    ("chrX", 57712300, 57719700), ("chrX", 58467900, 62522800), 
    ("chrX", 63129600, 63290600), ("chrX", 67311800, 67323800), 
    ("chrX", 67626800, 67632300), ("chrX", 68217300, 68230200), 
    ("chrX", 70600000, 70603800), ("chrX", 70640600, 70645000), 
    ("chrX", 70963600, 70964900), ("chrX", 71978800, 71980500), 
    ("chrX", 72489400, 72490800), ("chrX", 72743200, 73035800), 
    ("chrX", 73381000, 73387000), ("chrX", 73887000, 73891300), 
    ("chrX", 74660000, 74718100), ("chrX", 74789000, 74794000), 
    ("chrX", 74952200, 74995200), ("chrX", 78802400, 78804500), 
    ("chrX", 79765500, 79789600), ("chrX", 80534100, 80537000), 
    ("chrX", 82849700, 82859300), ("chrX", 83752100, 83756900), 
    ("chrX", 86046600, 86076600), ("chrX", 86395500, 86398100), 
    ("chrX", 86970000, 86975600), ("chrX", 87220500, 87222100), 
    ("chrX", 89060200, 89062700), ("chrX", 89202500, 89208400), 
    ("chrX", 91332900, 91336600), ("chrX", 93618000, 93633400), 
    ("chrX", 94863600, 94868300), ("chrX", 97509600, 97515000), 
    ("chrX", 100135800, 100141000), ("chrX", 100257100, 100261600), 
    ("chrX", 101471700, 101474900), ("chrX", 102188700, 102489200), 
    ("chrX", 103851800, 103897800), ("chrX", 106755500, 106769400), 
    ("chrX", 106813900, 106830900), ("chrX", 107515800, 107517200), 
    ("chrX", 109034800, 109069100), ("chrX", 109114900, 109119400), 
    ("chrX", 109520800, 109525700), ("chrX", 109985900, 109987300), 
    ("chrX", 110816700, 110833400), ("chrX", 111416100, 111418000), 
    ("chrX", 113141700, 113143600), ("chrX", 114701600, 114724300), 
    ("chrX", 115725600, 115889600), ("chrX", 116557600, 116595600), 
    ("chrX", 117874100, 117880000), ("chrX", 118009000, 118037800), 
    ("chrX", 118070900, 118072700), ("chrX", 121263700, 121268100), 
    ("chrX", 121299200, 121300600), ("chrX", 122528400, 122550000), 
    ("chrX", 124584300, 124588400), ("chrX", 125927600, 125937100), 
    ("chrX", 126463700, 126474200), ("chrX", 127116700, 127122600), 
    ("chrX", 127362200, 127368300), ("chrX", 128785000, 128788700), 
    ("chrX", 129337600, 129357900), ("chrX", 129388400, 129408400), 
    ("chrX", 130567700, 130572000), ("chrX", 131152200, 131157400), 
    ("chrX", 131378300, 131383300), ("chrX", 131664300, 131670000), 
    ("chrX", 132284600, 132320400), ("chrX", 133108600, 133116500), 
    ("chrX", 135718600, 135888700), ("chrX", 137074700, 137079100), 
    ("chrX", 137436600, 137439300), ("chrX", 138300600, 138302200), 
    ("chrX", 139437600, 139446800), ("chrX", 139621500, 139622800), 
    ("chrX", 140722400, 140726100), ("chrX", 141000400, 141108300), 
    ("chrX", 142478000, 142483800), ("chrX", 142892300, 142911600), 
    ("chrX", 143352000, 143356500), ("chrX", 144404500, 144475900), 
    ("chrX", 147281700, 147287100), ("chrX", 147653800, 147659900), 
    ("chrX", 148123500, 148129000), ("chrX", 148347100, 148378700),
    ("chrX", 149437900, 149441900), ("chrX", 150024800, 150026200), 
    ("chrX", 152173800, 152175100), ("chrX", 153251200, 153316400), 
    ("chrX", 154870000, 154890200), ("chrX", 154938900, 154945100), 
    ("chrX", 155299600, 155305100), ("chrX", 155454000, 155522000), 
    ("chrX", 155700400, 155727500), ("chrX", 155983500, 156040800),
    ("chrY", 4343800, 4345800), ("chrY", 10246200, 11041200),
    ("chrY", 11072100, 11335300), ("chrY", 11486600, 11757800),
    ("chrY", 26637300, 57227400)]

    # Check if the output BED file already exists; check column formatting. 
    if os.path.exists('annotation.bed'):
        try:
            print(f"\nAttempting to read regions from: 'annotation.bed'...")
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
        num_threads = min(os.cpu_count(), 4)
        min_coverage = 10 if min_coverage is None else int(min_coverage)

        mosdepth_args = ['mosdepth', '-n',
                        '-b', 'annotation.bed', 
                        '-t', str(num_threads), 
                        'tmp', sorted_bam]
        subprocess.run(mosdepth_args, check=True, stderr=subprocess.DEVNULL)

    if skip_chromosome:
        skipped = skip_chromosome.split(',')

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

def mosdepth_positions(sorted_bam, window_size=None, 
    skip_chromosome=None, chromosome=None):
    """
    Function that parses the provided sorted_bam file to determine mean
    coverage at the various nucleotide positions. Uses mosdepth to determine 
    reads at each nucleotide "window", 

    Args:
        sorted_bam:     PATH to BAM file
        window_size:    size of region to bin nucleotide positions into.
    
    (Optional)
        skip_chromosome:    comma separated list of chromosomes to skip 
                            from the input file.
        chromosome:         comma separated list of chromosomes to keep
                            from the input file.
    """
    valid_positions = []
    if os.path.exists('tmp.regions.bed.gz'):
        print(f"\nMean coverage summary for positions exists. Using for" 
            " analysis...")
        with gzip.open('tmp.regions.bed.gz', 'rt') as f:
            chromosomes = set()
            for line in f:
                chrom, start, end, cov_mean = line.strip().split()
                valid_positions.append((chrom, int(start), 
                    int(end), float(cov_mean)))
                chromosomes.add(chrom)
    else:
        num_threads = min(os.cpu_count(), 4)
        window = 1000 if window_size is None else int(window_size)
  
        mosdepth_args = ['mosdepth', '-n', 
                        '-b', str(window), 
                        '-t', str(num_threads), 
                        'tmp', sorted_bam]
        subprocess.run(mosdepth_args, check=True, stderr=subprocess.DEVNULL)

        with gzip.open('tmp.regions.bed.gz', 'rt') as f:
            chromosomes = set()
            for line in f:
                chrom, start, end, cov_mean = line.strip().split()
                valid_positions.append((chrom, int(start), 
                    int(end), float(cov_mean)))
                chromosomes.add(chrom)

    if chromosome: 
        kept = chromosome.split(',')
        valid_positions = ([pos for pos in valid_positions 
        if pos[0] in kept])

    if skip_chromosome: 
        skipped = skip_chromosome.split(',')
        valid_positions = ([pos for pos in valid_positions 
        if pos[0] not in skipped])

    if chromosome and skip_chromosome:
        common_elements = set(kept) & set(skipped)
        if common_elements:
            print("Warning: There are chromosomes listed in both the '-skip' "
                " and '-chrom' options. Remove from one option before "
                " continuing:", common_elements)

    return valid_positions

def plot_histogram(chrom, data):
    print(f"Plotting histogram for {chrom}.")
    start_positions = data['start_positions']
    end_positions = data['end_positions']
    cov_means = data['cov_means']
    bar_widths = [end - start for start, end in 
        zip(start_positions, end_positions)]
    # Change size and style.
    plt.ioff()
    plt.figure(figsize=(10, 6)); ax = plt.gca()
    plt.bar(start_positions, cov_means, width=bar_widths, linewidth=0.3, 
        color='lightskyblue', edgecolor='black')
    plt.xlim(left=0)

    # Despine plot.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Add titles, tickmarks, spacing.
    plt.title(f'{chrom} Mean Coverage')
    plt.xlabel('Nucleotide Position', fontsize=12)
    formatter = ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x))
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.xticks(fontsize=10)
    plt.ylabel('Mean Coverage', fontsize=12)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig(f'coverage_histogram_{chrom}.svg')
    plt.close()

def plot_coverage(valid_positions, window_size):
    """
    Reads the return from position based analysis to write output to file. 
    Args:
        valid_positions (list): List of tuples containing chromosome, start, 
                                end, and coverage values.
        output_dir (str):       Path to an output directory. 
    Returns:
        File:   SVG histogram plot for specific chromosomes.
    """    
    # chromosome_lengths for reference
    # ("chr1", 248956422), ("chr2", 242193529),
    # ("chr3", 198295559), ("chr4", 190214555), 
    # ("chr5", 181538259), ("chr6", 170805979), 
    # ("chr7", 159345973), ("chr8", 145138636),
    # ("chr9", 138394717), ("chr10", 133797422), 
    # ("chr11", 135086622), ("chr12", 133275309), 
    # ("chr13", 114364328), ("chr14", 107043718),
    # ("chr15", 101991189), ("chr16", 90338345), 
    # ("chr17", 83257441), ("chr18", 80373285), 
    # ("chr19", 58617616), ("chr20", 64444167),
    # ("chr21", 46709983), ("chr22", 50818468), 
    # ("chrX", 156040895), ("chrY", 57227415), 
    # ("RN7SK", 331), ("RN7SL", 288), 
    # ("RNU7", 63), ("RNY", 694), 
    # ("U13", 120), ("U14AB", 283), 
    # ("U17", 207), ("U3", 217), 
    # ("U8", 136), ("12S", 954), 
    # ("16S", 1559), ("18S", 12994), 
    # ("hs5S", 121), ("hssnRNA", 2058)

    # Iterate through valid_positions to populate coverage information.
    coverage_data = {}
    for entry in valid_positions:
        chrom, start, end, cov_mean = entry
        if chrom not in coverage_data:
            coverage_data[chrom] = {'start_positions': [], 
            'end_positions': [], 'cov_means': []}
        coverage_data[chrom]['start_positions'].append(start)
        coverage_data[chrom]['end_positions'].append(end)
        coverage_data[chrom]['cov_means'].append(cov_mean)

    try:
        nproc = subprocess.run(['nproc'], capture_output=True,
            text=True, check=True)
        thread_count = int(nproc.stdout.strip())
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        thread_count = 1

    print(f"\nUsing {thread_count} threads for plotting...")
    pool = multiprocessing.Pool(processes=thread_count)

    for chrom, data in coverage_data.items():
        if data['cov_means'] and any(mean > 0 for mean in data['cov_means']):
            pool.apply_async(plot_histogram, (chrom, data))
        else:
            print(f"No coverage data found for {chrom}.")

    pool.close()
    pool.join()

###########################################################################
# 3. Main, define accepted arguments. 
def parse_args():
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

NOTE: Arguments must be provided in the following order:

1. analysis_type:       Available options: 'gene', 'position'

                        Gene region-based uses mean coverage of reads for 
                        genomic coordinate regions to determine which genes 
                        are covered above a minimum coverage, writes 'genes' 
                        into a BED6 file for CRSSANT DG assembly. Requires
                        annotation file to be specified (see #4 [-a]).
                        
                        Position-based uses read coverage at nucleotide
                        positions of n=<window_size> to determine regions of 
                        high coverage (user defined) and plots a histogram
                        for read coverage profile.  

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

                        NOTE: Specific genomic regions or RNAs are masked and 
                        added back as a separate 'chromosome' in the PARIS
                        SHARC pipelines. These must manually be added as
                        separate lines to the annotation_file. Always check 
                        to see if the chromosomes are in the output BED6 files
                        generated using this script.

                        List of genes to be added as a separate 'chromosome':
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

-w, --window-size:      Defaults to [1,000] if no value is provided.
                        
                        Window size to check for positions from input file.  
                        Start with a larger value to minimize compute time.

-chrom, --chromosome:   Provides a list of comma separated chromosomes to 
                        keep when performing the histogram plotting.
                        Recommended for reducing processing time.

                            -chrom=chr1,chr2,chrM

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
                        
                            Initial arguments passed
                            Job timing
 
###########################################################################
"""),
    usage=
"""
    \npython3 %(prog)s analysis_type input output [options]

    Usage examples:
        
    - For gene region-based analysis: 
        %(prog)s gene input output -a=annotation_file \
[-min, -max, -skip, -stats, -k]\n

    - For position-based analysis: 
        %(prog)s position input output \
[-w, -max, -skip, -stats, -k]
""")

    # Main arguments to be provided in the command line.
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
        " Defaults to [1,000]")
    parser.add_argument('-chrom', '--chromosome', 
        help="Optional parameter to specify which chromosome to plot. "
        " Defaults to [None]")

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

    return parser.parse_args()

def main():
    ########################################################################  
    # Begin timing the job for the --statistics flag. 
    start_time = time.time()

    # Check for optional arguments to modify functions.
    min_coverage = args.min_coverage # Gene-based: Default 10
    window_size = args.window_size # Position-based: Default 20
    keep_files = args.keep_files # Default: False
    max_reads = args.max_reads # Default: None
    skip_chromosome = args.skip_chromosome # Default: None
    statistics = args.statistics # Default: None
    chromosome = args.chromosome # Default: None

    # Check for mandatory positional arguments.
    analysis_type = args.analysis_type # Must be 'gene' or 'position'
    input_file = args.input # Must be in unsorted BAM or SAM format
    output = args.output # Any string, no file extension necessary
    annotation = args.annotation # Annotation file in modified GTF format

    # Parse required arguments.
    if analysis_type != "gene" and analysis_type != "position": 
        print(f"Provided analysis type not supported."
            f" Use 'gene' or 'position'.")
        sys.exit()

    # Parse input file for analysis, set up the output.
    file_extension = os.path.splitext(input_file)[1].lower()
    if os.path.splitext(input_file)[1] != ".sam" and \
        os.path.splitext(input_file)[1] != ".bam":
        print(f"File type not supported."
            f" Check if the input file is in SAM or BAM format.")
        sys.exit()
    else:
        sorted_bam = sort_bam(input_file, max_reads)    

    if '.bed' in output:
        output = os.path.splitext(output)[0]

    # Gene-based analysis, output overlaps to BED6.
    if max_reads == True:
        if min_coverage is not None and float(min_coverage) < 10:
            print("WARNING: Using a low min_coverage value increases" 
            " analysis time. CTRL + C to cancel now if desired...")
            try:
                time.sleep(5)
                print("Continuing...")
            except KeyboardInterrupt:
                print("Operation canceled by user.")
                sys.exit(1)

    if analysis_type == "gene":
        gene_regions = collapse_gene_regions(annotation)
        covered_genes = mosdepth_regions(sorted_bam, 
            annotation, min_coverage, skip_chromosome)
        region_to_bed(covered_genes, gene_regions, output)

    # Coverage by nucleotide position, output to histograms.
    if analysis_type == "position":
        positions = mosdepth_positions(sorted_bam, 
            window_size, skip_chromosome, chromosome)
        plot_coverage(positions, window_size)
        output_directory = "coverage_histograms"
        os.makedirs(output_directory, exist_ok=True)
        for filename in os.listdir('.'):
            if filename.endswith('.svg'):
                shutil.move(filename, 
                    os.path.join(output_directory, filename))

    # End the job here and provide the statistics.
    end_time = time.time()
    print(f"\nBAMget.py completed on {timenow()}.\n")

    ########################################################################

    # Set up the statistics 
    if statistics == True:
        elapsed_time = "{:.2f}".format(end_time - start_time)
        if min_coverage and float(min_coverage) < 10:
            elapsed_time = "{:.2f}".format(float(elapsed_time) - 5.00)
        print(f"     Elapsed time: {elapsed_time}s")
        if analysis_type == "gene":
            count = sum("tRNA" not in gene for gene in covered_genes)
            print(f"  Genes collected: {count}")

    # Remove intermediate files unless -k argument is provided.
    if keep_files == True:
        print(f"\nWARNING: Keeping intermediate files between analysis types"
            f" may lead to incorrect coverage values. Remove and start over"
            f" if this was unintended.\n")
    else:
        os.remove('annotation.bed')
        os.remove('tmp.regions.bed.gz')
        os.remove('tmp.regions.bed.gz.csi')

    # Remove other output files from mosdepth.
    try:
        os.remove('tmp.mosdepth.global.dist.txt')
        os.remove('tmp.mosdepth.region.dist.txt')
        os.remove('tmp.mosdepth.summary.txt')
        
    #Ignore removal if the file does not exist.
    except FileNotFoundError:
        pass 

    ########################################################################

if __name__ == "__main__":
    args = parse_args()
    main()
    sys.exit()