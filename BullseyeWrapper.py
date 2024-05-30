#This is the main script to run the Bullseye pipeline
#Required dependencies and instructions for installing the proper conda environment can be found on https://github.com/mflamand/Bullseye
#Must install FastQC (conda install bioconda::fastqc) to environment if not found.
#
#This script was written as a wrapper for the Bullseye pipeline scripts and was written by Matt Tegowski, and is not sanctioned by the original writer of the pipeline
#
#Written for Python 3.7.12
#
#
# Outline
# 1) FastQC (Skip if not present)
# 2) Make STAR index is needed, this can be provided as an option. Make one by default if one doesn't exist.
# 3) Align to genome (Make organism of interest a required option?)
# 4) Perform feature counts for integration into DESeq2
# 5) parseBAM.pl
# 6) findEditSites.pl
# 7) RAC-filter
# 8) Metagene analysis (install R? Version?)
# 9) Final docs
# 	a) Unfiltered and filtered bed files
#	b) Doc with numbers of sites/RNAs, sites in each motif

# Important options
# Names of each input file
# Additional STAR options
# Different "modes" that allow you to skip sections of the analysis
# A no-index option, requires a path to index (GenomeDir)
# m6A filtering options
# Path to annotation files for metagene (Default is GenomeDir)

import argparse
import os
import sys
import gzip
import multiprocessing
import time
import builtins
import subprocess
from time import strftime
import re
from pathlib import Path
import asyncio
from multiprocessing import Process
import string
import random

sys.stdout = open("LogBullseye.txt", "w", buffering=1)
def print(text):
    builtins.print(text)
    os.fsync(sys.stdout)

parser = argparse.ArgumentParser(description='Wrapper for running the Bullseye pipeline for m6A identification from DART-seq experiments')
parser.add_argument('--fqR1Test', type=str, nargs='*', default=None, help="Input Read 1 Fastq files from APOBEC1-YTH-treated sample")
#parser.add_argument('--fqR1Test', type=argparse.FileType('r'), nargs='*', help="Input Read 1 Fastq files from APOBEC1-YTH-treated sample")
#f1testnames=parser.parse_args(['--fqR1Test'])
#_StoreAction(option_strings=['--fqR1Test'], dest='R1files', nargs=None, const=None, default=None, type=FileType('r'), choices=None, help=None, metavar=None)
parser.add_argument('--fqR2Test', type=str, nargs='*', default=None, help="Input Read 2 Fastq files from APOBEC1-YTH-treated sample")
parser.add_argument('--fqR1Control', type=str, nargs='*', default="Not supplied", help="Input Read 2 Fastq files from APOBEC1-YTHmut-treated sample")
parser.add_argument('--fqR2Control', type=str, nargs='*', default="Not supplied", help="Input Read 2 Fastq files from APOBEC1-YTHmut-treated sample")
parser.add_argument("--genomeDir", help="Directory containing the genome fasta file and refFlat annotation files")
parser.add_argument("--genomeFile", help="Path to and name of genome fasta file")
parser.add_argument("--STARgtf", help="GTF File for generating STAR index")
parser.add_argument("--noFastQC", help="Option to skip FastQC, even if mode is all", action='store_true')
parser.add_argument("--noSTARIndex", help="Skip generation of STAR index, even if mode is all", action='store_true')
parser.add_argument("--singleEnd", help="Use if not paired-end sequencing. Assumed data is paired-end by default", action='store_true')
parser.add_argument("--noControl", help="Use if there is no APOBEC1-YTHmut-treated control data", action='store_true')
parser.add_argument("--mode", nargs='?', type=str, help="all, bullseyeonly, findsites, RACfilter, metagene", default='all')
parser.add_argument("--indexDir", type=str, nargs='?', help="If STAR index already exists, specificy the path to its directory", default='GenomeDir/')
parser.add_argument("--memoryGB", type=int, help="Available memory in GB", default=48)
parser.add_argument("--workingDir", type=str, nargs='?', help="Base directory to use", default='./')
parser.add_argument("--matrixCov", type=int, help="Minimum reads of coverage for a nucleotide to be included in matrix generation", default="5")
parser.add_argument("--noAlign", help="Use if you supply fasta files, but already have aligned files in workdir/bamfiles", action='store_true')
parser.add_argument("--bamfiles", type=str, nargs='*', default="Not supplied", help="List of bamfiles (filenames only! Path is assumed to be workdir/bamfiles) for parsing into matrix, only needed if using bullseyeonly mode")

#Arguments for finding edit sites
parser.add_argument("--annotationFile", nargs='?', type=str, help="See Bullseye github for making flatREF annotation file. If not using genome, this file and path MUST be provided. It must not be provided if comparing to genome or knownSites", default="20")
parser.add_argument("--editType", nargs='?', type=str, help="Type of mutation being queried. Default is C2U, but any combination of _2_ nucleotide works, accepts U (treats as T) and I (treats as G)", default='C2U')
parser.add_argument("--minEdit", nargs='?', type=str, help="Minimum mutation percentage. Default is 10%", default="10")
parser.add_argument("--maxEdit", nargs='?', type=str, help="Maximum mutation percentage. Default is 95%", default="95")
parser.add_argument("--foldOverControl", nargs='?', type=str, help="Minimum fold-change in mutation percentage over control at any give site for calling. Default is 1.5", default="1.5")
parser.add_argument("--minCov", nargs='?', type=str, help="Minimum reads of coverage for site calling. If samples are sequenced deeply it is recommended to increase this number as high as reasonably possible.", default="20")
parser.add_argument("--minControlCov", nargs='?', type=str, help="Minimum reads of coverage in control sample for comparisons. Default is 10", default="10")
parser.add_argument("--cpuAvail", nargs='?', type=str, help="Number of available processing cores. Default is 1", default="1")
parser.add_argument("--minEdits", nargs='?', type=str, help="Minimum number of mutations observed at a given site for it to be called. It is highly recommended to never use 1. Default is 2. Higher than 2 is recommended is possible, but too high of a number may lead to false negatives for real sites with low editing.", default="2")
parser.add_argument("--Intron", help="Use to find sites within intronic sequences, otherwise only exons will be used.", action='store_true')
parser.add_argument("--UTRextend", help="Extend site calling 5kb from annotated end of 3'UTR.", action='store_true')
parser.add_argument("--fallback", help="Use this option if you want to consider sites with low/no coverage in control file.", action='store_true')
parser.add_argument("--genome", help="Do not use control matrix for comparison. Just use other thresholds set and compare to the genome instead.", action='store_true')
parser.add_argument("--score", help="score and filter reads based on probability of edit based on beta distribution curve. Sites with 10 fold higher probability in DART/control will be kept.", action='store_true')
parser.add_argument("--knownSites", nargs='?', type=str, help="Instead of finding all sites, only find sites contained within bed file provided with path. Must be bed-6", default='None')
parser.add_argument("--stranded", help="Use this option if you used a stranded sequencing kit.", action='store_true')
parser.add_argument("--filterBed", nargs='?', type=str, help="Do not include sites contained within bedfile. Must be bed-6. Provide path and filename.", default='None')
args = parser.parse_args()

global mode,indexdir,workdir,gendir,genfile,gtffile,testR1,testR2,controlR1,controlR2,r1testnames,r2testnames,r1controlnames,r2controlnames,r1allnames,r2allnames,edittype,mineditpct,maxeditpct,foldovercntrl,mincoverage,mincntrlcoverage,cpus,minedits,intron,utrext,fallback,genomecompare,bdist,knownsites,stranded,filtersites

mode = args.mode
indexdir = args.indexDir
workdir = args.workingDir
gendir = args.genomeDir
genfile = args.genomeFile
gtffile = args.STARgtf
memavail = args.memoryGB
matcov = args.matrixCov
bams = args.bamfiles

#Finding edit sites options
annotfile = args.annotationFile
edittype = args.editType
mineditpct = args.minEdit
maxeditpct = args.maxEdit
foldovercntrl = args.foldOverControl
mincoverage = args.minCov
mincntrlcoverage = args.minControlCov
cpus = args.cpuAvail
minedits = args.minEdits
intron = args.Intron
utrext = args.UTRextend
fallback = args.fallback
genomecompare = args.genome
bdist = args.score
knownsites = args.knownSites
stranded = args.stranded
filtersites = args.filterBed

if args.fqR1Control and args.fqR2Control is "Not supplied":
    print("WARNING: No control files supplied. Running without them. Will compare editing to genome.")
    testR1 = args.fqR1Test
    testR2 = args.fqR2Test
    r1testnames = [os.path.basename(file_path) for file_path in testR1]
    r2testnames = [os.path.basename(file_path) for file_path in testR2]
    r1allnames = r1testnames
    r2allnames = r2testnames
    r1allpathnames=testR1
    r2allpathnames=testR2
else:
    controlR1 = args.fqR1Control
    controlR2 = args.fqR2Control
    testR1 = args.fqR1Test
    testR2 = args.fqR2Test
    r1testnames = [os.path.basename(file_path) for file_path in testR1]
    r2testnames = [os.path.basename(file_path) for file_path in testR2]
    r1controlnames = [os.path.basename(file_path) for file_path in controlR1]
    r2controlnames = [os.path.basename(file_path) for file_path in controlR2]
    r1allnames = r1testnames + r1controlnames
    r2allnames = r2testnames + r2controlnames
    r1allpathnames=testR1 + controlR1
    r2allpathnames=testR2 + controlR2

print("Starting pipeline")
print("Test Read 1 Filenames:")
for file_name in r1testnames:
        print(file_name)
print("Test Read 2 Filenames:")
if args.singleEnd is not True:
    for file_name in r2testnames:
        print(file_name)
    else:
        print(" ")
print("Control Read 1 Filenames:")
if args.fqR1Control and args.fqR2Control is not "Not supplied":
    for file_name in r1controlnames:
        print(file_name)
    print("Control Read 2 Filenames:")
    for file_name in r2controlnames:
        print(file_name)
else:
    print(" ")
print("Location of fastq files: ")

for file_name in r1allpathnames:
    print(file_name)
    
print("*************************************************************************")
 
def indexing(gendir, genfile, gtffile):   
    if gendir and genfile and gtffile is not None:      
        indexcommand = 'STAR --runMode genomeGenerate --genomeFastaFiles ' + gendir + '/' + genfile + ' --sjdbGTFfile ' + gendir + '/' + gtffile + ' --genomeSAindexNbases 3'
        print(indexcommand)
        subprocess.call(indexcommand, shell=True)
    else:
        print("To generate the STAR index, please provide a directory (--genomeDir) containing a genome fasta file (--genomeFile) and a GTF (--STARgtf). If a STAR index exists, use the --noSTARIndex option")
    return("Done indexing")
 
def FastQC(r1testnames, r2testnames, r1controlnames, r2controlnames):
    print("Starting FastQC analysis")
    if args.singleEnd is True:
        singleoption = 1
    else:
        singleoption = 0
    if args.noControl is True:
        controloption = 1
    else:
        controloption = 0
    if singleoption > 0 and controloption > 0:
        for f in r1testnames:
            qccommand = 'fastqc ' + f
            print(qccommand)            
            subprocess.call(qccommand, shell=True)
    elif singleoption > 0 and controloption < 1:
        for f in r1testnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)
        for f in r1controlnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)
    elif singleoption < 1 and controloption > 0:
        for f in r1testnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)
        for f in r2testnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)
    elif singleoption < 1 and controloption < 1:   
        for f in r1testnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)
        for f in r2testnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)
        for f in r1controlnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)
        for f in r2controlnames:
            qccommand = 'fastqc ' + f
            print(qccommand)
            subprocess.call(qccommand, shell=True)

def Mapping(sem, r1allpathnames, r2allpathnames, indexdir, workdir, file1, file2):
    print("Using index directory " + indexdir)
    mapcommand = 'STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFilterMismatchNoverReadLmax 0.06 --outTmpDir ' + bamdir + ' --runThreadN 8 --genomeLoad LoadAndKeep --limitBAMsortRAM 20000000000 --genomeDir ' + indexdir + ' --outFileNamePrefix ' + workdir + '/bamfiles/' + stemname2 + ' --readFilesIn ' + file1 + ' ' + file2
    submitmapping = [subprocess.call(mapcommand, shell=True)]

def Parsing(r1allnames, workdir, matrixdir, matcov):
    parsecomm = 'perl ' + workdir + 'software/parseBAM.pl --input ' + bamfilein + ' --output ' + matrixdir + matstemname2 + '.matrix --verbose --removeDuplicates --minCoverage ' + str(matcov)
    print(parsecomm)
    subprocess.call(parsecomm, shell=True)
    
def FindSites(r1testnames, workdir, matrixin, annotfile, outbedname):
    if intron is True:
        intronoption = '--intron '
    else:
        intronoption = ''
    if utrext is True:
        utroption = 'extUTR:5000 '
    else:
        utroption = ''
    if fallback is True:
        fallbackoption = '--fallback '
    else:
        fallbackoption = ''
    if genomecompare is True:
        genomeoption = '--genome ' + genfile
    else:
        genomeoption = ''
    if bdist is True:
        bdistoption = '--score '
    else:
        bdistoption = ''
    if knownsites == "None":
        knownsitesopt = ''
    else:
        knownsitesopt = knownsites
    if stranded is True:
        strandopt = '--stranded '
    else:
        strandopt = ''
    if filtersites == "None":
        filteropt = ''
    else:
        filteropt = filtersites    
    sitescomm = 'perl ' + workdir + 'software/Find_edit_site.pl -a ' + annotfile + ' -d ' + matrixin + ' -c ' + controlmatrix + ' -e ' + edittype + ' --minEdit ' + mineditpct + ' --maxEdit ' + maxeditpct + ' -t ' + foldovercntrl + ' -o ' + outbedname + ' -v --EditedMinCoverage ' + mincoverage + ' --ControlMinCoverage ' + mincntrlcoverage + ' --bed6 --cpu ' + cpus + ' --mem ' + str(memavail) + ' -me ' + minedits + intronoption + utroption + fallbackoption + genomeoption + ' ' + bdistoption + knownsitesopt + ' ' + strandopt + filteropt
    print("find sites command \n")
    print(sitescomm)
    subprocess.call(sitescomm, shell=True)

if args.noSTARIndex is False:
    if mode == 'all':
        print("Generating STAR Index.")
        indexing(gendir, genfile, gtffile)
    else:
        print("Skipping Index Generation due to mode selection")
else:
    print("noSTARIndex selected, skipping genome index generation")

if args.noFastQC is False:
    if mode == 'all':
        FastQC(testR1, testR2, controlR1, controlR2)
    else:
        print("Skipping FastQC due to mode selection")
else:
    print("noFastQC selected, skipping fastqc analysis")

if args.noAlign is False:
    if mode == 'all':
        # Create a Process for each file
        if memavail <= 50:
            njobs = 1
        elif memavail > 50 and memavail <=100:
            njobs = 2
        elif memavail > 100 and memavail <=150:
            njobs = 3
        elif memavail > 150 and memavail <=200:
            njobs = 4
        elif memavail > 200 and memavail <=250:
            njobs = 5
        elif memavail > 250 and memavail <=300:
            njobs = 6
        elif memavail > 300 and memavail <=350:
            njobs = 7
        elif memavail > 350:
            njobs = 8
        else:
            print("Something went wrong.")
        print(str(njobs))
        print("Number of parallel jobs aligning at once is " + str(njobs) + " due to memory allocation. Available memory is " + str(memavail) + ". If it listed as 48, and you have more than 48GB avaialable, use the --memoryGB option to specify to allow more alignments to run at once. STAR is fast but memory intensive.")
        removedir = 'rm -r ' + workdir + '/tmp/'
        print(removedir)
        subprocess.call(removedir, shell = True)
        tmpdir = workdir + '/tmp'
        addtmp = 'mkdir -p ' + tmpdir
        subprocess.call(addtmp, shell = True)
        sem = asyncio.Semaphore(njobs)
        tottasks = len(r1allnames)
        processes = []
        for i in range(tottasks):
            sem.acquire()
            file1 = r1allpathnames[i-1]
            file2 = r2allpathnames[i-1]
            file3 = r1allnames[i-1]
            # for file1, file2, file3 in zip(r1allpathnames, r2allpathnames, r1allnames):
            stemname = os.path.splitext(os.path.basename(file3))[0]
            stemname2 = stemname.split(".fastq.gz")[0]
            bamdir = workdir + '/tmp/' + stemname2              
            p = Process(target=Mapping, args=(sem, r1allpathnames, r2allpathnames, indexdir, workdir, file1, file2))          
            p.start()
            processes.append(p)
            # Wait for all processes to complete
            for p in processes:
                p.join()
                print("All alignments have completed. \n")
    else:
        print("Skipping STAR alignment due to mode selection")
else:
    print("--noAlign selected, skipping alignment, but fastq files must be provided")

#In this section, nucleotide representation matrixes are made. They will be used in comparing %C2U values and calling sites in the next section
if mode == 'all':
    addmatrixdir = 'mkdir -p ' + workdir + '/matrix'
    subprocess.call(addmatrixdir, shell=True)
    matrixdir = workdir + '/matrix/'
    print("*************************************************************")
    print("Parsing bamfiles to generate nucleotide representation matrix \n")
    matrixprocesses = []
    for file in r1allnames:
        matstemname = os.path.splitext(os.path.basename(file))[0]
        matstemname2 = matstemname.split(".fastq.gz")[0]
        bamfilein = workdir + '/bamfiles/' + matstemname2 + 'Aligned.sortedByCoord.out.bam'
        print(bamfilein)
        ma = Process(target=Parsing, args=(r1allnames, workdir, matcov))
        matrixprocesses.append(ma)
        ma.start()
    # Wait for all processes to complete
    for ma in matrixprocesses:
        ma.join()
    print("All matrixes are completed. \n")
    print("Combining control matrixes. \n")
    controlr1basenames = [os.path.basename(path) for path in controlR1]
    controlr1basenames2 = [stemname.split("_")[0] for stemname in controlr1basenames]
    matrixsuffix = '.matrix'
    matrixdir = workdir + '/matrix/'
    controlmatrixlist = [matrixdir + s + matrixsuffix for s in controlr1basenames2]
    allr1controlnames = " ".join(controlmatrixlist)
    print(allr1controlnames)
    concatmats = 'zcat ' + allr1controlnames + ' > ' + matrixdir + 'tempcontrolmat.matrix'
    subprocess.call(concatmats, shell=True)
    sortmats = 'sort -k1,1 -k2,2n ' + matrixdir + 'tempcontrolmat.matrix -o ' + matrixdir + 'tempcontrolmat.matrix'
    subprocess.call(sortmats, shell=True)
    collapsemats = 'perl ' + workdir + 'software/collapse_matrix.pl -i '+ matrixdir + 'tempcontrolmat.matrix -o ' + matrixdir + 'combinedcontrolmatrixes.matrix'
    subprocess.call(collapsemats, shell=True)
    print("Finished combining control matrixes. \n")
elif mode == 'bullseyeonly':
    addmatrixdir = 'mkdir -p ' + workdir + '/matrix/'
    subprocess.call(addmatrixdir, shell=True)
    print("*************************************************************")
    print("Parsing bamfiles to generate nucleotide representation matrix \n")
    matrixprocesses = []
    for file in bams:
        matstemname = os.path.splitext(os.path.basename(file))[0]
        matstemname2 = matstemname.split(".fastq.gz")[0]
        bamfilein = matstemname2 + 'Aligned.sortedByCoord.out.bam'
        ma = Process(target=Parsing, args=(r1allnames, workdir, matcov))
        matrixprocesses.append(ma)
        ma.start()
    # Wait for all processes to complete
    for ma in matrixprocesses:
        ma.join()
    print("All matrixes are completed. \n")
    print("Combining control matrixes. \n")
    controlr1basenames = [os.path.basename(path) for path in controlR1]
    controlr1basenames2 = [stemname.split("_")[0] for stemname in controlr1basenames]
    matrixsuffix = '.matrix'
    matrixdir = workdir + '/matrix/'
    controlmatrixlist = [matrixdir + s + matrixsuffix for s in controlr1basenames2]
    allr1controlnames = " ".join(controlmatrixlist)
    print(allr1controlnames)
    concatmats = 'zcat ' + allr1controlnames + ' > ' + matrixdir + 'tempcontrolmat.matrix'
    subprocess.call(concatmats, shell=True)
    sortmats = 'sort -k1,1 -k2,2n ' + matrixdir + 'tempcontrolmat.matrix -o ' + matrixdir + 'tempcontrolmat.matrix'
    subprocess.call(sortmats, shell=True)
    collapsemats = 'perl ' + workdir + 'software/collapse_matrix.pl -i '+ matrixdir + 'tempcontrolmat.matrix -o ' + matrixdir + 'combinedcontrolmatrixes.matrix'
    subprocess.call(collapsemats, shell=True)
    print("Finished combining control matrixes. \n")
else:
    print("Skipping matrix generation due to mode selection")

if mode == 'all' or mode == 'bullseyeonly' or mode == 'findsites':
    addbeddir = 'mkdir -p ' + workdir + '/bedfiles'
    subprocess.call(addbeddir, shell=True)
    print("************************************************************* \n")
    print("Finding edit sites \n")
    bedprocesses = []
    for file in r1testnames:
        bedstemname = os.path.splitext(os.path.basename(file))[0]
        bedstemname2 = bedstemname.split("_")[0]
        matrixin = workdir + '/matrix/' + bedstemname2 + '.matrix.gz'
        controlmatrix = workdir + '/matrix/combinedcontrolmatrixes.matrix.gz'
        beddir = workdir + '/bedfiles/'
        outbedname = beddir + bedstemname2 + '.bed'
        print(matrixin)
        b = Process(target=FindSites, args=(r1testnames, workdir, matrixin, annotfile, outbedname))
        bedprocesses.append(b)
        b.start()
    for b in bedprocesses:
        b.join()
    print("Site finding complete. Bedfiles are in bedfiles directory \n") 
else:
    print("Not finding edit sites due to mode selection.")


			