#This is the main script to run the Bullseye pipeline
#Required dependencies and instructions for installing the proper conda environment can be found on https://github.com/mflamand/Bullseye
#Must install FastQC (conda install bioconda::fastqc) to environment if not found.
#
#This script was written as a wrapper for the Bullseye pipeline scripts and was written by Matt Tegowski, and is not sanctioned by the original writer of the Bullseye pipeline
#
#Written for Python 3.7.12

import argparse
import os
import sys
import gzip
import multiprocessing
import asyncio
import time
import builtins
import subprocess
from time import strftime
import re
from pathlib import Path
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
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
parser.add_argument("--FastQC", help="Option to perform FastQC, if mode is all", action='store_true')
parser.add_argument("--STARIndex", help="Generation of STAR index, if mode is all", action='store_true')
parser.add_argument("--singleEnd", help="Use if not paired-end sequencing. Assumed data is paired-end by default", action='store_true')
parser.add_argument("--noControl", help="Use if there is no APOBEC1-YTHmut-treated control data", action='store_true')
parser.add_argument("--mode", nargs='?', type=str, help="all, bullseyeonly, findsites, RACfilter, metagene", default='all')
parser.add_argument("--indexDir", type=str, nargs='?', help="If STAR index already exists, specificy the path to its directory", default='GenomeDir/')
parser.add_argument("--memoryGB", type=int, help="Available memory in GB", default=48)
parser.add_argument("--workingDir", type=str, nargs='?', help="Base directory to use", default='./')
parser.add_argument("--matrixCov", type=int, help="Minimum reads of coverage for a nucleotide to be included in matrix generation", default="5")
parser.add_argument("--noAlign", help="Use if you supply fasta files, but already have aligned files in workdir/bamfiles", action='store_true')
parser.add_argument("--bamfiles", type=str, nargs='*', default="Not supplied", help="List of bamfiles (filenames only! Path is assumed to be workdir/bamfiles) for parsing into matrix, only needed if using bullseyeonly mode")
parser.add_argument("--cpuAvail", nargs='?', type=str, help="Number of available processing cores. Default is 1", default="1")

#Arguments for finding edit sites
parser.add_argument("--annotationFile", nargs='?', type=str, help="See Bullseye github for making flatREF annotation file. If not using genome, this file and path MUST be provided. It must not be provided if comparing to genome or knownSites", default="20")
parser.add_argument("--editType", nargs='?', type=str, help="Type of mutation being queried. Default is C2U, but any combination of _2_ nucleotide works, accepts U (treats as T) and I (treats as G)", default='C2U')
parser.add_argument("--minEdit", nargs='?', type=str, help="Minimum mutation percentage. Default is 10%", default="10")
parser.add_argument("--maxEdit", nargs='?', type=str, help="Maximum mutation percentage. Default is 95%", default="95")
parser.add_argument("--foldOverControl", nargs='?', type=str, help="Minimum fold-change in mutation percentage over control at any give site for calling. Default is 1.5", default="1.5")
parser.add_argument("--minCov", nargs='?', type=str, help="Minimum reads of coverage for site calling. If samples are sequenced deeply it is recommended to increase this number as high as reasonably possible.", default="10")
parser.add_argument("--minControlCov", nargs='?', type=str, help="Minimum reads of coverage in control sample for comparisons. Default is 10", default="10")
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

global mode,indexdir,workdir,gendir,genfile,gtffile,testR1,testR2,controlR1,controlR2,r1testnames,r2testnames,r1controlnames,r2controlnames,r1allnames,r2allnames,edittype,mineditpct,maxeditpct,foldovercntrl,mincoverage,mincntrlcoverage,cpus,minedits,intron,utrext,fallback,genomecompare,bdist,knownsites,stranded,filtersites,matcov

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

def Mapping(mapcommand):  
    try:
        print(mapcommand)
        process = subprocess.Popen(mapcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        return process.returncode, stdout, stderr
    except Exception as e:
        return -1, None, str(e).encode()
    
def generateMapping(file1, file2, stem):
    mapcommand = ['STAR --runMode alignReads --outSAMtype BAM Unsorted --readFilesCommand zcat --outFilterMismatchNoverReadLmax 0.06 --outTmpDir ' + workdir + '/tmp/' + stem + ' --runThreadN 8 --genomeDir ' + gendir + ' --outFileNamePrefix ' + workdir + '/bamfiles/' + stem + ' --readFilesIn ' + file1 + ' ' + file2]
    return(mapcommand)

def main_Map():
    stemnames2 = [j.replace(".fastq.gz", "_") for j in r1allnames]
    
    commands = [generateMapping(file1, file2, stem) for file1, file2, stem in zip(r1allpathnames, r2allpathnames, stemnames2)]
    
    #Submit mapping commands concurrently
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(Mapping, cmd) for cmd in commands]

    #Wait for all commands to complete
    for future in concurrent.futures.as_completed(futures):
        returncode, stdout, stderr = future.result()
        print(f"Return code: {returncode}")
        if stderr:
                print(f"Error: {stderr.decode().strip()}")
        if stdout:
                print(f"Output: {stdout.decode().strip()}")
  
def Run_OpticalDups(dedupcommand):
    try:
        print(dedupcommand)
        process = subprocess.Popen(dedupcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        return process.returncode, stdout, stderr
    except Exception as e:
        return -1, None, str(e).encode()

async def generate_deduplication(file1, loop, executor):
    """Run a sequence of 4 jobs for a single file."""
    stemname2 = file1.replace(".fastq.gz", "_")
    
    # Define your commands here
    
    commands = [
        'samtools sort -n -o ' + workdir + '/bamfiles/' + stemname2 + '.name.bam ' + workdir + '/bamfiles/' + stemname2 + 'Aligned.out.bam',
        'samtools fixmate -r -m ' + workdir + '/bamfiles/' + stemname2 + '.name.bam ' + workdir + '/bamfiles/' + stemname2 + '.fix.bam',
        'samtools sort -o ' + workdir + '/bamfiles/' + stemname2 + '.possort.bam ' + workdir + '/bamfiles/' + stemname2 + '.fix.bam',
        'samtools markdup -r ' + workdir + '/bamfiles/' + stemname2 + '.possort.bam ' + workdir + '/bamfiles/' + stemname2 + '.duprm.bam'
    ]
    
    for dedupcommand in commands:
        returncode, stdout, stderr = await loop.run_in_executor(executor, Run_OpticalDups, dedupcommand)
        print(f"Command: {dedupcommand}")
        print(f"Return code: {returncode}")
        if stdout:
            print(f"Output: {stdout.decode().strip()}")
        if stderr:
            print(f"Error: {stderr.decode().strip()}")
    return(dedupcommand, stemname2)
    
async def dedupdeloop(r1allnames):
    loop = asyncio.get_running_loop()
    with ThreadPoolExecutor() as executor:
        tasks = [generate_deduplication(file1, loop, executor) for file1 in r1allnames]
        await asyncio.gather(*tasks)

def Parsing(r1allnames, matrixdir, stem, matcov):
    parsecomm = 'perl ' + workdir + 'software/parseBAM.pl --input ' + bamfilein + ' --output ' + matrixdir + stem + '.matrix --verbose --removeDuplicates --minCoverage ' + str(matcov)
    print(parsecomm)
    subprocess.call(parsecomm, shell=True)
    
def FindSites(r1testnames, stem, matrixin, annotfile, outbedname):
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

#Generate STAR index if desired
if args.STARIndex is True:
    if mode == 'all':
        print("Generating STAR Index.")
        indexing(gendir, genfile, gtffile)
    else:
        print("Skipping Index Generation due to mode selection")
else:
    print("noSTARIndex selected, skipping genome index generation")

#Perform FastQC if desired
if args.FastQC is True:
    if mode == 'all':
        FastQC(testR1, testR2, controlR1, controlR2)
    else:
        print("Skipping FastQC due to mode selection")
else:
    print("noFastQC selected, skipping fastqc analysis")

##################
#Align using STAR#
##################

if args.noAlign is False:
    if mode == 'all':
        removedir = 'rm -r ' + workdir + '/tmp/'
        print(removedir)
        subprocess.call(removedir, shell = True)
        tmpdir = workdir + '/tmp'
        addtmp = 'mkdir -p ' + tmpdir
        subprocess.call(addtmp, shell = True)
        genomeload = 'STAR --genomeLoad LoadAndExit --genomeDir ' + gendir
        print(genomeload)
        subprocess.call(genomeload, shell = True)
        main_Map()
        genomeremove = 'STAR --genomeLoad Remove --genomeDir ' + gendir
        print(genomeremove)
        subprocess.call(genomeremove, shell = True)
    else:
        print("Skipping STAR alignment due to mode selection")
else:
    print("--noAlign selected, skipping alignment, but fastq files must be provided")
    
###########################
#Remove optical duplicates#
###########################

if mode == 'all':
    asyncio.run(dedupdeloop(r1allnames))
else:
    print("Skipping optical deduplication due to mode selection")

###########################################################################################################################
#Make nucleotide representation matrixes. They will be used in comparing %C2U values and calling sites in the next section#
###########################################################################################################################

if mode == 'all' or mode == "bullseyeonly":
    addmatrixdir = 'mkdir -p ' + workdir + '/matrix'
    subprocess.call(addmatrixdir, shell=True)
    matrixdir = workdir + '/matrix/'
    print("*************************************************************")
    print("Parsing bamfiles to generate nucleotide representation matrix \n")
    print(str(matcov))
    matrixprocesses = []
    matstemname2 = [j.replace(".fastq.gz", "_") for j in r1allnames]
    for stem in matstemname2:        
        bamfilein = workdir + '/bamfiles/' + stem + '.duprm.bam'
        print(bamfilein)
        ma = Process(target=Parsing, args=(r1allnames, matrixdir, stem, matcov))
        matrixprocesses.append(ma)
        ma.start()
    # Wait for all processes to complete
    for ma in matrixprocesses:
        ma.join()
    print("All matrixes are completed. \n")
    print("Combining control matrixes. \n")
    controlstems = [j.replace(".fastq.gz", "_") for j in r1controlnames]
    matrixsuffix = '.matrix'
    matrixdir = workdir + '/matrix/'
    controlmatrixlist = [matrixdir + s + matrixsuffix for s in controlstems]
    allr1controlnames = " ".join(controlmatrixlist)
    concatmats = 'zcat ' + allr1controlnames + ' > ' + matrixdir + 'tempcontrolmat.matrix'
    print(concatmats)
    subprocess.call(concatmats, shell=True)
    sortmats = 'sort -k1,1 -k2,2n ' + matrixdir + 'tempcontrolmat.matrix -o ' + matrixdir + 'tempcontrolmat.matrix'
    print(sortmats)
    subprocess.call(sortmats, shell=True)
    collapsemats = 'perl ' + workdir + 'software/collapse_matrix.pl -i '+ matrixdir + 'tempcontrolmat.matrix -o ' + matrixdir + 'combinedcontrolmatrixes.matrix'
    print(collapsemats)
    subprocess.call(collapsemats, shell=True)
    print("Finished combining control matrixes. \n")
else:
    print("Skipping matrix generation due to mode selection")

#Here, edited sites will be identified
if mode == 'all' or mode == 'bullseyeonly' or mode == 'findsites':
    addbeddir = 'mkdir -p ' + workdir + '/bedfiles'
    subprocess.call(addbeddir, shell=True)
    print("************************************************************* \n")
    print("Finding edit sites \n")
    bedstemnames = [j.replace(".fastq.gz", "_") for j in r1testnames]
    bedprocesses = []
    for stem in bedstemnames:
        matrixin = workdir + '/matrix/' + stem + '.matrix.gz'
        controlmatrix = workdir + '/matrix/combinedcontrolmatrixes.matrix.gz'
        beddir = workdir + '/bedfiles/'
        outbedname = beddir + stem + '.bed'
        print(matrixin)
        b = Process(target=FindSites, args=(r1testnames, stem, matrixin, annotfile, outbedname))
        bedprocesses.append(b)
        b.start()
    for b in bedprocesses:
        b.join()
    print("Site finding complete. Bedfiles are in bedfiles directory \n") 
else:
    print("Not finding edit sites due to mode selection.")