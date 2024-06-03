# BullseyeWrapper
A python-based wrapper for identifying editing events using the Bullseye Pipeline. Please ensure enough memory is allocated for these jobs. STAR alignments are fast but memory intensive and not providing enough memory could slow the progress of the pipeline or lead to OOM errors. Commands below and .sh scripts are written for use within a linux system with SLURM management system. ***Fastq files must be provided as gzipped files with the suffix ".fastq.gz"***

***Default Behavior***

1. Provided fastq files, a STAR genome index, and refFlat annotation file, this will perform genome alignment with STAR, read deduplication (in progress), nucleotide representation matrix generation, and finding editing sites.

2. Optional arguments also allow for the running of FastQC on input fastq files and STAR index generation (See below)

***This is highly preliminary. The current version does not include reduction of UMI or optical duplicates. It may also be very buggy. If you run into any bugs, please raise an issue so it can be addressed***

This is wrapper for running Bullseye, the scripts of which are available on Github (https://github.com/mflamand/Bullseye). Please see the github for the Bullseye pipeline for detailed information on the pipeline and how it works.

# 1. Setup

It is highly recommended to run this wrapper within a conda environment designed for Bullseye. Detailed instructions for setting up this conda environment are also located at https://github.com/mflamand/Bullseye.

After the conda environment is setup, activate it
```
conda activate Bullseye
```
Then, you should navigate to the working directory (abbreviated here as '/path/to/workdir/') and make a directory called software and enter it
```
cd /path/to/workdir
mkdir software
cd software
```
Next, download the neccessary Bullseye scripts and the BullseyeWrapper.py
```
wget 'https://github.com/mflamand/Bullseye/blob/main/Code/Find_edit_site.pl'
wget 'https://github.com/mflamand/Bullseye/blob/main/Code/parseBAM.pl'
wget 'https://github.com/mflamand/Bullseye/blob/main/Code/scripts/collapse_matrix.pl'
wget 'https://github.com/tegowski/BullseyeWrapper/blob/main/BullseyeWrapper.py'
```

***Optional*** If you want to run FastQC, it will need to be installed within the Bullseye conda environment
```
conda install fastqc=0.12.1-0
```

# 2. Running BullseyeWrapper

The wrapper is simply run as one command from within the $workdir/software directory. See the bullseyewrapper.sh file for basic arguments.
```
sbatch bullseyewrapper.sh
```

# 3. Arguments for BullseyeWrapper

**Required Arguments**

--fqR1Test [Provide the path to and name of read 1 gzipped fastq files in which editing is to be measured. They **must** end with "fastq.gz"]

--fqR2Test [Provide the path to and name of read 2 gzipped fastq files in which editing is to be measured. They **must** end with "fastq.gz". This can be skipped only if the --singleEnd option is also used.]

--fqR1Control [Provide the path to and name of read 1 gzipped control (background, or YTHmut) fastq files to which editing is to be compared. They **must** end with "fastq.gz". This can be skipped only if the --noControl option is also used.]

--fqR2Control [Provide the path to and name of read 2 gzipped control (background, or YTHmut) fastq files to which editing is to be compared. They **must** end with "fastq.gz". This can be skipped only if the --noControl and/or --singleEnd options are also used.]

--mode [Mode to run the pipeline in. Can be "all", "bullseyeonly", or "findsites". "all" runs the whole pipeline (no genome generation or fastqc unless specified). "bullseyeonly" will run the matrix generation and find edit sites steps. "findsites" will only run the find edit sites step]

--workingDir [path to working directory for the project. Should be one level up from "software" directory containing scripts]


**Other important arguments**

--memoryGB [provide available memory in GB. Failing to set this may slow down the pipeline]

--cpuAvail [provide number of cpus available. Failing to set this may slow down the pipeline]


**Arguments for alignment with STAR**

--genomeDir [provide path to the STAR genome index]


**Arguments for matrix generation**

--matrixCov [Minimum read coverage for any nucleotide to be included in matrix. Default is 3.]


**Arguments for finding editing sites and setting thresholds**

*Required*

--annotationFile [Path to and file name of refFlat annotation file. See (https://github.com/mflamand/Bullseye) for detailed information on generating this file]

--editType [Type of mutation being queried. Default is C2U, but any combination of _2_ nucleotide works, accepts U (treats as T) and I (treats as G) default is C2U.]

*Optional/alternate arguments*

--minEdit [Minimum mutation percentage. Default is 10%]

--maxEdit [Maximum mutation percentage. Default is 95%]

--foldOverControl [Minimum fold-change in mutation percentage over control at any give site for calling. Default is 1.5]

--minCov [Minimum reads of coverage for site calling. **It is recommended to increase this number as high as reasonably possible.** Default is 10, which is low.]

--minControlCov [Minimum reads of coverage in control sample for comparisons. Default is 10.]

--minEdits [Minimum number of mutations observed at a given site for it to be called. It is highly recommended to never use 1. Default is 2. Higher than 2 is recommended is possible, but too high of a number may lead to false negatives for real sites with low editing. Default is 2.]

--fallback [Use this option if you want to consider sites with low/no coverage in control file. Must supply genomeFile if used.]

--genome [Do not use control matrix for comparison. Just use other thresholds set and compare to the genome instead. Provide argument alone, not with file. genome fasta file must be supplied with --genomeFile option.]

--score [score and filter reads based on probability of edit based on beta distribution curve. Sites with 10 fold higher probability in DART/control will be kept.]

--Intron [Use to find sites within intronic sequences, otherwise only exons will be used.]

--UTRextend [Extend site calling 5kb from annotated end of 3'UTR.]

--knownSites [Instead of finding all sites, only find sites contained within bed file provided with path. Must be bed-6. Must not use --annotationFile, --Intron, or --UTRextend options.]

--stranded [Use this option if you used a stranded sequencing kit.]

--filterBed [Do not include sites contained within bedfile. Must be bed-6. Provide path and filename.]


***Optional*** **Arguments for FastQC**
--FastQC [Use this argument to perform FastQC analysis on all input fastq files]

***Optional*** **Arguments for STAR Genome generation**
--STARIndex [Use this argument to generate the genome index]
--genomeDir [provide path to where the STAR genome index will be]
--genomeFile [provide path to and file name of genome fasta file]
--STARgtf [provide path to and file name of gtf file for generating index]
