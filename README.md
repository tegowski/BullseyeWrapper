# BullseyeWrapper
A python-based wrapper for identifying editing events using the Bullseye Pipeline. Please ensure enough memory is allocated for these jobs. STAR alignments are fast but memory intensive and not providing enough memory could slow the progress of the pipeline or lead to OOM errors. Commands below and .sh scripts are written for use within a linux system with SLURM management system. ***Fastq files must be provided as gzipped files with the suffix ".fastq.gz"***

***Default Behavior***

1. Provided fastq files, a STAR genome index, and refFlat annotation file, this will perform genome alignment with STAR, read deduplication (in progress), nucleotide representation matrix generation, and finding editing sites.

2. Optional arguments also allow for the running of FastQC on input fastq files and STAR index generation (See below)

***This is highly preliminary. The current version does not include reduction of UMI or optical duplicates. It may also be very buggy. If you run into any bugs, please raise an issue so it can be addressed***

This is wrapper for running Bullseye, the scripts of which are available on Github (https://github.com/mflamand/Bullseye).

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

