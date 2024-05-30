# BullseyeWrapper
A python-based wrapper for identifying editing events using the Bullseye Pipeline
***This is highly preliminary. The current version does not include reduction of UMI or optical duplicates. It may also be very buggy. If you run into any bugs, please raise an issue so it can be addressed***

This is wrapper for running Bullseye, the scripts of which are available on Github (https://github.com/mflamand/Bullseye).

1. Setup
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

