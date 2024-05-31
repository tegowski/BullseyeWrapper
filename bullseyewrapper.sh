#!/bin/bash
#SBATCH --mem=164G
#SBATCH -c 8
#SBATCH -p scavenger

fastqdir="/path/to/fastq/files"

python BullseyeWrapper.py \
	--mode all \
	--memoryGB 164 \
	--fqR1Test $fastqdir/APOYTH1_R1.fastq.gz $fastqdir/APOYTH2_R1.fastq.gz $fastqdir/APOYTH3_R1.fastq.gz $fastqdir/APOYTH4_R1.fastq.gz \
	--fqR2Test $fastqdir/APOYTH1_R2.fastq.gz $fastqdir/APOYTH2_R2.fastq.gz $fastqdir/APOYTH3_R2.fastq.gz $fastqdir/APOYTH4_R2.fastq.gz \
	--fqR1Control $fastqdir/APOYTHmut1_R1.fastq.gz $fastqdir/APOYTHmut2_R1.fastq.gz \
	--fqR2Control $fastqdir/APOYTHmut1_R2.fastq.gz $fastqdir/APOYTHmut2_R2.fastq.gz \
	--genomeDir /path/to/genome/index \
	--genomeFile /path/to/genome/species_releaseX.fasta \
	--STARgtf /path/to/gtf/species_releaseX.gtf \
	--workingDir /path/to/working/directory/ \
	--annotationFile /path/to/gtf/species_releaseX.refFlat \
	--editType C2U \
	--cpuAvail 8
