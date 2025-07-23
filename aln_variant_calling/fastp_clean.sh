#!/bin/bash
#SBATCH --job-name=clean
#SBATCH --mem=5gb
#SBATCH --partition=long
#SBATCH --cpus-per-task=2
echo 'cleaning' ${sample}
source activate fastp
sample=$1
fastp -i ${sample}_pass_1.fastq.gz -I ${sample}_pass_2.fastq.gz \
 -o ${sample}_R1_cleaned.fastq.gz -O ${sample}_R2_cleaned.fastq.gz

