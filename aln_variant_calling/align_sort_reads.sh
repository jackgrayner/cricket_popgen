#!/bin/bash
#SBATCH --job-name=aln_sort_rmdup
#SBATCH --export=ALL
#SBATCH -N1 --ntasks-per-node=16
#SBATCH --mem=16G
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jr228@st-andrews.ac.uk
#SBATCH --array=1-5
#SBATCH --output=aln_array_%A_%a.out

#this will run an array of high memory/CPU jobs each aligning a list of files 

while read line
do
	#extract sample name
	sample_name=$(echo "$line" | sed 's/_R1_cleaned.fastq.gz//')
	cat "${sample_name}" >> "samples_processed${SLURM_ARRAY_TASK_ID}"

	genome_index="/home/jrayner/scratch/PopGen/genomev3/TOC.asm.scaffold.fasta"  
	input_fastq1="./${sample_name}_R1_cleaned.fastq.gz"
	input_fastq2="./${sample_name}_R2_cleaned.fastq.gz"
	output_bam="${sample_name}_sorted.bam"

	#add read group info
	library_name="${sample_name}"
	platform="ILLUMINA"
	unit="unit1"
	rg_string="@RG\tID:${sample_name}\tSM:${sample_name}\tLB:${library_name}\tPL:${platform}\tPU:${unit}"

	#align reads and sort
	source activate bwa-mem2
	bwa mem -t 16 -R "$rg_string" "$genome_index" "$input_fastq1" "$input_fastq2" | samtools sort -o "$output_bam"
	samtools index "$output_bam"

	#mark duplicated and re-index
	source activate sambamba
	sambamba markdup -r "$output_bam" "${sample_name}.rmdup.bam"
	samtools index "${sample_name}.rmdup.bam"
	#rm "$output_bam"
done < "files${SLURM_ARRAY_TASK_ID}"
