#!/bin/bash
#SBATCH --job-name=hisat_stringtie
#SBATCH --export=ALL
#SBATCH --mem=32G #larger files will need more memory
#SBATCH --cpus-per-task=4
#SBATCH --partition=himem
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jr228@st-andrews.ac.uk
echo "processing ${1}"
echo $1

#clean reads
source activate fastp
fastp -i ./fastq_files/${1}_1.fastq.gz -I ./fastq_files/${1}_2.fastq.gz -o ./cleaned/${1}.R1.fastq.gz -O ./cleaned/${1}.R2.fastq.gz

#align to genome with hisat2
source activate hisat2 
hisat2 --rna-strandness RF --dta -x ~/scratch/Toc_genome_v3/TOC.asm.scaffold -1 ./fastq_files/${1}_1.fastq.gz -2 ./fastq_files/${1}_2.fastq.gz | samtools sort -o ./${1}.bam
samtools index $1.bam

#quantify expression of genes with stringtie
stringtie -B -p 4 -G ~/projects/uosa/Nathan_Bailey/Teleogryllus_oceanicus_v3_resources/03.protein_coding_genes/TOC.asm.scaffold.gene.gff3 -o ballgown/$1/stringtie.$1.gtf $1.bam

