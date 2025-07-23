#!/bin/bash
#SBATCH --job-name=tabix_vcf2fasta
#SBATCH --export=ALL
#SBATCH -N1 --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jrayner1@umd.edu

vcf_name=$1
fasta_file=~scratch/Toc_genome_v3/TOC.asm.scaffold.fasta
gff_file=/home/jrayner/scratch/Toc_genome_v3/flattened.gff3

./vcf2fasta/vcf2fasta.py --fasta ${fasta_file} \
	--vcf ${vcf_name} --gff ${gff_file} \
	--feat CDS --skip --out vcf2fasta_CDS


