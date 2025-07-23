#!/bin/bash
#SBATCH --job-name=mafft_fasttree
#SBATCH --export=ALL
#SBATCH -N1 --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --partition=long

source activate mafft
source activate seqtk

fasta_dir="./vcf2fasta_CDS/" 

# run mafft on each fasta file 
for fasta_file in ${fasta_dir}*.fas; do
    base_name=$(basename "$fasta_file" .fas)
    mafft --auto "$fasta_file" > "${fasta_dir}${base_name}_aligned.fasta"
done

#concatenate files
seqkit concat *aligned.fasta > combined_alignment.fasta

#make tree
source activate fasttree
fasttree -gtr -nt combined_alignment.fasta > tree_file 
