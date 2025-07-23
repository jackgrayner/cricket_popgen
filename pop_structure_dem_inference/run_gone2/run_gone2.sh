#!/bin/bash
#SBATCH --job-name=gone2
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=long
source activate plink
bed_file=$1
pop=$2
plink --bfile ${bed_file} --chr scaffold_5,scaffold_9_scaffold_14 \ 
	--snps-only --geno 0.1 --mind 0.1 --recode --keep ${pop}_ids --allow-extra-chr \ 
	--out ${pop}_filtered_chr5_9_14
./GONE2/gone2 -g 3 -r 1.1 -b 0.001 -s 100000 -o ${pop}_all -t 8 ${pop}_filtered_chr5_9_14.ped
