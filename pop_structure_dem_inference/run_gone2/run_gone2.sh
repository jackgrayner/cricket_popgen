#!/bin/bash
#SBATCH --job-name=gone2
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=long

pop=$1

#create pop bed file 
plink --bfile all_filtered_chr5_9_14 --chr scaffold_5 scaffold_9 scaffold_14 --thin-count 100000 --recode --keep ${pop} --allow-extra-chr --out ${pop}_filtered_chr5_9_14

#then, the full analysis (1M SNPs)
~/scratch/popgen/gone2/feb26/GONE2/gone2 -s 1000000 -g 0 -r 1.5 -o ./results_all/${pop}_r1.5_full -t 8 ${pop}_filtered_chr5_9_14.ped

#then, the subsampled analysis
for n in {1..20}
do
        plink --file ${pop}_filtered_chr5_9_14 --chr scaffold_5 scaffold_9 scaffold_14 --thin-count 100000 --thin-indiv-count 20 --recode --keep ${pop} --allow-extra-chr --out ./results_subset10/${pop}_filtered_chr5_9_14_n${n}
        ~/scratch/popgen/gone2/feb26/GONE2/gone2 -s 100000 -g 0 -r 1.5 -o ./results_subset10/${pop}_r1.5_n${n} -t 8 ./results_subset10/${pop}_filtered_chr5_9_14_n${n}.ped
done


