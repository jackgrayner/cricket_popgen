#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --export=ALL
#SBATCH -N1 --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --partition=long
#SBATCH --output=makebed_cwgwas.log

bed_name=$1

source activate GEMMA
gemma -bfile ${bed_name} -gk 1 -o ${bed_name}
gemma -bfile ${bed_name} -k ./output/${bed_name}.cXX.txt -lmm 2 -o ${bed_name}