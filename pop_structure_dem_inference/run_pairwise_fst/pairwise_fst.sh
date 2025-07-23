#!/bin/bash
#SBATCH --job-name=populations_fst
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jrayner1@umd.edu

source activate stacks
populations -V allpops_outgroup_filt_pruned_maxmiss0.9_scaff5_9_14.vcf.gz -M popmap.tsv -O ./populations --fstats
