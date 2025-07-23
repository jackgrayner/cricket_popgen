#!/bin/bash
#SBATCH --job-name=prunebed
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=long
#SBATCH --output=prunebed.log

source activate plink
vcf_name=$1

plink --allow-extra-chr --vcf ${vcf_name}.vcf.gz \
	--double-id --indep-pairwise 50 10 0.1 --memory 32000 --allow-no-sex \
	--out ${vcf_name}_outgroup_nodupes_pruned
	
plink --vcf ${vcf_name}.vcf.gz --memory 32000 \
	--double-id --allow-extra-chr  \
	--extract allpops_outgroup_nodupes_pruned.prune.in \
	--make-bed --out ${vcf_name}_pruned


plink --vcf ${vcf_name}.vcf.gz --memory 32000 --chr scaffold_5,scaffold_9,scaffold_14 \
	--double-id --allow-extra-chr  \
	--extract allpops_outgroup_nodupes_pruned.prune.in \
	--make-bed --out ${vcf_name}_chr5_9_14_pruned
