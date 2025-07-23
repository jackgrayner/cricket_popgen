#!/bin/bash
#SBATCH --job-name=vcffilt
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=long
#SBATCH --output=bcfmerge_filt.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jr228@st-andrews.ac.uk

vcf_name=$1

source activate bcftools
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ./${vcf_name}.vcf -O z -o ${vcf_name}.anno.vcf
bgzip ${vcf_name}.anno.vcf
bcftools index ${vcf_name}.anno.vcf.gz

source activate vcftools
vcftools --gzvcf ${vcf_name}.anno.vcf.gz --maf 0.05 --minDP 5 --maxDP 60 --minGQ 20 --minQ 30 --min-alleles 2 \
 --max-alleles 2 --max-missing 0.75 --recode --recode-INFO-all --stdout | bgzip -c > ${vcf_name}_filtered.vcf.gz

source activate plink
plink --vcf ${vcf_name}_filtered.vcf.gz \
 --double-id --memory 32000 --make-bed --allow-extra-chr --allow-no-sex \
 --out ${vcf_name}_filtered
