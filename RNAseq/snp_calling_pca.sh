#!/bin/bash
#SBATCH --job-name=bcftools   # Job name
#SBATCH --array=1-14                      # Array job for chromosomes 1 to 14
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=bcftools.out   # Standard output and error log
pwd; hostname; date

ref=~/scratch/Toc_genome_v3/TOC.asm.scaffold.fasta
bamlist=bamlist.txt

#use bcftools to call snps
source activate bcftools
bcftools mpileup -d 10000 -f $ref -r $chr -b $bamlist | bcftools call -mv -Oz -f GQ -o Toc.RNA.bcftools.vcf.gz

#perform basic filtering with vcftools
source activate vcftools
vcftools --gzvcf Toc.RNA.bcftools.vcf.gz --maf 0.05  --minGQ 20 --minQ 30 \ 
 --min-alleles 2 --max-alleles 2 --max-missing 0.75 --recode --recode-INFO-all \ 
 --stdout | bgzip -c > ./Toc.RNA.filtered.bcftools.vcf.gz
 
bcftools index Toc.all.RNA.bcftools.vcf.gz

#make bed file
plink --allow-extra-chr --vcf Toc.RNA.filtered.bcftools.vcf.gz \
 --double-id --make-bed --snps-only --allow-no-sex --set-missing-var-ids @:# \
 --out Toc.RNA.filtered

#run PCA for each chromosome
for chr1 in {1..14}
do
	chr="scaffold_${chr1}"
	plink --bfile Toc.RNA.filtered --double-id --allow-extra-chr --set-missing-var-ids @:# \
	--chr $chr --make-bed --pca --out "Toc_PCA_${chr}"
done
