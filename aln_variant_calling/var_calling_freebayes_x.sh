#!/bin/bash
#SBATCH --job-name=OACfreebayes_x
#SBATCH --cpus-per-task=16 
#SBATCH --mem=16gb
#SBATCH --partition=long
#SBATCH --output=freebayes1.3.6_x.log
pwd; hostname; date
source activate freebayes1.3.6

REF=~/scratch/PopGen/genomev3/TOC.asm.scaffold.fasta
var_vcf=~/scratch/PopGen/SNP_detection/alignments/freebayes/variants_chr1_to_14_variant.recode.vcf.gz

#-@ provides vcf to include variant calls for, -l specifies to ONLY include these sites in the output
freebayes -@ $var_vcf -l -f $REF -r scaffold_1 --report-genotype-likelihood-max --no-population-priors --use-best-n-alleles 4 \
 --hwe-priors-off --use-mapping-quality --theta 0.02 --haplotype-length -1 --genotype-qualities --ploidy 1 \
 --bam-list bam_files.txt > Variants_Chr1.vcf
