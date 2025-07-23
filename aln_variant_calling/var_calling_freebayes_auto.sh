#!/bin/bash
#SBATCH --job-name=freebayes_auto  
#SBATCH --cpus-per-task=16   
#SBATCH --mem=16gb                
#SBATCH --partition=long
#SBATCH --array=2-14
#SBATCH --output=freebayes1.3.6_array_%A_%a.log
pwd; hostname; date
source activate freebayes1.3.6

REF=~/scratch/PopGen/genomev3/TOC.asm.scaffold.fasta
var_vcf=~/scratch/PopGen/SNP_detection/alignments/freebayes/variants_chr1_to_14_variant.recode.vcf.gz
#-@ provides vcf to include variant calls for, -l specifies to ONLY include these sites in the output

freebayes -@ $var_vcf -l -f $REF -r "scaffold_${SLURM_ARRAY_TASK_ID}" --report-genotype-likelihood-max \
 --no-population-priors --use-best-n-alleles 4 \
 --hwe-priors-off --use-mapping-quality --theta 0.02 --haplotype-length -1 --genotype-qualities --ploidy 2 \
 --bam-list bam_files.txt > "Variants_Chr${SLURM_ARRAY_TASK_ID}.vcf"
