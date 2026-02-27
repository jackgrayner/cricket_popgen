
#first call variants in RAD-seq data
REF=TOC.asm.scaffold.fasta
var_vcf=./allpops_alloutgroup.vcf.gz
#-@ provides vcf to include variant calls for, -l specifies to ONLY include these sites in the output
source activate freebayes1.3.6
freebayes -@ $var_vcf -l -f $REF -r "scaffold_${SLURM_ARRAY_TASK_ID}" --report-genotype-likelihood-max \
 --no-population-priors --use-best-n-alleles 4 \
 --hwe-priors-off --use-mapping-quality --theta 0.02 --haplotype-length -1 --genotype-qualities --ploidy 2 \
 --bam-list bamlist.txt > "Variants_Chr${SLURM_ARRAY_TASK_ID}.vcf"

source activate bcftools
bcftools concat Variants_Chr1.vcf.gz Variants_Chr2.vcf.gz Variants_Chr3.vcf.gz Variants_Chr4.vcf.gz Variants_Chr5.vcf.gz Variants_Chr6.vcf.gz Variants_Chr7.vcf.gz Variants_Chr8.vcf.gz Varia>
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ./Variants_all_chr.vcf.gz -O z -o Variants_all_chr_anno.vcf.gz
bcftools index Variants_all_chr_anno.vcf.gz

source activate vcftools
vcftools --gzvcf Variants_all_chr_anno.vcf.gz  --maf 0.1 --min-alleles 2 --max-alleles 2 --max-missing 0.1 --minDP 5 --minGQ 20 --recode --recode-INFO-all --stdout | bgzip -c > Variants_all>

source activate bcftools
bcftools index Variants_all_chr_anno_filt.vcf.gz
bcftools merge Variants_all_chr_anno_filt.vcf.gz allpops_alloutgroup.vcf.gz -Oz -o RAD_WGS_merged.vcf.gz
bcftools index RAD_WGS_merged.vcf.gz

source activate vcftools
vcftools --gzvcf RAD_WGS_merged.vcf.gz  --maf 0.05 --min-alleles 2 --max-alleles 2 --max-missing 0.5 --minDP 5 --minGQ 20 --recode --recode-INFO-all --stdout | bgzip -c > RAD_WGS_merged_filt.vcf.gz

cd ./pca
source activate plink2
plink2 --allow-extra-chr --vcf ../RAD_WGS_merged_filt.vcf.gz \
 --double-id --geno 0.7 --max-alleles 2 --make-bed --maf 0.05 --snps-only \
 --out RAD_WGS_merged

for chr1 in {1..14}
do
        chr=scaffold_${chr1}
        plink2 --allow-extra-chr --remove exclude.txt --bfile RAD_WGS_merged \
                 --chr ${chr} --make-bed --mind 0.95 --snps-only --pca \
                 --out merged_rad_${chr}
done

cd ../winpca

for chr1 in {1..14}
do
	chr=scaffold_${chr1}
	end=$(grep -w ${chr} ~/scratch/Toc_genome_v3/TOC.asm.scaffold.fasta.fai | awk '{print $2}')
	ref_sample=$(head -n ${chr1} ref_samples.txt| tail -n 1 )
	~/scratch/RAD_Toc_Aus/genomev3/reads/aln/winpca/winpca/winpca pca -s keep_merged_samples.txt -p guide_samples -g ${ref_sample} -w 10000000 -i 100000 --np -v GT merged_${chr} ../RAD_WGS_merged_filt.vcf.gz $chr:1-$end
	~/scratch/RAD_Toc_Aus/genomev3/reads/aln/winpca/winpca/winpca chromplot -m pop_aus_HI.txt -g population merged_$chr $chr:1-$end -i 10
done
