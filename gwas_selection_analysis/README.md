# GWAS and Selection analyses

## association_test.sh
Takes a bed/fam file as input, calculates a relatedness matrix (GEMMA), and runs a linear mixed model (GEMMA) to test for genoetype-phenotype associations, accounting for relatedness between samples.

## genome_selection_analysis.R
Reads in a ped/map file and a csv file with the predictor variable (fly attack rate, i.e., selection), as well as a K value (num. of genetic groups in data) passed when submitting the job. The ped file is converted to LFMM format and missing genotypes are imputed with LEA's impute() function before running a latent factor mixed model (LFMM) to test for selection-associated variants.

## plotGWAS_selection_results.R
This script:
1. Reads in the results files from the above analyses, and identifies overlapping and nearby genes.
2. Plots genetic associations for Fw on Chr1 and Cw on Chr2, based on prior knowledge of their genetic architecture from mapping studies.
3. Plots the strength of association between fly attack rates and genotype across the whole genome and tests for overrepresentation of gene ontology categories (topGO).

Raw results files are too large to upload and will be made available upon manuscript submission.
