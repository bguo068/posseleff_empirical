#! /usr/bin/env bash

eval "$(/local/data/Malaria/Projects/Takala-Harrison/AFRIMS/miniconda3/bin/conda shell.bash hook)" 
conda activate snpeff

# test on biallelic snps
../annotate.py --vcf biallelic_ex.vcf.gz --out_prefix out_biallelic

# test on multiallelic snps
bcftools norm -m -any multiallelic_ex.vcf.gz -Oz -o multiallelic_norm.vcf.gz
../annotate.py --vcf multiallelic_norm.vcf.gz --out_prefix out_multiallelic