The test data include:

- genotype data: `test_data/all.passonly.subset_n200.vcf.gz`
  - this genotype file is a subset from Pf8 genotype
  - download script: `test_data/download_pf8_genotypes.py`
  - subset script: `test_data/subset_pf8_n200.sh`
  - sample list: `test_data/200_rand_good_samples.txt`
- files used to run the pipeline:
  - population map: `test_data/subpop.txt` and `test_data/200_rand_good_samples.txt`
  - vcf map: `test_data/subpop.txt` and `test_data/all.passonly.subset_n200.vcf.gz`
