The test data include:

- genotype data: `test_data/all.passonly.subset_n200.vcf.gz`
  - this genotype file is a subset from Pf8 genotype
  - download script: `test_data/download_pf8_genotypes.py`
  - subset script: `test_data/subset_pf8_n200.sh`
  - sample list: `test_data/200_rand_good_samples.txt`
- reduced genotype data used by the test run: `test_data/all.passonly.subset_n6.vcf.gz`
  - 4 clonal (Fws > 0.95) and 2 polyclonal (Fws < 0.8) samples picked from
    the Fws table of a full run on the 200-sample set; sites with no ALT
    allele among the 6 samples are dropped
  - subset script: `test_data/subset_n6.sh`
  - sample list: `test_data/6_selected_samples.txt`
- files used to run the pipeline:
  - population map: `test_data/subpop.txt` and `test_data/6_selected_samples.txt`
  - vcf map: `test_data/vcf.txt` and `test_data/all.passonly.subset_n6.vcf.gz`
