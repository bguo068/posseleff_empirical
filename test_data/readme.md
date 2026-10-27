The test data include:

- genotype data: `test_data/all.passonly.subset_n200.vcf.gz`
  - this genotype file is a subset from Pf8 genotype
  - download script: `test_data/download_pf8_genotypes.py`
  - subset script: `test_data/subset_pf8_n200.sh`
  - sample list: `test_data/200_rand_good_samples.txt`
- reduced genotype data used by the test run: `test_data/all.passonly.subset_n36.vcf.gz`
  - 30 clonal (Fws > 0.95) samples, 4 polyclonal samples with Fws 0.93-0.94
    (likely to carry a dominant clone) and 2 polyclonal samples with
    Fws < 0.8, picked from the Fws table of a full run on the 200-sample set;
    sites with no ALT allele among the 36 samples are dropped
  - 30 clonal samples are enough for the iHS scan in the IBD analyses
  - subset script: `test_data/subset_n36.sh`
  - sample list: `test_data/36_selected_samples.txt`
- files used to run the pipeline:
  - population map: `test_data/subpop.txt` and `test_data/36_selected_samples.txt`
  - vcf map: `test_data/vcf.txt` and `test_data/all.passonly.subset_n36.vcf.gz`
