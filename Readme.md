# `poseleff_empirical`

The code in this repository analyzes empirical WGS (Whole Genome Sequencing) data
to assess the effect of positive selection on IBD (Identity by Descent)-based
inference. This analysis process encompasses variant and sample filtering,
multiplicity of infection inference, deconvolution and imputation (for both
monoclonal and a subset of polyclonal samples), IBD calling, and IBD-based
inferences of effective population size (Ne) and population structure. These
functions are detailed and implemented across five modules, as described in the
"Pipeline outline" section below.

# System requirement and Software evironment

The pipeline has been tested on Linux Operation system and can be easily adapted to MacOS with
simple changes. Software dependencies and the version numbers are specified in the
'./env.yaml' Conda recipe. Additional depencies that are not available from Conda are
specified in the installation instruction below. The overall installation time
is about 5-30 minutes. The process is relatively long due to the
compliation of the R package 'moimix' and its dependencies.

1. Install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
2. Install [Conda](https://docs.conda.io/en/latest/miniconda.html)
3. Create a Conda evironment from a provided recipe: `conda env create -f ./env.ymal`
4. Activate conda environment: `conda activate empirical`
5. Install `moimix` package: Run `R` and in R console run:
```r
install.packages("BiocManager")
BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)
```
4. Modify (recombination rate) and install `hmmibd`:
```sh
    git clone https://github.com/glipsnort/hmmIBD.git; cd hmmIBD
    # checkout a specific version for reproducibility
    git checkout a2f796ef8122d7f6b983ae9ac4c6fba35afcd3aa
    sed -i -e 's/const double rec_rate = 7.4e-7/const double rec_rate = 6.67e-7/' hmmIBD.c
    # use the complier from the `empirical` conda environtment
    x86_64-conda_cos6-linux-gnu-gcc -o hmmIBD -O3 -lm -Wall hmmIBD.c
    cd ..
    cp hmmIBD/hmmIBD bin/hmmIBD
    rm -rf hmmIBD
```


# Pipeline outline

The entire pipeline is organized into five modules, which include:
- `01_filter_vcf.nf`: Filters VCF files by per-sample and per-site missingness
and minor allele frequency.
    - This involves two rounds of filtering to help retain more samples/sites.
- `02_calc_pop_fws.nf`: Calculates `Fws` statistics using the `moimix` R package.
    - Data is organized into populations (Pfv6 classification).
    - Within each population, `Fws` is inferred for individual samples.
- `03_phase_impute_monoclonal.nf`:
    - Filters samples and retains only monoclonal samples.
    - Deconvolutes each monoclonal sample (dEploid). 
    - Imputes phased genotype data within each population (Beagle).
    - Merges single-sample genotypes back to multiple-sample, phased & imputed
    VCF files.
- `04_phase_impute_polyclonal.nf`:
    - Uses a subset of phased & imputed monoclonal samples to build a
    deconvolution panel for polyclonal samples within each population.
    - Filters samples and retains only polyclonal samples.
    - Deconvolutes genotype data for each polyclonal sample using a
    population-specific panel (dEploid-IBD).
    - Filters polyclonal samples, retaining only those with a dominant clone/genome.
    - Combines phased polyclonal samples with the above monoclonal samples.
    - Imputes the combined genotype data (monoclonal/polyclonal samples) within
    each population (Beagle).
- `05_ibd_ne_ifm.nf main.nf`:
    - Calls IBD using `hmmIBD`.
    - Processes IBD to prepare input (with or without selection correction) for
    IBDNe analysis.
    - Processes IBD to prepare input  (with or without selection correction) for
    Infomap analysis.
    - Runs IBDNE to infer effective population size.
    - Runs Infomap to infer population structure.

The five modules will run sequentially if using `main.nf` as the the entry point.


## General notes:
- Thresholds can be configured in the top lines of each module file.
- The pipeline is designed to run in parallel for large datasets.
- Each module is independent of the others, provided that valid inputs are supplied.

# How to use the pipeline

After install the pipeline, first prepare the input files and then run the
nextflow pipeline with the input files.

## Input


### **1. Vcf file table**: 

The VCF file table is organized in a two-column format:
- Header: None
- Delimiter: Tab
- Column 1: Represents the genome interval
  - Format for genome interval: `[Chromosome]:[Start]-[End]`
- Column 2: Provides the full path to the joint-call VCF file corresponding to
the specified genome interval
  - It's essential that VCF files for different genome intervals include the same set of samples.
  - The VCF file should contain annotations for the hard filter/vqsr filter
  within the FILTER columns. Refer to the output from the
  [`snp_call_nf`](https://github.com/bguo068/snp_call_nf/tree/main) pipeline for
  more information.
  - For those interested in reproducing the VCF files via `snp_call_nf`
  pipeline, you can find the accession numbers of WGS data for the samples used
  in this study within the `samples` folder.


Example file `interval_vcf_table.tsv`
```
Pf3D7_01_v3:1-459121	    /full/path/to/filePf3D7_01_v3:1-459121.vqsrfilt.SNP.vcf.gz
Pf3D7_01_v3:459122-640851   /full/path/to/filePf3D7_01_v3:459122-640851.vqsrfilt.SNP.vcf.gz
Pf3D7_02_v3:1-448875	    /full/path/to/filePf3D7_02_v3:1-448875.vqsrfilt.SNP.vcf.gz
Pf3D7_02_v3:448876-947102   /full/path/to/filePf3D7_02_v3:448876-947102.vqsrfilt.SNP.vcf.gz
...
Pf3D7_14_v3:2121871-2706733 /full/path/to/filePf3D7_14_v3:2121871-2706733.vqsrfilt.SNP.vcf.gz
Pf3D7_14_v3:2706734-3291936 /full/path/to/filePf3D7_14_v3:2706734-3291936.vqsrfilt.SNP.vcf.gz
Pf3D7_14_v3:568025-1073306  /full/path/to/filePf3D7_14_v3:568025-1073306.vqsrfilt.SNP.vcf.tz
```
If the VCF file is not in 'gz' format, use `bgzip xx.vcf` to convert it to `*.vcf.gz` format.


### **2. Population table**
The population table is a two-column table
- Header: none
- Delimiter: tab
- The first column contains the population name
- The second column specifies the path to a file that contains all samples from the given population
    - In the list file, each sample name appears on a separate line
    - The sample names should be the same as in the input VCF files
    - All samples in the VCF files should be contained in one of the population sample list files

Example file `subpops_table.tsv`

```
CAF     subpops/CAF.txt
EAF	    subpops/EAF.txt
ESEA    subpops/ESEA.txt
Lab	    subpops/Lab.txt
OCE	    subpops/OCE.txt
SAM	    subpops/SAM.txt
SAS	    subpops/SAS.txt
WAF	    subpops/WAF.txt
WSEA    subpops/WSEA.txt
```
Example file `subpops/ESEA.txt`
```
Sample1
Sample2
Sample3
...
Samplen
```
The same format is shared maong all populaton sample files such as
`subpops/EAF.txt` (East Africa), `subpops/ESEA.txt` (Eastern Southeast Asia),
`subpops/Lab.txt` (Laboratory Strains), `subpops/OCE.txt` (Oceania),
`subpops/SAM.txt` (South America), `subpops/SAS.txt` (South Asia),
`subpops/WAF.txt` (West Africa), and `subpops/WSEA.txt` (Western Southeast
Asia).


## Run the pipeline
```sh
conda activate empirical
nextflow ./main.nf -profile -resume -profile sge \
    --input_vcf_fn ./interval_vcf_table.tsv \
    --input_subpop_fn ./subpops_table.tsv \
    --outdir results 
```

For large datasets, using a cluster such as SGE is recommended. An example `sge`
profile is provided in the `nextflow.config` file and should be adjusted to fit
your cluster system.

If run on a local computer, please remove the `-profile sge` option from the
above command.

## Test data

Test data for module 1-4 are not directly provided as they are extremely large but 
can be easily generated by running the `ena_test` test in the [`snp_call_nf` pipeline']
(https://github.com/bguo068/snp_call_nf).

A small dataset is provided, however, for module 5. The test of module 5 can be
run with the following command:
```
nextflow 05_ibd_ne_ifm.nf -profile sge --vcf 'test_data/*.vcf.gz'
```




## Check the Output
- The `results/01_filt_vcf/08_maf_bcf` folder contains the final result files
for module 1, i.e., VCF files filtered for samples or sites with high
missingness and sites with low minor allele frequency.
- The `results/02_pop_fws/04_moimix_fws` folder contains the final results for
module 2, i.e., a table of `Fws` statistics for each sample. Each population
will have a separate table.
- The
`results/04_phase_impute_polyclonal/05_phased_imputed_polyclonal_monoclonal`
folder contains the final result files for modules 3 and 4, i.e., combined VCF
files that have been deconvoluted and imputed. Each population will have a
separate VCF file.
- The `results/05_ibdanalysis/02_callibd/hmmibd` folder contains the IBD results
obtained using `hmmIBD`.
- The `results/05_ibdanalysis/{population}/ne_output` folder contains the
results of IBDNe (effective population size).
- The `results/05_ibdanalysis/{population}/ifm_output` folder contains the
results of Infomap (population structure).

# Citations

If you find this repository useful, please cite our preprint:
> Guo, B., Borda, V., Laboulaye, R., Spring, M. D., Wojnarski, M., Vesely, B.
A., Silva, J. C., Waters, N. C., O'Connor, T. D., & Takala-Harrison, S. (2023).
Strong Positive Selection Biases Identity-By-Descent-Based Inferences of Recent
Demography and Population Structure in Plasmodium falciparum. bioRxiv : the
preprint server for biology, 2023.07.14.549114.
https://doi.org/10.1101/2023.07.14.549114

Other citations:

- `Xir,s` statistics: 
> Henden, L., Lee, S., Mueller, I., Barry, A., & Bahlo, M. (2018).
Identity-by-descent analyses for measuring population dynamics and selection in
recombining pathogens. PLoS genetics, 14(5), e1007279.
https://doi.org/10.1371/journal.pgen.1007279

- `IBDNe`
> Browning, S. R., & Browning, B. L. (2015). Accurate Non-parametric Estimation
of Recent Effective Population Size from Segments of Identity by Descent.
American journal of human genetics, 97(3), 404–418.
https://doi.org/10.1016/j.ajhg.2015.07.012

- `Infomap` algorithm
> Rosvall, M., & Bergstrom, C. T. (2008). Maps of random walks on complex
networks reveal community structure. Proceedings of the National Academy of
Sciences of the United States of America, 105(4), 1118–1123.
https://doi.org/10.1073/pnas.0706851105

- `dEploid`:
> Zhu, S. J., Almagro-Garcia, J., & McVean, G. (2018). Deconvolution of multiple
infections in Plasmodium falciparum from high throughput sequencing data.
Bioinformatics (Oxford, England), 34(1), 9–15.
https://doi.org/10.1093/bioinformatics/btx530

- `dEploid-IBD`:
> Zhu SJ, Hendry JA, Almagro-Garcia J, Pearson RD, Amato R, Miles A, Weiss DJ,
Lucas TC, Nguyen M, Gething PW, Kwiatkowski D, McVean G; Pf3k Project. The
origins and relatedness structure of mixed infections vary with local prevalence
of P. falciparum malaria. Elife. 2019 Jul 12;8:e40845. doi: 10.7554/eLife.40845.
PMID: 31298657; PMCID: PMC6684230.

- `hmmIBD`
> Schaffner SF, Taylor AR, Wong W, Wirth DF, Neafsey DE. hmmIBD: software to infer
pairwise identity by descent between haploid genotypes. Malar J. 2018 May
15;17(1):196. doi: 10.1186/s12936-018-2349-7. PMID: 29764422; PMCID: PMC5952413.

- `Beagle`
> Browning, B. L., Zhou, Y., & Browning, S. R. (2018). A One-Penny Imputed
Genome from Next-Generation Reference Panels. American journal of human
genetics, 103(3), 338–348. https://doi.org/10.1016/j.ajhg.2018.07.015

- `moimix`
https://github.com/bahlolab/moimix


## Related Private Repository:

For internal use only, refer to the following link: https://github.com/bguo068/posseleff_empirical_arxiv.
