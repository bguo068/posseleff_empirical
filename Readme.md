# `poseleff_empirical`

Code for analyzing empirical WGS data (deconvoluted VCF) for the effect of
positive selection on IBD-based inference.

# Software evironment

1. Install Conda 
2. Install conda evironment: `conda env create -f ./env.ymal`
3. Activate conda environment: `conda activate empirical`
3. Install `moimix` package, Run `R` and in R console run:
`devtools::install_github("bahlolab/isoRelate")`
4. Modify and install `hmmibd`:
```sh
    git clone https://github.com/glipsnort/hmmIBD.git; cd hmmIBD
    git checkout a2f796ef8122d7f6b983ae9ac4c6fba35afcd3aa
    sed -i -e 's/const double rec_rate = 7.4e-7/const double rec_rate = 6.67e-7/' hmmIBD.c
    x86_64-conda_cos6-linux-gnu-gcc -o hmmIBD -O3 -lm -Wall hmmIBD.c
    cd ..
    cp hmmIBD/hmmIBD bin/hmmIBD
    rm -rf hmmIBD
```