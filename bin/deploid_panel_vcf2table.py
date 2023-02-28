#! /usr/bin/env python3

# %%
from typing import Optional
import allel
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="extact info from a homozgyous diploid VCF to a panel.txt file"
)
parser.add_argument("--vcf", type=str, required=True, help="input vcf file")
parser.add_argument("--panel", type=str, required=True, help="output vcf file")
parser.add_argument(
    "--samples",
    type=str,
    default=None,
    help="optional file contains sample list that will be kept in the output, if not specified all sample from the input vcf file is kept",
)
args = parser.parse_args()

# "/autofs/chib/toconnor_grp/bing/20210923_pf_malaria_gen_cambodia_ibd/test_deploid/panel.vcf.gz"
in_vcf = args.vcf
# "/autofs/chib/toconnor_grp/bing/20210923_pf_malaria_gen_cambodia_ibd/test_deploid/panel.txt"
out_panel = args.panel
sample_fn = args.samples

# %%
callset = allel.read_vcf(
    in_vcf,
    fields=[
        "variants/CHROM",
        "variants/POS",
        "variants/REF",
        "variants/ALT",
        "calldata/GT",
        "samples",
    ],
    alt_number=1,
)
# %%
CHROM = callset["variants/CHROM"]
POS = callset["variants/POS"]
# REF = callset['variants/REF']
# ALT = callset['variants/ALT']
GT = callset["calldata/GT"][:, :, 0]
SAMPLES = pd.Series(callset["samples"])

# remove the suffixed added by deploid process
all_samples = SAMPLES.str.replace("~.*$", "", regex=True)

# %%
if sample_fn:
    input_samples = pd.read_csv(sample_fn, header=None)[0].str.replace(
        "~.*$", "", regex=True
    )
else:
    input_samples = all_samples

sample_sel = all_samples.isin(input_samples)

# chr and pos
df = pd.DataFrame(
    {
        "CHROM": CHROM,
        "POS": POS,
    }
)

# gt
df_gt = pd.DataFrame(GT[:, sample_sel], columns=SAMPLES[sample_sel])

# conat horizontally
df = pd.concat([df, df_gt], axis=1)

df.to_csv(out_panel, index=None, sep="\t")
