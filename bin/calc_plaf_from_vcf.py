#! /usr/bin/env python3
import allel
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", type=str, required=True)
parser.add_argument("--out", type=str, required=True)
args = parser.parse_args()

vcf_fn = args.vcf
out_fn = args.out

callset = allel.read_vcf(
    vcf_fn, fields=["variants/CHROM", "variants/POS", "calldata/AD"], alt_number=1
)
AD = callset["calldata/AD"]
CHROM = callset["variants/CHROM"]
POS = callset["variants/POS"]

REF = AD[:, :, 0]
ALT = AD[:, :, 1]
PLAF = ALT.sum(axis=1) / (ALT.sum(axis=1) + REF.sum(axis=1))

df = pd.DataFrame({"CHROM": CHROM, "POS": POS, "PLAF": PLAF})

df.to_csv(out_fn, index=None, sep="\t")
