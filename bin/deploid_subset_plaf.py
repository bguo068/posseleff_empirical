#! /usr/bin/env python3
import allel
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="if plaf file does not have the same sites as vcf file, it will cause error" +
    " this code script is to make subset plaf file so that it ahs the same sites as input vcf")
parser.add_argument("--in_plaf", type=str, required=True,
                    help="input plaf file that have all the sites")
parser.add_argument("--in_vcf", type=str, required=True, help="input vcf file")
parser.add_argument("--out_plaf", type=str, required=True,
                    help="output plaf file that have a subset of sites")

args = parser.parse_args()

vcf = args.in_vcf
plaf = args.in_plaf
out = args.out_plaf


callset = allel.read_vcf(
    vcf, fields=['variants/CHROM', 'variants/POS'], alt_number=1)
CHROM = callset['variants/CHROM']
POS = callset['variants/POS']

df = pd.DataFrame({
    "CHROM": CHROM,
    "POS": POS,
})

df_plaf = pd.read_csv(plaf, sep='\t')

res = df_plaf.merge(df, how='right', on=['CHROM', 'POS'])

res.to_csv(out, index=None, sep='\t')
