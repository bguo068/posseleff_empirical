#! /usr/bin/env python3

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="extract the dominant haplotype from dEploid results. Output a tsv file")
parser.add_argument("--sample", type=str, help="sample name", required=True)
parser.add_argument(
    "--hap",
    type=str,
    help="input file .hap from dEploid",
    required=True)
parser.add_argument(
    "--log",
    type=str,
    help="input file .log from dEploid",
    required=True)

args = parser.parse_args()


with open(args.log, "r") as f:
    a = f.readlines()[-1].strip().split("\t")

prop = np.array([float(x.strip()) for x in a])
order = np.argsort(prop)
first = prop[order[-1]]
second = prop[order[-2]]


df = pd.read_csv(args.hap, sep='\t', index_col=[0, 1])
cols = list(df.columns)


which = cols[order[-1]].replace("h", "")
with open(f"which.txt", "w") as f:
    f.write(f"{args.sample}.{which}")

with open(f"new_sample_name.txt", "w") as f:
    f.write(f"{args.sample}~{first:0.2f}~{second:0.2f}")

if second/first < 1/3 and first > 0.7:
    df = df.iloc[:, order[-1]]
    df.rename(f"{args.sample}~{first:0.2f}~{second:0.2f}", inplace=True)
    df.reset_index().to_csv(
        f"{args.sample}~{first:0.2f}~{second:0.2f}.tsv",
        index=None,
        sep='\t')
    with open(f"{args.sample}~valid.txt", "w") as f:
        f.write('true')


else:
    df.reset_index()[['CHROM', 'POS']].to_csv(
        f"{args.sample}~{first:0.2f}~{second:0.2f}.tsv",
        index=None,
        sep='\t')
    with open(f"{args.sample}~valid.txt", "w") as f:
        f.write('false')
