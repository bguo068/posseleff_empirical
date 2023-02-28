#! /usr/bin/env python3
import pandas as pd
from pathlib import Path
import io
import argparse

from sklearn.metrics import top_k_accuracy_score

parser = argparse.ArgumentParser()
parser.add_argument('--bcftools_psc_stats_files', nargs='+', required=True, type=str)
parser.add_argument('--max_imiss', type=float, required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()

threshold = args.max_imiss
files= args.bcftools_psc_stats_files
out_fn = args.out

total_snps = 0
df_list = []
for p in files:
    with open(p) as f:
        psc_lines = []
        for line in f:
            if line.startswith("SN\t0\tnumber of SNPs"):
                nsnps = int(line.split('\t')[3])
                total_snps += nsnps
                print(nsnps)
            elif line.startswith('PSC'):
                psc_lines.append(line)
    df = pd.read_csv(io.StringIO(''.join(psc_lines)), sep='\t', header=None, usecols=[2, 13])
    df.columns = ['Sample', 'nMissing']
    df_list.append(df)

df_gw = pd.concat(df_list)
imiss = df_gw.groupby('Sample', sort=False)['nMissing'].sum() / total_snps
pd.Series(imiss[imiss<threshold].index).to_csv(out_fn, header=None, index=None)
