#! /usr/bin/env python3
import argparse
import pandas as pd
import allel

parser = argparse.ArgumentParser(
    description="Use the AD field of the single sample VCF to calculate high leverage sites which is used in dEploid -exclude opts")
parser.add_argument("--vcf", type=str, required=True,
                    help='input vcf file, single sample')
parser.add_argument(
    "--sites", type=str, required=True,
    help="output file containing sites to exclude")
args = parser.parse_args()

# vcf = '/autofs/chib/toconnor_grp/bing/20210923_pf_malaria_gen_cambodia_ibd/test_deploid/IVC17_OM_0661_clonal/IVC17_OM_0661.vcf.gz'
# exclude = '/autofs/chib/toconnor_grp/bing/20210923_pf_malaria_gen_cambodia_ibd/test_deploid/IVC17_OM_0661_clonal/exclude.txt'

vcf = args.vcf
exclude = args.sites

callset = allel.read_vcf(vcf, alt_number=1, fields=[
                         'variants/CHROM', 'variants/POS', 'calldata/AD'])

df = pd.DataFrame(
    {'CHROM': callset['variants/CHROM'],
     'POS': callset['variants/POS'],
     'AD0': callset['calldata/AD'][:, 0, 0],
     'AD1': callset['calldata/AD'][:, 0, 1]})
df['Coverage'] = df.AD0 + df.AD1
high_leverage = df.Coverage > df.Coverage.quantile(0.995)

df['Exclude'] = 0
df.loc[high_leverage, 'Exclude'] = 1

for i in df.index[high_leverage]:
    low = i - 10
    high = i + 10
    if low <= 0:
        low = 0
    if high > df.index[-1]:
        high = df.index[-1]

    df.loc[low:(high+1), 'Exclude'] += 1

df.loc[df.Exclude > 1, ['CHROM', 'POS']].to_csv(exclude, sep='\t', index=None)
