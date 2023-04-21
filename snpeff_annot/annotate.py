#! /usr/bin/env python3

import pandas as pd
from subprocess import run
from pathlib import Path
import argparse
import allel
import re

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--vcf", type=str, required=True, help="Input Vcf(.gz) file")
parser.add_argument(
    "--genome",
    type=str,
    default="PlasmoDB-44_Pfalciparum3D7",
    help='Prebuilt genome used by snpEff, see "setup.sh" for details',
)
parser.add_argument(
    "--genes",
    type=str,
    default=None,
    help="csv file (without header), col 1: gene symbol, col 2: gene id",
)
parser.add_argument("--out_prefix", type=str, default="out", help="output prefix")

test = False
if test:
    arg_lst = [
        "--vcf",
        "VCFs/afrims_2015_or_newer_biallelic.vcf.gz",
        "--out_prefix",
        "VCFs/annotated",
        "--genes",
        "genes.csv",
    ]
    args = parser.parse_args(arg_lst)
else:
    args = parser.parse_args()

vcf = args.vcf
genome = args.genome
out_prefix = args.out_prefix
out_vcf = f"{out_prefix}.vcf.gz"
out_stats = f"{out_prefix}_stats.html"
out_table = f"{out_prefix}.tsv"
out_table_filt = f"{out_prefix}_annotation_filt.tsv"
out_sites_filt = f"{out_prefix}_sites_filt.txt"
out_genotypes_filt = f"{out_prefix}_genotypes_filt.tsv"


# --- ----------run annotation -------------------------------
cmd = f""" snpEff -ud 0 {genome} {vcf} -s {out_stats} -o gatk | bgzip -c > {out_vcf} """
# print(cmd)
run(cmd, shell=True)

# ------------ parse annotation -------------------------------
# TODO: need to convert multi-allele site to biallelic sites
calldata = allel.read_vcf(
    out_vcf, fields=["CHROM", "POS", "REF", "ALT", "variants/EFF"], alt_number=1
)

# extract header
eff_header = [l for l in allel.read_vcf_headers(out_vcf).headers if "EFF" in l][0]
eff_header = re.split("'", eff_header)[1]
eff_header = re.split("\[", eff_header)[0]
eff_header = eff_header.replace("(", "|").replace(" ", "").replace("|", "\t").split()
eff_header

# extract annotation fileds
eff_df = (
    pd.Series(calldata["variants/EFF"])
    .str.replace("(", "|", regex=False)
    .str.replace(")", "", regex=False)
    .str.split("|", expand=True)
)
eff_df = eff_df.iloc[
    :, : (len(eff_header))
]  # remove the last two optional columns, error/warnings
eff_df.columns = eff_header

# variant info
variants_df = pd.DataFrame(
    {
        "Chrom": calldata["variants/CHROM"],
        "POS": calldata["variants/POS"],
        "REF": calldata["variants/REF"],
        "ALT": calldata["variants/ALT"],
    }
)

# combined table
res_df = pd.concat([variants_df, eff_df], axis=1)
res_df.to_csv(out_table, sep="\t", index=None)

# filtering
if args.genes is not None:
    # read gene ist
    genes = pd.read_csv(args.genes, names=["Gene_Symbol", "Gene_Name"])
    # filter annotation table by gene id
    res_filt_df = res_df.merge(genes, how="left", on="Gene_Name")
    res_filt_df = res_filt_df[lambda x: ~x.Gene_Symbol.isnull()]
    # Only keep non synonymous site
    res_filt_df = res_filt_df[lambda x: x.Effect == "NON_SYNONYMOUS_CODING"]
    res_filt_df.to_csv(out_table_filt, sep="\t", index=None)
    res_filt_df[["Chrom", "POS"]].to_csv(
        out_sites_filt, sep="\t", index=None, header=None
    )
    # filter vcf files by gene id and remove extra information
    cmd = f"""
    bcftools view -R {out_sites_filt} {vcf} | \
        bcftools annotate -x ^INFO/AN,^INFO/AC,^INFO/AF,^FMT/GT,^FMT/AD | \
        grep -v -e '^##' > {out_genotypes_filt}
    """
    run(cmd, shell=True)
    # read vcf as a table
    gt = pd.read_csv(out_genotypes_filt, sep="\t")
    # assign variants id to "{gene_symbol}:{aa_change}"
    df1 = gt.iloc[:, :2]
    df2 = res_filt_df[["Chrom", "POS", "Gene_Symbol", "Amino_Acid_change"]].copy()
    df2["ID"] = df2.Gene_Symbol + ":" + df2.Amino_Acid_change
    df2 = df2[["Chrom", "POS", "ID"]].rename(columns={"Chrom": "#CHROM"})
    IDs = df1.merge(df2, how="left", on=["#CHROM", "POS"]).ID
    gt["ID"] = IDs
    # out_put filterred genotypes information
    gt.to_csv(out_genotypes_filt, index=None, sep="\t")
