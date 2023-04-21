#! /usr/bin/env python3

import pandas as pd

fit_gt_annot = "./ESEA_5x70_annot_genotypes_filt.tsv"

df = pd.read_csv(fit_gt_annot, sep="\t")


df.set_index("ID", inplace=True)

df = df.iloc[:, 8:].transpose().copy()

df = df.replace("\|.*$", "", regex=True).astype(int)

df.columns.name = ""
df.index.name = "Sample"

cols = ["Sample"] + list(df.columns)

df = df.reset_index()[cols].copy()

df["Sample"] = df.Sample.str.replace("~.*$", "", regex=True)

df.to_parquet("./ESEA_5x70_annot_nonsyn_mutations.pq")
