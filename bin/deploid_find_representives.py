#! /usr/bin/env python3

import allel
import numpy as np
import pandas as pd
import argparse
import sys

# parsing arguments
parser = argparse.ArgumentParser(
    description="Find the n most diverse subset of samples (maxDists method)"
)
parser.add_argument(
    "--vcf",
    type=str,
    required=True,
    help="input vcf file with th AD field that contains contains all samples ",
)
parser.add_argument(
    "--clonal_samples",
    type=str,
    required=True,
    help="a file contains a list of sample id/names used in vcf file. They should be monclonal sampled as determined by Fws>-0.95",
)
parser.add_argument(
    "--n_samples", type=int, required=True, help="how many sample you want to return"
)
parser.add_argument(
    "--out",
    type=str,
    required=True,
    help="a file containing a subset of sample names, determined by maxDists method",
)

args = parser.parse_args()

vcf = args.vcf
clonal_samples_fn = args.clonal_samples
n_samples = args.n_samples
out = args.out

try:
    clonal_samples = pd.read_csv(clonal_samples_fn, header=None)[0].drop_duplicates()
except pd.errors.EmptyDataError:
    clonal_samples = pd.Series([], dtype=str)

if clonal_samples.size == 0:
    sys.exit(f"ERROR: no clonal samples listed in {clonal_samples_fn}; cannot build a panel")

callset = allel.read_vcf(
    vcf,
    alt_number=1,
    fields=["variants/CHROM", "variants/POS", "calldata/GT", "calldata/AD", "samples"],
    samples=clonal_samples.to_numpy(),  # only load data for selected samples
)

if callset is None or len(callset.get("samples", [])) == 0:
    sys.exit(f"ERROR: none of the clonal samples were found in {vcf}")

# check if vcf has few clonal sample than needed: use all of them as the panel

n_loaded = len(callset["samples"])
if n_loaded <= n_samples:
    print(
        f"WARNING: only {n_loaded} clonal sample(s) found in {vcf}, "
        f"fewer than or equal to --n_samples {n_samples}; using all of them as the panel",
        file=sys.stderr,
    )
    rep_names = callset["samples"]
    pd.DataFrame({"Sample": rep_names}).to_csv(out, sep="\t", header=None, index=None)
    sys.exit(0)


genotypes = allel.GenotypeChunkedArray(callset["calldata/GT"])

variants = allel.VariantChunkedTable(
    [callset["variants/CHROM"], callset["variants/POS"]], names=["CHROM", "POS"]
)

readdepths = allel.GenotypeChunkedArray(callset["calldata/AD"])


def cal_genotype_maf(genotypes):
    cnts = genotypes.count_alleles()
    af = cnts[:, 1] / cnts.sum(axis=1)
    maf = np.where(af <= 0.5, af, 1 - af)
    return maf


# Filter by minor allele frequency
maf_threshold = 0.05
maf = cal_genotype_maf(genotypes)
sel = maf >= maf_threshold
variants_maf = variants.compress(sel)
genotypes_maf = genotypes.subset(sel)
readdepths_maf = readdepths.subset(sel)

# ## LD Pruning
unlinked = allel.locate_unlinked(
    genotypes_maf.to_n_alt(), size=50, step=5, threshold=0.5
)
# print("remaining snps: ", np.count_nonzero(unlinked))
variants_unlinked = variants_maf.compress(unlinked)
genotypes_unlinked = genotypes_maf.subset(unlinked)
readdepths_unlinked = readdepths_maf.subset(unlinked)

# ## GT coded as the number of alt allele
gt = np.array(genotypes_unlinked.to_n_alt())


# # ## Impute missing data using KNNImputer
# gt = gt.T
# gt = KNNImputer().fit_transform(gt)
#
# # ## PCA via allel package
#
# PCs, model = allel.pca(gt.T, scaler="standard")


# %%
# fig, axs = plt.subplots(ncols=2)
# fig.set_size_inches(14, 6)
# ax = axs[0]
# ax.plot(PCs[:, 0], PCs[:, 1], ".")
# ax.set_xlabel("PC1")
# ax.set_ylabel("PC2")
#
# ax = axs[1]
# n = model.explained_variance_ratio_.size
# ax.bar(np.arange(n) + 1, model.explained_variance_ratio_)
# ax.set_ylabel("explained variance ratio")
# ax.set_xlabel("PCs")

# %% [markdown]
# ## Sample distance

# %%
wsaf = readdepths_unlinked[:, :, 1] / (
    readdepths_unlinked[:, :, 1] + readdepths_unlinked[:, :, 0]
)

wsaf[np.isnan(wsaf)] = 0

distance = np.matmul(wsaf.T, (1 - wsaf)) + np.matmul((1 - wsaf).T, wsaf)


def find_top_n_diverse_subset(distance, N=16):
    """Adapted from maxDists:

    Given a square matrix of pairwise distances, return indices of N objects with a maximal sum of pairwise distances.


    FROM https://rdrr.io/bioc/clstutils/src/R/taxTools.R"""

    m = np.copy(distance).astype(np.float32)

    assert m.shape[0] > N
    idx = []

    for i in range(N):
        if len(idx) > 0:
            row_sum = np.nansum(m[:, idx], axis=1)
            id_max = np.argmax(row_sum)
        else:
            #  starts with the object with the minimum median distance to all other objects.
            row_median = np.median(m, axis=1)
            id_max = np.argmin(row_median)

        idx.append(id_max)
        m[id_max, :] = np.nan

    return idx


representatives = find_top_n_diverse_subset(distance, n_samples)

rep_names = callset["samples"][representatives]

pd.DataFrame({"Sample": rep_names}).to_csv(out, sep="\t", header=None, index=None)
