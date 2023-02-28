#! /usr/bin/env python3

import allel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, leaves_list, to_tree, ward
from sklearn.cluster import AgglomerativeClustering
from sklearn.impute import KNNImputer
import argparse

# parsing arguments
parser = argparse.ArgumentParser(
    description="Using two method to find most n diverse subset of samples"
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
    "--out1",
    type=str,
    required=True,
    help="a file containing a subset of sample names, determined by hierachical method",
)
parser.add_argument(
    "--out2",
    type=str,
    required=True,
    help="a file containing a subset of sample names, determined by maxDists method",
)

args = parser.parse_args()

vcf = args.vcf
clonal_samples_fn = args.clonal_samples
n_samples = args.n_samples
out_1 = args.out1
out_2 = args.out2

clonal_samples = pd.read_csv(clonal_samples_fn, header=None)[0]

callset = allel.read_vcf(
    vcf,
    alt_number=1,
    fields=["variants/CHROM", "variants/POS", "calldata/GT", "calldata/AD", "samples"],
    samples=clonal_samples.to_numpy(),  # only load data for selected samples
)

# check if vcf has few clonal sample than needed

if clonal_samples.size <= n_samples:
    rep_names_1 = callset["samples"]
    rep_names_2 = callset["samples"]
    pd.DataFrame({"Sample": rep_names_1}).to_csv(
        out_1, sep="\t", header=None, index=None
    )
    pd.DataFrame({"Sample": rep_names_2}).to_csv(
        out_2, sep="\t", header=None, index=None
    )
    exit(0)


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

model = AgglomerativeClustering(
    n_clusters=None, distance_threshold=10000, linkage="average"
)

model = model.fit(distance)


def get_linkage_matrix(model):
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)
    return linkage_matrix


link_mat = get_linkage_matrix(model)
leaves = leaves_list(link_mat)
tree_root = to_tree(link_mat)


def get_map_from_old_id_to_reordered_id(leaves):
    tmp = pd.DataFrame(leaves).reset_index()
    tmp.columns = ["New", "Old"]
    tmp = tmp.sort_values(["Old"]).set_index("Old")
    return tmp["New"].to_numpy()


id_map = get_map_from_old_id_to_reordered_id(leaves)


def find_n_level_nodes(tree_root, n_levels):
    levels = [[tree_root]]
    while len(levels) <= n_levels:
        last_level = levels[-1]
        new_level = []
        for node in last_level:
            if node.is_leaf():
                new_level.append(node)
            else:
                new_level.append(node.get_left())
                new_level.append(node.get_right())
        levels.append(new_level)
    return levels[-1]


n_levels = np.log2(n_samples).astype(np.int32)
level_n = find_n_level_nodes(tree_root, n_levels)


def find_representatives_for_nodes(nodes):
    representatives = []
    for node in nodes:

        lst = node.pre_order()
        mid = len(lst) // 2
        representatives.append(lst[mid])
    return representatives


representatives = find_representatives_for_nodes(level_n)


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


representatives_2 = find_top_n_diverse_subset(distance, n_samples)


reordered = distance[:, leaves][leaves, :]
reordered.shape


# %%
# m = 0.5
# m2 = 0.2
# h = 3
# s = 10
#
# height = m + h + s + m + m2
# width = m + s + m
#
# fig = plt.gcf()
# fig.set_size_inches(width, height)
# ax_heatmap = fig.add_axes([m / width, m / height, s / width, s / height])
# ax_dendrogram = fig.add_axes([m / width, (m + s + m2) / height, s / width, h / height])
#
# dendrogram(
#     link_mat,
#     ax=ax_dendrogram,
#     link_color_func=lambda k: "red" if k in [node.id for node in level_4] else "black",
# )
# ax_dendrogram.get_xaxis().set_visible(False)
# ax_dendrogram.get_yaxis().set_visible(False)
#
#
# ax_heatmap.imshow(reordered)
#
# ax_heatmap.set_yticks(id_map[representatives])
# ax_heatmap.set_ylabel("Representatives\n clustering tree")
#
# ax_heatmap.set_xticks(id_map[representatives_2])
# ax_heatmap.set_xlabel("Representatives\n (maxDists)")


rep_names_1 = callset["samples"][representatives]


rep_names_2 = callset["samples"][representatives_2]


pd.DataFrame({"Sample": rep_names_1}).to_csv(out_1, sep="\t", header=None, index=None)
pd.DataFrame({"Sample": rep_names_2}).to_csv(out_2, sep="\t", header=None, index=None)
