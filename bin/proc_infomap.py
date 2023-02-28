#! /usr/bin/env python3

import ibdutils.utils.ibdutils as ibdutils
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--label", type=str, required=True)
args = pa.parse_args()

label = args.label


# read ibd
genome_Pf3D7_numeric = ibdutils.Genome.get_genome("Pf3D7")
ibd = ibdutils.IBD(genome=genome_Pf3D7_numeric, label="orig")
ibd.read_ibd(ibd_fn_lst=args.ibd_files)


# output files:
ofs_ifm_orig_ibd_pq = f"{label}_ifm_orig_ibd.pq"
ofs_ifm_rmpeaks_ibd_pq = f"{label}_ifm_rmpeaks_ibd.pq"


# remove highly relatedness samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)


# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()

ibd2 = ibd.duplicate("rmpeak")
ibd2.remove_peaks()

ibd._df.to_parquet(ofs_ifm_orig_ibd_pq)
ibd2._df.to_parquet(ofs_ifm_rmpeaks_ibd_pq)
