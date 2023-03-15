#! /usr/bin/env python3

import ibdutils.utils.ibdutils as ibdutils
import ibdutils.runner.ibdne as ibdne
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

ibd_jar_default = str(Path(__file__).parent / "ibdne.jar")

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--label", type=str, required=True)
pa.add_argument("--ibdne_jar", type=str, default=ibd_jar_default)
pa.add_argument("--ibdne_mincm", type=float, default=4)
pa.add_argument("--ibdne_minregion", type=float, default=10)
args = pa.parse_args()

assert args.ibdne_jar is not None

label = args.label
ibdne_mincm = args.ibdne_mincm
ibdne_minregion = args.ibdne_minregion

# create genome obj
genome_Pf3D7_numeric = ibdutils.Genome.get_genome("Pf3D7")

# create ibd obj
ibd = ibdutils.IBD(genome=genome_Pf3D7_numeric, label=f"{label}_orig")

# read ibd from files containing IBD segment info
ibd.read_ibd(ibd_fn_lst=args.ibd_files)

# ---- save origin IBD as parquet for distribution analyses
ofs_ibd_pq = f"{label}_ibddist_ibd.pq"
ibd._df.to_parquet(ofs_ibd_pq)


# ---- save ibd coverage (haploid) for all samples
ibd_all = ibd.duplicate()
ibd_all.filter_ibd_by_length(min_seg_cm=ibdne_mincm)
ibd_all.calc_ibd_cov()
ibd_all.find_peaks()
ibd_all._cov_df.to_parquet(f"{label}_orig_all.cov.pq")

# remove highly related samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)

# ---- save ibd coverage (haploid) for unrelated samples
ibd_unrel = ibd.duplicate()
ibd_unrel.filter_ibd_by_length(min_seg_cm=ibdne_mincm)
ibd_unrel.calc_ibd_cov()
ibd_unrel.find_peaks()
ibd_unrel._cov_df.to_parquet(f"{label}_orig_unrel.cov.pq")

# remove short segments
ibd.filter_ibd_by_length(min_seg_cm=ibdne_mincm)

#  remove ibd with tmrca < 1.5 (required according to IBDNe paper)
ibd.filter_ibd_by_time(min_tmrca=1.5)

# make diploid samples
# convert to heterzygous diploids
# Note: remove_hbd might not remove a lot segments as
# hdb only involves n/2 pairs of a total of n(n-1)/2 pairs (ratio: 1/(n-1)).
# When n is large, the difference is small
ibd.convert_to_heterozygous_diploids(remove_hbd=True)
ibd.flatten_diploid_ibd()

# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()

# ---- save ibd for ibdne call

# link ibdne.jar file
if not Path("ibdne.jar").exists():
    assert Path(args.ibdne_jar).exists()
    this = Path("ibdne.jar")
    target = Path(args.ibdne_jar).absolute()
    this.symlink_to(target)
    print(f"link {this} -> {target}")

# use nerunner (dry_run) to preare inputs and bash scripts
# for ibd before removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=6, mem_gb=20, dry_run=True)


# remove peaks

ibd2 = ibd.duplicate(f"{label}_rmpeaks")
ibd2.remove_peaks()


# ---- save ibd for ibdne call
# for ibd after removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd2,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=6, mem_gb=20, dry_run=True)


print(
    f"""
      output files:
        {ofs_ibd_pq}
        ibdne.jar
        {label}_orig.sh
        {label}_orig.map
        {label}_orig.ibd.gz
        {label}_orig.cov.pq
        {label}_rmpeaks.sh
        {label}_rmpeaks.map
        {label}_rmpeaks.ibd.gz
      """
)
