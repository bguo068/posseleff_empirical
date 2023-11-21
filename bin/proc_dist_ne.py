#! /usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import ibdutils.runner.ibdne as ibdne
import ibdutils.utils.ibdutils as ibdutils

ibd_jar_default = str(Path(__file__).parent / "ibdne.jar")

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--vcf_files", type=str, nargs=14, required=True)
pa.add_argument("--label", type=str, required=True)
pa.add_argument("--ibdne_jar", type=str, default=ibd_jar_default)
pa.add_argument("--ibdne_mincm", type=float, default=4)
pa.add_argument("--ibdne_minregion", type=float, default=10)
pa.add_argument("--ibdne_flatmeth", type=str, default="none")
pa.add_argument(
    "--peak_validate_meth", type=str, choices=["xirs", "ihs"], default="xirs"
)
pa.add_argument(
    "--ibdne_no_diploid_conversion",
    type=str,
    choices=["true", "false"],
    default="false",
)
args = pa.parse_args()
assert args.ibdne_jar is not None
if args.ibdne_flatmeth is None:
    args.ibdne_flatmeth == "none"
assert args.ibdne_flatmeth in ["none", "merge", "keep_hap_1_only"]

# output prefix string
label_str = "_".join(
    [
        str(x)
        for x in [
            args.label,
            args.ibdne_mincm,
            args.ibdne_minregion,
            args.ibdne_flatmeth,
        ]
    ]
)
ibdne_mincm = args.ibdne_mincm
ibdne_minregion = args.ibdne_minregion

# ------ ibd for all samples
# create genome obj
genome_Pf3D7_numeric = ibdutils.Genome.get_genome("Pf3D7")
# create ibd obj
ibd = ibdutils.IBD(genome=genome_Pf3D7_numeric, label=f"{label_str}_orig")
# read ibd from files each containing IBD segment info for a chromosome
# and remove sample name suffix (such as "~0.9-0.2") that was used
# for carrying haploid genome ratios
ibd.read_ibd(ibd_fn_lst=args.ibd_files, rm_sample_name_suffix=True)

# save origin IBD obj for distribution analyses
of_ibddist_ibdobj = f"{label_str}_ibddist.ibdobj.gz"
ibd.pickle_dump(of_ibddist_ibdobj)


# save IBD with coverage/peaks for all samples (diploid)
ibd_all = ibd.duplicate()
ibd_all.filter_ibd_by_length(min_seg_cm=ibdne_mincm)
# calculate XiR,s/iHS
if args.peak_validate_meth == "xirs":
    ibd_all.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
elif args.peak_validate_meth == "ihs":
    ibd_all.calc_ihs(vcf_fn_lst=args.vcf_files, min_maf=0.01, multitest_correction=None)
else:
    raise NotImplementedError(f"{args.peak_validate_meth} is not valid method")
# xirs_df = ibd_all.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)

if args.ibdne_no_diploid_conversion == "true":
    print("skip diploid conversion")
elif args.ibdne_no_diploid_conversion == "false":
    # convert to heterzygous diploids
    # Note: remove_hbd might not remove a lot segments as
    # hdb only involves n/2 pairs of a total of n(n-1)/2 pairs (ratio: 1/(n-1)).
    # When n is large, the difference is small
    ibd_all.convert_to_heterozygous_diploids(remove_hbd=True)

    if args.ibdne_flatmeth != "none":
        ibd_all.flatten_diploid_ibd(method=args.ibdne_flatmeth)
else:
    raise NotImplementedError(
        f"{args.ibdne_no_diploid_conversion} is not valid value is not a value value for option --ibdne_no_diploid_conversion"
    )

ibd_all.calc_ibd_cov()
ibd_all.find_peaks()
# Filter peaks by XiR,s/iHS
if args.peak_validate_meth == "xirs":
    xirs_df = ibd_all._xirs_df
    ibd_all.filter_peaks_by_xirs(xirs_df)
elif args.peak_validate_meth == "ihs":
    ibd_all.filter_peaks_by_ihs(min_ihs_hits=1, alpha=0.05)
else:
    raise NotImplementedError(f"{args.peak_validate_meth} is not valid method")
# ibd_all.filter_peaks_by_xirs(xirs_df)
ibd_all.pickle_dump(f"{label_str}_orig_all.ibdcov.ibdobj.gz")

# ------ ibd for samples
# remove highly related samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)

# save ibd coverage/peaks/ibdne for unrelated samples
ibd_unrel = ibd.duplicate(f"{label_str}_orig")
ibd_unrel.filter_ibd_by_length(min_seg_cm=ibdne_mincm)
# NOTE: xirs_df is calculated before converted to diploids
# calculate XiR,s/iHS
if args.peak_validate_meth == "xirs":
    ibd_unrel.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
elif args.peak_validate_meth == "ihs":
    ibd_unrel.calc_ihs(
        vcf_fn_lst=args.vcf_files, min_maf=0.01, multitest_correction=None
    )
else:
    raise NotImplementedError(f"{args.peak_validate_meth} is not valid method")
# xirs_df = ibd_all.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
xirs_df = ibd_unrel.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)


if args.ibdne_no_diploid_conversion == "true":
    print("skip diploid conversion")
elif args.ibdne_no_diploid_conversion == "false":
    # convert to heterzygous diploids
    # Note: remove_hbd might not remove a lot segments as
    # hdb only involves n/2 pairs of a total of n(n-1)/2 pairs (ratio: 1/(n-1)).
    # When n is large, the difference is small
    ibd_unrel.convert_to_heterozygous_diploids(remove_hbd=True)

    if args.ibdne_flatmeth != "none":
        ibd_unrel.flatten_diploid_ibd(method=args.ibdne_flatmeth)
else:
    raise NotImplementedError(
        f"{args.ibdne_no_diploid_conversion} is not valid value is not a value value for option --ibdne_no_diploid_conversion"
    )

ibd_unrel.calc_ibd_cov()
ibd_unrel.find_peaks()
# Filter peaks by XiR,s/iHS
if args.peak_validate_meth == "xirs":
    xirs_df = ibd_unrel._xirs_df
    ibd_unrel.filter_peaks_by_xirs(xirs_df)
elif args.peak_validate_meth == "ihs":
    ibd_unrel.filter_peaks_by_ihs(min_ihs_hits=1, alpha=0.05)
else:
    raise NotImplementedError(f"{args.peak_validate_meth} is not valid method")
# ibd_unrel.filter_peaks_by_xirs(xirs_df)
ibd_unrel.pickle_dump(f"{label_str}_orig.ibdne.ibdobj.gz")

# remove peaks and save object for Ne
ibd_unrel_rmpeaks = ibd_unrel.duplicate(f"{label_str}_rmpeaks")
ibd_unrel_rmpeaks.remove_peaks()
# cut and split. Note cut_and_split_ibd does not update self.
ibd_unrel_rmpeaks._df = ibd_unrel_rmpeaks.cut_and_split_ibd()
# save of OBJ copy
ibd_unrel_rmpeaks.pickle_dump(f"{label_str}_rmpeaks.ibdne.ibdobj.gz")


# save ibd for ibdne call

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
    ibd=ibd_unrel,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=6, mem_gb=20, dry_run=True)


# ---- save ibd for ibdne call
# for ibd after removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd_unrel_rmpeaks,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=6, mem_gb=20, dry_run=True)


print(
    f"""
      output files:
        ibdne.jar
        {label_str}_orig.sh
        {label_str}_orig.map
        {label_str}_orig.ibd.gz
        {label_str}_orig.cov.pq
        {label_str}_rmpeaks.sh
        {label_str}_rmpeaks.map
        {label_str}_rmpeaks.ibd.gz
        {label_str}_ibddist.ibdobj.gz
        {label_str}_orig_all.ibdcov.ibdobj.gz
        {label_str}_orig.ibdne.ibdobj.gz
        {label_str}_rmpeaks.ibdne.ibdobj.gz
    """
)
