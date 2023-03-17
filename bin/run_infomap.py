#! /usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import numpy as np
import pandas as pd
from ibdutils.utils.ibdutils import IBD


def parse_args():
    p = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("--ibd_obj", type=str, required=True)
    p.add_argument("--meta", type=str, default="")
    p.add_argument("--label", type=str, required=True)
    p.add_argument("--cut_mode", type=str, required=True)
    p.add_argument("--ntrails", type=int, default=1000)
    p.add_argument(
        "--transform", type=str, choices=["square", "cube", "none"], default="square"
    )
    p.add_argument("--ifm_mincm", type=float, default=2)
    p.add_argument("--ifm_mingwcm", type=float, default=5)
    args = p.parse_args()
    if args.transform == "none":
        args.transform = None

    return p.parse_args()


def run(args) -> pd.DataFrame:
    ibd = IBD.pickle_load(args.ibd_obj)

    # prepare meta dataframe
    if args.meta == "":
        # make meta data
        uniq_samples = ibd.get_samples_shared_ibd()
        meta = pd.DataFrame(
            {"Sample": uniq_samples, "Idx": np.arange(uniq_samples.size)}
        )
    else:
        meta = pd.read_csv(args.meta, sep="\t")

    mat = ibd.make_ibd_matrix(min_seg_cm=args.ifm_mincm, min_gw_ibd_cm=args.ifm_mingwcm)
    member_df = ibd.call_infomap_get_member_df(
        mat, meta, trials=args.ntrails, transform=args.transform
    )

    return member_df


def main():
    args = parse_args()
    member_df = run(args)

    ofs = f"{args.label}_{args.cut_mode}_member.pq"
    member_df.to_parquet(ofs)
    print(member_df)


main()
