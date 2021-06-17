"""
Opens the TAmap for each index.
Concatinates the maps and reports missing values.
Overwrites the all_TAmaps file.
"""
import argparse

import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment', type=str)
    parser.add_argument('--indexes', nargs='+', type=str)
    args = parser.parse_args()
    return args


def merge_maps(args):
    print(" * Merging TAmaps for {}".format(", ".join(args.indexes)))

    dfs = []
    for idx in args.indexes:
        fname = "data/{}/maps/{}_TAmaps.csv".format(args.experiment, idx)
        print(" * Loading TAmap : {}".format(fname))
        tamap = pd.read_csv(fname, delimiter=",")
        dfs.append(tamap)

    merged = pd.concat(dfs, axis=0)
    print(merged)

    # Save the comparison
    output_filename = "data/{}/maps/all_TAmaps.csv".format(args.experiment)
    print(" * Saved pairwise analysis to {}".format(output_filename))
    merged.to_csv(output_filename, header=True, index=False)


def main():
    print("\n * Running mergeMaps.py")
    args = get_args()
    merge_maps(args)


if __name__ == "__main__":
    main()
