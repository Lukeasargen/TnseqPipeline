"""

"""


import os
import re
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    # Init and setup
    parser.add_argument('--experiment', type=str)
    parser.add_argument('--index', type=str)
    parser.add_argument('--map', type=str)
    args = parser.parse_args()
    return args


def make_TAmap(args):
    # print(" * Creating a TAmap for {}.fasta and {}.gb".format(args.fasta, args.genbank))

    # Read the TAlist for the index

    # Read the map file

    # Add the mapped data to the TAlist
    # data/exp/maps/index_TAmap.tsv


    pass



def main():
    print("\n * Running readTAmap.py")
    args = get_args()

    make_TAlist(args)


if __name__ == "__main__":
    main()
