"""
Opens the alignments and the TAlist.
Bins the alignments for each TA site in forward and reverse direction.
Appends the counts to the TAlist.
"""

import os
import time
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from util import tamap_to_genehits
from util import total_count_norm


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment', type=str)
    parser.add_argument('--index', type=str)
    parser.add_argument('--map', type=str)
    args = parser.parse_args()
    return args


def make_TAmap(args):
    print(" * Creating a TAmap for {} and {}".format(args.index, args.map))
    t0 = time.time()

    talist_filename = "data/{}/references/{}_TAlist.csv".format(args.experiment, args.index)
    read_name, extension = os.path.splitext(args.map)
    mapped_filename = "data/{}/reads/processed/{}_trimmed_mapped".format(args.experiment, read_name)

    # Read the TAlist for the index into a pandas dataframe
    # Headers: Accession, Loci, Gene_ID, Locus_Tag, Start, End, Direction, TA_Site
    talist = pd.read_csv(talist_filename, delimiter=",")
    # print(talist)

    # Read the map file output by bowtie
    # The file has no headers so these headers below are added to the dataframe
    map_headers = ["Read_ID", "Direction", "Accession", "Location", "Sequence", "Quality", "Counts?", "??"]
    mapped = pd.read_csv(mapped_filename, delimiter="\t", header=None, names=map_headers)
    # print(mapped)

    # Convert the mapped locations strings to int, usually not necessary but it's here just in case
    mapped.iloc[:, 3] = pd.to_numeric(mapped.iloc[:, 3])

    seq_length = talist.iloc[-1]["End"]  # Last row, end of region
    # Create the bins for the couting
    bins = talist['TA_Site'].values
    total_ta_sites = len(bins)
    # Add the beginning and ending bins
    bins = np.append(bins, seq_length)
    bins = np.insert(bins, 0, 0)
    # print(bins[:4], bins[-4:])

    """
    This is how the TA alignments are binned for counting for 
    a genome with 4 TA sites.
              1           2           3           4
   -----------TA----------TA----------TA----------TA-----------
              | 1 forward | 2 forward | 3 forward | 4 forward |
    | 1 reverse | 2 reverse | 3 reverse | 4 reverse |

    Reasoning:
    Forward: Reads will be trimmed by a quality threshold and cropped.
    These changes shift the exact alignment. However the shift is
    always to a higher location. Therefore, it is reasonable to bin
    forwards reads in between the TA sites.
    Reverse: The same logic applies to reverse reads. However, the
    reads are shifted up by 2 so the TA site is within the bin.
    """

    # Forward Counts
    forward_locations = np.array(mapped["Location"][mapped["Direction"]=="+"].tolist())
    # print(forward_locations[:8])

    forward_hist, edges = np.histogram(forward_locations, bins=bins, density=False)
    # If the genome is CIRCULAR, add the first bin to the end
    forward_hist[-1] += forward_hist[1]
    # Now remove the first bin, skip it with [1:]
    forward_hist = forward_hist[1:]

    # Reverse Counts
    # Get the sequence length by finding string length
    reverse_lengths = np.array(mapped["Sequence"][mapped["Direction"]=="-"].str.len().tolist())
    # print(reverse_lengths[:8])
    reverse_locations = np.array(mapped["Location"][mapped["Direction"]=="-"].tolist())
    # print(reverse_locations[:8])
    # Alignment is also shifted by length+1. The +1 accounts for the 1 indexing of location.
    reverse_locations = reverse_locations + reverse_lengths + 1
    # print(reverse_locations[:8])

    reverse_hist, edges = np.histogram(reverse_locations, bins=bins, density=False)
    # If the genome is CIRCULAR, add the last bin to the beginning bin
    reverse_hist[0] += reverse_hist[-1]
    # Now remove the last bin, skip it with [:-1]
    reverse_hist = reverse_hist[:-1]

    # Save to map file
    # Mapped data for only this read
    tamap_filename = "data/{}/maps/{}_{}_TAmap.csv".format(args.experiment, args.index, read_name)
    print(" * Saving {}".format(tamap_filename))
    tamap = talist.copy()
    tamap[f"{read_name}_forward"] = pd.Series(forward_hist)
    tamap[f"{read_name}_reverse"] = pd.Series(reverse_hist)
    tamap[f"{read_name}_sum"] = tamap[f"{read_name}_forward"] + tamap[f"{read_name}_reverse"]
    tamap.to_csv(tamap_filename, header=True, index=False)
    duration = time.time() - t0

    # Plot all hits
    fig, ax = plt.subplots(2, figsize=[16,10])
    tamap.plot(x="TA_Site", y=f"{read_name}_forward", ax=ax[0], color='tab:green')
    ax[0].set_title("Forward TA Hits")
    tamap.plot(x="TA_Site", y=f"{read_name}_reverse", ax=ax[1], color='tab:red')
    ax[1].set_title("Reverse TA Hits")
    fig.tight_layout()
    plt.savefig("data/{}/maps/{}_{}_all_hits.png".format(args.experiment, args.index, read_name))

    # Load up the full map and merge
    merge_filename = "data/{}/maps/{}_TAmaps.csv".format(args.experiment, args.index)
    merge = pd.read_csv(merge_filename, delimiter=",")
    merged = pd.merge(merge, tamap)
    merged.to_csv(merge_filename, header=True, index=False)

    print(" * Creating GeneHits table...")
    genehits = tamap_to_genehits(merged)
    genehits_filename = "data/{}/maps/{}_Genehits.csv".format(args.experiment, args.index)
    print(" * Saving {}".format(genehits_filename))
    genehits.to_csv(genehits_filename, header=True, index=False)

    # Normalize the genehits table
    normed = total_count_norm(genehits)
    normed_filename = "data/{}/maps/{}_GenehitsNorm.csv".format(args.experiment, args.index)
    print(" * Saving {}".format(normed_filename))
    normed.to_csv(normed_filename, header=True, index=False)

    # Summarize Mapping
    # Remove all intergenic regions and count gene hits
    grouped = tamap.groupby("Gene_ID", dropna=True)
    gene_total_ta_sites = grouped["TA_Site"].count().sum(axis=0)
    genes_w_ta = len(grouped)

    # This data is a row for each gene with the sum of the counts
    hits_per_gene_forward = grouped[f"{read_name}_forward"].sum()
    hits_per_gene_reverse = grouped[f"{read_name}_reverse"].sum()

    # This is the number of TA sites within a gene that were hit
    genehits_fr = tamap[tamap['Gene_ID'].notna()]
    gene_ta_forward = genehits_fr[f"{read_name}_forward"].astype(bool).sum(axis=0)
    gene_ta_reverse = genehits_fr[f"{read_name}_reverse"].astype(bool).sum(axis=0)

    # These are integers of the number of genes hit
    hit_genes_forward = hits_per_gene_forward.astype(bool).sum(axis=0)
    hit_genes_reverse = hits_per_gene_reverse.astype(bool).sum(axis=0)
    hit_genes = (hits_per_gene_forward+hits_per_gene_reverse).astype(bool).sum(axis=0)


    stat_str = "Summary of {} mapping to {}.".format(read_name,  args.index)
    stat_str += "\nDuration : {:.3f} seconds.".format(duration)
    stat_str += "\nTotal Reads : {}".format(len(mapped.index))
    stat_str += "\nTotal TA Sites : {}".format(total_ta_sites)
    stat_str += "\nTotal Hit TA Sites : {}".format(np.count_nonzero(forward_hist+reverse_hist))
    stat_str += "\nPercentage Hit TA Sites: {:.2f}%".format(100*np.count_nonzero(forward_hist+reverse_hist)/total_ta_sites)
    stat_str += "\nTotal Gene TA Sites : {}".format(gene_total_ta_sites)
    stat_str += "\nTotal Genes with TA Sites : {}".format(genes_w_ta)
    stat_str += "\nHit Genes : {}".format(hit_genes)
    stat_str += "\nPercentage Hit Genes : {:.2f}%".format(100*hit_genes/genes_w_ta)

    # Forward
    stat_str += "\nForward Stats:"
    stat_str += "\nTotal Reads : {}".format(len(forward_locations))
    stat_str += "\nMax : {}".format(np.max(forward_hist))
    stat_str += "\nMean : {:.3f}".format(np.mean(forward_hist))
    nonzero_mean = np.true_divide(forward_hist.sum(), (forward_hist!=0).sum())
    stat_str += "\nNonzero Mean : {:.3f}".format(nonzero_mean)
    stat_str += "\nVariance : {:.3f}".format(np.var(forward_hist))
    stat_str += "\nStandard Deviation : {:.3f}".format(np.std(forward_hist))
    stat_str += "\nHit TA Sites : {}".format(np.count_nonzero(forward_hist))
    stat_str += "\nPercentage Hit TA Sites : {:.2f}%".format(100*np.count_nonzero(forward_hist)/total_ta_sites)
    stat_str += "\nHit Gene TA Sites : {}".format(gene_ta_forward)
    stat_str += "\nPercentage Hit Gene TA Sites : {:.2f}%".format(100*gene_ta_forward/gene_total_ta_sites)
    stat_str += "\nHit Genes : {}".format(hit_genes_forward)
    stat_str += "\nPercentage Hit Genes : {:.2f}%".format(100*hit_genes_forward/genes_w_ta)

    # Reverse
    stat_str += "\nReverse Stats:"
    stat_str += "\nTotal Reads : {}".format(len(reverse_locations))
    stat_str += "\nMax : {}".format(np.max(reverse_hist))
    stat_str += "\nMean : {:.3f}".format(np.mean(reverse_hist))
    nonzero_mean = np.true_divide(reverse_hist.sum(), (reverse_hist!=0).sum())
    stat_str += "\nNonzero Mean : {:.3f}".format(nonzero_mean)
    stat_str += "\nVariance : {:.3f}".format(np.var(reverse_hist))
    stat_str += "\nStandard Deviation : {:.3f}".format(np.std(reverse_hist))
    stat_str += "\nHit TA Sites : {}".format(np.count_nonzero(reverse_hist))
    stat_str += "\nPercentage Hit TA Sites : {:.2f}%".format(100*np.count_nonzero(reverse_hist)/total_ta_sites)
    stat_str += "\nHit Gene TA Sites : {}".format(gene_ta_reverse)
    stat_str += "\nPercentage Hit Gene TA Sites : {:.2f}%".format(100*gene_ta_reverse/gene_total_ta_sites)
    stat_str += "\nHit Genes : {}".format(hit_genes_reverse)
    stat_str += "\nPercentage Hit Genes : {:.2f}%".format(100*hit_genes_reverse/genes_w_ta)

    # Find top genes that cover % of all reads
    total = genehits[f"{read_name}_sum"].sum()
    temp = genehits.sort_values(by=f"{read_name}_sum", ascending=False).cumsum()
    for i in [0.1, 0.2, 0.5, 0.9, 0.95, 0.99, 0.999, 0.9999]:
        c =  (temp[f"{read_name}_sum"]<total*i).astype(bool).sum(axis=0)
        stat_str += "\n{:4.2f}% of reads accounted for by {:5.2f}% of genes ({})".format(i*100, 100*c/len(genehits), c)

    print(stat_str)

    stat_filename = "data/{}/maps/{}_{}_stats.txt".format(args.experiment, args.index, read_name)
    with open(stat_filename, "w") as f:
        f.write(stat_str)


def main():
    print("\n * Running readTAmap.py")
    args = get_args()
    make_TAmap(args)


if __name__ == "__main__":
    main()
