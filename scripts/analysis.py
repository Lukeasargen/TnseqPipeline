import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from util import tamap_to_genehits, gc_content
from util import total_count_norm, quantile_norm, length_norm
from zinb_glm import zinb_glm_llr_test


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment', type=str)
    parser.add_argument('--index', type=str)
    parser.add_argument('--controls', nargs='+', default=[], type=str)
    parser.add_argument('--samples', nargs='+', default=[], type=str)
    # parser.add_argument('--store_true', default=False, action='store_true')
    args = parser.parse_args()
    return args


def col_stats(table, columns):
    for n in columns:
        print( "\nColumn Stats : {}".format(n) )
        print( "Min={:.4f}. Max={:.4f}.".format(table[n].min(), table[n].max()) )
        print( "Median={:.4f}.".format(table[n].median()) )
        print( "Mean={:.3f}. nz mean={:.3f}".format(table[n].mean(), table[table[n]!=0][n].mean()) )
    print()


def pairwise_comparison(args):
    """ Compare the controls and samples.
        Support replicates. Input as lists.
    """
    # Pairwise Comparisons
    # Supports biological replicates
    print("Controls :", args.controls)
    print("Samples :", args.samples)

    # Append "_sum" to match the column names in the TAmap and Genehits
    controls = [c+"_sum" for c in args.controls]
    samples = [c+"_sum" for c in args.samples]

    # Load the tamap 
    tamap_filename = "data/{}/maps/{}_TAmaps.csv".format(args.experiment, args.index)
    print(" * Loading TAmap : {}".format(tamap_filename))
    tamap = pd.read_csv(tamap_filename, delimiter=",")

    # TODO : normalize the TAmap, TA counts are used in the ZINB model

    # Compress the data into a genehits table
    genehits = tamap_to_genehits(tamap)

    print(" * Normalizing...")
    # genehits = length_norm(genehits, length="Gene_Length")
    genehits = total_count_norm(genehits)
    # genehits = quantile_norm(genehits, q=0.75)

    # Compute the gc content
    # print(" * Calculating GC content...")
    # fasta_filename = r"data/demo/references/14028s_chromosome.fasta"
    # genehits = gc_content(genehits, fasta_filename)

    # Combine the replicates
    print(" * Combining replicates...")
    # For now, simply mean the reads together
    # Use mean in case there are different number of replicates for each condition
    genehits["Control_Hits"] = genehits[controls].mean(axis=1)
    genehits["Sample_Hits"] = genehits[samples].mean(axis=1)
    col_stats(genehits, columns=["Control_Hits", "Sample_Hits"])

    # TODO : renorm after combining replicates?
    # genehits = total_count_norm(genehits, columns=["Control_Hits", "Sample_Hits"])

    # Pairwise analysis below

    # TODO : move this into the tamap_to_genehits
    print(" * Calculating insertion density per gene...")
    # Need remove intergenic and group by gene
    temp = tamap[tamap['Gene_ID'].notna()].copy()  # Remove intergenic
    temp["Control"] = temp[controls].sum(axis=1)  # Combine control replicates
    temp["Sample"] = temp[samples].sum(axis=1)  # Combine sample replicates
    grouped = temp.groupby("Gene_ID", as_index=False).nunique()  # Builtin method to get unique
    # Unique_Insertions : Unique TA hits / TA Sites
    # Diversity = Unique_Insertions / TA_Counts
    # Subtract 1 bc it counts the 0 counts as unique but we don't care about this
    genehits["Control_Unique_Insertions"] = grouped["Control"]-1
    genehits["Sample_Unique_Insertions"] = grouped["Sample"]-1
    genehits["Control_Diversity"] = genehits["Control_Unique_Insertions"] / genehits["TA_Count"]
    genehits["Sample_Diversity"] = genehits["Sample_Unique_Insertions"] / genehits["TA_Count"]
    # Don't need these anymore
    del temp
    del grouped
    col_stats(genehits, columns=["Control_Diversity", "Sample_Diversity"])

    # NOTE : Below starts the actual statistical analysis, everything above was just building the table

    # Count thresholding
    print(" * Trimming by minimum count...")
    min_count = 0
    trimmed = genehits[ (genehits[["Control_Hits","Sample_Hits"]] > min_count).all(axis=1) ].copy()
    removed = genehits[ (genehits[["Control_Hits","Sample_Hits"]] <= min_count).any(axis=1) ].copy()
    print("Threshold={}. {}({:.2f}%) Genes Removed. {} Genes Remaining.".format(min_count, len(removed), 100*len(removed)/len(genehits), len(trimmed)))

    # Save genes that are removed
    removed_filename = "data/{}/analysis/removed.csv".format(args.experiment)
    print(" * Saved removed genes to {}".format(removed_filename))
    removed.to_csv(removed_filename, header=True, index=False)

    # Ratio, Log2FC, LinearDiff
    print(" * Calculating fold change...")
    trimmed["Mean"] = (trimmed["Sample_Hits"] + trimmed["Control_Hits"]) / 2.0
    smoothing = 0.0
    trimmed["Ratio"] = (trimmed["Sample_Hits"] + smoothing) / (trimmed["Control_Hits"] + smoothing)
    trimmed["Log2FC"] = np.log2(trimmed["Ratio"])
    trimmed["LinearDiff"] = trimmed["Sample_Hits"] - trimmed["Control_Hits"]
    col_stats(trimmed, columns=["Log2FC", "LinearDiff"])

    # TODO : Fitness?
    print(" * Calculating fitness...")
    # trimmed["Sample_Fitness"] = np.log10( 1 + (trimmed["Sample_Hits"]*12e6/trimmed["Sample_Hits"].sum()) / (1000*trimmed["Gene_Length"]) )
    trimmed["Sample_Fitness"] = np.log10( 1000*(1+trimmed["Sample_Hits"])/trimmed["Gene_Length"] )

    # TODO : statistical genehits here, get p-value
    print(" * Calculating statistical significance...")
    # data = np.array(trimmed[["Control_Hits", "Sample_Hits"]])
    # conditions = [0, 1]

    data = np.array(trimmed[controls+samples])
    conditions = [0]*len(controls) + [1]*len(samples)    
    print("data.shape :", data.shape)
    print("conditions :", conditions)

    # pvalues = zinb_glm_llr_test(data, conditions)
    # sig_bool = pvalues<0.05
    # num_genes = len(pvalues)
    # print("Genes :", num_genes)
    # sig_genes = np.sum(sig_bool!=(pvalues==0))
    # print("Significant genes : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/num_genes))
    # print("Nan pvalues : {}".format(np.sum(np.isnan(pvalues))))
    # print("Zero pvalues : {}".format(np.sum(pvalues==0)))
    # trimmed[["P_Value"]] = pvalues
    # trimmed[["log10"]] = -np.log10(pvalues)
    # print(trimmed.sort_values(by="P_Value", ascending=False))

    # Save the comparison
    pairwise_filename = "data/{}/analysis/pairwise.csv".format(args.experiment)
    print(" * Saved pairwise analysis to {}".format(pairwise_filename))
    trimmed.to_csv(pairwise_filename, header=True, index=False)

    # print(trimmed.sort_values(by="Sample_Fitness", ascending=False)[["Gene_ID", "Sample_Fitness"]])
    # exit()

    # Plotting below
    print(" * Calculating generating plots...")
    print("Stop plotting by pressing Ctrl+C in the terminal.")

    # These are scatter plots for now
    combine_plots = [ # x, y suffix, s, xlog, ylog
        # ["Gene_Length", "Unique_Insertions", 2, True, False],
        # ["Gene_Length", "Diversity", 2, True, False],
        # ["Gene_Length", "Hits", 1, True, True],
        # ["TA_Count", "Unique_Insertions", 2, True, False],
        # ["TA_Count", "Hits", 1, True, True],
        # ["Start", "Diversity", 8, False, False],
        # ["Start", "Hits", None, False, False],
        # ["GC", "Hits", 4, False, True],
    ]
    single_plots = [ # x, y, s, xlog, ylog
        # ["Gene_Length", "TA_Count", 2, True, True],
        # ["Start", "Log2FC", None, False, False],
        # ["Start", "LinearDiff", None, False, False],
        # ["Log2FC", "LinearDiff", None, False, False],
        # ["Control_Hits", "Sample_Hits", None, True, True],
        # ["Log2FC", "log10", 1, False, False],
    ]

    for x, y, s, xlog, ylog in combine_plots:
        print("Plotting x={}. y={}.".format(x, y))
        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(111)
        ax.scatter(x=trimmed[x], y=trimmed[f"Control_{y}"], s=s, color="tab:green", label="Control")
        ax.scatter(x=trimmed[x], y=trimmed[f"Sample_{y}"], s=s, color="tab:red", label="Sample")
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if xlog: ax.set_xscale('log')
        if ylog: ax.set_yscale('log')
        plt.legend()
        fig.tight_layout()
        plt.savefig(f"data/{args.experiment}/analysis/{x}_vs_{y}.png")

    for x, y, s, xlog, ylog in single_plots:
        print("Plotting x={}. y={}.".format(x, y))
        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(111)
        ax.scatter(x=trimmed[x], y=trimmed[y], s=s, color="tab:green")
        # plt.hlines(trimmed[y].median(), xmin=trimmed[x].min(), xmax=trimmed[x].max(), color="tab:red", label="Median={:.2f}".format(trimmed[y].median()))
        # plt.hlines(trimmed[y].mean(), xmin=trimmed[x].min(), xmax=trimmed[x].max(), color="tab:blue", label="Mean={:.2f}".format(trimmed[y].mean()))
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if xlog: ax.set_xscale('log')
        if ylog: ax.set_yscale('log')
        fig.tight_layout()
        plt.savefig(f"data/{args.experiment}/analysis/{x}_vs_{y}.png")


    # Always make the MA plot
    fig = plt.figure(figsize=[12, 8])
    A = 0.5 * ( np.log2(trimmed["Sample_Hits"]) + np.log2(trimmed["Control_Hits"]) )
    M = trimmed["Log2FC"]
    plt.scatter(x=A, y=M, s=10, color="tab:green")
    plt.hlines(M.median(), xmin=A.min(), xmax=A.max(), color="tab:red", label="Median={:.2f}".format(M.median()))
    plt.hlines(M.mean(), xmin=A.min(), xmax=A.max(), color="tab:blue", label="Mean={:.2f}".format(M.mean()))
    plt.legend()
    plt.xlabel("A = 1/2 * ( log2(sample) + log2(control) )")
    plt.ylabel("M = log2(sample) - log2(control)")
    fig.tight_layout()
    plt.savefig(f"data/{args.experiment}/analysis/MA_Plot.png")


if __name__ == "__main__":
    args = get_args()
    pairwise_comparison(args)
