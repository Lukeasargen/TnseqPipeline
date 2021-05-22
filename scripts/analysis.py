import time
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from util import tamap_to_genehits, column_stats
from util import total_count_norm, quantile_norm, gene_length_norm
from util import ttr_norm
from util import bh_procedure
from util import time_to_string

from zinb_glm import zinb_glm_llr_test


# TODO : arg help message
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment', type=str, required=True)
    parser.add_argument('--index', type=str, required=True)
    parser.add_argument('--controls', nargs='+', default=[], type=str, required=True)
    parser.add_argument('--samples', nargs='+', default=[], type=str, required=True)
    parser.add_argument('--debug', default=False, action='store_true')
    parser.add_argument('--plot', default=False, action='store_true')
    
    # lowest number of hits per gene, removes zeros to eliminate log errors
    parser.add_argument('--min_count', default=1, type=int)

    # lowest amount of observed TA site with hits
    parser.add_argument('--min_inserts', default=2, type=int)
    
    # smooth the ratio for small counts ratio=(sample+smoothing)/(control+smoothing)
    parser.add_argument('--smoothing', default=1, type=float)

    # calculates the GC content of each gene
    parser.add_argument('--gc', default=False, action='store_true')

    args = parser.parse_args()
    return args



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
    test_columns = controls+samples

    # Load the tamap 
    tamap_filename = "data/{}/maps/{}_TAmaps.csv".format(args.experiment, args.index)
    print(" * Loading TAmap : {}".format(tamap_filename))
    tamap = pd.read_csv(tamap_filename, delimiter=",")

    print(" * Normalizing...")
    if args.debug: print("\nStats before norm:"); column_stats(tamap, columns=test_columns)
    # tamap = gene_length_norm(tamap, columns=test_columns, debug=args.debug)
    # tamap = total_count_norm(tamap, columns=test_columns, debug=args.debug)
    # tamap = quantile_norm(tamap, q=0.75, columns=test_columns, debug=args.debug)
    tamap = ttr_norm(tamap, trim=0.1, columns=test_columns, debug=args.debug)
    if args.debug: print("\nStats after norm:"); column_stats(tamap, columns=test_columns)

    # exit()

    print(" * Compressing TAmap into Genehits table...")
    fasta_filename = "data/{}/references/{}.fasta".format(args.experiment, args.index) if args.gc else None
    if args.debug: print(f"GC content fasta filename : {fasta_filename}")
    genehits = tamap_to_genehits(tamap, fasta_filename)

    # Combine the replicates, average the counts
    print(" * Combining replicates...")
    genehits["Control_Hits"] = genehits[controls].mean(axis=1)
    genehits["Sample_Hits"] = genehits[samples].mean(axis=1)
    column_stats(genehits, columns=["Control_Hits", "Sample_Hits"])

    # Pairwise analysis below

    # TODO : can this be moved into tamap_to_genehits
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
    column_stats(genehits, columns=["Control_Diversity", "Sample_Diversity"])

    # NOTE : Below starts the actual statistical analysis, everything above was just building the table

    # Count and Observation thresholding
    print(" * Trimming by minimum counts and unique insertions...")
    keep = (genehits[["Control_Hits","Sample_Hits"]] >= args.min_count).all(axis=1) & (genehits[["Control_Unique_Insertions","Sample_Unique_Insertions"]] >= args.min_inserts).all(axis=1)
    trimmed = genehits[keep].copy()
    removed = genehits[~keep].copy()
    print("Thresholds: min_count={}. min_inserts={}.".format(args.min_count, args.min_inserts))
    print("{}({:.2f}%) Genes Removed. {} Genes Remaining.".format(len(removed), 100*len(removed)/len(genehits), len(trimmed)))

    if args.debug: column_stats(trimmed, columns=["Control_Hits", "Sample_Hits", "Control_Unique_Insertions", "Sample_Unique_Insertions"])

    # Save genes that are removed
    removed_filename = "data/{}/analysis/removed.csv".format(args.experiment)
    print(" * Saved removed genes to {}".format(removed_filename))
    removed.to_csv(removed_filename, header=True, index=False)

    # Ratio, Log2FC, LinearDiff
    print(" * Calculating fold change...")
    trimmed["Mean"] = (trimmed["Sample_Hits"] + trimmed["Control_Hits"]) / 2.0
    trimmed["Ratio"] = (trimmed["Sample_Hits"] + args.smoothing) / (trimmed["Control_Hits"] + args.smoothing)
    trimmed["Log2FC"] = np.log2(trimmed["Ratio"])
    trimmed["LinearDiff"] = trimmed["Sample_Hits"] - trimmed["Control_Hits"]
    column_stats(trimmed, columns=["Log2FC"])

    # Possible Metrics
    trimmed["Diversity_Ratio"] = (trimmed["Sample_Diversity"]) / (trimmed["Control_Diversity"])
    trimmed["Log2_Diversity_Ratio"] = np.log2(trimmed["Diversity_Ratio"])
    # print("TRIMMED diversity")
    # column_stats(trimmed, columns=["Control_Diversity", "Sample_Diversity"])
    # column_stats(trimmed, columns=["Diversity_Ratio", "Log2_Diversity_Ratio"])


    # TODO : Fitness?
    # print(" * Calculating fitness...")
    # trimmed["Sample_Fitness"] = np.log10( 1 + (trimmed["Sample_Hits"]*12e6/trimmed["Sample_Hits"].sum()) / (1000*trimmed["Gene_Length"]) )
    # trimmed["Sample_Fitness"] = np.log10( 1000*(1+trimmed["Sample_Hits"])/trimmed["Gene_Length"] )


    print(" * Calculating statistical significance...")
    print("Stop early by pressing Ctrl+C in the terminal.")
    t0 = time.time()  # Start time
    trimmed["P_Value"] = np.nan
    trimmed["P_Sig"] = False
    c = 0  # genes counter, this is different than the index value
    for i in trimmed.index:
        try:
            # Helpful time estimate
            if (c+1) % 10 == 0:
                duration = time.time()-t0
                remaining = duration/(c+1) * (len(trimmed)-c+1)
                print("gene {}/{}. {:.1f} genes/second. elapsed={}. remaining={}.".format(c+1, len(trimmed), (c+1)/duration, time_to_string(duration), time_to_string(remaining)), end="\r")
            c += 1
            # # gene_name is used to index the full TAmap 
            # # size is used to get the length of the condition array
            gene_name, size = trimmed.loc[i][["Gene_ID", "TA_Count"]]
            if args.debug: print("gene_name :", gene_name)
            df = tamap[tamap['Gene_ID']==gene_name]
            gene_data = np.array(df[test_columns]).T.reshape(-1)
            conditions = np.array([0]*size*len(controls) + [1]*size*len(samples))
            pvalue = zinb_glm_llr_test(gene_data, conditions, dist="nb", debug=args.debug)
            trimmed.loc[i, "P_Value"] = pvalue
        except KeyboardInterrupt:
            break
    
    print()  # ^time estimate ended with return character so this prints a newline
    
    trimmed["P_Sig"] = np.logical_and(trimmed["P_Value"]<0.05, trimmed["P_Value"]!=0)
    trimmed["Q_Value"] = bh_procedure(np.nan_to_num(trimmed["P_Value"]))
    trimmed["Q_Sig"] = np.logical_and(trimmed["Q_Value"]<0.05, trimmed["Q_Value"]!=0)

    sig_genes = trimmed["P_Sig"].sum()
    print("Significant p-values : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/len(trimmed)))
    print("Nan pvalues : {}".format(np.sum(np.isnan(trimmed["P_Value"]))))
    print("Zero pvalues : {}".format(np.sum(trimmed["P_Value"]==0)))
    sig_genes = trimmed["Q_Sig"].sum()
    print("Significant q-values : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/len(trimmed)))
    print("Zero qvalues : {}".format(np.sum(trimmed["Q_Value"]==0)))

    cutoff = 1e-4
    trimmed["Log10Q"] = -np.log10(trimmed[trimmed["Q_Value"]>cutoff]["Q_Value"])
    trimmed["Log10P"] = -np.log10(trimmed[trimmed["P_Value"]>cutoff]["P_Value"])
    # print(trimmed.sort_values(by="P_Value", ascending=False))

    # TODO : create booleans for filtering against the metrics
    # abs(log2fc) > 1, pvalue<0.05, qvalue<0.05

    # Save the comparison
    pairwise_filename = "data/{}/analysis/pairwise.csv".format(args.experiment)
    print(" * Saved pairwise analysis to {}".format(pairwise_filename))
    trimmed.to_csv(pairwise_filename, header=True, index=False)

    # print(trimmed.sort_values(by="Sample_Fitness", ascending=False)[["Gene_ID", "Sample_Fitness"]])
    # exit()

    # Plotting below
    # These are scatter plots for now
    combine_plots = [ # x, y suffix, s, xlog, ylog
        ["Gene_Length", "Diversity", 2, True, False],
        ["Gene_Length", "Hits", 1, True, True],
        ["Start", "Diversity", 8, False, False],
        ["Start", "Hits", None, False, False],
    ]
    single_plots = [ # x, y, s, xlog, ylog
        ["Gene_Length", "TA_Count", 2, True, True],
        ["Control_Hits", "Sample_Hits", 4, True, True],
        ["Start", "Log2FC", None, False, False],
        ["Start", "LinearDiff", None, False, False],
        ["Log2FC", "LinearDiff", None, False, False],
        ["Log2FC", "Log2_Diversity_Ratio", None, False, False],
        ["Log2FC", "Log10P", 8, False, False],
        ["Log2FC", "Log10Q", 8, False, False],
    ]
    if args.gc:
        combine_plots.append(["GC", "Hits", 4, False, True])
        single_plots.append(["Start", "GC", None, False, False])

    if args.plot:
        print(" * Calculating generating plots...")
        for x, y, s, xlog, ylog in combine_plots:
            print("Plotting x={} y={}".format(x, y))
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
            print("Plotting x={} y={}".format(x, y))
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
