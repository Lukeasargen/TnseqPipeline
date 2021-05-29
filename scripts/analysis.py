import os
import sys
import time
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon

from util import tamap_to_genehits, column_stats
from util import total_count_norm, quantile_norm, gene_length_norm
from util import ttr_norm
from util import bh_procedure
from util import time_to_string

from zinb_glm import zinb_glm_llr_test


def get_args():
    parser = argparse.ArgumentParser(description="Pairwise Comparison (Supports Replicates).")
    parser.add_argument('--experiment', type=str, required=True,
        help="Experiment folder name.")
    parser.add_argument('--index', type=str, required=True,
        help="Index name.")
    parser.add_argument('--controls', nargs='+', type=str, required=True,
        help="List read names without the filetype and separated by a space.")
    parser.add_argument('--samples', nargs='+', type=str, required=True,
        help="List read names without the filetype and separated by a space.")
    parser.add_argument('--output', type=str, default="default",
        help="Output name. Analysis outputs to a folder with this name. default=default.")

    parser.add_argument('--debug', default=False, action='store_true',
        help="Boolean flag that outputs all my debugging messages. default=False.")
    parser.add_argument('--plot', default=False, action='store_true',
        help="Boolean flag that automatically makes a few plots of the data. default=False.")
    parser.add_argument('--gc', default=False, action='store_true',
        help="Boolean flag that calculates the GC content of each gene. Not used in any test, but it makes some plots and gets saved in the output. default=False.")

    parser.add_argument('--pooling', type=str, default="sum",
        choices=["sum", "average"],
        help="String argument. Sum or average the hits PER GENE to get a merged value for the expression at the gene level. default=sum.")
    parser.add_argument('--min_count', default=1, type=int,
        help="Integer argument. Threshold for lowest number of hits PER GENE after pooling. Removes genes with low pooled hits. These genes are not tested for significance and saved in a separate output table. default=1.")
    parser.add_argument('--min_inserts', default=2, type=int,
        help="Integer argument. Threshold for lowest number of insertion sites with hits BY GENE. Removes genes with low TA sites or low hit diversity (diversity=unique hit sites/total gene sites). These genes are not tested for significance and saved in a separate output table. default=2.")
    parser.add_argument('--smoothing', default=1, type=float,
        help="Float argument. Smooth the ratio for small counts. ratio=(sample+smoothing)/(control+smoothing). default=1.")
    parser.add_argument('--alpha', default=0.05, type=float,
        help="Float argument. Significance level. default=0.05.")
    parser.add_argument('--insert_weighting', default=False, action='store_true',
        help="Boolean flag that scales PER GENE based on unique inserts. It increases the Formula is new_hits=old_hits*(unique_inserts/average_unique_inserts). default=False.")

    args = parser.parse_args()
    return args


def pairwise_comparison(args):
    """ Compare the controls and samples.
        Support replicates. Input as lists.
    """

    output_folder = "data/{}/analysis/{}".format(args.experiment, args.output)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    command = "python3 " + " ".join(sys.argv)
    command_filename = "{}/run_command.txt".format(output_folder)
    with open(command_filename, "w") as f:
        f.write(command)

    # Supports biological replicates
    print("Controls :", args.controls)
    print("Samples :", args.samples)

    # Append "_sum" to match the column names in the TAmap and Genehits
    # this matches the step done in readTAmap.py, where only the "_sum" is merged into the TAmap
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
    tamap = ttr_norm(tamap, trim=0.05, columns=test_columns, debug=args.debug)
    if args.debug: print("\nStats after norm:"); column_stats(tamap, columns=test_columns)

    print(" * Compressing TAmap into Genehits table...")
    fasta_filename = "data/{}/references/{}.fasta".format(args.experiment, args.index) if args.gc else None
    if args.debug: print(f"GC content fasta filename : {fasta_filename}")
    genehits = tamap_to_genehits(tamap, fasta_filename=fasta_filename, pooling=args.pooling)

    # Combine the replicates, average the counts
    print(" * Combining replicates...")
    genehits["Control_Hits"] = genehits[controls].mean(axis=1)
    genehits["Sample_Hits"] = genehits[samples].mean(axis=1)
    column_stats(genehits, columns=["Control_Hits", "Sample_Hits"])

    # NOTE : Pairwise analysis below

    # TODO : can this be moved into tamap_to_genehits
    print(" * Calculating insertion density per gene...")
    # Need remove intergenic and group by gene
    temp = tamap[tamap['Gene_ID'].notna()].copy()  # Remove intergenic
    temp["Control"] = temp[controls].mean(axis=1)  # Combine control replicates
    temp["Sample"] = temp[samples].mean(axis=1)  # Combine sample replicates
    grouped = temp.groupby("Gene_ID", as_index=False).nunique()  # Builtin method to get unique
    # Unique_Insertions : Unique TA hits / TA Sites
    # Diversity = Unique_Insertions / TA_Counts
    # Subtract 1 bc it counts the 0 counts as unique but we don't care about this, genes with 0 hits should have 0 unique insertions
    genehits["Control_Unique_Insertions"] = grouped["Control"]-1
    genehits["Sample_Unique_Insertions"] = grouped["Sample"]-1
    genehits["Control_Diversity"] = genehits["Control_Unique_Insertions"] / genehits["TA_Count"]
    genehits["Sample_Diversity"] = genehits["Sample_Unique_Insertions"] / genehits["TA_Count"]
    # Don't need these anymore
    del temp
    del grouped
    column_stats(genehits, columns=["Control_Diversity", "Sample_Diversity"])

    # TSAS Weighting
    if args.insert_weighting:
        """ https://github.com/srimam/TSAS """
        print(" * Calculating insertion weighting (TSAS by Saheed Imam)...")
        if args.debug: print("\nStats before insertion weighting:"); column_stats(genehits, columns=["Control_Hits","Sample_Hits"])
        avg_unique = (genehits["Control_Unique_Insertions"].mean()+genehits["Sample_Unique_Insertions"].mean()) / 2.0
        genehits["Control_Hits"] = genehits["Control_Hits"] * (genehits["Control_Unique_Insertions"]/avg_unique)
        genehits["Sample_Hits"] = genehits["Sample_Hits"] * (genehits["Sample_Unique_Insertions"]/avg_unique)
        if args.debug: print("\nStats after insertion weighting:"); column_stats(genehits, columns=["Control_Hits","Sample_Hits"])

    # NOTE : Below starts the actual statistical analysis, everything above was just building the table

    print(" * Calculating differential statistics...")
    # Reads differences
    genehits["Ratio_Reads"] = (genehits["Sample_Hits"] + args.smoothing) / (genehits["Control_Hits"] + args.smoothing)
    genehits["Log2FC_Reads"] = np.log2(genehits["Ratio_Reads"])
    genehits["LinearDiff_Reads"] = genehits["Sample_Hits"] - genehits["Control_Hits"]
    # Inserts differences
    genehits["Ratio_Inserts"] = (genehits["Sample_Unique_Insertions"] + args.smoothing) / (genehits["Control_Unique_Insertions"] + args.smoothing)
    genehits["Log2FC_Inserts"] = np.log2(genehits["Ratio_Inserts"])
    genehits["LinearDiff_Inserts"] = genehits["Sample_Unique_Insertions"] - genehits["Control_Unique_Insertions"]


    # Fitness
    print(" * Calculating fitness...")
    # genehits["Sample_Fitness"] = np.log10( 1 + (genehits["Sample_Hits"]*12e6/genehits["Sample_Hits"].sum()) / (1000*genehits["Gene_Length"]) )
    # genehits["Sample_Fitness"] = np.log10( (1+genehits["Sample_Hits"])/genehits["Gene_Length"] )
    # genehits["Sample_Fitness"] = (np.log(genehits["Sample_Diversity"]) - np.log(genehits["Control_Diversity"])) / (np.log(1-genehits["Sample_Diversity"]) - np.log(1-genehits["Control_Diversity"]))
    # column_stats(genehits, columns=["Sample_Fitness"])


    # Survival Index
    print(" * Calculating survival index...")
    dval_ctl = (genehits["Control_Hits"]*genehits["Gene_Length"].sum()) / (genehits["Gene_Length"]*genehits["Control_Hits"].sum())
    dval_exp = (genehits["Sample_Hits"]*genehits["Gene_Length"].sum()) / (genehits["Gene_Length"]*genehits["Sample_Hits"].sum())
    si = dval_exp/dval_ctl
    genehits["SI"] = si
    genehits["Log2SI"] = np.log2(genehits["SI"])

    # Count and Insertion thresholding
    print(" * Trimming away low saturated genes...")

    # First, find "possibly essential" genes
    # TODO : citation. Tn-seq; high-throughput parallel sequencing for fitness and genetic interaction studies in microorganisms
    possibly = (genehits[["Control_Unique_Insertions","Sample_Unique_Insertions"]] < 3).any(axis=1) & (genehits["Gene_Length"] >= 400)
    possibly_filename = "{}/possibly_essential.csv".format(output_folder)
    print("{}/{}({:.2f}%) possibly essential genes.".format( possibly.sum(), len(genehits), 100*possibly.sum()/len(genehits) ) )

    # Save the possibly essential genes
    print(" * Saved possibly essential genes to {}".format(possibly_filename))
    genehits[possibly].to_csv(possibly_filename, header=True, index=False)

    # This bool saves whether the gene has counts in all groups
    hit_bool = ~(genehits[["Control_Unique_Insertions", "Sample_Unique_Insertions"]] == 0).any(axis=1)
    num_hit = hit_bool.sum()
    print("{}/{}({:.2f}%) had no hits in one group. {}/{}({:.2f}%) had at least one hit in all groups.".format( len(genehits)-num_hit, len(genehits), 100*(len(genehits)-num_hit)/len(genehits), num_hit, len(genehits), 100*num_hit/len(genehits) ))

    # Save genes that had no insertions in at least one group
    no_hit_filename = "{}/no_hits.csv".format(output_folder)
    print(" * Saved genes with no hits to {}".format(no_hit_filename))
    genehits[~hit_bool].to_csv(no_hit_filename, header=True, index=False)

    # This bool saves the user defined thresholds
    keep = (genehits[["Control_Hits","Sample_Hits"]] >= args.min_count).all(axis=1) & (genehits[["Control_Unique_Insertions","Sample_Unique_Insertions"]] >= args.min_inserts).all(axis=1)
    # Separate the genes that will be tested from the rest
    trimmed = genehits[keep&hit_bool].copy()
    removed = genehits[~keep&hit_bool].copy()
    print(" * Thresholds: min_count={}. min_inserts={}.".format(args.min_count, args.min_inserts))
    print("{}/{}({:.2f}%) more genes removed by threshold. {}/{}({:.2f}%) genes remaining.".format( len(removed), len(genehits), 100*len(removed)/len(genehits), len(trimmed), len(genehits), 100*len(trimmed)/len(genehits) ))
    if args.debug: column_stats(trimmed, columns=["Control_Hits", "Sample_Hits", "Control_Unique_Insertions", "Sample_Unique_Insertions"])

    # Save genes that are removed
    removed_filename = "{}/removed.csv".format(output_folder)
    print(" * Saved removed genes to {}".format(removed_filename))
    removed.to_csv(removed_filename, header=True, index=False)


    print(" * Calculating statistical significance...")
    print("Stop early by pressing Ctrl+C in the terminal.")
    t0 = time.time()  # Start time
    trimmed["P_Value"] = np.nan
    trimmed["P_Sig"] = False
    c = 0  # genes counter, this is different than the index value
    for i in trimmed.index:
        try:
            # Helpful time estimate
            c += 1
            if c%10 == 0:
                duration = time.time()-t0
                remaining = duration/c * (len(trimmed)-c)
                print("gene {}/{}. {:.1f} genes/second. elapsed={}. remaining={}.".format(c, len(trimmed), c/duration, time_to_string(duration), time_to_string(remaining)), end="\r")
            # # gene_name is used to index the full TAmap 
            # # size is used to get the length of the condition array
            gene_name, size = trimmed.loc[i][["Gene_ID", "TA_Count"]]
            df = tamap[tamap["Gene_ID"]==gene_name]
            if args.debug: print("\nStatistical Test for: ", gene_name)

            # gene_data = np.array(df[test_columns]).T.reshape(-1)
            # conditions = np.array([0]*size*len(controls) + [1]*size*len(samples))
            # pvalue = zinb_glm_llr_test(gene_data, conditions, dist="nb", debug=args.debug)
            
            data1 = np.array(df[controls].mean(axis=1))  # Combine control replicates
            data2 = np.array(df[samples].mean(axis=1))  # Combine sample replicates
            # u_stat, pvalue = mannwhitneyu(data1, data2)
            t_stat, pvalue = ttest_ind(data1, data2)
            # t_stat, pvalue = wilcoxon(data1, data2)

            trimmed.loc[i, "P_Value"] = pvalue

            if args.debug: print(trimmed.loc[i])
            if args.debug and c > 7: break        
        except KeyboardInterrupt:
            break

    # ^time estimate ended with return character so this prints a newline
    duration = time.time()-t0
    remaining = duration/c * (len(trimmed)-c)
    print("gene {}/{}. {:.1f} genes/second. elapsed={}. remaining={}.".format(c, len(trimmed), c/duration, time_to_string(duration), time_to_string(remaining)))


    trimmed["P_Sig"] = np.logical_and(trimmed["P_Value"]<args.alpha, trimmed["P_Value"]!=0)
    trimmed["Log10P"] = -np.log10(trimmed["P_Value"])
    sig_genes = trimmed["P_Sig"].sum()
    print("Significant p-values : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/len(trimmed)))
    print("Genes not tested : {}".format(np.sum(np.isnan(trimmed["P_Value"]))))
    print("Test failures : {}".format(np.sum(trimmed["P_Value"]==0)))
   
    print(" * Adjusting p-values for multiple test...")
    trimmed["Q_Value"] = bh_procedure(np.nan_to_num(trimmed["P_Value"]))
    trimmed["Q_Sig"] = np.logical_and(trimmed["Q_Value"]<args.alpha, trimmed["Q_Value"]!=0)
    trimmed["Log10Q"] = -np.log10(trimmed["Q_Value"])
    sig_genes = trimmed["Q_Sig"].sum()
    print("Significant q-values : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/len(trimmed)))

    # Save the comparison
    pairwise_filename = "{}/pairwise.csv".format(output_folder)
    print(" * Saved pairwise analysis to {}".format(pairwise_filename))
    trimmed.to_csv(pairwise_filename, header=True, index=False)

    # print(trimmed.sort_values(by="P_Value", ascending=False))
    # print(trimmed.sort_values(by="Sample_Fitness", ascending=False)[["Gene_ID", "Sample_Fitness"]])
    # exit()

    # TODO : move plotting to separate script, import functions to use here
    # Plotting below
    # These are scatter plots for now
    combine_plots = [ # x, y suffix, s, xlog, ylog
        ["Gene_Length", "Diversity", 2, True, False],
        ["Gene_Length", "Hits", 1, True, True],
        ["Start", "Diversity", 8, False, False],
        ["Start", "Hits", None, False, True],
    ]
    single_plots = [ # x, y, s, xlog, ylog
        ["Gene_Length", "TA_Count", 2, True, True],
        ["Control_Hits", "Sample_Hits", 4, True, True],
        # Reads differences
        ["Start", "Log2FC_Reads", None, False, False],
        ["Start", "LinearDiff_Reads", None, False, False],
        ["Log2FC_Reads", "LinearDiff_Reads", None, False, False],
        # Inserts differences
        ["Start", "Log2FC_Inserts", None, False, False],
        ["Start", "LinearDiff_Inserts", None, False, False],
        ["Log2FC_Inserts", "LinearDiff_Inserts", None, False, False],

        ["Log2FC_Reads", "Log2SI", None, False, False],
    ]
    if args.gc:
        combine_plots.append(["GC", "Hits", 4, False, True])
        single_plots.append(["Start", "GC", None, False, False])

    colors = {False:'tab:green', True:'tab:red'}
    p_sig_colors = trimmed["P_Sig"].map(colors)

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
            plt.savefig(f"{output_folder}/{x}_vs_{y}.png")

        for x, y, s, xlog, ylog in single_plots:
            print("Plotting x={} y={}".format(x, y))
            fig = plt.figure(figsize=[16, 8])
            ax = fig.add_subplot(111)
            ax.scatter(x=trimmed[x], y=trimmed[y], s=s, color=p_sig_colors)
            # plt.hlines(trimmed[y].median(), xmin=trimmed[x].min(), xmax=trimmed[x].max(), color="tab:red", label="Median={:.2f}".format(trimmed[y].median()))
            # plt.hlines(trimmed[y].mean(), xmin=trimmed[x].min(), xmax=trimmed[x].max(), color="tab:blue", label="Mean={:.2f}".format(trimmed[y].mean()))
            ax.set_xlabel(x)
            ax.set_ylabel(y)
            if xlog: ax.set_xscale('log')
            if ylog: ax.set_yscale('log')
            fig.tight_layout()
            plt.savefig(f"{output_folder}/{x}_vs_{y}.png")

        hist_plots = ["Log2SI", "Log2FC_Reads", "Log2FC_Inserts"]
        for col in hist_plots:
            print("Plotting col={}".format(col))
            fig = plt.figure(figsize=[12, 8])
            ax = fig.add_subplot(111)
            trimmed["Log2SI"].plot.hist(bins=200)
            plt.xlabel(f"{col}")
            plt.savefig(f"{output_folder}/{col}_hist.png")

        # Make the MA plot
        print("Plotting MA plot")
        fig = plt.figure(figsize=[12, 8])
        A = 0.5 * ( np.log2(trimmed["Sample_Hits"]) + np.log2(trimmed["Control_Hits"]) )
        M = trimmed["Log2FC_Reads"]
        plt.scatter(x=A, y=M, s=10, color=p_sig_colors)
        plt.hlines(M.median(), xmin=A.min(), xmax=A.max(), color="tab:red", label="Median={:.2f}".format(M.median()))
        plt.hlines(M.mean(), xmin=A.min(), xmax=A.max(), color="tab:blue", label="Mean={:.2f}".format(M.mean()))
        plt.hlines([-1, 1], xmin=A.min(), xmax=A.max(), color="tab:blue", linestyle='dashed')
        plt.legend(loc="upper left")
        plt.xlabel("A = 1/2 * ( log2(sample) + log2(control) )")
        plt.ylabel("M = log2(sample) - log2(control)")
        fig.tight_layout()
        plt.savefig(f"{output_folder}/MA_Plot.png")


        # Make the volcano plot
        print("Plotting Volcano plot")
        fig = plt.figure(figsize=[12, 8])
        X = trimmed["Log2FC_Reads"]
        Y = trimmed["Log10P"]
        plt.scatter(x=X, y=Y, s=8, color=p_sig_colors)
        plt.hlines(-np.log10(args.alpha), xmin=X.min(), xmax=X.max(), color="tab:blue", linestyle='dashed')
        plt.vlines([-1, 1], ymin=Y.min(), ymax=Y.max(), color="tab:blue", linestyle='dashed')
        plt.xlabel("log2 fold Change")
        plt.ylabel("-log10(p-value)")
        fig.tight_layout()
        plt.savefig(f"{output_folder}/Volcano_Plot.png")
        plt.ylim(Y.min(), min(4, Y.max()))
        plt.savefig(f"{output_folder}/Volcano_Plot_Trim.png")



if __name__ == "__main__":
    args = get_args()
    pairwise_comparison(args)
