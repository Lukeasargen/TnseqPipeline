import os
import sys
import time
import argparse
from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon

from plotting import pairwise_plots
from util import tamap_to_genehits, column_stats
from util import total_count_norm, quantile_norm, gene_length_norm, ttr_norm, nzmean_norm
from util import exclude_sites_tamap
from util import bh_procedure
from util import time_to_string
from zinb_glm import zinb_glm_llr


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
    parser.add_argument('--norm', type=str, default="ttr",
        choices=["total", "quantile", "ttr", "nzmean"],
        help="String argument. Choose the normalization between total count, quantile, ttr (Total Trimmed Reads), nzmean (Non-Zero Mean). default=ttr.")
    parser.add_argument('--quantile', default=0.75, type=float,
        help="Float argument. Quantile used to normalize. default=0.75.")
    parser.add_argument('--ttr', default=0.05, type=float,
        help="Float argument. Percentage of the highest and lowest values which are excluded before calculating the mean. default=0.05.")
    parser.add_argument('--strand', type=str, default="both",
        choices=["both", "forward", "reverse"],
        help="String argument. Specify strand for analysis. default=both.")
    parser.add_argument('--stat', type=str, default=None,
        choices=["mannu", "ttest", "wilcoxon", "zinb"],
        help="String argument. Choose the statistical test between mannu (Mann Whitney U test), ttest (T test), wilcoxon, zinb (Zero-Inflated Binomial Regression). default=None.")
    parser.add_argument('--alpha', default=0.05, type=float,
        help="Float argument. Significance level. default=0.05.")
    parser.add_argument('--min_count', default=1, type=int,
        help="Integer argument. Threshold for lowest number of insertions PER GENE after pooling. Removes genes with insertions. These genes are not tested for significance and saved in a separate output table. default=1.")
    parser.add_argument('--min_inserts', default=2, type=int,
        help="Integer argument. Threshold for lowest number of insertion sites with hits BY GENE. Removes genes with low hit diversity (unique insertion sites). These genes are not tested for significance and saved in a separate output table. default=2.")
    parser.add_argument('--min_sites', default=0, type=int,
        help="Integer argument. Threshold for lowest number of insertion sites with hits BY GENE. Removes genes with low TA sites (TA_Count < min_sites). These genes are not tested for significance and saved in a separate output table. default=0.")
    parser.add_argument('--pooling', type=str, default="sum",
        choices=["sum", "average"],
        help="String argument. Sum or average the hits PER GENE to get a merged value for the expression at the gene level. default=sum.")
    parser.add_argument('--smoothing', default=1, type=float,
        help="Float argument. Smooth the ratio for small counts. ratio=(sample+smoothing)/(control+smoothing). default=1.")
    parser.add_argument('--expansion', default=250, type=float,
        help="Float argument. Expansion factor used in the fitness formula described by Opijnen (Nature 2009). default=250.")
    parser.add_argument('--insert_weighting', default=False, action='store_true',
        help="Boolean flag that scales PER GENE based on unique inserts. The Formula is new_hits=old_hits*(unique_inserts/average_unique_inserts). default=False.")
    parser.add_argument('--length_norm', default=False, action='store_true',
        help="Boolean flag that scales PER GENE based on gene length. default=False.")
    parser.add_argument('--gc', default=False, action='store_true',
        help="Boolean flag that calculates the GC content of each gene. Not used in any test, but it makes some plots and gets saved in the output. default=False.")
    parser.add_argument('--ef', "--exclude_first", default=0, type=float, dest="exclude_first",
        help="Float argument. Exclude insertions in the first X percent of the gene. default=0.")
    parser.add_argument('--el', "--exclude_last", default=0, type=float, dest="exclude_last",
        help="Float argument. Exclude insertions in the last X percent of the gene. default=0.")

    args = parser.parse_args()

    if args.expansion <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid expansion factor." % args.expansion)

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

    # Append "_"+args.strand to match the column names in the TAmap and Genehits
    controls = [c+"_"+str(args.strand) for c in args.controls]
    samples = [c+"_"+str(args.strand) for c in args.samples]
    test_columns = controls+samples

    # Load the tamap 
    tamap_filename = "data/{}/maps/{}_TAmaps.csv".format(args.experiment, args.index)
    print(" * Loading TAmap : {}".format(tamap_filename))
    tamap = pd.read_csv(tamap_filename, delimiter=",")

    print(" * Normalizing...")
    if args.debug: print("\nStats before norm:"); column_stats(tamap, columns=test_columns)

    if args.length_norm:
        tamap = gene_length_norm(tamap, columns=test_columns, debug=args.debug)

    norms = {
        "total": partial(total_count_norm, tamap, columns=test_columns, debug=args.debug),
        "quantile": partial(quantile_norm, tamap, q=args.quantile, columns=test_columns, debug=args.debug),
        "ttr": partial(ttr_norm, tamap, trim=args.ttr, columns=test_columns, debug=args.debug),
        "nzmean": partial(nzmean_norm, tamap, columns=test_columns, debug=args.debug),
    }
    tamap = norms[args.norm]()
    
    if args.debug: print("\nStats after norm:"); column_stats(tamap, columns=test_columns)

    print(" * Compressing TAmap into Genehits table...")
    fasta_filename = "data/{}/references/{}.fasta".format(args.experiment, args.index) if args.gc else None
    if args.debug: print(f"GC content fasta filename : {fasta_filename}")

    # Removes sites that are in the first or last percentage of a gene
    tamap = exclude_sites_tamap(tamap, exclude_first=args.exclude_first, exclude_last=args.exclude_last)

    print(" * Merging into genehits table...")
    # NOTE : This function removes intergenic regions
    # This might mess up the normalization totals, but this is likely fine
    # since the normalization above accounts for sequencing depth of the entire
    # library and not just for the gene regions.
    genehits = tamap_to_genehits(tamap, fasta_filename=fasta_filename, pooling=args.pooling)

    # Combine the replicates, average the counts
    print(" * Combining replicates...")
    tamap["Control_Hits"] = tamap[controls].mean(axis=1)
    tamap["Sample_Hits"] = tamap[samples].mean(axis=1)
    genehits["Control_Hits"] = genehits[controls].mean(axis=1)
    genehits["Sample_Hits"] = genehits[samples].mean(axis=1)
    column_stats(genehits, columns=["Control_Hits", "Sample_Hits"])

    # NOTE : Pairwise analysis below

    # TODO : can this be moved into tamap_to_genehits
    print(" * Calculating insertion density per gene...")
    # Need remove intergenic and group by gene
    temp = tamap[tamap['Gene_ID'].notna()].copy()  # Remove intergenic
    temp["Control"] = temp[controls].mean(axis=1).astype(bool)  # Combine control replicates
    temp["Sample"] = temp[samples].mean(axis=1).astype(bool)  # Combine sample replicates
    grouped = temp.groupby("Gene_ID", as_index=False).sum()

    # Unique_Insertions : Unique TA hits / TA Sites
    # Diversity = Unique_Insertions / TA_Counts
    genehits["Control_Unique_Insertions"] = grouped["Control"]
    genehits["Sample_Unique_Insertions"] = grouped["Sample"]
    genehits["Control_Diversity"] = genehits["Control_Unique_Insertions"] / genehits["TA_Count"]
    genehits["Sample_Diversity"] = genehits["Sample_Unique_Insertions"] / genehits["TA_Count"]
    # Don't need these anymore
    del temp
    del grouped
    column_stats(genehits, columns=["Control_Diversity", "Sample_Diversity"])

    # TSAS Weighting
    if args.insert_weighting:
        """ https://github.com/srimam/TSAS """
        print(" * Calculating insertion weighting (Idea from : TSAS by Saheed Imam)...")
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


    print(" * Calculating fitness...")
    """
    Idea from : "Tn-seq; high-throughput parallel sequencing for fitness and genetic interaction studies in microorganisms"
    As far as I can tell, expansion factor (args.expansion) is an arbitrary value.
    I explored it's effects on the fitness formula in desmos.
    Formula = fitness = ln(ef*(sf)/(cf)) / ln(ef*(1-sf)/(1-cf))
    where sf=sample frequency(t2), cf=control frequency(t1), and ef=expansion factor.
    Changing expansion factor appears to change the steepness of the slope
    around the neutral value. Larger expansion, smaller slope.
    """
    # Copy the columns so the original table is not changed
    df = pd.DataFrame()
    df[["Control_Hits", "Sample_Hits"]] = genehits[["Control_Hits", "Sample_Hits"]]
    # Re adjust the table to with the same behavior as vanOpijnenLab "correction factors"
    # "Total Counts" effectively makes the columns have the same sum, which is identical to the vanOpijnenLab procedure
    df = total_count_norm(df, columns=["Control_Hits", "Sample_Hits"])
    # I'm pretty sure frequency is just the site counts over total conts
    df["contrl_freq"] = df["Control_Hits"] / genehits["Control_Hits"].sum()
    df["sample_freq"] = df["Sample_Hits"] / genehits["Sample_Hits"].sum()
    # This is the formula from Opijnen (Nature 2009) Online Methods
    df["Sample_Fitness"] = np.log(args.expansion*(df["sample_freq"])/(df["contrl_freq"])) / np.log(args.expansion*(1-df["sample_freq"])/(1-df["contrl_freq"]))
    genehits["Sample_Fitness"] = df["Sample_Fitness"]
    genehits["Log2Fitness"] = np.log2(genehits["Sample_Fitness"])

    # Survival Index
    print(" * Calculating survival index...")
    # Idea from : "Genetic basis of persister tolerance to aminoglycosides (2015) Shan...Lewis"
    dval_ctl = (genehits["Control_Hits"]*genehits["Gene_Length"].sum()) / (genehits["Gene_Length"]*genehits["Control_Hits"].sum())
    dval_exp = (genehits["Sample_Hits"]*genehits["Gene_Length"].sum()) / (genehits["Gene_Length"]*genehits["Sample_Hits"].sum())
    si = dval_exp/dval_ctl
    # This is equivalent to above
    # si = (genehits["Sample_Hits"]*genehits["Control_Hits"].sum())/(genehits["Control_Hits"]*genehits["Sample_Hits"].sum())
    genehits["Survival_Index"] = si
    # This will throw an error if there are 0 sample hits, but whatever
    # numpy is good at handling errors, log2(0)=-inf
    genehits["Log2SI"] = np.log2(genehits["Survival_Index"])

    # First, find "possibly essential" genes
    print(" * Finding possibly essential genes...")
    # Idea from: "Tn-seq; high-throughput parallel sequencing for fitness and genetic interaction studies in microorganisms"
    possibly = (genehits[["Sample_Unique_Insertions"]] < 3).any(axis=1) & (genehits["Gene_Length"] >= 400)
    possibly_filename = "{}/possibly_essential.csv".format(output_folder)
    print("{}/{}({:.2f}%) possibly essential genes.".format( possibly.sum(), len(genehits), 100*possibly.sum()/len(genehits) ) )

    # Save the possibly essential genes
    print(" * Saved possibly essential genes to {}".format(possibly_filename))
    genehits[possibly].to_csv(possibly_filename, header=True, index=False)

    # Count and Insertion thresholding
    print(" * Trimming away low saturated genes by thresholds...")

    # This bool saves whether the gene has counts in all groups
    hit_bool = ~(genehits[["Control_Unique_Insertions", "Sample_Unique_Insertions"]] == 0).any(axis=1)
    num_hit = hit_bool.sum()
    print("{}/{}({:.2f}%) had no insertions. {}/{}({:.2f}%) had at least one insertion.".format( len(genehits)-num_hit, len(genehits), 100*(len(genehits)-num_hit)/len(genehits), num_hit, len(genehits), 100*num_hit/len(genehits) ))

    # Save genes that had no insertions in at least one group
    no_hit_filename = "{}/no_hits.csv".format(output_folder)
    print(" * Saved genes with no hits to {}".format(no_hit_filename))
    genehits[~hit_bool].to_csv(no_hit_filename, header=True, index=False)

    # Create boolean columns for user defined thresholds
    keep_count = (genehits[["Control_Hits","Sample_Hits"]] >= args.min_count).all(axis=1)
    keep_inserts = (genehits[["Control_Unique_Insertions","Sample_Unique_Insertions"]] >= args.min_inserts).all(axis=1)
    keep_sites = (genehits["TA_Count"] >= args.min_sites)
    keep = keep_count & keep_inserts & keep_sites
    # Separate the genes that will be tested from the rest
    trimmed = genehits[keep&hit_bool]
    removed = genehits[~keep&hit_bool]
    print(" * Thresholds: min_count={}. min_inserts={}. min_sites={}.".format(args.min_count, args.min_inserts, args.min_sites))
    print("{}/{}({:.2f}%) more genes removed by threshold. {}/{}({:.2f}%) genes remaining.".format( len(removed), len(genehits), 100*len(removed)/len(genehits), len(trimmed), len(genehits), 100*len(trimmed)/len(genehits) ))
    if args.debug: column_stats(trimmed, columns=["Control_Hits", "Sample_Hits", "Control_Unique_Insertions", "Sample_Unique_Insertions", "Control_Diversity", "Sample_Diversity"])

    # Save genes that are removed
    removed_filename = "{}/removed.csv".format(output_folder)
    print(" * Saved removed genes to {}".format(removed_filename))
    removed.to_csv(removed_filename, header=True, index=False)

    # Values per column for ZINB-GLM offsets
    nzmean = np.array( genehits[test_columns].replace(0, np.NaN).mean() )
    diversity = np.array( tamap[test_columns].astype(bool).sum()/len(tamap) )

    trimmed["P_Value"] = np.nan
    trimmed["P_Sig"] = False
    if args.stat:
        print(" * Calculating statistical significance...")
        print("Stop early by pressing Ctrl+C in the terminal.")
        t0 = time.time()  # Start time
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

                if args.stat=="zinb":
                    gene_data = np.array(df[test_columns]).T.reshape(-1)
                    conditions = np.array([0]*size*len(controls) + [1]*size*len(samples))
                    pvalue = zinb_glm_llr(gene_data, conditions, nzmean, diversity, dist="nb", rescale=0, debug=args.debug)
                else:
                    data1 = np.array(df[controls].mean(axis=1))  # Combine control replicates
                    data2 = np.array(df[samples].mean(axis=1))  # Combine sample replicates
                    if args.stat=="mannu":
                        u_stat, pvalue = mannwhitneyu(data1, data2)
                    elif args.stat=="ttest":
                        t_stat, pvalue = ttest_ind(data1, data2)
                    elif args.stat=="wilcoxon":
                        t_stat, pvalue = wilcoxon(data1, data2)

                if pvalue==0:
                    print(f"\n{gene_name} failed.")
                    print(trimmed.loc[i])
                trimmed.loc[i, "P_Value"] = pvalue

                if args.debug and c > 5: break
            except KeyboardInterrupt:
                break

        # ^time estimate ended with return character so this prints a newline
        duration = time.time()-t0
        remaining = duration/c * (len(trimmed)-c)
        print("gene {}/{}. {:.1f} genes/second. elapsed={}. remaining={}.".format(c, len(trimmed), c/duration, time_to_string(duration), time_to_string(remaining)))

    # Make a boolean column for the P-value and a negative log10 for the volcano plot
    trimmed["P_Sig"] = np.logical_and(trimmed["P_Value"]<args.alpha, trimmed["P_Value"]!=0)
    trimmed["Log10P"] = -np.log10(trimmed["P_Value"])
    sig_genes = trimmed["P_Sig"].sum()
    print("Significant p-values : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/len(trimmed)))
    print("Genes not tested : {}".format(np.sum(np.isnan(trimmed["P_Value"]))))
    fails = np.sum(trimmed["P_Value"]==0)
    print("Test failures : {} ({:.2f}%)".format(fails, 100*fails/len(trimmed)))
   
    # The same as above but for adjusted Q-values
    print(" * Adjusting p-values for multiple test...")
    qvalues, new_alpha = bh_procedure(np.nan_to_num(trimmed["P_Value"]))
    trimmed["Q_Value"] = qvalues
    trimmed["Q_Sig"] = np.logical_and(trimmed["Q_Value"]<args.alpha, trimmed["Q_Value"]!=0)
    trimmed["Log10Q"] = -np.log10(trimmed["Q_Value"])
    sig_genes = trimmed["Q_Sig"].sum()
    print("Significant q-values : {} ({:.2f}%)".format(sig_genes, 100*sig_genes/len(trimmed)))

    # Save the comparison
    pairwise_filename = "{}/pairwise.csv".format(output_folder)
    print(" * Saved pairwise analysis to {}".format(pairwise_filename))
    trimmed.to_csv(pairwise_filename, header=True, index=False)

    if args.plot:
        print(" * Generating plots...")
        pairwise_plots(trimmed, output_folder, args.alpha)


if __name__ == "__main__":
    args = get_args()
    pairwise_comparison(args)
