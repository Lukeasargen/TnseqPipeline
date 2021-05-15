
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from util import tamap_to_genehits
from util import total_count_norm, quantile_norm
from util import length_norm

# Pairwise Comparisons
# Supports biological replicates
controls = [
    # "c25k_sum",
    # "c50k_sum",
    # "c100k_sum",
    # "c200k_sum",
    # "c400k_sum",
    # "c800k_sum",
    # "c1600k_sum",
    "c3200k_sum",
]
samples = [
    # "s25k_sum",
    # "s50k_sum",
    # "s100k_sum",
    # "s200k_sum",
    # "s400k_sum",
    # "s800k_sum",
    # "s1600k_sum",
    "s3200k_sum",
]

print("Controls :", controls)
print("Samples :", samples)

# TODO : save all the inputs and outputs to txt, also save input cmd for cli

# Load the tamap 
tamap_filename = r"data/demo/maps/14028c_TAmaps.csv"
tamap = pd.read_csv(tamap_filename, delimiter=",")
fasta_filename = r"data/demo/references/14028s_chromosome.fasta"
genehits = tamap_to_genehits(tamap)

# print(genehits.sort_values(by="Gene_Length", ascending=False)[["Gene_ID", "Gene_Length"]])
# exit()

# Normalization here
# genehits = length_norm(genehits, length="Gene_Length")
genehits = total_count_norm(genehits)
# genehits = quantile_norm(genehits, q=0.75)


# TODO : PCA here
# compare biological replicates
# Check for outliers here
# Cook's distance cutoff- outliers between samples


# For now, simply mean the reads together
genehits["Control_Hits"] = genehits[controls].mean(axis=1)
genehits["Sample_Hits"] = genehits[samples].mean(axis=1)

# TODO : renorm after combining replicates?
# genehits = total_count_norm(genehits, columns=["Control_Hits", "Sample_Hits"])


# Pairwise analysis below
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


# NOTE : Below starts the actual statistical analysis, everything above was just building the table

# Count thresholding
min_count = 0
trimmed = genehits[ (genehits[["Control_Hits","Sample_Hits"]] > min_count).all(axis=1) ].copy()

# Save which gene were cut by min threshold
removed = len(genehits)-len(trimmed)
print("Count thresholding={}. {} Genes Removed. {:.2f}% dropped. {} Genes Remaining.".format(min_count, removed, 100*removed/len(genehits), len(trimmed)))

# Ratio, Log2FC, LinearDiff
trimmed["Mean"] = (trimmed["Sample_Hits"] + trimmed["Control_Hits"]) / 2.0
smoothing = 0.0
trimmed["Ratio"] = (trimmed["Sample_Hits"] + smoothing) / (trimmed["Control_Hits"] + smoothing)
trimmed["Log2FC"] = np.log2(trimmed["Ratio"])
trimmed["LinearDiff"] = trimmed["Sample_Hits"] - trimmed["Control_Hits"]

# TODO : Fitness?
# trimmed["Sample_Fitness"] = np.log10( 1 + (trimmed["Sample_Hits"]*12e6/trimmed["Sample_Hits"].sum()) / (1000*trimmed["Gene_Length"]) )
trimmed["Sample_Fitness"] = np.log10( 1000*(1+trimmed["Sample_Hits"])/trimmed["Gene_Length"] )

# TODO : statistical genehits here, get p-value

# print(trimmed.sort_values(by="Sample_Fitness", ascending=False)[["Gene_ID", "Sample_Fitness"]])
# exit()

# Save the comparison
trimmed_filename = "data/demo/output.csv"
trimmed.to_csv(trimmed_filename, header=True, index=False)


# Column stats
for n in ["Control_Hits", "Sample_Hits", "Log2FC",
        "Control_Diversity", "Sample_Diversity", "Gene_Length"]:
    print( "\nColumn Stats : {}".format(n) )
    print( "Max = {}".format(trimmed[n].max()) )
    print( "Min = {}".format(trimmed[n].min()) )
    print( "Median = {}".format(trimmed[n].median()) )
    print( "Mean = {}".format(trimmed[n].mean()) )


# Plotting below

# TODO : In future
# Gene length axis : usually 0 to 5000

# These are scatter plots for now
combine_plots = [ # x, y suffix, s, xlog, ylog
    ["Gene_Length", "Unique_Insertions", 2, True, False],
    ["Gene_Length", "Diversity", 2, True, True],
    ["Gene_Length", "Hits", 1, True, True],
    ["TA_Count", "Unique_Insertions", 2, True, False],
    ["TA_Count", "Hits", 1, True, True],
    ["Start", "Diversity", 8, False, False],
    ["Start", "Hits", None, False, False],
]
single_plots = [ # x, y, s, xlog, ylog
    ["Gene_Length", "TA_Count", 2, True, True],
    ["Start", "Log2FC", None, False, False],
    ["Start", "LinearDiff", None, False, False],
    ["Log2FC", "LinearDiff", None, False, False],
    ["Control_Hits", "Sample_Hits", None, True, True],
]

if type(genehits.loc[0, "GC"])==float:
    combine_plots.append(["GC", "Hits", 4, False, True])

for x, y, s, xlog, ylog in combine_plots:
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
    plt.savefig(f"data/demo/{x}_vs_{y}.png")

for x, y, s, xlog, ylog in single_plots:
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
    plt.savefig(f"data/demo/{x}_vs_{y}.png")


# MA plot
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
plt.savefig("data/demo/MA_Plot.png")



