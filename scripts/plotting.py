import numpy as np
import matplotlib.pyplot as plt


def pairwise_plots(table, output_folder, alpha=0.05):
    """ Makes a bunch scatter plots and histograms """
    combine_plots = [ # x, y suffix, s, xlog, ylog
        ["Gene_Length", "Diversity", 2, True, False],
        ["Gene_Length", "Hits", 2, True, True],
        ["TA_Count", "Diversity", 4, True, False],
        ["TA_Count", "Hits", 4, True, True],
        ["Start", "Diversity", 8, False, False],
        ["Start", "Hits", None, False, True],
        ["GC", "Diversity", 4, False, True],
        ["GC", "Hits", 4, False, True],
        ["Log2FC_Reads", "Diversity", 4, False, False],
    ]
    single_plots = [ # x, y, s, xlog, ylog
        ["Gene_Length", "TA_Count", 6, True, True],
        ["Control_Hits", "Sample_Hits", 6, True, True],
        ["Start", "GC", None, False, False],
        # Reads differences
        ["Start", "Log2FC_Reads", None, False, False],
        ["Start", "LinearDiff_Reads", None, False, False],
        ["Log2FC_Reads", "LinearDiff_Reads", None, False, False],
        # Inserts differences
        ["Start", "Log2FC_Inserts", None, False, False],
        ["Start", "LinearDiff_Inserts", None, False, False],
        ["Log2FC_Inserts", "LinearDiff_Inserts", None, False, False],
        # Survival Index
        ["Start", "Log2SI", None, False, False],
        ["Log2FC_Reads", "Log2SI", None, False, False],
        # Fitness
        ["Start", "Sample_Fitness", None, False, False],
        ["Log2FC_Reads", "Sample_Fitness", None, False, False],
    ]
    hist_plots = [ # column, bins
        ["Log2SI", 100],
        ["Log2FC_Reads", 100],
        ["Log2FC_Inserts", 100],
        ["Sample_Fitness", 100],
        ["P_Value", 50],
        ["Q_Value", 50],
    ]

    # print(table.columns)
    # f = table.boxplot(column=['Control_Hits', 'Sample_Hits'], showfliers=False)
    # plt.savefig("temp.png")
    # exit()

    colors = {False:'tab:green', True:'tab:red'}
    p_sig_colors = table["P_Sig"].map(colors)

    for x, y, s, xlog, ylog in combine_plots:
        if not all(i in table.columns for i in [x, f"Control_{y}", f"Sample_{y}"]):
            continue
        print("Plotting Combined: x={} y={}".format(x, y))
        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(111)
        ax.scatter(x=table[x], y=table[f"Control_{y}"], s=s, alpha=0.7, color="tab:green", label="Control")
        ax.scatter(x=table[x], y=table[f"Sample_{y}"], s=s, alpha=0.7, color="tab:red", label="Sample")
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if xlog: ax.set_xscale('log')
        if ylog: ax.set_yscale('log')
        plt.legend()
        fig.tight_layout()
        plt.savefig(f"{output_folder}/{x}_vs_{y}.png", bbox_inches="tight")
        plt.close(fig)

    for x, y, s, xlog, ylog in single_plots:
        if not all(i in table.columns for i in [x,y]):
            continue
        print("Plotting: x={} y={}".format(x, y))
        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(111)
        ax.scatter(x=table[x], y=table[y], s=s, alpha=0.7, color=p_sig_colors)
        # plt.hlines(table[y].median(), xmin=table[x].min(), xmax=table[x].max(), color="tab:red", label="Median={:.2f}".format(table[y].median()))
        # plt.hlines(table[y].mean(), xmin=table[x].min(), xmax=table[x].max(), color="tab:blue", label="Mean={:.2f}".format(table[y].mean()))
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if xlog: ax.set_xscale('log')
        if ylog: ax.set_yscale('log')
        fig.tight_layout()
        plt.savefig(f"{output_folder}/{x}_vs_{y}.png", bbox_inches="tight")
        plt.close(fig)

    for col, bins in hist_plots:
        if col not in table.columns:
            continue
        print("Plotting Histogram: col={}".format(col))
        fig = plt.figure(figsize=[12, 8])
        ax = fig.add_subplot(111)
        df = table[col].replace([np.inf, -np.inf], np.nan, inplace=False)
        df.plot.hist(bins=bins)
        plt.xlabel(f"{col}")
        plt.savefig(f"{output_folder}/{col}_hist.png", bbox_inches="tight")
        plt.close(fig)


    # Make the MA plot
    print("Plotting MA plot")
    fig = plt.figure(figsize=[12, 8])
    A = 0.5 * ( np.log2(table["Sample_Hits"]) + np.log2(table["Control_Hits"]) )
    M = table["Log2FC_Reads"]
    plt.scatter(x=A, y=M, s=10, alpha=0.9, color=p_sig_colors)
    plt.hlines(M.median(), xmin=A.min(), xmax=A.max(), color="tab:red", label="Median={:.2f}".format(M.median()))
    plt.hlines(M.mean(), xmin=A.min(), xmax=A.max(), color="tab:blue", label="Mean={:.2f}".format(M.mean()))
    plt.hlines([-1, 1], xmin=A.min(), xmax=A.max(), color="tab:blue", linestyle='dashed')
    plt.legend(loc="upper left")
    plt.xlabel("A = 1/2 * ( log2(sample) + log2(control) )")
    plt.ylabel("M = log2(sample) - log2(control)")
    fig.tight_layout()
    plt.savefig(f"{output_folder}/MA_Plot.png", bbox_inches="tight")
    plt.close(fig)


    # Make the volcano plot
    print("Plotting Volcano plot")
    for col in ["Log10P", "Log10Q"]:
        if col not in table.columns:
            continue
        fig = plt.figure(figsize=[12, 8])
        X = table["Log2FC_Reads"]
        Y = table[col]
        plt.scatter(x=X, y=Y, s=8, alpha=0.9, color=p_sig_colors)
        plt.hlines(-np.log10(alpha), xmin=X.min(), xmax=X.max(), color="tab:blue", linestyle='dashed')
        plt.vlines([-1, 1], ymin=Y.min(), ymax=Y.max(), color="tab:blue", linestyle='dashed')
        plt.xlabel("log2 fold Change")
        plt.ylabel(f"-log10({col})")
        fig.tight_layout()
        plt.ylim(0, min(30, Y.max()))
        plt.savefig(f"{output_folder}/Volcano_Plot_{col}.png", bbox_inches="tight")
        plt.ylim(0, min(4, Y.max()))
        plt.savefig(f"{output_folder}/Volcano_Plot_{col}_Trim.png", bbox_inches="tight")
        plt.close(fig)

