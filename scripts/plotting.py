import numpy as np
import matplotlib.pyplot as plt


def pairwise_plots(table, output_folder, alpha=0.05):
    """ Makes a bunch scatter plots and histograms """
    combine_plots = [ # x, y suffix, s, xlog, ylog
        ["Gene_Length", "Diversity", 2, True, False],
        ["Gene_Length", "Hits", 2, True, True],
        ["Start", "Diversity", 8, False, False],
        ["Start", "Hits", None, False, True],
        ["GC", "Hits", 4, False, True]
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
        ["Start", "Log2Fitness", None, False, False],
        ["Log2FC_Reads", "Log2Fitness", None, False, False],
    ]
    hist_plots = [ # column, bins
        ["Log2SI", 100],
        ["Log2FC_Reads", 100],
        ["Log2FC_Inserts", 100],
        ["Sample_Fitness", 100],
        ["P_Value", 20],
        ["Q_Value", 20],
    ]

    colors = {False:'tab:green', True:'tab:red'}
    p_sig_colors = table["P_Sig"].map(colors)

    for x, y, s, xlog, ylog in combine_plots:
        print("Plotting x={} y={}".format(x, y))
        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(111)
        ax.scatter(x=table[x], y=table[f"Control_{y}"], s=s, color="tab:green", label="Control")
        ax.scatter(x=table[x], y=table[f"Sample_{y}"], s=s, color="tab:red", label="Sample")
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if xlog: ax.set_xscale('log')
        if ylog: ax.set_yscale('log')
        plt.legend()
        fig.tight_layout()
        plt.savefig(f"{output_folder}/{x}_vs_{y}.png")
        plt.close(fig)

    for x, y, s, xlog, ylog in single_plots:
        print("Plotting x={} y={}".format(x, y))
        fig = plt.figure(figsize=[16, 8])
        ax = fig.add_subplot(111)
        ax.scatter(x=table[x], y=table[y], s=s, color=p_sig_colors)
        # plt.hlines(table[y].median(), xmin=table[x].min(), xmax=table[x].max(), color="tab:red", label="Median={:.2f}".format(table[y].median()))
        # plt.hlines(table[y].mean(), xmin=table[x].min(), xmax=table[x].max(), color="tab:blue", label="Mean={:.2f}".format(table[y].mean()))
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        if xlog: ax.set_xscale('log')
        if ylog: ax.set_yscale('log')
        fig.tight_layout()
        plt.savefig(f"{output_folder}/{x}_vs_{y}.png")
        plt.close(fig)

    for col, bins in hist_plots:
        print("Plotting col={}".format(col))
        fig = plt.figure(figsize=[12, 8])
        ax = fig.add_subplot(111)
        df = table[col].replace([np.inf, -np.inf], np.nan, inplace=False)
        df.plot.hist(bins=bins)
        plt.xlabel(f"{col}")
        plt.savefig(f"{output_folder}/{col}_hist.png")
        plt.close(fig)


    # Make the MA plot
    print("Plotting MA plot")
    fig = plt.figure(figsize=[12, 8])
    A = 0.5 * ( np.log2(table["Sample_Hits"]) + np.log2(table["Control_Hits"]) )
    M = table["Log2FC_Reads"]
    plt.scatter(x=A, y=M, s=10, color=p_sig_colors)
    plt.hlines(M.median(), xmin=A.min(), xmax=A.max(), color="tab:red", label="Median={:.2f}".format(M.median()))
    plt.hlines(M.mean(), xmin=A.min(), xmax=A.max(), color="tab:blue", label="Mean={:.2f}".format(M.mean()))
    plt.hlines([-1, 1], xmin=A.min(), xmax=A.max(), color="tab:blue", linestyle='dashed')
    plt.legend(loc="upper left")
    plt.xlabel("A = 1/2 * ( log2(sample) + log2(control) )")
    plt.ylabel("M = log2(sample) - log2(control)")
    fig.tight_layout()
    plt.savefig(f"{output_folder}/MA_Plot.png")
    plt.close(fig)


    # Make the volcano plot
    print("Plotting Volcano plot")
    fig = plt.figure(figsize=[12, 8])
    X = table["Log2FC_Reads"]
    Y = table["Log10P"]
    plt.scatter(x=X, y=Y, s=8, color=p_sig_colors)
    plt.hlines(-np.log10(alpha), xmin=X.min(), xmax=X.max(), color="tab:blue", linestyle='dashed')
    plt.vlines([-1, 1], ymin=Y.min(), ymax=Y.max(), color="tab:blue", linestyle='dashed')
    plt.xlabel("log2 fold Change")
    plt.ylabel("-log10(p-value)")
    fig.tight_layout()
    plt.ylim(max(0, Y.min()), Y.max())
    plt.savefig(f"{output_folder}/Volcano_Plot.png")
    plt.ylim(max(0, Y.min()), min(4, Y.max()))
    plt.savefig(f"{output_folder}/Volcano_Plot_Trim.png")
    plt.close(fig)

