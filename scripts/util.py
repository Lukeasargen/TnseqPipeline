import numpy as np
from scipy.stats import trim_mean


def time_to_string(t):
    if t > 3600: return "{:.2f} hours".format(t/3600)
    if t > 60: return "{:.2f} minutes".format(t/60)
    else: return "{:.2f} seconds".format(t)


def column_stats(table, columns):
    if type(columns)==str: columns=[columns]
    for n in columns:
        print( "\nColumn Stats : {}".format(n) )
        print( "Sum={:.4f}.".format(table[n].sum()) )
        print( "Min={:.4f}. Max={:.4f}. Median={:.4f}.".format(table[n].min(), table[n].max(), table[n].median()) )
        print( "Mean={:.3f}. Non-Zero Mean={:.3f}".format(table[n].mean(), table[table[n]!=0][n].mean()) )
    print()


def read_fasta(filename, ret_name=False):
    # This line reads every single line from the fastas file and removes the endline character
    unedited = [line.rstrip('\n') for line in open(filename, 'r')]  # encoding="utf-8"
    # Combine all the lines into one long string
    # This is indexed from 1: because the first line is not part of the sequence
    fullseq = "".join(unedited[1:])
    if ret_name:
        # We get name from the first line, it's the first string before space but without the first character
        # It's wrapped as string just in case
        name = unedited[0].split()[0][1:]
        return fullseq, name

    return fullseq


def get_read_columns(table):
    # Only read columns have one of these substrings
    subs = ["_both", "_forward", "_reverse"]
    # ret = list(filter(lambda x: any(s in x for s in subs), table.columns))
    ret = [c for c in table.columns if any(s in c for s in subs)]
    return ret


def exclude_sites_tamap(tamap, exclude_first=0, exclude_last=0):
    """ Remove a a percentage of the reads from the start and end of a region.
        exclude_first and exclude_last are percentages as decimals
        Idea from : "Tn-seq; high-throughput parallel sequencing for fitness and genetic interaction studies in microorganisms"
    """
    tamap = tamap.copy()
    map_names = get_read_columns(tamap)
    # We want to keep the data inside the exclude_first and
    # exclude_last region and remove the rest.  The simplest way
    # I could come up with was to set everything outside the region to zero.

    # This section calculates the start and end of the regions we want to keep
    # 1+End-Start is the gene length and exlcued
    if exclude_first>0 and exclude_first<1.0:
        new_start = tamap["Start"] + (1+tamap["End"]-tamap["Start"])*exclude_first
    else:
        new_start = tamap["Start"]
    if exclude_first>0 and exclude_first<1.0:
        new_end = tamap["End"] - (1+tamap["End"]-tamap["Start"])*exclude_last
    else:
        new_end = tamap["End"]
    # This is a boolean column with all the rows which are set to 0
    # the vertical bar | is the python bitwise "or" operator
    remove = (tamap["TA_Site"] < new_start) | (tamap["TA_Site"] > new_end)
    # Set all the read colums for the excluded regions to 0
    tamap.loc[remove, map_names] = 0
    return tamap


def tamap_to_genehits(tamap, fasta_filename=None, pooling="sum"):
    """ Covert a TAmap table into a GeneHits table """
    tamap = tamap.copy()
    # Remove the intergenic regions
    genemap = tamap[tamap['Gene_ID'].notna()]
    map_names = get_read_columns(genemap)
    # Get other gene data into a genehits df
    grouped = genemap.groupby("Gene_ID", as_index=False)
    genehits = grouped.agg({"Start": "first", "End": 'first', "Direction": "first"})
    genehits["TA_Count"] = grouped["TA_Site"].count()["TA_Site"]
    genehits["Gene_Length"] = 1+genehits["End"] - genehits["Start"]
    # Compute GC content with the genome fasta file
    genehits["GC"] = None  # Blank
    if fasta_filename:
        full_seq = read_fasta(fasta_filename)
        for i, row in genehits.iterrows():
            seq = full_seq[row["Start"]:row["End"]+1].lower()
            gc = seq.count('g') + seq.count('c')
            genehits.loc[i, "GC"] = (gc)/len(seq)
    # Add the forward and reverse reads
    if pooling=="sum":
          # Get the total number of hits per gene TA sites
        gene_hits_sum = grouped.sum()
    if pooling=="average":
        """ Used in "Genome-Wide Fitness and Genetic Interactions Determined by Tn-seq, a High-Throughput Massively Parallel Sequencing Method for Microorganisms"
            https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696536/
            Mentioned under "Anticipated Results"
        """
        # Get the average of all TA sites per gene
        gene_hits_sum = grouped.mean()
        
    for name in map_names:
        genehits[name] = gene_hits_sum[name]
    return genehits


""" Normalization Operations
    Applied to TAmap, but they should work on Genehits tables.
"""

def apply_multiply_factors(table, factors):
    """ table is a has the count data
        factors is a dictionary with the correct column names
    """
    for k,v in factors.items():
        table[k] = v*table[k]
    return table


def total_count_norm(genehits, columns=None, debug=False):
    """ Normalize genehits using the total reads. """
    temp = genehits.copy()
    if columns==None:
        columns = get_read_columns(temp)
    multiply_factor = {}
    for name in columns:
        multiply_factor.update({name: temp[name].sum()})
    max_total = max(multiply_factor.values())
    multiply_factor = {k:max_total/v for k,v in multiply_factor.items()}
    temp = apply_multiply_factors(temp, multiply_factor)
    if debug:
        print("\ntotal_count_norm")
        print("max_total :", max_total)
        print("multiply_factor :", multiply_factor)
    return temp


def quantile_norm(genehits, q=0.5, columns=None, debug=False):
    """ Normalize genehits using the q'th quantile of the non-zero values."""
    temp = genehits.copy()
    if columns==None:
        columns = get_read_columns(temp)
    multiply_factor = {}
    for name in columns:
        multiply_factor.update({name: temp[temp[name]!=0][name].quantile(q=q)})
    max_total = max(multiply_factor.values())
    multiply_factor = {k:max_total/v for k,v in multiply_factor.items()}
    temp = apply_multiply_factors(temp, multiply_factor)
    if debug:
        print("\nquantile_norm q={}".format(q))
        print("max_total :", max_total)
        print("multiply_factor :", multiply_factor)
    return temp


def ttr_norm(genehits, trim=0.05, columns=None, debug=False):
    """ Normalize genehits using the .
        Assumes most genes are not differentially expressed.
        Does not consider gene length or library size.
    """
    temp = genehits.copy()
    if columns==None:
        columns = get_read_columns(temp)
    multiply_factor = {}
    for name in columns:
        a = np.mean(temp[name][temp[name]>0], axis=0)
        b = trim_mean(temp[name][temp[name]>0], trim)
        multiply_factor.update({name: a*b})
    max_total = max(multiply_factor.values())
    multiply_factor = {k:max_total/v for k,v in multiply_factor.items()}
    temp = apply_multiply_factors(temp, multiply_factor)
    if debug:
        print("\nttr_norm trim={}".format(trim))
        print("max_total :", max_total)
        print("multiply_factor :", multiply_factor)
    return temp


def nzmean_norm(genehits, columns=None, debug=False):
    """ Normalize genehits using the .
        Assumes most genes are not differentially expressed.
        Does not consider gene length or library size.
    """
    temp = genehits.copy()
    if columns==None:
        columns = get_read_columns(temp)
    multiply_factor = {}
    for name in columns:
        nzmean = genehits[name].replace(0, np.NaN).mean()
        multiply_factor.update({name: nzmean})
    max_total = max(multiply_factor.values())
    multiply_factor = {k:max_total/v for k,v in multiply_factor.items()}
    temp = apply_multiply_factors(temp, multiply_factor)
    if debug:
        print("\nnzmean_norm trim")
        print("max_total :", max_total)
        print("multiply_factor :", multiply_factor)
    return temp


def gene_length_norm(genehits, columns=None, debug=False):
    """ Normalize genehits using the a length column. """
    temp = genehits.copy()
    if columns==None:
        columns = get_read_columns(temp)
    # Normalize by gene length
    lengths = 1+temp["End"]-temp["Start"]
    norm_length = lengths.mean()
    for name in columns:
        temp[name] = (lengths/norm_length)*temp[name]
    if debug:
        print("\ngene_length_norm")
        print("columns :", columns)
        print("norm_length :", norm_length)
    return temp


""" Stats Functions """

def bh_procedure(pvalues, alpha=0.05):
    """ Benjamini-Hochberg Procedure.
        Formula from this video:
            "False Discovery Rates, FDR, clearly explained"
            by StatQuest with Josh Starmer
            https://www.youtube.com/watch?v=K8LQSvtjcEo
    """
    pvalues = np.array(pvalues)
    m = len(pvalues)

    # First step is to sort the pvalues
    sort_idx = np.argsort(pvalues)
    qvalues = np.take(pvalues, sort_idx)

    new_alpha = alpha
    # Then step through the list from largest to smallest
    # Largest qvalue is the same, so use m-1
    # The rest is a simple formula that takes the minimum
    for i in reversed(range(m-1)):
        new_q = qvalues[i] * m/(i+1)
        qvalues[i] = min(new_q, qvalues[i+1])
        if qvalues[i] < alpha and qvalues[i]!=0 and new_alpha==alpha:
            new_alpha = pvalues[sort_idx[i]]

    # Double argsort trick. https://www.berkayantmen.com/rank.html
    out = np.empty_like(qvalues)
    out[sort_idx] = qvalues
    # unsort_idx = np.argsort(sort_idx)
    # qvalues = np.take(qvalues, unsort_idx)
    return out, new_alpha

""" Misc """

def calc_sample_fitness(table, control, sample, expansion):
    """
    Idea from : "Tn-seq; high-throughput parallel sequencing for fitness and genetic interaction studies in microorganisms"
    As far as I can tell, expansion factor (args.expansion) is an arbitrary value.
    I explored it's effects on the fitness formula in desmos.
    Formula = fitness = ln(ef*(sf)/(cf)) / ln(ef*(1-sf)/(1-cf))
    where sf=sample frequency(t2), cf=control frequency(t1), and ef=expansion factor.
    Changing expansion factor appears to change the steepness of the slope
    around the neutral value. Larger expansion, smaller slope.
    """
    # Re adjust the table to with the same behavior as Opijnen "correction factors"
    # "Total Counts" effectively makes the columns have the same sum, which is identical to the Opijnen procedure
    df = total_count_norm(table, columns=[control, sample])
    # I'm pretty sure frequency is just the site counts over total conts
    df["contrl_freq"] = (df[control]) / df[control].sum()
    df["sample_freq"] = (df[sample]) / df[sample].sum()
    # This is the formula from Opijnen (Nature 2009) Online Methods
    num = expansion*(df["sample_freq"])/(df["contrl_freq"])
    dem = expansion*(1-df["sample_freq"])/(1-df["contrl_freq"])
    fitness = np.log(num, out=np.zeros_like(num), where=(num!=0)) / np.log(dem, out=np.zeros_like(dem), where=(dem!=0))

    # df["Sample_Fitness"] = fitness
    # temp = df.sort_values(by="Sample_Fitness", ascending=True)
    # for idx in temp.index[:2]:
    #     print("idx :", idx)
    #     print(df.iloc[idx])

    return fitness


def calc_survival_index(table, control, sample):
    # Idea from : "Genetic basis of persister tolerance to aminoglycosides (2015) Shan...Lewis"
    # dval_ctl = (genehits["Control_Hits"]*genehits["Gene_Length"].sum()) / (genehits["Gene_Length"]*genehits["Control_Hits"].sum())
    # dval_exp = (genehits["Sample_Hits"]*genehits["Gene_Length"].sum()) / (genehits["Gene_Length"]*genehits["Sample_Hits"].sum())
    # si = dval_exp/dval_ctl
    # This is equivalent to above
    si = (table[sample]*table[control].sum())/(table[control]*table[sample].sum())
    return si

