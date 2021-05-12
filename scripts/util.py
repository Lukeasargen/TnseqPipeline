


def normalize_genehits(genehits, total=True, length=None, length_first=True):
    """ Normalize a genehits dataframe.
        total: bool, True norms by total number of reads, False 
            skips this step. default is True.
        length: str, None does not normalize by length, "Gene_Length" 
            and "TA_Count" normalize by their respective column, default
            is None
        length_first: bool, True performs length first, False is second,
            default is True

        returns: normalized dataframe
    """
    # Note : all normalization is done by multipling by greater than 1
    # This makes sure that you can threshold by minimum count of 1 as a lower bound

    map_names = genehits.columns[6:]

    # Normalize by gene length
    if length and length_first:
        max_length = max(genehits[length])
        for name in map_names:
            genehits[name] = (max_length/genehits[length])*genehits[name]

    # Normalize for total reads
    if total:
        totals = {}
        for name in map_names:
            totals.update({name: genehits[name].sum()})
        max_total = max(totals.values())
        for k,v in totals.items():
            genehits[k] = (max_total/v)*genehits[k]

    # Normalize by gene length
    if length and not length_first:
        max_length = max(genehits[length])
        for name in map_names:
            genehits[name] = (max_length/genehits[length])*genehits[name]

    # Normalization is done
    return genehits




