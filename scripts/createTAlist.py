import os
import re
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    # Init and setup
    parser.add_argument('--experiment', type=str)
    parser.add_argument('--fasta', type=str)
    parser.add_argument('--genbank', type=str)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()
    return args


def make_TAlist(args):
    print(" * Creating a TAlist for {}.fasta and {}.gb".format(args.fasta, args.genbank))

    # We know these files exist if the bowtie step worked
    fasta_filename = "data/{}/references/{}.fasta".format(args.experiment, args.fasta)
    gb_filename = "data/{}/references/{}.gb".format(args.experiment, args.genbank)
    
    # This is where the TAlist is output
    output_filename = "data/{}/references/{}_TAlist.tsv".format(args.experiment, args.output)
    print("Output Location:", output_filename)
    if os.path.exists(output_filename):
        print("* TAlist already exists")
        # TODO : what is the behavior here?

    # This line reads every single line from the fastas file and removes the endline character
    unedited = [line.rstrip('\n') for line in open(fasta_filename, 'r')]  # encoding="utf-8"

    # Combine all the lines into one long string
    # This is indexed from 1: because the first line is not part of the sequence
    fullseq = "".join(unedited[1:])
    print("Full Sequence Length =", len(fullseq))

    # However, we do need the name of the genenome
    # We get that from the first line, it's the first string before space but without the first character
    # It's wrapped as string just in case
    geneome_name = str(unedited[0].split()[0][1:])
    print("Geneome name:", geneome_name)

    # Now, we find every single TA site from the fullseq string
    # startswith finds the index before the T so that's why idx+1 is needed
    # Add 1 if you want to know what the A location is
    ta_sites = [idx+1 for idx in range(len(fullseq)) if fullseq.startswith("TA", idx)]
    print("Total TA Sites =", len(ta_sites))
    # This is the first TA
    # The -1 is because arrays are 0 indexed
    print("First TA is at {}: {}".format(ta_sites[0], fullseq[ta_sites[0]-1:ta_sites[0]+1]))
    # Now we know were every TA site is


    # Open genbank and removes the endline character
    unedited = [line.rstrip('\n') for line in open(gb_filename, 'r')]  # encoding="utf-8"

    # Output starts with headers
    output_array = ["Genome Name\tLoci\tGene ID\tLocus Tag\tStart\tEnd\tDirection\tTA Site"]

    # Here we get the gene feature idexes for the next for loop
    features_idxs = [idx for idx in range(len(unedited)) if unedited[idx].startswith("     gene")]
    print("Total Genes:", len(features_idxs))

    ta_idx = 0  # step through all the ta sites
    prev_gene_end = 0  # used to marked the first intergenic region could start at 0

    # Go through each gene in genbank, get the features, make the TA list
    # for readablity, anythin that starts with "gene_" will be in the TA list
    for i in range(len(features_idxs)):

        # 'i' also is the intergenic location minus 1
        # TODO : what is the template for naming loci?
        gene_loci = "G_{}".format(i+1)  # Starts at 1
        ig_loci = "IG_{}".format(i+1)

        # 'j' is the line number in the genbank
        # 'j' will iterate as the feature is parsed
        j = features_idxs[i]

        # First, we parse the start, end, and direction of the gene
        loc_line = unedited[j].split()[1]
        loci_nums = re.findall(r'\d+', loc_line)
        gene_start, gene_end = int(min(loci_nums)), int(max(loci_nums))
        gene_direction = "-" if "complement" in loc_line else "+"

        # Find lines with identation to get feature data
        j += 1  # go to next line after the gene loci
        identation = 21  # feature lines have this many spaces
        # unedited[j] is a string of the j'th line
        # this while loop checks if the identation is the same, a new feature will be false so the loop breaks
        gene_id = ""  # Default case
        gene_lucas_tag = ""  # Default case
        while unedited[j].lstrip() == unedited[j][identation:]:
            line = unedited[j].strip()  # remove spaces on both sides
            if line[1:6] == "gene=":
                gene_id = line[7:-1]
            if line[1:11] == "locus_tag=":
                gene_lucas_tag = line[12:-1]
            j += 1  # go to the next lin

        # Now we have all the information about the gene
        # Using the loci, we can get the TA sites for the intergenic region and the gene

        # Intergenic region first
        # from prev_gene_end to gene_start
        # use less than bc the gene_start is part of the gene and not intergenic
        while ta_sites[ta_idx] < gene_start:
            out = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(geneome_name, ig_loci, "", "", prev_gene_end, gene_start, "", ta_sites[ta_idx])            
            # print(out)
            output_array.append(out)
            ta_idx += 1  # next site

        # Now the gene
        # from gene_start to gene_end
        # use less thna or equal bc the gene_end is part of the gene
        while ta_sites[ta_idx] <= gene_end:
            out = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(geneome_name, gene_loci, gene_id, gene_lucas_tag, gene_start, gene_end, gene_direction, ta_sites[ta_idx])
            # print(out)
            output_array.append(out)
            ta_idx += 1  # next site

        # if i > 4:  # leave early for debug
        #     break

        prev_gene_end = gene_end  # start for the next intergenic region
        # END OF THE GENBANK LOOP

    # Final intergenic region goes to the end to the sequence
    # prev_gene_end = len(fullseq)
    ig_loci = "IG_{}".format(i+2)  # the last was +1, so now +2
    while ta_sites[ta_idx] < len(fullseq):
        out = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(geneome_name, ig_loci, "", "", prev_gene_end, gene_start, "", ta_sites[ta_idx])            
        # print(out)
        output_array.append(out)
        ta_idx += 1  # next site
        if ta_idx >= len(ta_sites):
            break

    # Save everything
    with open(output_filename, "w") as filehandle:
        for line in output_array:
            filehandle.writelines("%s\n" % line)

    print(" * Saved TAlist to {}".format(output_filename))


def main():
    print("\n * Running createTAlist.py")
    args = get_args()

    make_TAlist(args)


if __name__ == "__main__":
    main()
