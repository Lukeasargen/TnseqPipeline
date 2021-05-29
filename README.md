
## Overview
- Setup Workspace
- Stage 1 - Process the References
- Stage 2 - Process the Reads
- Stage 3 - Analyze


# Setup Workspace

First, complete the install instructions on [docs/install.md](docs/install.md).

Everything should now be installed in your terminal you are using to run the analysis. Then, open a terminal and change directory to the location you want as your workspace. 

Run this command to clone this repository:

```
git clone https://github.com/Lukeasargen/TnseqPipeline.git
```

Install additional python requirements:
```
python3 -m pip install numpy==1.17.4 pandas==1.2.4 matplotlib==3.1.2 scipy==1.3.3
```

## All Done Setup

The project folder will look like this:
```
  TnseqPipeline/
  ├── data
  │   ├── blank - empty folders with the correct structure
  │   ├── demo
  │   |   ├── adapters - removed before alignment, fasta or fa files
  │   |   ├── analysis - genehits tables are output here
  │   |   ├── indexes - mapping output by bowtie-build
  │   |   ├── maps - TA maps
  │   |   ├── reads - reads straight from the machine, likely fastq files
  │   |   |   └── processed - trimming and alignment outputs go here
  │   |   └── references - genomes and TAlist
  │   └── Your experiments ...
  ├── docs
  │   ├── readme
  │   ├── images and screenshots
  │   └── instructions
  ├── scripts
  │   ├── R scripts
  │   ├── Python scripts
  │   └── Shell scripts
  └── tools
      └── Trimmomatic-0.36
```

# Starting a New Experiment

1. Make a copy of the `blank` folder found in `data` and rename. Choose a short name with no spaces.

```
cp -a data/blank /data/your_experiment_name
```
2. You need a reference genome and genebank. Put these in your new experiment folder in the `references` folder. The genome should be a *.fa or *.fasta and the genebank is a *.gb.

3. Adapters are not required. If you have adapters, put them in the `adapters` folder. These are *.fa or *.fasta files. These are used by Trimmomatic to trim these bases from the raw reads. _If you what to use a different trimming tool, you can put the output files in `reads/processed` and `reads1.sh` will try and bowtie align them. However, to avoid issues, you can just do the alignment manually following the steps for Stage 2 detailed below._

4. Finally, put your sequenced reads in the `reads` folder. These fastq files with extensions *.fq or *.fastq.

# Stage 1 - Process the Reference

Inputs:
- Sequenced Reference - fasta file
- Genbank file - gb file

Outputs:
- Indexes - used by Bowtie for alignment
- TAmap - csv marking all TA sites with indicators for genome and gene. This file has TA sites as rows and read data down the columns. 

## *Stage 1 only needs to be done once for a reference.*

## Create the Bowtie index and TAlist

<!-- Bowtie 1 and Bowtie 2 output different files. Be consistent on which version you use in stage 1 and stage 2. -->
This is the template for creating the correct reference files:
```
./scripts/reference1.sh -e experiment -f fasta -g gb -o index_name
```

Here is the command filled out to run the demo:
```
./scripts/reference1.sh -e demo -f 14028s_chromosome.fasta -g 14028s_genome.gb -o 14028c
```
<!-- # Bowtie 2
./scripts/reference2.sh -e demo -f 14028s_chromosome -g 14028s_genome -o 14028c -->

The indexes are output into the indexes folder of your experiment. These files will only be used by bowtie during the alignment step of processing reads. The name of these indexes, specified after -o, is used as the index name for the following steps.

## Notice: The TAlist indicates the location using 1 indexing. This means the first base is labeled 1. For instance, if the gene start is 10, there will be 9 bases before it.


# Stage 2 - Process the Reads

Inputs:
- Reads - fastq files
- Adapters - fa files
- Indexes - from stage 1
- TA map - from stage 1

Outputs:
- Trimmed Reads - no adpaters and trimmed
- Aligned Reads - mapped reads to the genome
- TA maps - total alignments for every TA site

<!-- Use the same version of Bowtie as stage 1. -->
This is the template for processing reads:
```
./scripts/reads1.sh -e experiment_name -i index_name -a adapters -r reads
```

This is a command filled out for the demo. Notice the reads are listed with the filetype extension and separated by a comma.
```
./scripts/reads1.sh -e demo -i 14028c -a PolyC_Adapter.fa -r c25k.fastq,c50k.fastq,s25k.fastq,s50k.fastq
```

The `reads#.sh` scripts are most useful if you want to use default settings for trimming and aligning. If you need to change these parameters, this is the breakdown of each process called in the read scripts.

## Notice: If you want to run bowtie with different options, you have to change the options in reads1.sh and remove the files from the `reads/processed` folder.


## 1. Trim and Crop - Trimmomatic

```
java -jar tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 data/$EXPERIMENT_NAME/reads/${READ_NAME}.fastq data/$EXPERIMENT_NAME/reads/processed/${READ_NAME}_trimmed ILLUMINACLIP:data/$EXPERIMENT_NAME/adapters/${ADPATER_NAME}:2:30:10 LEADING:20 TRAILING:20 MINLEN:48 CROP:20
```

The best way to understand these arguments is to read the manual: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

TRAILING and MINLEN are dependent on the experiment. From the manual:
- LEADING: Cut bases off the start of a read, if below a threshold quality
- TRAILING: Cut bases off the end of a read, if below a threshold quality
- MINLEN: Drop the read if it is below a specified length
- CROP: Cut the read to a specified length by removing bases from the end

These defaults were chosen to work with 50bp reads from an Illumina MiSeq System. The combination of the LEADING, TRAILING, and MINLEN removes all reads that start or end with more than 2 low quality bases. Even with this strict threshold, only around 5% of reads are dropped. Then the reads are cropped to 20 bp because bowtie works well with short reads (usually 5-18% of alignments fail).


## 2. Align - Bowtie

```
bowtie -t -v 3 -a -m 1 --best --strata data/$EXPERIMENT_NAME/indexes/$INDEX data/$EXPERIMENT_NAME/reads/processed/${READ_NAME}_trimmed data/$EXPERIMENT_NAME/reads/processed/${READ_NAME}_trimmed_mapped
```

The bowtie manual is a longer read http://bowtie-bio.sourceforge.net/manual.shtml

Here are the default arguments used:

|||
|-|-|
| -t | Prints the duration of the alignment |
| -v 3 | v mode alignments only check mismatches and ignores quality. There will be at most 3 mismatches |
| -a | finds all mapping alignments (doesn’t stop at first hit) |
| -m 1 | if there are multiple alignments, only report the best one |
| --best | reports in best-to-worst order |
| --strata | reports only the best -m matches |

## 3. Map - Python

```
python3 scripts/readTAmap.py --experiment=$EXPERIMENT_NAME --index=$INDEX --map=${READ_NAME}
```
This script is documented in place. Nearly every line has a comment explaining what it does.

# Stage 3 - Analysis

For now, [analysis.py](scripts/analysis.py) only performs pairwise test for essential genes. It does support multiple biological replicates. The first operation is a normalization (total counts). Then, the replicates are used in 2 ways: 1) each replicate's raw counts are used to calculate conditional statistical significance and 2) the replicates are merged by averaging to calculate ratios and other pairwise metrics.



Here are the essential arguments for getting started:

Display help message
```
python3 scripts/analysis.py -h
```

Single reads. This is the minimum required inputs to run a pairwise comparison.
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k
```

Biological replicates, separated by spaces, no filetype
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c25k c50k --samples s25k s50k
```

--output = Create a folder in analysis with this name. This example outputs to a folder called "group1v3".
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k --output=group1v3
```

There are several other arguments. This is the help message.

```
usage: analysis.py [-h] --experiment EXPERIMENT --index INDEX --controls
                   CONTROLS [CONTROLS ...] --samples SAMPLES [SAMPLES ...]
                   [--output OUTPUT] [--debug] [--plot] [--gc]
                   [--pooling {sum,average}] [--min_count MIN_COUNT]
                   [--min_inserts MIN_INSERTS] [--smoothing SMOOTHING]
                   [--alpha ALPHA] [--insert_weighting]

Pairwise Comparison (Supports Replicates).

optional arguments:
  -h, --help            show this help message and exit
  --experiment EXPERIMENT
                        Experiment folder name.
  --index INDEX         Index name.
  --controls CONTROLS [CONTROLS ...]
                        List read names without the filetype and separated by
                        a space.
  --samples SAMPLES [SAMPLES ...]
                        List read names without the filetype and separated by
                        a space.
  --output OUTPUT       Output name. Analysis outputs to a folder with this
                        name. default=default.
  --debug               Boolean flag that outputs all my debugging messages.
                        default=False.
  --plot                Boolean flag that automatically makes a few plots of
                        the data. default=False.
  --gc                  Boolean flag that calculates the GC content of each
                        gene. Not used in any test, but it makes some plots
                        and gets saved in the output. default=False.
  --pooling {sum,average}
                        String argument. Sum or average the hits PER GENE to
                        get a merged value for the expression at the gene
                        level. default=sum.
  --min_count MIN_COUNT
                        Integer argument. Threshold for lowest number of hits
                        PER GENE after pooling. Removes genes with low pooled
                        hits. These genes are not tested for significance and
                        saved in a separate output table. default=1.
  --min_inserts MIN_INSERTS
                        Integer argument. Threshold for lowest number of
                        insertion sites with hits BY GENE. Removes genes with
                        low TA sites or low TA diversity (diversity=unique hit
                        sites/total gene sites). These genes are not tested
                        for significance and saved in a separate output table.
                        default=2.
  --smoothing SMOOTHING
                        Float argument. Smooth the ratio for small counts.
                        ratio=(sample+smoothing)/(control+smoothing).
                        default=1.
  --alpha ALPHA         Float argument. Significance level. default=0.05.
  --insert_weighting    Boolean flag that scales PER GENE based on unique
                        inserts. It increases the Formula is new_hits=old_hits
                        *(unique_inserts/average_unique_inserts).
                        default=False.
```
