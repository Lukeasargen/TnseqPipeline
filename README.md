
## Overview
- Setup Workspace
- Stage 1 - Process the References
- Stage 2 - Process the Reads
- Stage 3 - Analyze


# Setup Workspace

First, complete the install instructions on [docs/install.md](docs/install.md).

Then, open a terminal and change directory to the location you want as your workspace. Run this command to clone this repository:

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
```
./scripts/reference1.sh -e experiment -f fasta -g gb -o out

# Bowtie 1
./scripts/reference1.sh -e demo -f 14028s_chromosome.fasta -g 14028s_genome.gb -o 14028c
```
<!-- # Bowtie 2
./scripts/reference2.sh -e demo -f 14028s_chromosome -g 14028s_genome -o 14028c -->

The indexes are output into the indexes folder of your experiment. These files will only be used by bowtie during the alignment step of processing reads.

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
```
./scripts/reads1.sh -e experiment_name -i index_name -a adapters -r reads

# Bowtie 1
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
| -v 3 | v mode alignments only check mismatches and ignores quality. There will be at most 3 mismatches. |
| -a | Report all alignments. |
| -m 1 | If there are multiple alignments, only report the best one. |
| --best --strata | Best reports in best-to-worst order. Strata reports reads with the least number mismatches. |

## 3. Map - Python

```
python3 scripts/readTAmap.py --experiment=$EXPERIMENT_NAME --index=$INDEX --map=${READ_NAME}
```
This script is documented in place. Nearly every line has a comment explaining what it does.

# Stage 3 - Analysis

For now, [analysis.py](scripts/analysis.py) only performs pairwise test for essential genes. It does support multiple biological replicates. The first operation is a normalization (total counts). Then, the replicates are used in 2 ways: 1) each replicate's raw counts are used to calculate conditional statistical significance and 2) the replicates are merged by averaging to calculate ratios and other pairwise metrics.

There are several arguments. Below are examples of each argument.

Display help message
```
python3 scripts/analysis.py -h
```

Single reads
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k
```

Biological replicates
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c25k c50k --samples s25k s50k
```

--output = Create a folder in analysis with this name, output here (default="default")
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k --output=group1v3
```

--plot = Create plots of some metrics
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k --plot
```

--gc = Calculates the GC content of each gene 
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k --gc
```

--min_count = only use genes that have more counts than this (default=1)
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k --min_count=1
```

--min_inserts = only use genes that have more insertion sites than this (default=2)
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k --min_inserts=2
```

--smoothing = smooths the ratio for small counts (default=1)
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c50k --samples s50k --smoothing=1
```

You can combine the arguments

Control replicates, GC content, minimum of 10 hits per gene, plot
```
python3 scripts/analysis.py --experiment demo --index 14028c --controls c25k c50k --samples s50k --gc --min_count=10 --plot
```
