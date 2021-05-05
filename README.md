
## Overview
- Setup Workspace
- Stage 1 - Process the References
- Stage 2 - Process the Reads
- Stage 3 - Analyze
- Stage 4 - Plot


# Setup Workspace

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
      ├── Bowtie
      └── Trimmomatic-0.36
```

# Stage 1 - Process the Reference

Inputs:
- Sequenced Reference - fasta file
- Genbank file - gb file

Outputs:
- Indexes - used by Bowtie for alignment
- TAlist - table seperated file marking all TA sites with indicators for genome and gene

Stage 1 only needs to be done once for a reference. If you add another reference, the same scripts have an optional argument to append the existing TAlists into one merged TAlist.

## Create the Bowtie index and TAlist

Bowtie 1 and Bowtie 2 output different files. Be consistent on which version you use in stage 1 and stage 2.
```
./scripts/reference1.sh -e experiment -f fasta -g gb -o out -m

# Bowtie 1
./scripts/reference1.sh -e demo -f 14028s_chromosome -g 14028s_genome -o 14028c

# Bowtie 2
./scripts/reference2.sh -e demo -f 14028s_chromosome -g 14028s_genome -o 14028c
```
The indexes are output into the indexes folder of your experiment. These files will only be used by bowtie during the alignment step of processing reads.

## Notice: The TAlist indicates the location using 1 indexing. This means the first base is labeled 1. For instance, if the gene start is 10, there will be 9 bases before it.

# Stage 2 - Process the Reads

Inputs:
- Reads - fastq files
- Adapters - fa files
- Indexes - from stage 1
- TA map - from stage 1

Outputs:
- Trimmed Reads - removed adpaters
- Cropped Reads - removed excessively long reads
- Aligned Reads - mapped reads to the genome
- TA maps - total alignments for every TA site

Use the save version of Bowtie as stage 1.
```
./scripts/reads1.sh -e experiment_name -i index_name -a adapters -r reads

# Bowtie 1
./scripts/reads1.sh -e demo -i 14028c -a PolyC_Adapter.fa -r pt1_S1_L001_R1_001.fastq,cf1_S3_L001_R1_001.fastq

# Bowtie 2

```

The `reads#.sh` scripts are most useful if you want to use default settings for trimming and aligning. If you need to change these parameters, this is the breakdown of each process called in the read scripts.

## 1. Trim

## 2. Crop

## 3. Align

## 4. Map


# Stage 3 - Analysis


