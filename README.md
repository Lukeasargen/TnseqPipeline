
## Overview
1. Setup Workspace
2. Stage 1 - Process the Genome
3. Stage 2 - Process the Reads
4. Repeat 4 for all reads
5. Stage 3 - Analyze
6. Stage 4 - Plot


# Setup Workspace

The project folder will look like this:
```
  TnseqPipeline/
  ├── data
  │   ├── demo
  │   |   ├── adapters - used for trimming the reads before alignment
  │   |   ├── analysis - TA sequence hits are output here
  │   |   ├── indexes - mapping output by bowtie-build
  │   |   ├── maps - TA maps
  │   |   ├── reads - tnseq reads and temporary files
  │   |   └── references - genomes and genelists
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

# Stage 1 - Process the Genome

Inputs:
- Sequenced Genome - fasta, fastq, or fa file

Outputs:
- Indexes - used by Bowtie for alignment
- TA map - master file for all following alignment and counting steps

Stage 1 only needs to be done once for a genome. The

## Create the bowtie index
```
bowtie-build data/demo/reference/14028s_chromosome.fasta data/demo/indexes/14028c
bowtie2-build data/demo/reference/14028s_chromosome.fasta data/demo/indexes/14028c
```

## Create the Blank map
```

```

# Stage 2 - Process the Reads

Inputs:
- Reads - fastq files
- Adapters - fa files
- Indexes - from stage 1
- TA map - from stage 1

Outputs:
- Trimmed Reads
- Cropped Reads
- TA maps

If you're feeling lucky, these scripts automate each step. Run these with the names of your genome, reads, and adapters.

```
./scripts/reads1.sh experiment_name genome_name read_name adapter_name

# Bowtie1
./scripts/reads1.sh demo 14208c pt1_S1_L001_R1_001 PolyC_Adapter
./scripts/reads1.sh demo 14208c cf1_S3_L001_R1_001 PolyC_Adapter

# Bowtie2
./scripts/reads2.sh demo 14208c pt1_S1_L001_R1_001 PolyC_Adapter
./scripts/reads2.sh demo 14208c cf1_S3_L001_R1_001 PolyC_Adapter
```

The read scripts are most useful if you want to use my default settings for trimming and aligning. If you need to change these parameters, this is the breakdown of each process called in the read scripts.

## 1. Trim

## 2. Crop

## 2. Align


# Stage 3 - Analysis





