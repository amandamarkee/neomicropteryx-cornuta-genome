## _STOP: Please note, I will be following methods outlined in [this repository](https://github.com/amandamarkee/actias-luna-genome/blob/main/annotation-notes-ACTIVE.md)_ ##
Goal: to replicate the BRAKER2 annotation pipeline in linked repository above for the Neomicropteryx cornuta genome, beginning after the RepeatModeler/RepeatMasker step.

## Feature Annotation – Building the Gene Model in BRAKER2 ##
![image](https://github.com/amandamarkee/neomicropteryx-cornuta-genome/assets/56971761/3860a858-ad95-4a2f-9629-76e4dd25d5b0)

Goal: to replicate the BRAKER3 pipeline for this genome annotation, beginning after the RepeatModeler/RepeatMasker step.

## **5/23/2023; Annotation – BRAKER 2 Setup and Input Files**

Resources:

- Running BRAKER with proteins of any evolutionary distance: https://github.com/Gaius-Augustus/BRAKER#braker-with-protein-data
- Running BRAKER with long-read IsoSeq data: https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md
- Output files: https://github.com/Gaius-Augustus/BRAKER#output-of-braker

All BRAKER2 annotation will be done in the following directory:
```
/blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2
```

Input files for BRAKER2:

1. soft-masked genome (fasta)
2. protein database (fasta)
3. long-read isoseq (fastq)

Directories for each input file:
```
soft-masked genome
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/Neomicropteryx_cornuta/Neomicropteryx_cornuta_softmasked.fasta

protein database
/blue/kawahara/yimingweng/LepidoPhylo_Project/OrthoDBv11/Arthropoda.fasta

long-read isoseq
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/IsoSeq/Neomicropteryx_cornuta/SRR14882579.fastq.gz
```

## Step 1 Feature Annotation – Running with protein sequences


## Step 2 Feature Annotation – Running with IsoSeq (long-read) data


## Step 3 Gene Model Evaluation – Running BUSCO on individual gene models


## Step 4 Combining Gene Models – Running TSEBRA transcript collapser 


## Step 5 Combined Model Evaluation – Running BUSCO on combined gene model
