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
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/IsoSeq/Neomicropteryx_cornuta/Neomicropteryx_cornuta_softmasked.fasta

protein database
/blue/kawahara/yimingweng/LepidoPhylo_Project/OrthoDBv11/Arthropoda.fasta

long-read isoseq
/blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/IsoSeq/Neomicropteryx_cornuta/SRR14882579.fastq.gz
```

## Step 1 Feature Annotation – Running with protein sequences
## (a) Run [ProtHint](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) to create protein gff file

Here, I will run ProtHint on the following file containing sequence data for arthropod protein sequences found within the OrthoDBv11 database:
```
sbatch -J Nc_ProtHint prothint.sh /blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/IsoSeq/Neomicropteryx_cornuta/Neomicropteryx_cornuta_softmasked.fasta /blue/kawahara/yimingweng/LepidoPhylo_Project/OrthoDBv11/Arthropoda.fasta
```

Script for running ProtHint using a SLURM submission:
```
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

dates;hostname;pwd

genome=${1}
protein=${2}

module load prothint/2.6.0
module load genemark_es/4.69

prothint.py --threads ${SLURM_CPUS_ON_NODE:-1} ${genome} ${protein}
```

## (b) Run BRAKER2 with protein evidence

Now that the ProtHint .gff file is completed, we can use this as protein evidence for runnning BRAKER2. Use the following command to execute the code below it.

```
sbatch -J Nc_prot_braker2 Nc_braker2_protein.sh /blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/IsoSeq/Neomicropteryx_cornuta/Neomicropteryx_cornuta_softmasked.fasta /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/prothint/prothint_augustus.gff Neomicropteryx_cornuta
```
```
#!/bin/bash
#SBATCH --job-name=%j_Nc_braker2_prot
#SBATCH --output=%j_Nc_braker2_prot.log
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

genome=${1}
protein_gff=${2}
species=${3}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/amanda.markee/neomicropteryx_annotation/Augustus/config \
--genome=${genome} --species ${species} --hints=${protein_gff} --softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
```

## Step 2 Feature Annotation – Running with IsoSeq (long-read) data


## Step 3 Gene Model Evaluation – Running BUSCO on individual gene models


## Step 4 Combining Gene Models – Running TSEBRA transcript collapser 


## Step 5 Combined Model Evaluation – Running BUSCO on combined gene model
