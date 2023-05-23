# _Neomicropteryx cornuta_ Genome Project
Genome annotation workflow for the mandibulate moth _Neomicropteryx cornuta_ using protein and long-read IsoSeq evidence. 

## Background
- _Neomicropteryx cornuta_ is a species of mandibulate archaic moths (Lepidoptera: Micropterigidae) belonging to a lineage thought to have split from all other Lepidoptera more than 300 Ma.
- The genome for this species was recently published by Li et al. (2021) using long-read sequencing methods; the MAKER annotation pipeline was used with protein evidence alone.
- The focus of this repository will be to provide the udpated BRAKER2 workflow for annotating a high-quality reference genome for this species, with long-read PacBio IsoSeq evidence along with protein evidence.

## General Workflow
This section describes the methods I used for annotating the _Neomicropteryx cornuta_ genome with the BRAKER2 pipeline with PacBio High Fidelity (HiFi) reads. 
Note: I will be using input files generated by Xuankun Li, including the soft-masked genome fasta file from RepeatModeler/RepeatMasker, an arthropoda protein database fasta file from OrthoDB v.11, and the IsoSeq long-read RNA fastq file.

### Step 1 Feature Annotation – Running with protein sequences


### Step 2 Feature Annotation – Running with IsoSeq (long-read) data


### Step 3 Gene Model Evaluation – Running BUSCO on individual gene models


### Step 4 Combining Gene Models – Running TSEBRA transcript collapser 


### Step 5 Combined Gene Model Evaluation – Running BUSCO on combined gene model


## References

Xuankun Li and others, First Annotated Genome of a Mandibulate Moth, Neomicropteryx cornuta, Generated Using PacBio HiFi Sequencing, Genome Biology and Evolution, Volume 13, Issue 10, October 2021, evab229, https://doi.org/10.1093/gbe/evab229
