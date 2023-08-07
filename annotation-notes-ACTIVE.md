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

soft-masked genome w/trucated sample name
/blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/Ncor.softmasked_tk.fasta

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

Now that the ProtHint .gff file is completed, we can use this as protein evidence for runnning BRAKER2. Use the following command to execute the code below it. Note: Here I use a truncated
version of the soft-masked genome, with white spaces removed. This file contains shortened IDs to match the ProtHint augustus gff output. IDs must match identically or braker will not run properly.

```
sbatch -J Nc_prot_braker2 Nc_braker2_protein.sh /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/Ncor.softmasked_tk.fasta /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/prothint/prothint_augustus.gff Neomicropteryx_cornuta
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


These commands were used for removing all whitespaces in the softmasked genome. I first copy over YiMing's softmasked genome to my directory using this command:
```
cp /blue/kawahara/yimingweng/LepidoPhylo_Project/annotations/IsoSeq/Neomicropteryx_cornuta/Neomicropteryx_cornuta_softmasked.fasta /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot
```

Then use the following awk command to remove whitespaces:
```
sed 's, ,_,g'  Neomicropteryx_cornuta_softmasked.fasta >nospace_Neomicropteryx_cornuta_softmasked.fasta
```

## Step 2 Feature Annotation – Running with IsoSeq (long-read) data
Important directories for running braker on long-read IsoSeq data:
```
# soft-masked genome truncated and white-space removed
/blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/Ncor.softmasked_tk.fasta

# iso-seq fastq, unzipped
/blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_longread/SRR14882579.fastq
```

To run braker on long-read data, I am following instructions from [my annotation](https://github.com/amandamarkee/actias-luna-genome/blob/main/annotation-notes-ACTIVE.md#step-3-of-braker2-feature-annotation--running-with-isoseq-long-read-data) for _Actias luna_ which primarily takes steps from [this github's modified pipeline](https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md#braker2) for running braker2 with long-read IsoSeq evidence. However, I use this modified pipeline up until the mapping step. For mapping and collapsing, I use [PacBio bulk IsoSeq workflow](https://isoseq.how/classification/workflow.html) using the pbmm and collapse funcitons to ensure script compatibility.

First, I run pbmm (PacBio minimap) to map the IsoSeq transcripts back to the _Neomicropteryx cornuta_ genome using the following script:
```
#!/bin/bash
#SBATCH --job-name=pbmm2
#SBATCH -o pb_minimap.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 24:00:00
#SBATCH -c 32

module load pbmm2
module load isoseq3

pbmm2 align --preset ISOSEQ --sort /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_longread/SRR14882579.fastq \
/blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/Ncor.softmasked_tk.fasta \
nc_isoseq_mapped.bam
```

After mapping, I use the PacBio IsoSeq collapse function, which serves the same purpose as Cupcake, and collapses redundant IsoForms:
```
#!/bin/bash
#SBATCH --job-name=isoseq_collapse
#SBATCH -o isoseq_collapse.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 24:00:00
#SBATCH -c 32

module load isoseq3

isoseq3 collapse nc_isoseq_mapped.bam nc_isoseq_collapse.gff
```

After collapsing, I return to the modified BRAKER2 pipeline to conduct GeneMarkS-T predictions. I first run the Augustus stringtie2fa.py script within the BRAKER2 pipeline using the masked genome, and output collapsed .gff file from the previous isoseq collapse step:
```
# stringtie2fa.py -g genome.fa -f isoseq.collapsed.gff -o isoseq_collapsed.fa
module load python
python stringtie2fa.py -g /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/Ncor.softmasked_tk.fasta -f nc_isoseq_collapse.gff -o nc_isoseq_collapsed.fa
```

stringtie2fa.py:
```
#!/usr/bin/env python3

# Author: Katharina J. Hoff
# E-Mail: katharina.hoff@uni-greifswald.de
# Last modified on November 8th 2021
#
# This Python script extracts exon features from a GTF file, excises
# corresponding sequence windows from a genome FASTA file, stitches the
# codingseq parts together, makes reverse complement
# if required
# Output file is:
#    * file with mRNA squences in FASTA format
# Beware: the script assumes that the gtf input file is sorted by coordinates!
# This script is also compatible with cupcake gtf format

try:
    import argparse
except ImportError:
    raise ImportError(
        'Failed to import argparse. Try installing with \"pip3 install argparse\"')

try:
    import re
except ImportError:
    raise ImportError(
        'Failed to import argparse. Try installing with \"pip3 install re\"')

try:
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    raise ImportError(
        'Failed to import biophython modules. Try installing with \"pip3 install biopython\"')


parser = argparse.ArgumentParser(
    description='Generate *.codingseq and *.aa FASTA-format files from genes \
                 in a GTF-file produced by AUGUSTUS auxprogs tool joingenes \
                 and a corresponding genomic FASTA-file.')
parser.add_argument('-g', '--genome', required=True,
                    type=str, help='genome sequence file (FASTA-format)')
parser.add_argument('-o', '--out', required=True, type=str,
                    help="name stem pf output file with coding sequences and \
                    protein sequences (FASTA-format); will be extended by \
                    .codingseq/.aa")
parser.add_argument('-p', '--print_format_examples', required=False, action='store_true',
                    help="Print gtf input format examples, do not perform analysis")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-f', '--gtf',
                    type=str, help='file with CDS coordinates (GTF-format)')
args = parser.parse_args()


if args.print_format_examples:
    print('This script requires an annotation of a transcript with exon ' +
        'features in GTF format. We here provide a compatible ' +
        'example. ' +
        'This script will only process the exon lines. The input data may contain other feature lines that ' +
        'are ignored.')
    print('\nGTF format example:\n')
    print('1\tStringTie\ttranscript\t1555206\t1572018\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; cov "5.737374"; FPKM "5.261884"; TPM "18.775906";\n' +
         '1\tStringTie\texon\t1555206\t1555441\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "1"; cov "6.080509";\n' +
         '1\tStringTie\texon\t1565008\t1565346\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "2"; cov "5.917404";\n' +
         '1\tStringTie\texon\t1571901\t1572018\t1000\t-\t.\tgene_id "STRG.16"; transcript_id "STRG.16.1"; exon_number "3"; cov "4.533898";\n')
    print('\nThis script has successfully been tested with GTF format produced by Stringtie2.')
    exit(0)

# output file names:
mrnaFile = args.out + ".mrna"

# Read GTF file exon entries for transcripts
tx2seq = {}
tx2str = {}
mrna = {}

if args.gtf:
    try:
        with open(args.gtf, "r") as gtf_handle:
            for line in gtf_handle:
                if re.match(
                        r"\S+\t[^\t]+\texon\t\d+\t\d+\t\S+\t\S\t\.\t.*transcript_id \"(\S+)\"", line):
                    #print("I attempt to store exon info")
                    seq_id, st, en, stx, tx_id = re.match(
                        r"(\S+)\t[^\t]+\texon\t(\d+)\t(\d+)\t\S+\t(\S)\t\.\t.*transcript_id \"(\S+)\"", line).groups()
                    if seq_id not in mrna:
                        mrna[seq_id] = {}
                    if tx_id not in mrna[seq_id]:
                        mrna[seq_id][tx_id] = []
                    mrna[seq_id][tx_id].append(
                        {'start': int(st), 'end': int(en), 'strand': stx})
                    if not tx_id in tx2seq:
                        tx2seq[tx_id] = seq_id
                        tx2str[tx_id] = stx
    except IOError:
        print("Error: Failed to open file " + args.gtf + "!")
        exit(1)
else:
    print("Error: No annotation file in GTF format was provided!")
    exit(1)

# Read genome file (single FASTA entries are held in memory, only), extract
# exon sequence windows
seq_len = {}
mrnaseq = {}
try:
    with open(args.genome, "r") as genome_handle:
        for record in SeqIO.parse(genome_handle, "fasta"):
            seq_len[record.id] = len(record.seq)
            if record.id in mrna:
                for tx in mrna[record.id]:
                    #print("I do something for tx")
                    if tx not in mrnaseq:
                        if mrna[record.id][tx][0]['strand'] == '.':
                            descr = tx + ' strand_unknown'
                        else:
                            descr = tx
                        mrnaseq[tx] = SeqRecord(Seq(""), id=tx, description=descr)
                    nExons = len(mrna[record.id][tx])
                    for i in range(0, nExons):
                        mrna_line = mrna[record.id][tx][i]
                        mrnaseq[tx].seq += record.seq[mrna_line['start'] - 1:mrna_line['end']]
                        if i == (nExons - 1) and mrna_line['strand'] == '-':
                            mrnaseq[tx].seq = mrnaseq[tx].seq.reverse_complement()
except IOError:
    print("Error: Failed to open file " + args.genome + "!")
    exit(1)

# Print mRNA sequences to file
try:
    with open(mrnaFile, "w") as mrna_handle:
        for tx_id, seq_rec in mrnaseq.items():
            SeqIO.write(seq_rec, mrna_handle, "fasta")
except IOError:
    print("Error: Failed to open file " + mrnaFile + "!")
    exit(1)
```

After running stringtie2fa.py, I ran GeneMarkST to us the IsoSeq fata to train the gene model:
```
module load genemark_s/t-3.10.001
gmst.pl --strand direct nc_isoseq_collapsed.fa.mrna --output gmst.out --format GFF
```

Lastly, I use the GeneMarkS-T coordinates and the long-read transcripts to create a gene set in GTF format. Note, the gmst2globalCoords.py script is only in the long_read BRAKER documentation, so you must export this script following the [installation instructions](https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md#installation). The output file should be gmst.global.gtf:
```
gmst2globalCoords.py -t nc_isoseq_collapse.gff -p gmst.out -o gmst.global.gtf -g /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/Ncor.softmasked_tk.fasta
```

I want to see how the gene model will look by combining all three outputs (protein evidence, RNA-seq evidence, and IsoSeq long read evidence), as well as the IsoSeq long read evidence on it's own. To assess combined model quality, I use TSEBRA to combine the gene models and then run BUSCO on the resulting .gtf file (converted to .aa) to assess quality. Please see the "Evaluate gene models" section for these results.

## Step 3 Gene Model Evaluation – Running BUSCO on individual gene models
## (a) From Arthropod protein evidence 
```
sbatch Nc_prot_model_busco.sh
```
```
#!/bin/bash
#SBATCH --job-name=Nc_prot_busco
#SBATCH -o Nc_prot_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="/blue/kawahara/amanda.markee/Aluna_genome/aluna_assembly/hifiasm/BUSCO/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/config"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/braker/augustus.hints.aa \
 -o ./Nc_prot_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12
```
```
	***** Results: *****

	C:77.9%[S:72.5%,D:5.4%],F:1.9%,M:20.2%,n:5286	   
	4116	Complete BUSCOs (C)			   
	3831	Complete and single-copy BUSCOs (S)	   
	285	Complete and duplicated BUSCOs (D)	   
	101	Fragmented BUSCOs (F)			   
	1069	Missing BUSCOs (M)			   
	5286	Total BUSCO groups searched	
```

Eek. Results for just protein from arthropod database are poor. Let's retry using the Endopterygota database instead of Lepidotpera specific. I changed the script as shown below, and reran.

```
#!/bin/bash
#SBATCH --job-name=Nc_prot_busco
#SBATCH -o Nc_prot_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=amanda.markee@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="/blue/kawahara/amanda.markee/Aluna_genome/aluna_assembly/hifiasm/BUSCO/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/amanda.markee/Aluna_genome/aluna_annotation/braker2/Augustus/config"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/amanda.markee/neomicropteryx_annotation/braker2/braker_prot/braker/augustus.hints.aa \
 -o ./Nc_prot_busco_out \
 -l /data/reference/busco/v5/lineages/endopterygota_odb10 \
 -m protein -c 12
```

```
	***** Results: *****

	C:90.4%[S:85.6%,D:4.8%],F:3.3%,M:6.3%,n:2124	   
	1921	Complete BUSCOs (C)			   
	1819	Complete and single-copy BUSCOs (S)	   
	102	Complete and duplicated BUSCOs (D)	   
	71	Fragmented BUSCOs (F)			   
	132	Missing BUSCOs (M)			   
	2124	Total BUSCO groups searched
```
Much better!

## Step 4 Combining Gene Models – Running TSEBRA transcript collapser 


## Step 5 Combined Model Evaluation – Running BUSCO on combined gene model
