# GCS-seq

GCS-seq is a ChIP-seq based approach that identifies genome-wide gyrase cleavage sites (GCSs) in *E. coli* popluation grown to staionary-phase. In this project, we performed GCS-seq with a panel of differnt fluorquinolone (FQs) antibiotics that are known to have differential killing capacity, aiming to identify the association between cleavage patterns (e.g., number, strength, location) and bacterial persister levels. 

This repository contains scripts that used to identify genome-wide GCSs and the downstream analysis performed related to the manuscript *Genome-wide mapping of fluoroquinolone-stabilized DNA gyrase cleavage sites displays drug specific effects that correlate with bacterial persistence* For full project description, please refer to our manuscript (in submission). 


## File Description

### venn_diagram.R

Draw Venn diagram of identified GCSs from levofloxacin (LEVO), moxifloxacin (MOXI), norfloxacin (NOR), ciprofloxacin (CIP),and gemifloxacin (GEMI) treatment.

**Input**: folder contains the identified GCSs

### utils_data.R

Script with supporting functions

### transcription_sets_association.R

Script visualizes GCSs within different sets of regions with different transcription levels (upstream, downstream, 5' region, and 3' region). 

**Input**: 1. Annotation file: GCF_000005845.2_ASM584v2_genomic_insert_hierexon.gffread.gtf 2. folder contains the identified GCSs 3. RNA expression data (raw read counts)


### pssm_construction.py
This script construct motifs for MOXI, NOR, CIP, GEMI, and LEVO treated samples

**Input**: folder contains the identified GCSs

### plot.R
Plot GCS number/cleavage strength vs persistence level

see **data** folder for numbers used in the scripts

### motif_scoring.py 
The script is used to scan log-odds score of a given seqence or acrross the genome to predict cleavage probability based on obtained motif 

**Input**: 1. folder contains the identified GCSs 2. genome fasta file


### gcs_distribution.R
This script plot the GCS distribution across the chromosome and perform statistical tests

**Input**: 1. folder contains the identified GCSs



### coverage_plot.R
This scripts takes in the raw read count of FLAG and FLAGless control strains, zoom into 2 GCSs (Mu, mobA) and 2 control sites (nuoN, trkH) to plot the coverage depth.

**Input** raw read count data

see **data** folder for numbers used in the scripts


### PCA.R

perform PCA analysis on GCS distribution

**Input**: merged csv containing the cleavage strength for each replicate across 5 FQ treatment 

### GCS_clustering.R

generate GCS hierarchical clustering heat map

**Input**: merged csv containing the cleavage strength for each replicate across 5 FQ treatment 


### GCS_calling.py 

This scripts is used for the identification of GCSs.

**Input** the folder contains the coverage data from NGS
**Output** files containing GCSs

### DEseq.R
This script performs differential gene expression analysis and saves the results

**Input** 1. annotation file (GCF_000005845.2_ASM584v2_genomic_insert_hierexon.gffread.gtf) 2. GCS file 3. gene expression files (raw counts)


