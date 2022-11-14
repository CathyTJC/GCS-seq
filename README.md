# GCS-seq

GCS-seq is a ChIP-seq based approach that identifies genome-wide gyrase cleavage sites (GCSs) in *E. coli* popluation grown to staionary-phase. In this project, we performed GCS-seq with a panel of differnt fluorquinolone (FQs) antibiotics that are known to have differential killing capacity, aiming to identify the association between cleavage patterns (e.g., number, cleavage strength, location) and bacterial persister levels. 

This repository contains scripts that used to identify genome-wide GCSs and the downstream analysis performed related to the manuscript *Genome-wide mapping of fluoroquinolone-stabilized DNA gyrase cleavage sites displays drug specific effects that correlate with bacterial persistence* 

For full project description, please refer to our manuscript (in revision and available on https://www.biorxiv.org/content/10.1101/2022.10.27.514060v1). 

Raw sequencing data were deposited with GEO under accession number GSE206610.



## File Description

### coverage_read.py
This scripts read in coverage file (processed from Galaxy, available upon request), preprocess, and prepare for GCS calling 

### GCS_calling.py 

This scripts is used for the identification of GCSs.

**Input** the folder contains the coverage data from NGS that preprocessed with coverage_read.py

**Output** files containing GCSs (GCS calling output is avaialbe in folder **GCS_calling**)

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

### GCS_sort_plot.R

Sort GCS based on the cleavage strengths and save the file of sites and strength for motif construction; Further, retain the top sites which are identified as real GCSs for motif construction 

Also heatmap and barplots showng the  top strengths sites for each FQ treatment 

Next step in pipeline for motif construction: pssm_construction.py


### DEseq.R
This script performs differential gene expression analysis and saves the results

**Input** 1. annotation file (GCF_000005845.2_ASM584v2_genomic_insert_hierexon.gffread.gtf) 2. GCS file 3. gene expression files (raw counts)

### GO_prepare_updated.R

This script Add GCS annotation/ prepare for functional enrichment analysis

### GO_cutoff_plot.R

Ranked GCS plot based on cleavage strength

### GO_enrichment_plot.R

Plot GO enrichment results


