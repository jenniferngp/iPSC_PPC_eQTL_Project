# iPSC-PPC eQTL Repository

[![DOI](https://zenodo.org/badge/666282159.svg)](https://zenodo.org/badge/latestdoi/666282159)

# Analysis
This folder contains scripts that were developed in the study to process and analyze RNA-seq and eQTL data.

0. `rna_pipeline.pbs`: shell commands to process bulk RNA-seq data
0. `identity.pbs`: shell commands to match bulk RNA-seq sample to the correct donor using WGS calls 
0. `liftover.R`: R script showing how to perform LiftOver from hg38 to hg19 for GTEx v8 eQTLs
1. `map_eqtls.R`: R script showing how to map eQTLs. Developed by Matteo D'Antonio, PhD (now Assistant Professor at UCSD)
2. `run_coloc.R `: R scripts showing how to perform colocalization between eQTL-eQTL and eQTL-GWAS. P-values, MAF, and N were used as input
3. `finemap.R`: R sciript showing how to fine-map each eQTL association using p-values, MAF, and N as input into coloc
4. `eqtl_networks.R`: script showing how eQTL modules were identified using iGraph with colocalized eQTL pairs as input. Community clustering was performed by chromosome to increase sensitivity.
5. `LD_analysis.ipynb`: jupyter notebook showing how LD analysis was performed between eQTL pairs
6. `process_gwas_coloc.R`: R script showing how GWAS colocalizations was processed and summarized
7. `check_numbers.ipynb`: R jupytern notebook showing how the numbers in the manuscript were calculated using the Supplementary Tables.

# Figures
Contains scripts to generate figures in the study



