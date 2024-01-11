# iPSC-PPC eQTL Repository

[![DOI](https://zenodo.org/badge/666282159.svg)](https://zenodo.org/badge/latestdoi/666282159)

# Analysis
This folder contains scripts that were developed in the study. 

0. `rna_pipeline.pbs`: shell commands to process bulk RNA-seq data
0. `identity.pbs`: shell commands to match bulk RNA-seq sample to the correct donor using WGS calls 
0. `liftover.R`: R script showing how to perform LiftOver from hg38 to hg19 for GTEx v8 eQTLs
1. `map_eqtls.R`: R script showing how to map eQTLs. Developed by Matteo D'Antonio et al., Nature Communications, 2023.
2. `run_coloc.R `: R scripts showing how to perform colocalization between eQTL-eQTL and eQTL-GWAS. P-values, MAF, and N were used as input
3. `finemap.R`: R sciript showing how to fine-map each eQTL association using p-values, MAF, and N as input into coloc
4. `eqtl_networks.R`: script showing how eQTL modules were identified using iGraph with colocalized eQTL pairs as input. Community clustering was performed by chromosome to increase sensitivity.
5. `LD_analysis.ipynb`: jupyter notebook showing how LD analysis was performed between eQTL pairs
6. `process_gwas_coloc.R`: R script showing how GWAS colocalizations was processed and summarized
7. `check_numbers.ipynb`: R jupytern notebook showing how the numbers in the manuscript were calculated using the Supplementary Tables.


# Figures
Contains R jupyter notebooks to generate figures in the study

# Supplemental Tables and Other Data
Available on Figshare (https://figshare.com/projects/Large-scale_eQTL_analysis_of_iPSC-PPC/156987)

# Reference
Nguyen JP, Arthur TD, Fujita K, Salgado BM, Donovan MK, Matsui H, Kim JH, D’Antonio-Chronowska A, D’Antonio M, Frazer KA. eQTL mapping in fetal-like pancreatic progenitor cells reveals early developmental insights into diabetes risk. Nature Communications. 2023 Oct 30;14(1):6928.

D’Antonio M, Nguyen JP, Arthur TD, Matsui H, D’Antonio-Chronowska A, Frazer KA. Fine mapping spatiotemporal mechanisms of genetic variants underlying cardiac traits and disease. Nature communications. 2023 Feb 28;14(1):1132.
