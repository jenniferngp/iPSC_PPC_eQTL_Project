# iPSC-PPC eQTL Repository

[![DOI](https://zenodo.org/badge/666282159.svg)](https://zenodo.org/badge/latestdoi/666282159)

This repository contains code and scripts for the analyses conducted in the manuscript. 

## 0. rna_pipeline.pbs
Commands for processing bulk RNA-seq data 

## 0. identity.pbs
Commands for running sample identity

## 0. liftover.R
This script shows how to perform LiftOver from hg38 -> hg19 for GTEx v8 eQTLs.

## 1. Map eQTLs
This pipeline was developed by Matteo D'Antonio, PhD.
We are currently organizing our eQTL pipeline and will publish as soon as we can. 

## 2. Colocalization of eQTLs and GWAS
See `run_coloc.R `. This script shows how colocalization was run in this study. We used p-values, MAF, and N as input to identify shared eQTL and GWAS signals.

## 3. Finemapping eQTL signals
See `finemap.R`. We fine-mapped each eQTL signal using coloc using p-values, MAF, and N as input. 

## 4. eQTL networks
See `eqtl_networks.R`. This script shows how to identify eQTL modules using colocalized eQTL pairs as input. We performed the community clustering by chromosome to increase sensitivity.

## 5. LD between iPSC-PPC, adult islet, and adult panc eQTLs
See `LD_analysis.ipynb`. This notebook for LD analysis between eQTL pairs and how it was used to update the eQTL annotations.

## 6. Process GWAS colocalization results
See `process_gwas_coloc.R`. This script shows how GWAS colocalization data was processed.

## 7. check_numbers.ipynb
This script shows how the numbers in the manuscript were calculated using the supplementary tables.





