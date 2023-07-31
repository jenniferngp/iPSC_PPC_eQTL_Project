# iPSC-PPC eQTL Repository

This repository contains code and scripts for the analyses conducted in the manuscript. 

### check_numbers.ipynb
- This script shows how the numbers in the manuscript were calculated using the supplementary tables.

### run_coloc.R 
- This script shows how colocalization was run in this study. We used p-values, MAF, and N as input to identify shared eQTL and GWAS signals.

### liftover.R
- This script shows how to perform LiftOver from hg38 -> hg19 for GTEx v8 eQTLs.

### finemap.R 
- This script shows how to fine-map eQTL and GWAS signals using p-values, MAF, and N as input.

### eqtl_networks.R 
- This script shows how to identify eQTL modules using colocalized eQTL pairs as input. We performed the community clustering by chromosome to increase sensitivity.

### LD_analysis.ipynb 
- This notebook contains commands for how to run LD between pairs of eQTLs and how it was used to update the eQTL annotations.

### run_eqtl.R (pending)
- We are currently organizing our eQTL pipeline and will publish as soon as we can. 
