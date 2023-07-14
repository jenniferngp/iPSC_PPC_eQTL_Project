setwd("/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/functions.R" )
source("scripts/input_data.R")

suppressMessages(library(rtracklayer  ))
suppressMessages(library(GenomicRanges))
suppressMessages(library(coloc        ))

option_list = list(make_option("--taskid"  , type = "integer"  , default = 0 , help = "SGE task ID"))

opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)
taskid     = opt$taskid

qtls2run = fread("eqtls.txt", data.table = F)
transcript_id = qtls2run[taskid,]$transcript_id
filename = qtls2run[taskid,]$filename
N = qtls2run[taskid,]$N 

qtl = fread(filename, data.table = F)

results = finemap.abf(list(snp = qtl$id, 
					   pvalues = qtl$pval, 
					   N = N, 
					   MAF = qtl$maf, 
					   type = "quant", 
					   beta = qtl$beta, 
					   varbeta = (qtl$se**2)/N), 
					   p1 = 1 / nrow(qtl))

out_filename = paste("finemap", paste(transcript_id, "robj", sep = "."), sep = "/")
save(results, file = out_filename)
message(paste("Saved:", out_filename))