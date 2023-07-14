setwd("/frazer01/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/input_data.R")
source("scripts/functions.R" )
source("scripts/coloc_functions.R")
source("scripts/3.4.coloc_adult/functions.R")

suppressMessages(library(coloc))

option_list = list(make_option("--taskid", type = "integer", help = "SGE task id"))

opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)
taskid     = opt$taskid

n_tissues = list("iPSC-PPC" = 107, "Islets" = 420, "Pancreas" = 305)

all_pairs = fread("all_pairs.txt", data.table = F)

filename1 = all_pairs[taskid,]$filename1
filename2 = all_pairs[taskid,]$filename2

qtl_data1 = fread(filename1, data.table = F) %>% select(chrom, pos, ref, alt, pval, beta, se, maf) %>% mutate(id = paste(chrom, pos, sep = "_"))
qtl_data2 = fread(filename2, data.table = F) %>% select(chrom, pos, ref, alt, pval, beta, se, maf) %>% mutate(id = paste(chrom, pos, sep = "_"))

tiss1 = ifelse(qtl_id1 %like% "Islets", "Islets", ifelse(qtl_id1 %like% "Pancreas", "Pancreas", "iPSC-PPC"))
tiss2 = ifelse(qtl_id2 %like% "Islets", "Islets", ifelse(qtl_id2 %like% "Pancreas", "Pancreas", "iPSC-PPC"))

N1 = n_tissues[[tiss1]]
N2 = n_tissues[[tiss2]]

tocoloc = merge(qtl_data1, qtl_data2, by = "id", suffixes = 1:2)

data1 = tocoloc[tocoloc$ref1 == tocoloc$ref2 & tocoloc$alt1 == tocoloc$alt2,]
data2 = tocoloc[tocoloc$ref1 == tocoloc$alt2 & tocoloc$alt1 == tocoloc$ref2,]

data2$beta2 = -data2$beta2

tocoloc        = rbind(data1, data2)
tocoloc$snp_id = paste("VAR", tocoloc$chrom1, tocoloc$pos1, tocoloc$ref1, tocoloc$alt1, sep = "_")

tocoloc = tocoloc %>% mutate(maf1 = as.numeric(maf1), maf2 = as.numeric(maf2))
tocoloc = tocoloc[tocoloc$maf1 > 0 & 
                  tocoloc$maf1 < 1 & 
                  tocoloc$maf2 > 0 & 
                  tocoloc$maf2 < 1 &
                  !is.na(tocoloc$maf1) & 
                  !is.na(tocoloc$maf2) &
                  duplicated(tocoloc$snp_id) == F,]

coloc = suppressWarnings(coloc.abf(dataset1 = list(snp = tocoloc$snp_id, pvalues = tocoloc$pval1, N = N1, type = "quant", MAF = tocoloc$maf1),
                                   dataset2 = list(snp = tocoloc$snp_id, pvalues = tocoloc$pval2, N = N2, type = "quant", MAF = tocoloc$maf2)))


outfile = paste("pipeline/3.4.coloc_adult/allpairs/data", paste(taskid1, "robj", sep = "."), sep = "/")
save(gene_list, file = outfile)

message(paste("Saved:", outfile))