setwd("/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/functions.R" )
source("scripts/input_data.R")

suppressMessages(library(rtracklayer  ))
suppressMessages(library(GenomicRanges))

option_list = list(make_option("--taskid"  , type = "integer"  , default = 0 , help = "SGE task ID"))

opt_parser = OptionParser(option_list = option_list)
opt        = parse_args(opt_parser)
taskid     = opt$taskid

liftover  = function(bed)
{
    suppressMessages(library(rtracklayer  ))
    suppressMessages(library(GenomicRanges))
    
    path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch   = import.chain(path)
    
    gr38 = makeGRangesFromDataFrame(bed, T)
    gr19 = liftOver(gr38, ch)
    return(gr19) 
}

qtls2run = fread("eqtls.txt", data.table = F)
transcript_id = qtls2run[taskid,]$transcript_id
filename = qtls2run[taskid,]$filename

qtl_data = fread(filename, data.table = F)
qtl_data$chr = as.numeric(unlist(lapply(qtl_data$id, function(x) { unlist(strsplit(x, "_"))[1] })))
qtl_data$pos = as.numeric(unlist(lapply(qtl_data$id, function(x) { unlist(strsplit(x, "_"))[2] })))

tolift = data.frame(chr = qtl_data$chr, start = qtl_data$pos, end = qtl_data$pos, id = qtl_data$id)
gr19   = data.frame(liftover(tolift)) %>% select(end, id) %>% dplyr::rename(hg19_pos = end)

qtl_data = merge(qtl_data, gr19, by = "id") 
qtl_data$hg38_pos = qtl_data$pos
qtl_data$pos = qtl_data$hg19_pos

qtl_data$chr = gsub("chr", "", qtl_data$chr)
qtl_data$id = paste("VAR", qtl_data$chr, qtl_data$pos, qtl_data$ref, qtl_data$alt, sep = "_")

colnames(qtl_data)[which(colnames(qtl_data) == "slope")] = "beta"
colnames(qtl_data)[which(colnames(qtl_data) == "slope_se")] = "se"
colnames(qtl_data)[which(colnames(qtl_data) == "pval_nominal")] = "pval"

qtl_data$maf = ifelse(qtl_data$af > 0.5, 1-qtl_data$af, qtl_data$af)

out_filename = paste(transcript_id, "hg19", sep = ".")
fwrite(qtl_data, out_filename, row.names = F, sep = "\t")

message(paste("Saved:", out_filename))