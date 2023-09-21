setwd("/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/functions.R" )
source("scripts/input_data.R")

manifest = fread("pipeline/3.3.gwas/meta/manifest_to_use_April2022.txt", data.table = F)
files = list.files("pipeline/3.3.gwas/coloc", full.names = T)

print(Sys.time())
out = as.data.frame(rbindlist(lapply(files, function(f)
{
    load(f)
    as.data.frame(rbindlist(lapply(names(outlist), function(trait_id)
    {
        as.data.frame(rbindlist(lapply(names(outlist[[trait_id]]), function(type)
        {
            this = data.frame(t(data.frame(outlist[[trait_id]][[type]]$summary )))
            results = outlist[[trait_id]][[type]]$results 
            this = this %>% mutate(trait_id = trait_id, 
                                   type = unlist(strsplit(type, " "))[2], 
                                   transcript_id = unlist(strsplit(f, "[.]"))[4],
                                   topsnp_pp = max(results$SNP.PP.H4),
                                   topsnp = paste(results[results$SNP.PP.H4 == max(results$SNP.PP.H4),]$snp, collapse = ";"),
                                   topsnp_gwas_pval = paste(results[results$SNP.PP.H4 == max(results$SNP.PP.H4),]$pvalues.df2, collapse = ";"),
                                   topsnp_eqtl_pval = paste(results[results$SNP.PP.H4 == max(results$SNP.PP.H4),]$pvalues.df1, collapse = ";"),
                                   max_pp = max(this[1,c(2:6)]),
                                   likely_hyp = colnames(this)[2:6][which(this[1,c(2:6)] == max(this[1,c(2:6)]))])
            return(this)
        })))
    })))
})))
