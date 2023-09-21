setwd("/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/functions.R" )
source("scripts/input_data.R")

# manifest containing trait ids
manifest = fread("pipeline/3.3.gwas/meta/manifest_to_use_April2022.txt", data.table = F)

# list all files from gwas coloc
files = list.files("pipeline/3.3.gwas/coloc", full.names = T)

print(Sys.time())

# interate through each file (corresponds to each eGene))
out = as.data.frame(rbindlist(lapply(files, function(f)
{
    load(f)

    # iterate through each trait
    as.data.frame(rbindlist(lapply(names(outlist), function(trait_id)
    {

        # iterate through each conditional eqtl
        as.data.frame(rbindlist(lapply(names(outlist[[trait_id]]), function(type)
        {

            # get pp for each model
            this = data.frame(t(data.frame(outlist[[trait_id]][[type]]$summary )))

            # get finemapping results
            results = outlist[[trait_id]][[type]]$results 

            # save
            this = this %>% mutate(trait_id = trait_id, 
                                   type = unlist(strsplit(type, " "))[2], # i.e., discovery order
                                   transcript_id = unlist(strsplit(f, "[.]"))[4],
                                   topsnp_pp = max(results$SNP.PP.H4), # pp of putative causal variant
                                   topsnp = paste(results[results$SNP.PP.H4 == max(results$SNP.PP.H4),]$snp, collapse = ";"), # lead putative causal variant
                                   topsnp_gwas_pval = paste(results[results$SNP.PP.H4 == max(results$SNP.PP.H4),]$pvalues.df2, collapse = ";"), 
                                   topsnp_eqtl_pval = paste(results[results$SNP.PP.H4 == max(results$SNP.PP.H4),]$pvalues.df1, collapse = ";"), 
                                   max_pp = max(this[1,c(2:6)]), # maximum pp 
                                   likely_hyp = colnames(this)[2:6][which(this[1,c(2:6)] == max(this[1,c(2:6)]))] # most likely model
                                  )
            return(this)
        })))
    })))
})))
