setwd("/projects/PPC/analysis/ppc_eqtls")

# Load libraries
source("scripts/packages.R"  )
source("scripts/input_data.R")
source("scripts/functions.R" )
source("scripts/3.4.coloc_adult/functions.R")
source("scripts/eqtl_network_functions.R")

suppressWarnings(library(igraph))
suppressWarnings(library(lsa))

suppressWarnings(packageVersion("igraph"))

set.seed(5366)

make_network = function(analysis, graph_input, num_snps, chr, res)
{

	tograph = graph_input[[analysis]] %>% 
				filter(chrom == chr) %>% 
                filter(qtl_id.1 != qtl_id.2) %>% 
                select(qtl_id.1, qtl_id.2, tiss.1, tiss.2, PP.H4.abf, nsnps) %>%
               	filter(nsnps >= num_snps)

    g = graph_from_data_frame(tograph, directed = F)
    c = cluster_leiden(g, objective_function = "modularity", n_iterations = 500, resolution = res)
    
    mod = modularity(g, c$membership)
    
    return(list("graph" = g, "community" = c, "input" = tograph, "modularity" = mod))
}

# coloc_list contains all colocalization results for gene and isoform analysis
graphlist = list()
graphlist[["gene"]] = lapply(c(1:22), function(x) { make_network(analysis = "gene", graph_input = coloc_list, num_snps = 500, x, res = 3) })
graphlist[["isof"]] = lapply(c(1:22), function(x) { make_network(analysis = "isof", graph_input = coloc_list, num_snps = 500, x, res = 3) })

names(graphlist[["gene"]]) = paste0("chr", c(1:22))
names(graphlist[["isof"]]) = paste0("chr", c(1:22))

graphlist[["gene"]]$summary = as.data.frame(rbindlist(lapply(c(1:22), function(x) { summarize_communities_gene("gene", graphlist[["gene"]][[paste0("chr", x)]]$community) %>% mutate(locus_id = paste(x, locus_id, sep = "_")) %>% rename(module_id = locus_id) })))
graphlist[["isof"]]$summary = as.data.frame(rbindlist(lapply(c(1:22), function(x) { summarize_communities_isof("isof", graphlist[["isof"]][[paste0("chr", x)]]$community) %>% mutate(locus_id = paste(x, locus_id, sep = "_")) %>% rename(module_id = locus_id) })))

graphlist[["gene"]]$summary = classify(graphlist[["gene"]]$summary)
graphlist[["isof"]]$summary = classify(graphlist[["isof"]]$summary)

out_filename = "graph.robj"

save(graphlist, file = out_filename)


