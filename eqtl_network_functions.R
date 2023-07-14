summarize_communities_gene = function(analysis, c)
{
    geneinfo = fread("pipeline/1.2.expression/gene_isoform_info_notation.txt", data.table = F)
    out = as.data.frame(rbindlist(lapply(c(1:length(c)), function(id)
    {
        genes = c[[id]]
        ppc_ids = genes[which(!genes %like% "Islets" & !genes %like% "Pancreas")]
        islet_ids = genes[which(genes %like% "Islets")]
        panc_ids = genes[which(genes %like% "Pancreas")]

        if (analysis == "gene")
        {
            ppc_genes  = unlist(lapply(ppc_ids   , function(x) { unlist(strsplit(x, " "))[2] })) 
            islet_genes = unlist(lapply(islet_ids, function(x) { unlist(strsplit(x, " "))[2] }))
            panc_genes  = unlist(lapply(panc_ids , function(x) { unlist(strsplit(x, " "))[2] }))
        } else
        {
            ppc_genes = unlist(lapply(ppc_ids, function(x) { unlist(strsplit(x, " "))[2] }))
            islet_genes = unlist(lapply(islet_ids, function(x) { gsub("Islets ", "", unlist(strsplit(x, "[.]"))[1]) }))
            panc_genes = unlist(lapply(panc_ids, function(x) { simplify_id(unlist(strsplit(x, ":"))[5]) }))
        }

        data.frame(analysis              = analysis,
                   locus_id              = id,
                   signals               = paste(genes, collapse = ","),
                   number_signals        = length(unique(genes)),
                   number_ppc_signals    = length(unique(ppc_ids)),
                   number_islet_signals  = length(unique(islet_ids)),
                   number_panc_signals   = length(unique(panc_ids)),
                   number_genes          = length(unique(c(ppc_genes, islet_genes, panc_genes))),
                   number_ppc_genes      = length(unique(ppc_genes)),
                   number_islet_genes    = length(unique(islet_genes)),
                   number_panc_genes     = length(unique(panc_genes)),
                   same_genes.all        = length(unique(Reduce(intersect, list(ppc_genes,islet_genes,panc_genes)))),
                   same_genes.ppc_islet  = length(unique(Reduce(intersect, list(ppc_genes,islet_genes)))),
                   same_genes.ppc_panc   = length(unique(Reduce(intersect, list(ppc_genes,panc_genes)))),
                   diff_genes.ppc.islet_panc = length(unique(setdiff(ppc_genes, c(islet_genes, panc_genes)))),
                   diff_genes.islet_panc.ppc = length(unique(setdiff(c(islet_genes, panc_genes), ppc_genes))),
                   diff_genes.ppc.islet  = length(unique(setdiff(ppc_genes, islet_genes))),
                   diff_genes.islet.ppc  = length(unique(setdiff(islet_genes, ppc_genes))),
                   diff_genes.ppc.panc   = length(unique(setdiff(ppc_genes, panc_genes))),
                   diff_genes.panc.ppc   = length(unique(setdiff(panc_genes, ppc_genes)))
                  )
    })))

    return(out)
}

summarize_communities_isof = function(analysis, c)
{
  out = as.data.frame(rbindlist(lapply(c(1:length(c)), function(id)
    {
        signals = c[[id]]

        gene_list = list()

        for (s in signals)
        {
            if (s %like% "Pancreas")
            {
                gene_list[["Pancreas"]] = c(gene_list[["Pancreas"]], simplify_id(unlist(strsplit(s, ":"))[5]))
            } else if (s %like% "Islets")
            {
                gene_list[["Islets"]] = c(gene_list[["Islets"]], simplify_id(unlist(strsplit(s, " "))[2]))
            } else
            {
                gene_list[["iPSC-PPC"]] = simplify_id(geneinfo[geneinfo$transcript_id %like% unlist(strsplit(s, " "))[2],]$gene_id)
            }
        }

        row = data.frame(analysis = "isof", module_id = paste("AS", id, sep = "_"), signals = a[row,]$signals)
        
        row$number_signals = length(signals)
        row$number_ppc_signals = length(signals[which(!signals %like% "Pancreas" & !signals %like% "Islets")])
        row$number_islet_signals = length(signals[which(signals %like% "Islets")])
        row$number_panc_signals = length(signals[which(signals %like% "Pancreas")])
        row$number_genes = length(unique(c(gene_list[["Pancreas"]], gene_list[["Islets"]], gene_list[["iPSC-PPC"]])))
        row$number_ppc_genes = length(unique(gene_list[["iPSC-PPC"]]))
        row$number_islet_genes = length(unique(gene_list[["Islets"]]))
        row$number_panc_genes = length(unique(gene_list[["Pancreas"]]))
        row$same_genes.all = length(intersect(gene_list[["Pancreas"]], intersect(gene_list[["Islets"]], gene_list[["iPSC-PPC"]])))
        row$same_genes.ppc_islet = length(intersect(gene_list[["iPSC-PPC"]], gene_list[["Islets"]]))
        row$same_genes.ppc_panc = length(intersect(gene_list[["iPSC-PPC"]], gene_list[["Pancreas"]]))
        row$diff_genes.ppc.islet_panc = length(setdiff(gene_list[["iPSC-PPC"]], c(gene_list[["Islets"]], gene_list[["Pancreas"]])))
        row$diff_genes.islet_panc.ppc = length(setdiff(c(gene_list[["Islets"]], gene_list[["Pancreas"]]), gene_list[["iPSC-PPC"]]))
        row$diff_genes.ppc.islet = length(setdiff(gene_list[["iPSC-PPC"]], gene_list[["Islets"]]))
        row$diff_genes.islet.ppc = length(setdiff(gene_list[["Islets"]], gene_list[["iPSC-PPC"]]))
        row$diff_genes.ppc.panc = length(setdiff(gene_list[["iPSC-PPC"]], gene_list[["Pancreas"]]))
        row$diff_genes.panc.ppc = length(setdiff(gene_list[["Pancreas"]], gene_list[["iPSC-PPC"]]))

        return(row)
    })))

    return(out)
}

classify = function(summary)
{
    summary$islet_label = "NA"
    summary$islet_label = ifelse(summary$number_ppc_signals != 0 & summary$number_islet_signals != 0 & summary$same_genes.ppc_islet != 0 & summary$diff_genes.ppc.islet == 0 & summary$diff_genes.islet.ppc == 0, "same", summary$islet_label)
    summary$islet_label = ifelse(summary$number_ppc_signals != 0 & summary$number_islet_signals != 0 & summary$same_genes.ppc_islet == 0 & (summary$diff_genes.ppc.islet != 0 | summary$diff_genes.islet.ppc != 0), "diff", summary$islet_label)
    summary$islet_label = ifelse(summary$number_ppc_signals != 0 & summary$number_islet_signals != 0 & summary$same_genes.ppc_islet != 0 & (summary$diff_genes.ppc.islet != 0 | summary$diff_genes.islet.ppc != 0), "partial", summary$islet_label)
    summary$islet_label = ifelse(summary$number_ppc_signals != 0 & summary$number_islet_signals == 0 & summary$number_panc_signals == 0, "ppc-unique", summary$islet_label)
    summary$islet_label = ifelse(summary$number_ppc_signals != 0 & summary$number_islet_signals == 0 & summary$number_panc_signals != 0, "zero", summary$islet_label)
    summary$islet_label = ifelse(summary$number_ppc_signals == 0 & summary$number_islet_signals != 0 & summary$number_panc_signals == 0, "islet-unique", summary$islet_label)
    summary$islet_label = ifelse(summary$number_ppc_signals == 0 & summary$number_islet_signals == 0 & summary$number_panc_signals != 0, "panc-unique", summary$islet_label)
    summary$islet_label = ifelse(summary$number_ppc_signals == 0 & summary$number_islet_signals != 0 & summary$number_panc_signals != 0, "adult-unique", summary$islet_label)

    summary$panc_label = "NA"
    summary$panc_label = ifelse(summary$number_ppc_signals != 0 & summary$number_panc_signals != 0 & summary$same_genes.ppc_panc != 0 &  summary$diff_genes.ppc.panc == 0 & summary$diff_genes.panc.ppc == 0 , "same", summary$panc_label)
    summary$panc_label = ifelse(summary$number_ppc_signals != 0 & summary$number_panc_signals != 0 & summary$same_genes.ppc_panc == 0 & (summary$diff_genes.ppc.panc != 0 | summary$diff_genes.panc.ppc != 0), "diff", summary$panc_label)
    summary$panc_label = ifelse(summary$number_ppc_signals != 0 & summary$number_panc_signals != 0 & summary$same_genes.ppc_panc != 0 & (summary$diff_genes.ppc.panc != 0 | summary$diff_genes.panc.ppc != 0), "partial", summary$panc_label)
    summary$panc_label = ifelse(summary$number_ppc_signals != 0 & summary$number_islet_signals == 0 & summary$number_panc_signals == 0, "ppc-unique", summary$panc_label)
    summary$panc_label = ifelse(summary$number_ppc_signals != 0 & summary$number_panc_signals == 0 & summary$number_islet_signals != 0, "zero", summary$panc_label)
    summary$panc_label = ifelse(summary$number_ppc_signals == 0 & summary$number_panc_signals != 0 & summary$number_islet_signals == 0, "panc-unique", summary$panc_label)
    summary$panc_label = ifelse(summary$number_ppc_signals == 0 & summary$number_panc_signals == 0 & summary$number_islet_signals != 0, "islet-unique", summary$panc_label)
    summary$panc_label = ifelse(summary$number_ppc_signals == 0 & summary$number_panc_signals != 0 & summary$number_islet_signals != 0, "adult-unique", summary$panc_label)

    summary$final_class = ifelse(summary$islet_label == "ppc-unique"   & summary$panc_label == "ppc-unique"  , "ppc-unique", "NA")
    summary$final_class = ifelse(summary$islet_label == "islet-unique" & summary$panc_label == "islet-unique", "islet-unique", summary$final_class)
    summary$final_class = ifelse(summary$islet_label == "panc-unique"  & summary$panc_label == "panc-unique" , "panc-unique", summary$final_class)
    summary$final_class = ifelse(summary$islet_label == "adult-unique"  & summary$panc_label == "adult-unique", "adult-unique", summary$final_class)

    summary$final_class = ifelse(summary$final_class == "NA", paste(summary$islet_label, summary$panc_label), summary$final_class)

    return(summary)
}

filter_H3_modules = function(summary, coloc_summary)
{
  h3_dist = as.data.frame(rbindlist(lapply(c(1:nrow(summary)), function(row)
  {
      signals = summary[row,]$signals
      signals = unlist(strsplit(signals, ","))

      this = coloc_summary %>% 
          filter(qtl_id.1 %in% signals & qtl_id.2 %in% signals & max_pp >= 0.8)  %>% 
          distinct() %>% 
          mutate(pair_id = paste(qtl_id.1, qtl_id.2))

      this2 = this %>% select(qtl_id.1, qtl_id.2)

      for (i in 1:nrow(this2)) { this2[i, ] = suppressWarnings(sort(this2[i, ])) }

      this2 = this2[!duplicated(this2),] %>% mutate(pair_id = paste(qtl_id.1, qtl_id.2))

      this = merge(this2, this[,c("pair_id", "likely_hyp", "max_pp")], by = "pair_id")

      out = coloc_summary[row,c("analysis", "module_id", "signals", "number_signals")]

      out$PP.H0.abf = nrow(this[this$likely_hyp == "PP.H0.abf",])
      out$PP.H1.abf = nrow(this[this$likely_hyp == "PP.H1.abf",])
      out$PP.H2.abf = nrow(this[this$likely_hyp == "PP.H2.abf",])
      out$PP.H3.abf = nrow(this[this$likely_hyp == "PP.H3.abf",])
      out$PP.H4.abf = nrow(this[this$likely_hyp == "PP.H4.abf",])
          
  })))

  h3_dist$poss_pairs = unlist(lapply(c(1:nrow(h3_dist)), function(x) { ncol(combn(unlist(strsplit(h3_dist[x,]$signals, ",")), 2)) }))
  h3_dist$pct_h3 = h3_dist$PP.H3.abf / h3_dist$poss_pairs
  h3_dist$pct_h4 = h3_dist$PP.H4.abf / h3_dist$poss_pairs
  h3_dist$ratio = h3_dist$PP.H4.abf / h3_dist$PP.H3.abf
  h3_dist$ratio = ifelse(is.infinite(h3_dist$ratio), 1, h3_dist$ratio)

  h3_dist.1 = h3_dist %>% filter(PP.H3.abf != 0) %>% filter(pct_h4 >= 0.3) %>% filter(ratio >= 2) %>% arrange(ratio)
  h3_dist.2 = h3_dist %>% filter(PP.H3.abf == 0 & pct_h4 >= 0.3)

  keep = rbind(h3_dist.1, h3_dist.2)

  return(list(dist = h3_dist, keep = keep))
}

    