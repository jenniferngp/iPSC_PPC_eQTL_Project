     ca#message(paste("0:", Sys.info()["nodename"]))

setwd("/projects/PPC/analysis/ppc_eqtls")

source("scripts/packages.R"  )
source("scripts/input_data.R")
source("scripts/functions.R" )

suppressPackageStartupMessages(library(rhdf5))

option_list   = list(make_option("--taskid"          , type="integer"  , default=0      , help="SGE task ID"       , metavar="character"),
					 make_option("--name"            , type="character", default="test1", help="output folder name", metavar="character"),
					 make_option("--gene_file"       , type="character", default=""     , help="gene IDs file"     , metavar="character"),
					 make_option("--interaction_file", type="character", default=""     , help="file containing covariates for interaction test", metavar="character"),
					 make_option("--analysis"        , type="character", default="gene" , help="gene or isoform"   , metavar="character"),
                     make_option("--run_analysis_int", type="character", default=""     , help="run int or not"    , metavar="character")
					) 
					
opt_parser        = OptionParser(option_list=option_list)
opt               = parse_args(opt_parser)
taskid            = opt$taskid
name              = opt$name
gene_file         = opt$gene_file
interaction_file  = opt$interaction_file
analysis          = opt$analysis
run_analysis_int  = opt$run_analysis_int

#message(paste("1:", Sys.info()["nodename"]))

tryCatch(
{

        # taskid            = 2
        # name              = "test_pca.20"
        # gene_file         = "/frazer01/projects/PPC/analysis/ppc_eqtls/pipeline/3.2.eqtls/genes_to_test.txt"
        # interaction_file  = NULL
        # analysis          = "gene"
        # run_analysis_int  = FALSE
    
        #message(paste("2:", Sys.info()["nodename"]))

        if(analysis == "gene")
        {
            geneinfo   = fread("pipeline/1.2.expression/gene_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
            exp_folder = paste(getwd(), "pipeline", "1.2.expression", "tpm_gene", sep = "/")
            gt_folder  = paste(getwd(), "pipeline", "1.3.genotype"  , "tpm_gene", sep = "/")
            out_folder = paste(getwd(), "pipeline", "3.2.eqtls"     , "tpm_gene", sep = "/")
        }
        if(analysis == "isoform")
        {
            geneinfo   = fread("pipeline/1.2.expression/isoform_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
            exp_folder = paste(getwd(), "pipeline", "1.2.expression", "use_isoform", sep = "/")
            gt_folder  = paste(getwd(), "pipeline", "1.3.genotype"  , "use_isoform", sep = "/")
            out_folder = paste(getwd(), "pipeline", "3.2.eqtls"     , "use_isoform", sep = "/")
        }

        if (gene_file == ""){gene_ids = geneinfo$transcript_id}
        if (gene_file != ""){gene_ids = readLines(gene_file)}

        rownames(geneinfo) = geneinfo$transcript_id
        geneinfo           = geneinfo[gene_ids,]

        chromsizes  = read.table("/frazer01/reference/public/hg19/hg19.size.txt", sep = "\t", header = FALSE, col.names = c("chrom", "size"))[,1:2]
        gene_id     = geneinfo[taskid, "transcript_id"]
        message(paste(taskid, gene_id))

        tmp_folder  = paste("~/scratch", paste(name, analysis, gene_id, sep = "."), sep = "/")

        system(paste("rm -rf", tmp_folder))

        out_folder  = paste(getwd()   , "pipeline/3.2.eqtls", "eqtls_by_gene", paste(name, analysis, sep = "."), sep = "/")

        dir.create(out_folder, showWarnings = FALSE)
        dir.create(tmp_folder, showWarnings = FALSE)

        run_analysis_int = FALSE
        run_analysis_qtl = TRUE

        if(file.exists(paste(out_folder, paste("fdr", gene_id, "txt", sep = "."), sep = "/")) == T)
        {
            fdr_out = fread(paste(out_folder, paste("fdr", gene_id, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE, data.table = FALSE)
            totest  = fdr_out[fdr_out$fdr <= 0.05,]

            run_analysis_qtl = FALSE

        }

        if(file.exists(paste(out_folder, paste("qtl", gene_id, "txt", sep = "."), sep = "/")) == F)
        {
            run_analysis_qtl = TRUE
        }

        vars0             = readLines(paste(getwd(), "pipeline/3.2.eqtls", "vars", paste(name, analysis, "txt", sep = "."), sep = "/"))
        covariates        = add_rownames(fread    ("pipeline/3.1.covariates/covariates.txt", sep = "\t", header = TRUE , data.table = FALSE))
        metadata          = add_rownames(fread    ("pipeline/3.1.covariates/metadata.txt"  , sep = "\t", header = TRUE , data.table = FALSE))
        kinship           = add_rownames(fread    ("pipeline/3.1.covariates/kinship.txt"   , sep = "\t", header = TRUE , data.table = FALSE))

        expdata              = add_rownames(fread(paste(exp_folder, paste(           gene_id, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE , data.table = FALSE))
        gtinfo               =              fread(paste(gt_folder , paste("gt_info", gene_id, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE , data.table = FALSE)
        gtdata               = add_rownames(fread(paste(gt_folder , paste("gt_data", gene_id, "txt", sep = "."), sep = "/"), sep = "\t", header = TRUE , data.table = FALSE))

        gt_file              = paste(tmp_folder, paste("TMP", gene_id, "gt"     , "txt", sep = "."), sep = "/")
        exp_file             = paste(tmp_folder, paste("TMP", gene_id, "exp"    , "csv", sep = "."), sep = "/")
        cov_file             = paste(tmp_folder, paste("TMP", gene_id, "cov"    , "csv", sep = "."), sep = "/")
        kin_file             = paste(tmp_folder, paste("TMP", gene_id, "kin"    , "csv", sep = "."), sep = "/")
        bed_file             = paste(tmp_folder, paste("TMP", gene_id, "gt"            , sep = "."), sep = "/")
        h5_file              = paste(tmp_folder, paste("TMP", gene_id, "gt"     , "h5" , sep = "."), sep = "/")
        qtl_file_tmp         = paste(tmp_folder, paste("TMP", gene_id, "qtl"    , "csv", sep = "."), sep = "/")
        geneloc_file         = paste(tmp_folder, paste("TMP", gene_id, "geneloc", "txt", sep = "."), sep = "/")
        snploc_file          = paste(tmp_folder, paste("TMP", gene_id, "snploc" , "txt", sep = "."), sep = "/")
        qtl_file             = paste(out_folder, paste("qtl", gene_id,            "txt", sep = "."), sep = "/")
        fdr_file             = paste(out_folder, paste("fdr", gene_id,            "txt", sep = "."), sep = "/")
        #int_file             = paste(out_folder, paste("int", gene_id,            "txt", sep = "."), sep = "/")

        covariates           = add_rownames(merge(covariates, metadata, by = "row.names"))
        expdata              = as.data.frame(t(as.matrix(expdata)))["norm",]
        covariates           = covariates[colnames(expdata),]
        covariates           = covariates[covariates$wgs_id %in% colnames(gtdata),]
        covariates$sex       = ifelse(covariates$sex == "M", 1, 0)
        gtdata               = gtdata[,covariates$wgs_id]
        colnames(gtdata)     = rownames(covariates)
        expdata              = expdata[,rownames(covariates)]
        kinship_all          = kinship
        kinship              = kinship[covariates$wgs_id, covariates$wgs_id]
        colnames(kinship)    = rownames(covariates)
        rownames(kinship)    = rownames(covariates)
        rownames(gtinfo )    = gtinfo$id
        covariates_all       = covariates
        covariates           = covariates[,vars0]
        rownames(expdata)    = c("trait")
        gtdata               = as.matrix(gtdata)
        geneloc              = geneinfo[geneinfo$transcript_id == gene_id, c("gene_id", "chrom", "start", "end")]
        geneloc$gene_id      = gene_id
        snploc               = gtinfo[,c("id", "chrom", "pos")]
        snploc$chrom         = paste("chr", gsub("chr", "", snploc$chrom), sep = "")
        chrom                = geneloc$chrom

        if(length(gtdata[is.na(gtdata) == TRUE]) > 0){gtdata[is.na(gtdata) == TRUE] = 0}
        gtdata = as.data.frame(gtdata)

        fwrite(gtdata        , gt_file       , sep = "\t", row.names = TRUE , col.names = TRUE)
        fwrite(expdata       , exp_file      , sep = "," , row.names = TRUE , col.names = TRUE)
        fwrite(covariates    , cov_file      , sep = "," , row.names = TRUE , col.names = TRUE)
        fwrite(kinship       , kin_file      , sep = "," , row.names = TRUE , col.names = TRUE)
        fwrite(geneloc       , geneloc_file  , sep = "\t", row.names = FALSE, col.names = TRUE)
        fwrite(snploc        , snploc_file   , sep = "\t", row.names = FALSE, col.names = TRUE)
    
        if((run_analysis_qtl == TRUE)|(run_analysis_int == TRUE))
        {
            if(run_analysis_qtl == TRUE)
            {
                print("here")
                
                qtl_0_list     = run_eqtl(gene_id, analysis, tmp_folder, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, chrom)
                qtl            = qtl_0_list$qtl
                fdr            = qtl_0_list$lead
                lead_var       = fdr[1, "id"]
                lead_vars      = list()
                qtl_list       = list()
                fdr_list       = list()
                qtl$type       = 0
                fdr$type       = 0
                lead_vars[[1]] = lead_var
                qtl_list [[1]] = qtl
                fdr_list [[1]] = fdr

                if(fdr$fdr <= 0.05)
                {
                    for(type in 1:5)
                    {
                        covariates[,paste("gt", type, sep = "")] = as.numeric(gtdata[lead_var, rownames(covariates)])

                        fwrite(covariates, cov_file, sep = "," , row.names = TRUE , col.names = TRUE)
                        qtl_1_list            = run_eqtl(gene_id, analysis, tmp_folder, expdata, gtinfo, gtdata, h5_file, cov_file, kin_file, qtl_file_tmp, geneloc_file, snploc_file, gt_file, chrom)
                        qtl                   = qtl_1_list$qtl
                        fdr                   = qtl_1_list$lead
                        lead_var              = fdr[1, "id"]
                        qtl$type              = type
                        fdr$type              = type
                        lead_vars[[type + 1]] = lead_var
                        qtl_list [[type + 1]] = qtl
                        fdr_list [[type + 1]] = fdr

                        if(fdr$fdr > 0.05) {break}
                    }
                }

                qtl_out = as.data.frame(rbindlist(qtl_list), stringsAsFactors = FALSE)
                fdr_out = as.data.frame(rbindlist(fdr_list), stringsAsFactors = FALSE)
                totest  = fdr_out[fdr_out$fdr <= 0.05,]

                fwrite(qtl_out, qtl_file, sep = "\t", row.names = FALSE, col.names = TRUE)
                fwrite(fdr_out, fdr_file, sep = "\t", row.names = FALSE, col.names = TRUE)

                message(paste(qtl_file, "written"))
            } else {
                message(paste(taskid, qtl_file, "exists"))
            }

            system(paste("rm -r", tmp_folder))
        } else {
            message(paste(paste(out_folder, paste("qtl", gene_id, "txt", sep = "."), sep = "/"), "exists"))
        }
    
    system(paste("rm -r", tmp_folder))
}, 
    error = function(cond)
    {
        geneinfo   = fread("pipeline/1.2.expression/gene_info.txt", sep = "\t", header = TRUE , data.table = FALSE)
        rownames(geneinfo) = geneinfo$transcript_id
        gene_ids   = readLines(gene_file)
        geneinfo   = geneinfo[gene_ids,]
        message(paste("error:", cond, taskid, gene_id, Sys.info()["nodename"]))
        #message(paste("Error:", cond))
        system(paste("rm -r", tmp_folder))
    }
)

