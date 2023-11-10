#!/usr/bin/env Rscript
library("optparse")

option_list <-list(
  optparse::make_option(c("-i", "--idx"), type="integer", default=1,
                        help="PBS_ARRAY_INDEX", metavar="character"),
  optparse::make_option( c("-n", "--ncpus"), type="integer", default=1,
                         help="Number of CPUs to use.", metavar="character"),
  optparse::make_option( c("-b", "--batches"), type="integer", default=7015,
                         help="Number of total batches.", metavar="character")
);
opt_parser <- optparse::OptionParser(option_list=option_list)
opt <- optparse::parse_args(opt_parser)
root <- "/rds/general/project/neurogenomics-lab/ephemeral/rare_disease_min_genes1"


library(MultiEWCE)
library(data.table)
data.table::setDTthreads(threads = opt$ncpus)

#### Load CTD ####
ctd <- MultiEWCE::load_example_ctd(file = "ctd_DescartesHuman.rds")
annotLevel <- 2 
#### Get phenotype gene lists ####
gene_data <- HPOExplorer::load_phenotype_to_genes() 
## Subset to only those found in the CTD
gene_data <- gene_data[gene_symbol %in% rownames(ctd[[annotLevel]]$mean_exp)]
## Subset to only those with >=4 genes
gene_data[,n_gene:=(length(unique(gene_symbol))),by="hpo_id"]
gene_data <- gene_data[n_gene<4,]

#### Split HPO IDs into N chunks ####
ids <- unique(gene_data$hpo_id)
chunks <- split(ids, cut(seq_along(ids),opt$batches,labels = FALSE))
list_names <- chunks[[opt$idx]]
## Define save dir
save_dir <- gsub(":","",
                 file.path(root,
                           paste0("subjob.",paste(list_names,collapse = "_"))
                           )
                 )

#### Optimise memory usage ####
## Reduce the CTD to save memory: 1.3Gb --> 70.1Mb
lvls <- seq(length(ctd))
lvls <- lvls[lvls != annotLevel]
ctd[lvls] <- NULL
annotLevel <- 1 ## important! redefine after removing all other levels
## Reduce the gene_data to save memory: 44Mb --> 66Kb
gene_data <- gene_data[hpo_id %in% list_names,]
## Ensure memory is returned
remove(chunks)
gc()


#### Run enrichment analyses ##### 
all_results <- MultiEWCE::gen_results(
  ctd = ctd,
  list_name_column = "hpo_id",
  list_names = list_names,
  gene_data = gene_data,
  annotLevel = annotLevel,
  reps = 100000,
  cores = 1,#opt$ncpus,
  min_genes = 1,
  force_new = FALSE,
  save_dir = save_dir)

