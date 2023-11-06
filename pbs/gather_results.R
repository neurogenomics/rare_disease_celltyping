#### Find all results files ####
ephemeral <- "/rds/general/project/neurogenomics-lab/ephemeral/rare_disease/"
ephemeral_logs <- file.path(dirname(ephemeral),"rare_disease.pbs_output")
f <- list.files(ephemeral, 
                full.names = TRUE)
f2 <- sapply(f, function(x){list.files(x,pattern = ".rds", full.names = TRUE)})
# f3 <- sapply(f2, function(x){file.rename(x,file.path(dirname(x),"gen_results.rds"))})

length(unlist(f2)) 
dat <- (
  lapply(unlist(f2), readRDS) |> 
    data.table::rbindlist(use.names = TRUE, idcol = "batch")
)[,batch:=basename(batch)][,q:=stats::p.adjust(p,method = "bonf")]

#### Find missing IDs ####
ctd <- MultiEWCE::load_example_ctd(file = "ctd_DescartesHuman.rds")
annotLevel <- 2
gene_data <- HPOExplorer::load_phenotype_to_genes()
gene_data <- gene_data[gene_symbol %in% rownames(ctd[[annotLevel]]$mean_exp)]
gene_data[,n_gene:=(length(unique(gene_symbol))),by="hpo_id"]
gene_data <- gene_data[n_gene>=4,][,batch_id:=.GRP, by="hpo_id"]
missing_dat <- gene_data[!hpo_id %in% unique(dat$hpo_id)]
length(unique(missing_dat$hpo_id)) 
head(sort(unique(missing_dat$hpo_name)),100) 

#### Inspect PBS logs ####
## Peak at the most recent pbs output logs 
## to figure out why some of the jobs failed.
logs <- system(paste("ls",ephemeral_logs,"-Artlsh | tail -10"), intern = TRUE)
logs_df <- data.table::fread(text = logs)[,path:=paste0(file.path(ephemeral_logs,V10))]
logs_df[,type:=ifelse(grepl("pbs.e",V10,fixed = TRUE),"e","o")]
out <- lapply(logs_df[type=="e",]$path,function(x){
  message("\n~~~~~~~~~~",basename(x),"~~~~~~~~~~")
  readLines(x)|>cat(sep = "\n")
})

#### Compute % of significant tests ####
nrow(dat[q<0.05,]) / nrow(dat)*100
#### Compute % of significant phenotypes ####
length(unique(dat[q<0.05,]$hpo_id)) / length(unique(dat$hpo_id))*100
#### Compute % of significant cell types ####
length(unique(dat[q<0.05,]$CellType)) / length(unique(dat$CellType))*100
