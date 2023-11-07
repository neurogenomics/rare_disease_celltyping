#### Set paths ####
ephemeral <- "/rds/general/project/neurogenomics-lab/ephemeral/rare_disease1/"
ephemeral_logs <- file.path(dirname(ephemeral),"rare_disease.pbs_output")

#### Find all results files ####
f <- list.files(ephemeral, 
                full.names = TRUE)
f2 <- sapply(f, function(x){list.files(x,pattern = ".rds", full.names = TRUE)})
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
head(sort(unique(missing_dat$hpo_name)),10) 

#### Inspect PBS logs ####
## Peak at the most recent PBS output logs 
## to figure out why some of the batches/subjobs failed.
## Get the most recent jobID
## Would normally use `qstat -xf` to get this info, 
## but our HPC is not configured to enable this.
logs <- system(paste("ls",ephemeral_logs,"-Artlsh | tail -10"), intern = TRUE)
jobID <- gsub("[a-z]","",rev(strsplit(tail(logs,1),"\\.")[[1]])[2])
logs_df <- data.table::data.table(batch_id=unique(missing_dat$batch_id))[,
  path:=file.path(ephemeral_logs,
                 paste("rare_disease_celltyping.pbs",
                       paste0("e",jobID),batch_id,sep=".")
  )
][file.exists(path),]
## Print the logs for each batch
out <- lapply(logs_df$path,function(x){
  message("\n~~~~~~~~~~",basename(x),"~~~~~~~~~~")
  readLines(x)|>cat(sep = "\n")
})

#### Compute % of significant tests ####
nrow(dat[q<0.05,]) / nrow(dat)*100
#### Compute % of significant phenotypes ####
length(unique(dat[q<0.05,]$hpo_id)) / length(unique(dat$hpo_id))*100
#### Compute % of significant cell types ####
length(unique(dat[q<0.05,]$CellType)) / length(unique(dat$CellType))*100
