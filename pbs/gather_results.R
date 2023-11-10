#### Set paths ####
ephemeral <- "/rds/general/project/neurogenomics-lab/ephemeral/rare_disease/"
# ephemeral <- "/rds/general/project/neurogenomics-lab/ephemeral/rare_disease_min_genes1/"
ephemeral_logs <- file.path(dirname(ephemeral),"rare_disease.pbs_output")

#### Find all results files ####
f <- list.files(ephemeral, 
                full.names = TRUE)
f2 <- sapply(f, function(x){list.files(x,pattern = ".rds", full.names = TRUE)})
length(unlist(f2)) 
DAT <- (
  lapply(unlist(f2), readRDS) |> 
    data.table::rbindlist(use.names = TRUE, idcol = "batch")
)[,batch:=basename(batch)]
saveRDS(DAT,here::here("pbs/rare_disease_min_genes4.rds"))
extra_dat <- readRDS(here::here("pbs/rare_disease_min_genes1.rds")) 
length(unique(extra_dat[q<0.05]$hpo_id))/length(unique(extra_dat$hpo_id))
length(unique(DAT[q<0.05]$hpo_id))/length(unique(DAT$hpo_id))

dat <- rbind(
  DAT,
  extra_dat,
  fill=TRUE
) |>
  merge(unique(gene_data[,c("hpo_id","n_gene")]), on="hpo_id")
dat <- dat[n_gene>0,]
# table(dat$n_gene)

## Apply MTC
dat[,q:=stats::p.adjust(p,method = "fdr")]
## Report n significant results
message(dplyr::n_distinct(dat[q<0.05]$hpo_id)," significant phenotypes.")

## Add metadata
dat[,hpo_name:=HPOExplorer::harmonise_phenotypes(hpo_id)]
dat <- HPOExplorer::add_ancestor(dat)
## Add gene number
dat <- dat[unique(gene_data[,c("hpo_id","n_gene")]),on="hpo_id"]
## Get top celltype association per phenotype
top_dat <- dplyr::group_by(dat[q<.05],hpo_name)|>
  dplyr::slice_max(fold_change,n=1) |>
  data.table::data.table()
View(top_dat)
saveRDS(dat, here::here("pbs",paste0(basename(ephemeral),".rds")))

#### Find missing IDs ####
ctd <- MultiEWCE::load_example_ctd(file = "ctd_DescartesHuman.rds")
annotLevel <- 2
gene_data <- HPOExplorer::load_phenotype_to_genes()
gene_data <- HPOExplorer::add_ont_lvl(gene_data)
gene_data <- gene_data[gene_symbol %in% rownames(ctd[[annotLevel]]$mean_exp)]
gene_data[,n_gene:=(length(unique(gene_symbol))),by="hpo_id"]
gene_data <- gene_data[n_gene>=4,][,batch_id:=.GRP, by="hpo_id"]
missing_dat <- gene_data[!hpo_id %in% unique(dat$hpo_id)]
## Summarise missing data
length(unique(missing_dat$hpo_id)) 
head(sort(unique(missing_dat$hpo_name)),10) 
summary(unique(missing_dat[,c("hpo_id","n_gene")])$n_gene)
summary(unique(missing_dat[,c("hpo_id","ontLvl")])$ontLvl)
# length(unique(gene_data[n_gene>=1428]$hpo_id))
saveRDS(missing_dat,here::here("pbs/missing_dat.rds"))

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
