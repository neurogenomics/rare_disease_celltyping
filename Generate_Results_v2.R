# Results Gen V2
# Generates EWCE results on HPO gene lists in parallel

# USER INPUT ###################################################################
phenotype_to_genes = read.delim("data/phenotype_to_genes.txt",
                                skip = 1,
                                header=FALSE)
colnames(phenotype_to_genes) = c("ID", "Phenotype", "EntrezID", "Gene",
                                 "Additional", "Source", "LinkID")

Phenotypes <- unique(phenotype_to_genes$Phenotype)
background <- unique(phenotype_to_genes$Gene)
results_dir <- "results_test/"
overwrite_past_analysis = FALSE
reps = 10
annotLevel = 1
genelistSpecies = "human"
sctSpecies = "human"
cores = 1
MergeResults <- TRUE

# Dependencies #################################################################
if (!require(HPOEWCE)) {
  #install.packages("rd_package/RareDiseaseEWCE", repos=NULL, type="source")
  devtools::install_github("ovrhuman/HPOEWCE")
  library(HPOEWCE)
}
if (!require(MultiEWCE)) {
  #install.packages("rd_package/RareDiseaseEWCE", repos=NULL, type="source")
  devtools::install_github("ovrhuman/MultiEWCE")
  library(MultiEWCE)
}
if (!require(HPOExplorer)) {
  #install.packages("rd_package/RareDiseaseEWCE", repos=NULL, type="source")
  devtools::install_github("ovrhuman/HPOExplorer")
  library(HPOExplorer)
}

# Run ##########################################################################
# remove gene lists with not enought genes
ctd_genes = rownames(ctd[[1]]$specificity_quantiles)
validLists = c()
for (p in Phenotypes) {
  if (sum(unique(HPOExplorer::get_gene_list(p,phenotype_to_genes)) %in% ctd_genes) >= 4) {
    validLists = append(validLists, p)
  }
}
Phenotypes = validLists
rm(validLists)


# Create results directory and remove already analysied phenos from input
if (!file.exists(results_dir)) {
  dir.create(results_dir) # <- create results dir
}

if (!overwrite_past_analysis) { # <- remove finished phenos
  Phenotypes <- MultiEWCE::get_unfinished_list_names(Phenotypes,results_dir)
}


# Create list of gene lists, indexed by phenotype
GeneLists <- list()
for (p in Phenotypes) {
  GeneLists[[p]] <- HPOExplorer::get_gene_list(p,phenotype_to_genes)
}


# Run analysis

MultiEWCE::ewce_para(list_names = Phenotypes,
          gene_lists = GeneLists,
          results_directory = results_dir,
          ctd_file = ctd,
          background_genes = background,
          bootstrap_reps = reps,
          annotation_Level= annotLevel,
          genes_Species = genelistSpecies,
          ctd_Species = sctSpecies,
          cores = cores)


# Merge results ################################################################
if (MergeResults) {
  results_final <- MultiEWCE::merge_results(results_dir = results_dir,
                                            description_col = "phenotype")
}
