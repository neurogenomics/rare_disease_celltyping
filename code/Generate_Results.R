# Results gen v3 (using MultiEWCE package)

## Dependencies ################################################################
# Install MultiEWCE package with devtools::install_github("ovrhuman/MultiEWCE")
# and HPOExplorer with devtools::install_github("ovrhuman/HPOExplorer")

## Load data ###################################################################
phenotype_to_genes <- HPOExplorer::load_phenotype_to_genes("data/phenotype_to_genes.txt")
ctd <- readRDS("data/CTD_Descartes_withplot.rds")

## Run analysis ################################################################
all_results <- MultiEWCE::gen_results(ctd,
                                      phenotype_to_genes,
                                      list_names = unique(phenotype_to_genes$Phenotype),
                                      background_genes = unique(phenotype_to_genes$Gene),
                                      list_name_column = "Phenotype",
                                      gene_column = "Gene",
                                      results_dir = "results_test_191021",
                                      overwrite_past_analysis = FALSE,
                                      reps = 100000,
                                      annotLevel = 1,
                                      genelistSpecies = "human",
                                      cores = 8,
                                      MergeResults = TRUE)

saveRDS(all_results, paste0("HPO_EWCE_Results ",stringr::str_replace_all(Sys.time(),":",""),".rds"))
