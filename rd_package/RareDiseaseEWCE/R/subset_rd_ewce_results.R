#' Subset RD EWCE Results and add HPO Id column
#'
#' This subsets  the Rare disease EWCE results by cell type, q threshold and fold change.
#' It calls other functions which add HPO Id column and it removes NA or invalid
#' HPO Ids. These checks are necessary for the interactive plots on the web app
#' to avoid errors, but not needed if manually plotting.
#'
#' @param phenotype_to_genes The list of HPO terms with their assocaited gene lists taken from HPO website
#' @param all_results_merged The dataframe of RD EWCE Results
#' @param hpo The HPO ontology data object
#' @param cell_type A string representing the cell type of interest.
#' @param q_threshold The q threshold. The subset of results will have a q lower than this
#' @param fold_threshold The fold change threshold. The subest of results will have a fold change greater than this.
#'
#' @returns A data frame of results taken from the main data frame of results
#' @export
subset_phenos = function(phenotype_to_genes, all_results_merged, hpo,cell_type = "Neurons", q_threshold =0.0005, fold_threshold = 1) {
  phenos = get_cell_ontology(cell_type,all_results_merged,q_threshold = q_threshold, fold_threshold = fold_threshold, phenotype_to_genes, hpo)
  phenos = phenos[!is.na(phenos$HPO_term_Id) & phenos$HPO_term_valid,]
  return (phenos)
}
