#' Create interactive network plot start to finish
#'
#' This puts all the functions together from gettig the subest of results to creating
#' the final interactive plot.
#'
#' @param phenotype_to_genes The phenotype gene lists taken from the HPO
#' @param all_results_merged The RD EWCE Results
#' @param hpo The HPO ontology data object
#' @param disease_descriptions The dataframe of all disease descriptions in the HPO
#' @param cell_type The cell type of interest to be plotted
#' @param q_threshold The q value threshold for the subset of results to be plotted
#' @param fold_threshold The minimum fold change in specific expression for the subest of results to be plotted
#'
#' @returns A interactive network plot of the selected subset of results from RD EWCE analysis
#' @export
ggnetwork_plot_full = function(phenotype_to_genes=phenotype_to_genes,
                               all_results_merged=all_results_merged, hpo,
                               disease_descriptions,cell_type = "Neurons",
                               q_threshold =0.0005, fold_threshold = 1){
  phenos = subset_phenos(phenotype_to_genes, all_results_merged, hpo,cell_type =cell_type, q_threshold =q_threshold, fold_threshold = fold_threshold)
  adjacency = adjacency_matrix(unique(phenos$HPO_term_Id), hpo)
  phenoNet = make_network_object(phenos,adjacency,hpo)
  network_plot = ggnetwork_plot(phenoNet, phenos,disease_descriptions)
  return(network_plot)
}
