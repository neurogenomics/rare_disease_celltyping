#' Subset results and add HPO term id
#'
#' This is called by the subset_phenos function to add the HPO term id column and
#' subset the EWCE results specific to a particular cell.
#'
#' @param cell The cell type of interest <string>
#' @param results EWCE Results <data.frame>
#' @param q_threshold The q value threshold of significance
#' @param fold_threshold The fold change threshold
#' @param phenotype_to_genes The HPO Ids with associated gene lists downloaded from HPO website
#' @param hpo The HPO Ontology data object
#' @returns A data frame of the selected subset of RD EWCE results with HPO ID column added.
#' @export
get_cell_ontology = function(cell, results, q_threshold, fold_threshold, phenotype_to_genes,hpo){
  signif_cell_data = results[results$CellType == cell & results$q <= q_threshold & results$fold_change >= fold_threshold,]
  signif_cell_data = add_hpo_termid_col(signif_cell_data, phenotype_to_genes , hpo)
  return (signif_cell_data)
}
