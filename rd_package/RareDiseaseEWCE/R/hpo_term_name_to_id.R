#' Get HPO Id from phenotype name.
#'
#' I have done this more efficiently elsewhere using the hpo data object.
#' May be worth replacing, or just add the HPO Id to all datapoints in the results permanently.
#' Alternative method: \code{hpo$id[match(term_name, hpo$name)]}
#' This function is called by the add_hpo_termid_col function, which is called by the get_cell_ontology
#' function when selecting a subset of the data and then adding a HPO id column.
#' @param phenotype Phenotype name from the HPO <string>
#' @param phenotype_to_genes The hpo terms with gene list annotations data frame from hpo website
#'
#' @returns The HPO Id <string>
#'
#' @export
get_hpo_termID = function(phenotype, phenotype_to_genes){
  return(phenotype_to_genes$ID[phenotype_to_genes$Phenotype == phenotype][1])
}
#' Get HPO Id from phenotype name.
#'
#' I have done this more efficiently elsewhere using the hpo data object.
#' May be worth replacing, or just add the HPO Id to all datapoints in the results permanently.
#' Alternative method: \code{hpo$id[match(term_name, hpo$name)]}
#' This function is called by the add_hpo_termid_col function, which is called by the get_cell_ontology
#' function when selecting a subset of the data and then adding a HPO id column.
#' @param phenotype Phenotype name from the HPO <string>
#' @param phenotype_to_genes The hpo terms with gene list annotations data frame from hpo website
#'
#' @returns The HPO Id <string>
#'
#' @export
get_hpo_termID = function(phenotype, phenotype_to_genes){
  return(phenotype_to_genes$ID[phenotype_to_genes$Phenotype == phenotype][1])
}
