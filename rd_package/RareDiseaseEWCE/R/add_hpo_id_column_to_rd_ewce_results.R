#' Add HPO term Id column to dataframe.
#'
#' This adds the HPO term id column to the subest of ewce results data to be plotted
#' in the cell select app. It also checks if it is a valid HPO term id to pevent error and adds
#' a boolean column where TRUE if term is valid. If the HPO Id is not correct, it caused
#' an error in the ontologyPlot package
#' @param cells The dataframe of subset of RD EWCE results to be plotted in the cell select app.
#' @param phenotype_to_genes The hpo terms with gene list annotations data frame from hpo website
#' @param hpo The HPO ontology data object
#' @examples
#' \dontrun{
#' add_hpo_termid_col(cells,phenotype_to_genes,hpo)}
#' @returns The subset of ewce result data frame with a HPO Id column added.
#' @export
add_hpo_termid_col = function(cells, phenotype_to_genes, hpo) {
  HPOtermID = c()
  ValidTerm = c()
  for (p in cells$list){
    termid = get_hpo_termID(p, phenotype_to_genes)
    ValidTerm = append(ValidTerm,(termid %in% hpo$id))
    HPOtermID = append(HPOtermID, termid)
  }
  cells$HPO_term_Id =HPOtermID
  cells$HPO_term_valid = ValidTerm
  return(cells)
}
