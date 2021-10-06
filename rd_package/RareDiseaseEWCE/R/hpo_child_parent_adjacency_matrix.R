#' Create adjacency matrix of HPO child-parent relationships
#'
#' This is needed for plots created using ggnetwork as it can coerce an
#' adjacency matrix into a directed graph object, and assigns
#' each datapoint an x,y coordinate. Maybe shouldn't use a for loop for this.
#' It also may be possible to use a hash table for ggnetwork, which may be more
#' efficent for the web app
#'
#' @param pheno_ids a character vector of HPO Ids
#' @param hpo ontology object (available in ontologyIndex package)
#' @examples
#' \dontrun{
#' adjacency_matrix(c("HP:000001","HP:000002"),hpo)
#' }
#' @returns adjacency matrix with HPO Ids for col and row names.
#' If adjacency[i,j] == 1 then phenotype[i] is a parent of phenotype[j]
#' @export

adjacency_matrix <- function(pheno_ids, hpo) {
  HPO_id = unique(pheno_ids)
  size = length(HPO_id)
  adjacency = data.frame(matrix(nrow = size, ncol = size))
  rownames(adjacency) = HPO_id
  colnames(adjacency) = HPO_id
  adjacency[is.na(adjacency)] = 0
  for (id in HPO_id) {
    children = hpo$children[id][[1]]
    adjacency[id, children] = 1
  }
  return(adjacency[HPO_id,HPO_id])
}
