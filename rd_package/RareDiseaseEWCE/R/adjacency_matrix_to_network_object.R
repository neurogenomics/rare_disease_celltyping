#' Make network object
#'
#' This uses the network package to coerce the adjacency matrix into a
#' network object. It also adds the fold change, label, and relative ontology level
#' parameters to each node in the network.
#'
#' @param phenos The subset of the results to be plotted
#' @param adjacency The adjacency matrix of all HPO terms
#' @param hpo The HPO ontology data object
#' @examples
#' \dontrun{make_network_object(phenos,adjacency,hpo)}
#' @returns A ggnetowrk graph/ network object of a subset of the RD EWCE results.
#' @export
make_network_object = function(phenos, adjacency, hpo) {
  ValidTerms = phenos$HPO_term_Id
  phenoAdj = adjacency[ValidTerms,ValidTerms]
  # MAKE NETWORK OBJECT
  phenoNet = network(phenoAdj, directed = TRUE)
  # To add a another value to the nodes do this
  phenoNet %v% "fold" = as.numeric(phenos$fold_change)
  phenoNet %v% "label" = as.character(phenos$list)
  phenoNet %v% "heirarchy" = as.numeric(get_heirarchy(phenoAdj,hpo,reverse=TRUE)+1)
  # (%v% is for setting vertex attributes, %e% is for edge, %n% is for network attributes)
  return(phenoNet)
}
