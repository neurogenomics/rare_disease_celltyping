#' Get component relative ontology level of all terms within a subset of HPO
#'
#' This calls the \code{find_parent} function on all phenotypes in the subset of
#' the HPO to be plotted. The subest chosen when creating the phenoAdj from the main
#' adjacency matrix of all phenotypes. So, the phenotypes to be plotted can be
#' found in the row and column names of phenoAdj.
#'
#' @param phenoAdj A adjacency matrix of phenotypes where 1 represents i is parent of j
#' and 0 represents that i is not a parent of j. It is a subset of the main phenotype adjacency matrix
#' @param hpo The HPO ontology data object
#' @param reverse A boolean, if TRUE it will reverse the ontology level numbers so that
#' the parent terms are larger than the child terms.
#' @returns A named vector of relative ontology level, where names are HPO Ids and
#' value is relative ontology level.
#' @export
get_heirarchy <- function (phenoAdj,hpo,reverse=TRUE) {
  heirarchy = c()
  for (p in rownames(phenoAdj)) {
    heirarchy[p] = find_parent(p,phenoAdj,hpo)
  }
  if (reverse) {
    heirarchy = max(heirarchy) - heirarchy
  }
  return (heirarchy)
}
