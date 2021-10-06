#' Get absolute ontology level
#'
#' This gets the absolute ontology level of a term (without consideration for the
#' particular subset of the data you are looking at, as in the find_parent function).
#'
#' @param hpo The HPO ontology data object
#' @param term_id HPO term ID <string>
#' @example \dontrun{get_ont_level(hpo,"HP:0000003")}
#' @return returns the ontology level <numeric>
#' @export
get_ont_level = function(hpo,term_id) {
  children = unique(setdiff(unlist(hpo$children[term_id]), term_id))
  if (length(children) == 0) {
    return(0)
  } else {
    return(1 + get_ont_level(hpo,children)) #<- recursion..
  }
}
