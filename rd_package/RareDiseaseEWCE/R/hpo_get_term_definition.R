#' Get HPO term definition
#'
#' This gets the disease description from a data frame of disease descriptions.
#' The rows names of the data frame are the HPO ID and the description column is
#' called "description". It also adds new lines to the description so that the
#' hover box in the web app does not get too wide. This is done by calling the
#' \code{newlines_to_definition} function.
#'
#' @param ontologyId The HPO Id of the term (string)
#' @param disease_descriptions A data frame of disease descriptions corresponding to each HPO Id
#'
#' @return The disease description with new lines added.
#'
#' @examples
#' \dontrun{hpo_get_term_definition("HP:123456", disease_descriptions)}
#'
#' @export
hpo_get_term_definition <- function(ontologyId, disease_descriptions) {
  definition = disease_descriptions[ontologyId,"description"]
  definition = newlines_to_definition(definition)
  return (definition)
}
