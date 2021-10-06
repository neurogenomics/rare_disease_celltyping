#' Add new lines to disease description
#'
#' Adds new lines to the description so that hover boxes dont get too wide.
#'
#' @param definition A disease description string
#' @param line_length A integer representing the desired words per line.
#'
#' @returns The disease description with newline symbols added every nth word.
#'
#' @examples
#' \dontrun{newlines_to_definition(disease_description, 10)}
#' @export
newlines_to_definition <- function(definition, line_length = 10) {
  definition = strsplit(definition, split = " ")[[1]]
  if (length(definition) > line_length) {
    remainder = length(definition) %% line_length
    n_new_lines = floor((length(definition)/line_length))
    new_line_index = seq(line_length,(n_new_lines*line_length),line_length)
    definition[new_line_index] = paste0("\n", definition[new_line_index])
  }
  definition = paste(definition,collapse = " ")
  return(definition)
}
