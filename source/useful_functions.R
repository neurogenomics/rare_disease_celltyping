# Useful functions #############################################################


#' Load RDA file and assign to specific variable
#'
#' Can be useful, but sort of pontles now as we have decided to only use \code{.rds}
#'
#' @param file the path to a .rda file <string>
#' @example RDA_assign_load("data/results.rda")
#' @returns the data contained in the .rda file
#' @export
RDA_assign_load <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  return(tmp[[ls(tmp)[1]]])
}


