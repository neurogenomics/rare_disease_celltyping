# GENERATE WRITEUP FUNCTION

#' Render writeup
#'
#' This takes Jais .Rmd file and passes in the desired parameters for the writeup.
#' It was intended so that the full analysis could be done from generating the CTD
#' to creating the writuep. For my project report I ended up just creating a new
#' writeup. It is parameterised in a similar way so i could create a similar function,
#' however I now think its unlikely to want to generate different versions of the
#' writeup using different subsets of the results without having to rewrite significant
#' chunks of text. So maybe its not really worth it.
#'
#' The params for this function mostly are variables passed to the writeup.Rmd file.
#' changing them would vary the parts of the data shown in plots etc.
#' @returns NULL but saves the final writeup in a desired path
#' @export
render_writeup = function (markdown_file, writeup_name,
                           title, author, date,results, phenotype2genes,
                           keyword_1,keyword_1_search_terms,
                           keyword_2,keyword_2_search_terms,
                           keyword_3,keyword_3_search_terms,
                           parent_term_example = parent_term_example,
                           parent_term_example2 = parent_term_example2,
                           expected_terms = expected_terms,
                           cell_of_interest = cell_of_interest,
                           related_patterns = related_patterns,
                           enriched_cell = enriched_cell,
                           enriched_cell_related_cells = enriched_cell_related_cells,
                           not_enriched_cell = not_enriched_cell,
                           not_enriched_cell_related_cells = not_enriched_cell_related_cells,
                           working_directory) {


  if (!file.exists(paste0(working_directory,"data/hp.obo"))){
    download.file("http://purl.obolibrary.org/obo/hp.obo",
                  paste0(working_directory,"data/hp.obo"), mode = "wb")
  }
  rmarkdown::render(markdown_file, params = list(
    title = title,
    author = author,
    date = date,
    results_data =paste0(working_directory,results),
    phenotype_to_genes_txt_file = paste0(working_directory,phenotype2genes),
    keyword_1 = keyword_1,
    keyword_1_search_terms = keyword_1_search_terms,
    keyword_2 = keyword_2,
    keyword_2_search_terms = keyword_2_search_terms,
    keyword_3 = keyword_3,
    keyword_3_search_terms = keyword_3_search_terms,
    parent_term_example = parent_term_example,
    parent_term_example2 = parent_term_example2,
    expected_terms = expected_terms,
    cell_of_interest = cell_of_interest,
    related_patterns = related_patterns,
    enriched_cell = enriched_cell,
    enriched_cell_related_cells = enriched_cell_related_cells,
    not_enriched_cell = not_enriched_cell,
    not_enriched_cell_related_cells = not_enriched_cell_related_cells),
  output_file = paste0(working_directory,writeup_name))#   <-, ".pdf" ?
}
