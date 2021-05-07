# GENERATE WRITEUP FUNCTION

render_writeup = function (markdown_file, writeup_name,
                           title, author, date,results, phenotype2genes,
                           keyword_1,keyword_1_search_terms,
                           keyword_2,keyword_2_search_terms,keyword_3,keyword_3_search_terms,
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
    keyword_3_search_terms = keyword_3_search_terms),
  output_file = paste0(working_directory,writeup_name))#   <-, ".pdf" ?
}
