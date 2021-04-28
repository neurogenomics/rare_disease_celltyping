# GENERATE WRITEUP FUNCTION

render_writeup = function (markdown_file, writeup_name,
                           title, author, date,
                             results, phenotype2genes, working_directory){
  if (!file.exists(paste0(working_directory,"/hp.obo"))){
    download.file("http://purl.obolibrary.org/obo/hp.obo",
                  paste0(working_directory,"hp.obo"), mode = "wb")
  }
  rmarkdown::render(markdown_file, params = list(
    title = title,
    author = author,
    date = date,
    results_data = results,
    phenotype_to_genes_txt_file = phenotype2genes),
  output_file = paste0(writeup_name ))#, ".pdf"))
}
