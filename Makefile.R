# MAKEFILE #

source("source/CTD_Gen.R")
source("source/Results_Gen.R")
source("source/Writeup_Gen.R")

# Inputs
title = "Rare Disease EWCE Tabula Muris"
date = "06/05/2021"
author = "author"
# Choose search terms for section: "Searching for enrichments in phenotypes.."
keyword_1 = "heart"
keyword_1_search_terms = c("heart","atria","aorta")
keyword_2 = "circulating"
keyword_2_search_terms = c("circulating","circulate","circulation")
keyword_3 = "glucose"
keyword_3_search_terms = c("gluc","glucose")
generate_ctd = FALSE
generate_results = FALSE


# File locations
ctd_file_name = "data/CTD_tm_l1l2_nz.rda" # <- will be created if generate_ctd = TRUE
ewce_output_all_phenotypes = "results/EWCE_output_all_phenotypes"
merged_results_rda_filename = "results/all_results_merged_fixednames.rda" # <- "all_results_merged_fixednames.rda" (use these results to test the final knit)
unmerged_results_rda_filename = "results/resultsTest_unmerged.rda"
phenotype_to_genes_txt_file = "data/phenotype_to_genes.txt" # don't include Makefile/ before filename
r_markdown_file = "Writeup.Rmd"
working_directory = paste0(getwd() ,"/")

# Generate CTD
if (generate_ctd){
  gen_ctd(output_ctd_file_name = ctd_file_name)
}

# Generate Results
if(generate_results) {
  gen_results(CTD_file = ctd_file_name,
              mc.cores = 1,
              reps = 100,
              output_path = ewce_output_all_phenotypes,
              phenotype_to_genes_txt_file = paste0(working_directory,phenotype_to_genes_txt_file),
              output_merged_rda_filename = paste0(working_directory,merged_results_rda_filename),
              output_unmerged_rda_filename = paste0(working_directory,unmerged_results_rda_filename))
}

# Render Writeup
render_writeup(markdown_file = r_markdown_file,
               writeup_name =  paste0(str_replace_all(title, " ","_"),"_",str_replace_all(date, "/","")),
               title = title, author = author, date = date,
               results = merged_results_rda_filename,
               phenotype2genes = phenotype_to_genes_txt_file,
               keyword_1 = keyword_1,
               keyword_1_search_terms = keyword_1_search_terms,
               keyword_2 = keyword_2,
               keyword_2_search_terms = keyword_2_search_terms,
               keyword_3 = keyword_3,
               keyword_3_search_terms = keyword_3_search_terms,
               working_directory = working_directory)






## testing
#load("results/all_results_merged_fixednames.rda")
#genedata = read.delim(phenotype_to_genes_txt_file, skip =1) #HPO annotation
#colnames(genedata) = c("ID", "Phenotype", "EntrezID", "Gene", "Additional", "Source", "LinkID")

