# MAKEFILE #

source("Makefile/source/CTD_Gen.R")
source("Makefile/source/Results_Gen.R")
source("Makefile/source/Writeup_Gen.R")

# User Inputs
output_ctd_file_name = "Makefile/ctd_tm_l1l2_nz.rda"
merged_results_rda_filename = "resultsTest_merged.rda"
unmerged_results_rda_filename = "resultsTest_unmerged.rda"
phenotype_to_genes_txt_file = "phenotype_to_genes.txt" # don't include Makefile/ before filename
r_markdown_file = "Makefile/Writeup.Rmd"
makefile_directory = "Makefile"
title = "Rare Disease EWCE"
date = "26/04/2021"
author = "author"

# Create working directory of markdownfile
working_directory = paste0(getwd() ,"/", makefile_directory,"/")

# Generate CTD
gen_ctd(output_ctd_file_name = output_ctd_file_name)

# Generate Results (note ive redcuded the dataset in the gen_results function to make it run fast)
gen_results(CTD_file = output_ctd_file_name,
            mc.cores = 1,
            reps = 100,
            output_path = "Makefile/Output",
            phenotype_to_genes_txt_file = paste0(working_directory,phenotype_to_genes_txt_file),
            output_merged_rda_filename = paste0(working_directory,merged_results_rda_filename),
            output_unmerged_rda_filename = paste0(working_directory,unmerged_results_rda_filename))

# Render Writeup
render_writeup(r_markdown_file, writeup_name = "Test_Writeup",
               title = title, author = author, date = date,
               results = merged_results_rda_filename,
               phenotype2genes = phenotype_to_genes_txt_file,
               working_directory = working_directory)

