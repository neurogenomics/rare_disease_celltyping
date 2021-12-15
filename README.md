# rare_disease_ewce
  - Code for generating the CTD files is stored within CTD_code
  - The actual CTD files are in the lab dropbox
  - Around 31st October 2017, I wrote a script ([process_hpo_to_genelists.r] which read in [TermsToUse.csv], which contained a list of brain relevant traits, and extracted those gene sets. At some point, I then modified this file to create [process_hpo_to_genelists-wholebody.r]. This no longer limited itself to brain relevant traits. It saved all the gene lists to the [GeneLists-WholeBody] folder which is also included in the repo.

# The raw data on genes->phenotypes
You can find the link from here... it changes every once in a while...
https://hpo.jax.org/app/download/annotation

Currently the correct file is: http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt

# Paper 
* [Paper Rmd](https://github.com/neurogenomics/rare_disease_ewce/blob/master/paper.rmd): This is a fork of the masters project writeup.
* [Masters Rmd](https://github.com/neurogenomics/rare_disease_ewce/blob/master/report_BobbyGordonSmith_masters_RD_project_writeup.Rmd): This is the original masters writeup.
* [Masters Knitted](https://github.com/neurogenomics/rare_disease_ewce/blob/master/report%20(1).pdf): Knitted versoin of the masters writeup.
* [Grouped figures Rmd](https://github.com/neurogenomics/rare_disease_ewce/blob/master/groupedPlots.rmd): These are the plots for the paper grouped into larger figures. This is still under development and there is a markdown [here](https://github.com/neurogenomics/rare_disease_ewce/blob/master/figure_group_planning.md) for planning the groupings. 
* [Grouped figures knitted](https://github.com/neurogenomics/rare_disease_ewce/blob/master/groupedPlots.pdf): This is the leatest knitted version of grouped figures.
