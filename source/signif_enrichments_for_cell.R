
# Plot ontology associated with cell type - issue 10 ###########################
# This uses ontologyX packages rather than ggnetwork. The figures are more printable but not interactive

library(ontologyIndex)
library(ontologyPlot)
data(hpo)

# Try this data instead? -
# hpo = get_OBO("data/hp.obo", propagate_relationships = "is_a", extract_tags = "minimal")

# (not needed?)
phenotype_to_genes = read.delim("data/phenotype_to_genes.txt", skip = 1, header=FALSE)
colnames(phenotype_to_genes) = c("ID", "Phenotype", "EntrezID", "Gene",
                       "Additional", "Source", "LinkID")

# load ewce results (not needed?)
load("data/Descartes_All_Results_extras.rda")




## Functions ####################################################################

#' Subset RD EWCE results data by cell type, fold change and q value
#'
#' I have written similar functions in other scripts in the source, but I think
#' this is the one I refernce most often. May be worth checking through for redundancies.
#' directory. This also adds a HPO term id column for the subset though.
#' @param cell The cell type of interest <string>
#' @param q_threshold The q value threshold of significance
#' @param fold_threshold The fold change threshold
#' @param phenotype_to_genes The HPO Ids with associated gene lists downloaded from HPO website
#' @param hpo The HPO Ontology data object
#' @returns A data frame of the selected subset of RD EWCE results with HPO ID column added.
#' @export
get_cell_ontology = function(cell, results, q_threshold, fold_threshold, phenotype_to_genes,hpo){
  phenotype_to_genes = read.delim("data/phenotype_to_genes.txt", skip = 1, header=FALSE)
  colnames(phenotype_to_genes) = c("ID", "Phenotype", "EntrezID", "Gene",
                                   "Additional", "Source", "LinkID")
  signif_cell_data = results[results$CellType == cell & results$q <= q_threshold & results$fold_change >= fold_threshold,]
  signif_cell_data = add_hpo_termid_col(signif_cell_data, phenotype_to_genes , hpo)
  return (signif_cell_data)
}


#' Get HPO Id from phenotype name.
#'
#' I have done this more efficiently elsewhere using the hpo data object.
#' May be worth replacing, or just add the HPO Id to all datapoints in the results permanently.
#' Alternative method: \code{hpo$id[match(term_name, hpo$name)]}
#' This function is called by the add_hpo_termid_col function, which is called by the get_cell_ontology
#' function when selecting a subset of the data and then adding a HPO id column.
#' @param phenotype Phenotype name from the HPO <string>
#' @param phenotype_to_genes The hpo terms with gene list annotations data frame from hpo website
#'
#' @returns The HPO Id <string>
#'
#' @export
get_hpo_termID = function(phenotype, phenotype_to_genes){
  return(phenotype_to_genes$ID[phenotype_to_genes$Phenotype == phenotype][1])
}


#' Add HPO term Id column to dataframe.
#'
#' This adds the HPO term id column to the subest of ewce results data to be plotted
#' in the cell select app. It also checks if it is a valid HPO term id to pevent error and adds
#' a boolean column where TRUE if term is valid. If the HPO Id is not correct, it caused
#' an error in the ontologyPlot package
#' @param cells The dataframe of subset of RD EWCE results to be plotted in the cell select app.
#' @param phenotype_to_genes The hpo terms with gene list annotations data frame from hpo website
#' @param hpo The HPO ontology data object
#' @returns The subset of ewce result data frame with a HPO Id column added.
#' @export
add_hpo_termid_col = function(cells, phenotype_to_genes, hpo) {
  HPOtermID = c()
  ValidTerm = c()
  for (p in cells$list){
    termid = get_hpo_termID(p, phenotype_to_genes)
    ValidTerm = append(ValidTerm,(termid %in% hpo$id))
    HPOtermID = append(HPOtermID, termid)
  }
  cells$HPO_term_Id =HPOtermID
  cells$HPO_term_valid = ValidTerm
  return(cells)
}





# Some inputs for testing the function (probalby delete now .. )
results=all_results_merged
cell = "Bladder cells"
heatmapped_value = "q"
q_threshold = 0.05
fold_threshold =1



#' Plot graph of phenotypes in the RD EWCE Results subsetted by Cell type
#'
#' This subsets the results, using other functions from this script, by cell type,
#' q value, and fold change. It then plots a graph/ network of the phenotypes and
#' colors the nodes by their fold change or q value (see the \code{heatmapped_value} param).
#' The ontologyPlot package did not have prebuilt options to create the heatmap, so the
#' colours are assigned manually to each phenotype in this function. There must be a more
#' efficent way to do this. Atleas it may be good to replace the for loop with one of the
#' apply functions.
#' @param results The RD EWCE results dataframe
#' @param cell The cell type of interest <string>
#' @heatmapped_value "q", "fold change" or "p" <string>
#' @param q_threshold The q value threshold of significance
#' @param fold_threshold The fold change threshold
#' @param phenotype_to_genes The HPO Ids with associated gene lists downloaded from HPO website
#' @param hpo The HPO Ontology data object
#' @returns A ontologyPlot plot of the network of phenotypes in a subset of RD EWCE Results
one_cell_ontology_plot_heatmap = function(results, cell = "Bladder cells", heatmapped_value = "q",
                                          q_threshold, fold_threshold, phenotype_to_genes, hpo){
  #' heatmapped_value = "q", "fold_change", or "p". In other words, any continuous variable from the all_cell_ontology to be mapped on to the heatmap colors
  #' reverse_heatmap - reverse reccomeneded for q or p values, so the lowest "most significant" value is red
  cells = get_cell_ontology(cell, results, q_threshold, fold_threshold, phenotype_to_genes, hpo)
  cells = cells[cells$HPO_term_valid,]

  if (heatmapped_value == "q"){
    values = cells$q
  } else if (heatmapped_value == "p"){
    values = cells$p
  } else if (heatmapped_value == "fold change"){
    values = cells$fold_change
  } else {
    print("invalid heatmapped_value, enter 'p', 'q', or 'fold change'"); return(0)
    }
  names(values) = cells$HPO_term_Id
  values = sort(values)

  # creating the list of colors for the heatmap
  heatpallette = heat.colors(length(values))
  heat = c()
  prev = 0
  index = 1
  next_index = 1
  # this was necessary so equal values have the same color mapped to them
  for (v in values){
    if (v == prev) {
      heat = append(heat, heatpallette[index])
      next_index = next_index + 1
    } else if (v > prev){
      index = index + next_index
      next_index = 1
      prev = v
      heat = append(heat, heatpallette[index])
    }
  }
  if (heatmapped_value == "fold change") {
    heat = rev(heat)}
  # return (onto_plot(hpo, terms=names(values), fillcolor = heat, label = character_ID))
  # could use the above to add lables like ** for significance?
  return (onto_plot(hpo,terms=names(values), fillcolor = heat, shape="rect"))
}






## Testing ######################################################################

# plot = one_cell_ontology_plot_heatmap(all_results_merged, "Bladder cells",
#                                       heatmapped_value = "q", q_threshold = 0.05,
#                                       fold_threshold = 1,
#                                       phenotype_to_genes = phenotype_to_genes,
#                                       hpo = hpo)


#png("test.png")
#plot
#dev.off()



