
# Plot ontology associated with cell type - issue 10 ###########################


library(ontologyIndex)
library(ontologyPlot)
data(hpo)
phenotype_to_genes = read.delim("data/phenotype_to_genes.txt", skip = 1, header=FALSE)


# load ewce results
load("results/all_results_merged_fixednames.rda")




## Functions ####################################################################


get_cell_ontology = function(cell, results, q_threshold, fold_threshold, phenotype_to_genes,hpo){
  signif_cell_data = results[results$CellType == cell & results$q <= q_threshold & results$fold_change >= fold_threshold,]
  signif_cell_data = add_hpo_termid_col(signif_cell_data, phenotype_to_genes , hpo)
  return (signif_cell_data)
}

get_all_cell_ontology = function(results, q_threshold, fold_threshold, phenotype_to_genes, hpo){
  all_cell_ontologies = list()
  for (c in unique(results$CellType)) {
    all_cell_ontologies[[c]] =get_cell_ontology(c, results, q_threshold, fold_threshold, phenotype_to_genes, hpo)
  }
  return (all_cell_ontologies)
}

get_hpo_termID = function(phenotype, phenotype_to_genes){
  return(phenotype_to_genes$V1[phenotype_to_genes$V2 == phenotype][1])
}

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


cell_ontology_plot_heatmap = function(all_cell_ontologies , cell = "Bladder cells", heatmapped_value = "q", reverse_heatmap = TRUE){
  #' heatmapped_value = "q", "fold_change", or "p". In other words, any continuous variable from the all_cell_ontology to be mapped on to the heatmap colors
  #' reverse_heatmap - reverse reccomeneded for q or p values, so the lowest "most significant" value is red
  cells = all_cell_ontologies[[cell]];
  cells = cells[cells$HPO_term_valid,]

  if (heatmapped_value == "q"){values = cells$q}
  else if (heatmapped_value == "p"){values = cells$p}
  else if (heatmapped_value == "fold_change"){values = cells$fold_change}
  else {print("invalid heatmapped_value, enter p q or fold_change"); return(0)}
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

  if (reverse_heatmap) {
    heat = rev(heat)}
  # return (onto_plot(hpo, terms=names(values), fillcolor = heat, label = character_ID))
  # could use the above to add lables like ** for significance?
  return (onto_plot(hpo, terms=names(values), fillcolor = heat))
}




## Testing ######################################################################

all_cell_ontologies = get_all_cell_ontology(all_results_merged, 0.0005, -10, phenotype_to_genes, hpo)
# note: make a version where you dont need to genrate the ontologies for all cell types to do the plot


png("plot_test.png", width = 2600,height = 2600, units = "px")
cell_ontology_plot_heatmap(all_cell_ontologies, cell = "T cells", heatmapped_value = "fold_change", reverse_heatmap = FALSE)
dev.off()
