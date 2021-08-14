
# Plot ontology associated with cell type - issue 10 ###########################


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


get_cell_ontology = function(cell, results, q_threshold, fold_threshold, phenotype_to_genes,hpo){
  phenotype_to_genes = read.delim("data/phenotype_to_genes.txt", skip = 1, header=FALSE)
  colnames(phenotype_to_genes) = c("ID", "Phenotype", "EntrezID", "Gene",
                                   "Additional", "Source", "LinkID")
  signif_cell_data = results[results$CellType == cell & results$q <= q_threshold & results$fold_change >= fold_threshold,]
  signif_cell_data = add_hpo_termid_col(signif_cell_data, phenotype_to_genes , hpo)
  return (signif_cell_data)
}



get_hpo_termID = function(phenotype, phenotype_to_genes){
  return(phenotype_to_genes$ID[phenotype_to_genes$Phenotype == phenotype][1])
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



results=all_results_merged
cell = "Bladder cells"
heatmapped_value = "q"
q_threshold = 0.05
fold_threshold =1




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



