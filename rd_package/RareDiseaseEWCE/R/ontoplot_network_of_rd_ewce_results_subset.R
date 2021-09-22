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
#' @param heatmapped_value "q", "fold change" or "p" <string>
#' @param q_threshold The q value threshold of significance
#' @param fold_threshold The fold change threshold
#' @param phenotype_to_genes The HPO Ids with associated gene lists downloaded from HPO website
#' @param hpo The HPO Ontology data object
#' @returns A ontologyPlot plot of the network of phenotypes in a subset of RD EWCE Results
#' @export
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
  heatpallette = grDevices::heat.colors(length(values))
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
  return (ontologyPlot::onto_plot(hpo,terms=names(values), fillcolor = heat, shape="rect"))
}
