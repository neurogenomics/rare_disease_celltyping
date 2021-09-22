#' Dataframe of number of significant RD EWCE Phenotypes per cell type
#'
#' This creates the dataframe of n signif enrichments for each cell type from the
#' RD EWCE results. It is used in the plot_phenos_per_cell function. It also orders the
#' results by number of phenotypes for better plot.
#'
#' @param all_results_merged Data frame of RD EWCE results
#' @param cell_mappings tissue-cell mappings for colouring results by tissue type in the plot
#' @param fold Fold threshold for deciding which results count as significant
#' @param q_val The q value threshold for significant results
#' @param dataset "Descartes" or "other" string. Probably not needed anymore
#'
#' @returns A data frame of enrichment counts per cell
#'
#' @export
dataframe_phenos_per_cell = function (all_results_merged, cell_mappings, fold = 1, q_val = 0.005, dataset = "Descartes") {
  phenos_per_cell = data.frame()

  # count n phenotypes per cell
  for (c in unique(all_results_merged$CellType)){
    cur_cell = c
    n_phenos = length(all_results_merged[all_results_merged$CellType==c & all_results_merged$q < q_val & all_results_merged$fold_change > fold, "list"])
    phenos_per_cell = rbind(phenos_per_cell, data.frame("Cell"=cur_cell,"n_phenos"=n_phenos))
  }
  n_pheno_text = phenos_per_cell
  phenos_per_cell$Tissue = rep(NA, length(phenos_per_cell$Cell))

  phenos_per_cell$Cell = stats::reorder(phenos_per_cell$Cell, phenos_per_cell$n_phenos)

  return(phenos_per_cell)
}
