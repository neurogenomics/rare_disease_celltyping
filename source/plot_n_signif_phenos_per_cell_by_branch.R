# PLOT N SIGNIF PHENOS PER CELL BY HPO BRANCH


# Dependencies
library(ggplot2)
library(cowplot)
library(ontologyIndex)
data(hpo)
phenotype_to_genes = read.delim("data/phenotype_to_genes.txt", skip = 1, header=FALSE)
colnames(phenotype_to_genes) = c("ID", "Phenotype", "EntrezID", "Gene",
                                 "Additional", "Source", "LinkID")


load("data/Descartes_All_Results.rda")


#' Plot number of significant enrichments per cell within specific HPO branch
#'
#' This plots how many significant enrichments there were per cell from within
#' specific branches of the HPO. It also runs hypergeometric tests to see if the
#' number of significant enrichments for a particular cell, within a particular branch
#' is higher than expected by chance. This can tell you if a branch as a whole is
#' commonly associated with a particular cell type.
#'
#' This function could be split up into more smaller functions and possibly include
#' options for including the hypergeometric test.
#'
#' @param all_results_merged The RD EWCE results data frame
#' @plot_branches A character vector of the HPO branches to be plotted (using name not Id, maybe change this)
#' @q_threshold The maximum q value accepted as a significant result.
#' @returns A plot of n enrichments per cell type, faceted by branch. Astrics are
#' used to represent significance in the hypergeometric test.
#' @export

plot_n_signif_phenos_per_cell_by_branch <- function(all_results_merged,
                                                    plot_branches = c("Abnormality of the nervous system",
                                                                      "Abnormality of the cardiovascular system",
                                                                      "Abnormality of the immune system"),
                                                    q_threshold = 0.05) {

  # Get main hpo branches (child nodes of PHenotypic abnormality)

  hpo_branches = hpo$children["HP:0000118"]
  hpo_branches_names = c()
  for (b in hpo_branches){
    hpo_branches_names = c(hpo_branches_names, hpo$name[b])
  }
  main_branch_df = data.frame("hpo_id"=hpo_branches,"phenotype"=hpo_branches_names)

  descendants = list ()
  for (b in hpo_branches[[1]] ) {
    print(b)
    print(hpo$name[b])
    descendants[[hpo$name[b]]] = get_descendants(hpo,b)
  }
  # Label all results based on branch
  if (! "branch" %in% colnames(all_results_merged)) {
    phenotypes = unique(phenotype_to_genes$Phenotype)
    all_results_merged$branch = NA
    pheno_branches = c()
    for (p in phenotypes){
      print(p)
      index = match(p,hpo$name)
      if (!is.na(index)){
        id = hpo$id[[index]]
        for (i in seq(1,length(descendants))) {
          if (id %in% descendants[[i]]) {
            pheno_branches[p] = names(descendants[i])
          }
        }
      }
    }

    all_results_merged$branch = NA
    remaining = length(unique(all_results_merged$list))
    cat("Asigning branch to phenotypes\nN remaining:\n")
    for (p in unique(all_results_merged$list)) {
      cat("\r",remaining,"      ")
      remaining = remaining - 1
      all_results_merged$branch[all_results_merged$list == p] = pheno_branches[p]
    }
  }
  # plot it
  signif_results = all_results_merged[all_results_merged$q <= q_threshold & !is.na(all_results_merged$branch), ]
  n_signif_per_cell_by_branch = data.frame()
  for (b in unique(plot_branches)) {
    for (c in unique(all_results_merged$CellType)){
      n = length(signif_results[signif_results$CellType == c & signif_results$branch == b,]$list)
      n_signif_per_cell_by_branch = rbind(n_signif_per_cell_by_branch,
                                          data.frame("branch" = b,"CellType"=c,"n_signif" = n))
    }
  }


  # Run hypergeometric tests
  # To show which cell types within a branch have more enriched phenotypes than expected
  n_signif_per_cell_by_branch$hypergeo_p = NA
  total_signif = length(signif_results$q)

  for (i in seq(1, length(n_signif_per_cell_by_branch$CellType))) {
    cell = n_signif_per_cell_by_branch$CellType[i]
    branch = n_signif_per_cell_by_branch$branch[i]
    group1 = sum(n_signif_per_cell_by_branch$n_signif[n_signif_per_cell_by_branch$CellType == cell])
    group2 = sum(n_signif_per_cell_by_branch$n_signif[n_signif_per_cell_by_branch$branch == branch])
    overlap = n_signif_per_cell_by_branch$n_signif[i]
    n_signif_per_cell_by_branch$hypergeo_p[i] = phyper(overlap - 1, group2, total_signif - group2, group1, lower.tail=FALSE)
  }
  n_signif_per_cell_by_branch$hypergeo_q = p.adjust(n_signif_per_cell_by_branch$hypergeo_p,method = "BH")

  n_signif_per_cell_by_branch$signif_asterics = ""
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.05] = "*"
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.005] = "**"
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.0005] = "***"
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.00005] = "****"

  # Ordering cell types by dendrogram grouping
  cell_order = factor(gsub("_"," ",ctd[[1]]$plotting$cell_ordering))
  n_signif_per_cell_by_branch$cell_order = match(n_signif_per_cell_by_branch$CellType, cell_order)
  n_signif_per_cell_by_branch$CellType = reorder(n_signif_per_cell_by_branch$CellType, n_signif_per_cell_by_branch$cell_order)

  facet_branch_plt <- ggplot(n_signif_per_cell_by_branch[n_signif_per_cell_by_branch$branch %in% plot_branches,] , aes(x=CellType,y=n_signif,fill=branch)) +
    geom_col() +
    geom_text(mapping= aes(label = signif_asterics, y = n_signif + 5))+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle = 90, hjust=1),legend.position="none") +
    #scale_x_discrete(breaks =cell_order)+
    facet_wrap(~branch,ncol=1)

  return(facet_branch_plt)
}
