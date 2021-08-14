
# Number of significant phenotypes per cell

cell_mappings = read.csv("data/DescartesHuman_celltype_mapping.csv")
cell_mappings$level1 = gsub("_"," ",cell_mappings$level1)


plot_phenos_per_cell1 = function (all_results_merged, cell_mappings, fold = 1, q_val = 0.005, dataset = "Descartes") {
  phenos_per_cell = data.frame()

  # count n phenotypes per cell
  for (c in unique(all_results_merged$CellType)){
    cur_cell = c
    n_phenos = length(all_results_merged[all_results_merged$CellType==c & all_results_merged$q < q_val & all_results_merged$fold_change > fold, "list"])
    phenos_per_cell = rbind(phenos_per_cell, data.frame("Cell"=cur_cell,"n_phenos"=n_phenos))
  }
  n_pheno_text = phenos_per_cell
  phenos_per_cell$Tissue = rep(NA, length(phenos_per_cell$Cell))

  if (dataset == "Descartes") {
    phenos_per_cell_tissue = data.frame()
    for (i in seq(1:length(phenos_per_cell$Cell))) {
      cur = cell_mappings[cell_mappings$level1 == paste(phenos_per_cell$Cell[i]), ]$tissue
      for (t in unique(cur)) {
        if (! is.na(t)){
          phenos_per_cell$Tissue[i] = t
          phenos_per_cell_tissue = rbind(phenos_per_cell_tissue, phenos_per_cell[i,])
        }
      }
    }

    phenos_per_cell = phenos_per_cell_tissue
  }
  phenos_per_cell$Cell = reorder(phenos_per_cell$Cell, phenos_per_cell$n_phenos)

  if (dataset == "Descartes") {
    for (i in seq(1,length(phenos_per_cell$Cell))) {
      cur = phenos_per_cell[phenos_per_cell$Cell == phenos_per_cell$Cell[i],]
      phenos_per_cell$n_phenos[i] = phenos_per_cell$n_phenos[i] / length(unique(cur$Tissue))
    }
  }

  explan_plt <- ggplot(phenos_per_cell, aes(n_phenos ,  as.factor(Cell))) +
    #geom_point(size = 12, color = "black") +

    scale_x_continuous("",expand = c(0,0),  limits = c(0, max(n_pheno_text$n_phenos)+20))+#, position = "top") + #
    theme_classic()+
    scale_y_discrete(phenos_per_cell$Cell)+
    labs(title = paste0("Significantly enriched phenotypes per cell (q < ", q_val, ", fold change > ", fold,")" )) +
    ylab("Cell Type")

  if (dataset == "Descartes") {
    explan_plt = explan_plt + geom_col(aes(fill = Tissue),size = 2,position = "stack")
  } else {
    explan_plt = explan_plt + geom_segment(aes(xend = 0, yend = Cell), size = 2)
  }
  # explan_plt = explan_plt +
  #   geom_vline(xintercept = mean(n_phenos), color = "grey40", linetype = 3, size = 2, alpha = 0.8) +
  #   annotate("text", x = mean(n_phenos)-7,size =5, y = 4,color = "grey40", label = "Mean",angle=90)

  explan_plt= explan_plt + geom_text(data=n_pheno_text,aes(label = n_phenos, x = n_phenos+ 6, y = Cell), size = 4, color = "black") +
    ylab("Cell Type") +#+ xlab("Cell Type")
    theme(axis.title.y = element_blank())
  # adding a vertical line to show mean
  return (explan_plt)
}





################################################################################

# get the phenos per cell dataframe
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

  phenos_per_cell$Cell = reorder(phenos_per_cell$Cell, phenos_per_cell$n_phenos)

 return(phenos_per_cell)
}

