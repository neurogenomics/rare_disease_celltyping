#' Modification to ewce.plot function to automatically align the dendrogram
#'
#' @param total_res Results to be plotted ?
#' @param mtc_method Multiple comparison adjustment method
#' @param ctd Cell type data object <list>
#' @examples
#' \dontrun{ewce.plot(results_dataframe,"bonferroni",CTD)}
#' @returns A bar chart with dendrogram of EWCE results in each cell type.
#' @export
ewce.plot <-function (total_res, mtc_method = "bonferroni", ctd = NULL)
{
  if (!mtc_method %in% c("holm", "hochberg", "hommel",
                         "bonferroni", "BH", "BY", "fdr",
                         "none")) {
    stop("ERROR: Invalid mtc_method argument. Please see '?p.adjust' for valid methods.")
  }
  multiList = TRUE
  if (is.null(total_res$list)) {
    multiList = FALSE
  }
  make_dendro = FALSE
  if (!is.null(ctd)) {
    make_dendro = TRUE
    cells.in.ctd <- function(ctdIN, cells) {
      if (sum(!cells %in% colnames(ctdIN$specificity) ==
              0)) {
        return(1)
      }
      else {
        return(0)
      }
    }
    if (length(ctd[[1]]$plotting) > 0) {
      annotLevel = which(unlist(lapply(ctd, FUN = cells.in.ctd,
                                       cells = as.character(total_res$CellType))) ==
                           1)
      if (length(annotLevel) == 0) {
        stop("All of the cells within total_res should come from a single annotation layer of the CTD")
      }
    }
    if (length(ctd[[annotLevel]]$plotting) > 0) {
      total_res$CellType = factor(total_res$CellType, levels = ctd[[annotLevel]]$plotting$cell_ordering)
    }
  }
  total_res$q = stats::p.adjust(total_res$p, method = mtc_method)
  ast_q = rep("", dim(total_res)[1])
  ast_q[total_res$q < 0.05] = "*"
  total_res$ast_q = ast_q
  total_res$sd_from_mean[total_res$sd_from_mean < 0] = 0
  graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
    theme(panel.grid.major = element_line(size = 0.5, color = "grey"),
          axis.line = element_line(size = 0.7, color = "black"),
          text = element_text(size = 14), axis.title.y = element_text(vjust = 0.6))
  upperLim = max(abs(total_res$sd_from_mean))
  total_res$y_ast = total_res$sd_from_mean * 1.05
  total_res$abs_sd = abs(total_res$sd_from_mean)
  if ("Direction" %in% colnames(total_res)) {
    the_plot = ggplot(total_res) + geom_bar(aes_string(x = "CellType",
                                                       y = "abs_sd", fill = "Direction"), position = "dodge",
                                            stat = "identity") + graph_theme
  }
  else {
    the_plot = ggplot(total_res, mapping = aes(x=CellType)) + geom_bar(aes_string(x = "CellType",
                                                       y = "abs_sd"), fill = "red", stat = "identity") +
      graph_theme + theme(legend.position = "none")
  }
  the_plot = the_plot + theme(plot.margin = unit(c(1, 0, 0,
                                                   0), "mm"), axis.text.x = element_text(angle = 90,
                                                                                         hjust = 1,vjust=0.2)) + theme(panel.border = element_rect(colour = "black",
                                                                                                                                         fill = NA, size = 1)) + xlab("") + theme(strip.text.y = element_text(angle = 0)) +
    coord_cartesian(ylim = c(0, 1.1 * upperLim)) + ylab("Std.Devs. from the mean") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  the_plot = the_plot + scale_y_continuous(breaks = c(0, ceiling(upperLim *
                                                                   0.66))) + geom_text(aes_string(label = "ast_q",
                                                                                                  x = "CellType", y = "y_ast"), size = 10)
    #scale_x_discrete(breaks = ctd[[annotLevel]]$plotting$cell_ordering)
  # CHANGE 2

  if (multiList) {
    the_plot = the_plot + facet_grid("list ~ .", scales = "free",
                                     space = "free_x")
  }
  output = list()
  output$plain = the_plot
  if (make_dendro) {
    the_dendrogram = ctd[[annotLevel]]$plotting$ggdendro_horizontal +
      theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm")) +
    #  scale_x_discrete(breaks = total_res$CellType)
    scale_x_discrete(breaks = ctd[[annotLevel]]$plotting$cell_ordering)
    # CHANGE ^: added scale_x_discrete to set the mapping of the dendro to the x axis scale

    combined_plot = cowplot::plot_grid(the_dendrogram, the_plot,axis = "lr",
                                       align = "v", ncol = 1,
                                       rel_heights = c(0.3, 1))
    # CHANGE ^: align argument to "v" and rel_heights to c(0.3,1) to make dend and barchart closer

    output$withDendro = combined_plot
  }
  return(output)
}
