#' Plot RD EWCE results subset as interactive network plot
#'
#' This coerces the network object into a plot with ggnetwork (which assigns
#' x, y coordinates to each node in the network). The hover box with results and
#' disease description for each node is also added here in the
#' \code{pheno_ggnetwork$hover} column. Once the x and y coordiantes have been
#' added, it can be plot using ggplot2.
#' @param phenoNet The network object created using \code{create_network_object}
#' @param phenos The subset of results to be plotted (data frame)
#' @param disease_descriptions The data frame of all disease descriptions, This is
#' where new lines are added using the hpo_term_definition_list function.
#' @returns A interactive plot of the network of phenotypes in the selected subset of results.
#' @import ggplot2
#' @export
ggnetwork_plot <- function(phenoNet,phenos,disease_descriptions) {
  term_definitions = hpo_term_definition_list(phenos$HPO_term_Id,disease_descriptions)
  pheno_ggnetwork = ggnetwork::ggnetwork(phenoNet, arrow.gap=0)
  pheno_ggnetwork$hover = paste(pheno_ggnetwork$label,
                                "\nId:",phenos$HPO_term_Id,
                                "\nFold:",phenos$fold_change,
                                "\nq:",phenos$q,
                                "\nDefinition:",term_definitions)

  network_plot <-  ggplot(pheno_ggnetwork, aes(x=x,y=y,xend=xend,yend=yend,text=hover)) +
    ggnetwork::geom_edges(color = "darkgray")+
    geom_point(aes(colour = fold, size = heirarchy)) +  #, text= hover)) +

    geom_text(aes(label = label), color = "black",alpha = 0.7) +
    scale_colour_gradient2(low = "white", mid = "yellow", high = "red") +
    scale_size(trans = "exp") +
    guides(size = FALSE)+
    labs(colour="Fold")+
    theme_blank() #, tooltip = "hover") put ggplotly back before ggplot above
  return(network_plot)
}
