#devtools::install_github("briatte/ggnet")

library(ggnetwork)
library(ggnet)
library(Rgraphviz)
library(sna)
library(network)
library(ontologyIndex)
source("source/signif_enrichments_for_cell.R")

data(hpo)
phenotype_to_genes = read.delim("data/phenotype_to_genes.txt", skip = 1, header=FALSE)
colnames(phenotype_to_genes) = c("ID", "Phenotype", "EntrezID", "Gene",
                               "Additional", "Source", "LinkID")
disease_descriptions = readRDS("data/disease_descriptions.Rda")
rownames(disease_descriptions) = disease_descriptions$HPO_id

# HPO API GET REQUEST - disease description

require(httr)
require(jsonlite)


#' Get HPO term definition
#'
#' This gets the disease description from a data frame of disease descriptions.
#' The rows names of the data frame are the HPO ID and the description column is
#' called "description". It also adds new lines to the description so that the
#' hover box in the web app does not get too wide. This is done by calling the
#' \code{newlines_to_definition} function.
#'
#' @param ontologyId The HPO Id of the term (string)
#' @param disease_descriptions A data frame of disease descriptions corresponding to each HPO Id
#'
#' @return The disease description with new lines added.
#'
#' @examples
#' hpo_get_term_definition("HP:123456", disease_descriptions)
#'
#' @export
hpo_get_term_definition <- function(ontologyId, disease_descriptions) {
   definition = disease_descriptions[ontologyId,"description"]
   definition = newlines_to_definition(definition)
   return (definition)
}



#' Get HPO term description for a list of terms
#'
#' Need to redo this without a for loop (use \code{lapply} or something). This just
#' applys the \code{hpo_get_term_definition} function to a character vector of terms.
#'
#' @param ontologyId_list A character vector of HPO Ids
#' @param disease_descriptions A data frame of disease descriptions for all HPO Id
#'
#' @retun A named vector of disease descriptions, with HPO Id as names and descriptions
#' as values.
#'
#' @examples
#' hpo_term_definition_list(HPO_terms_char_vector, Disease_description_df)
#'
#' @export
hpo_term_definition_list <- function(ontologyId_list, disease_descriptions) {
   term_details <- c()
   for (term in ontologyId_list) {
      term_details[term] <- hpo_get_term_definition(term, disease_descriptions)
   }
   return (term_details)
}



#' Add new lines to disease description
#'
#' Adds new lines to the description so that hover boxes dont get too wide.
#'
#' @param definition A disease description string
#' @param line_length A integer representing the desired words per line.
#'
#' @returns The disease description with newline symbols added every nth word.
#'
#' @examples
#' newlines_to_definition(disease_description, 10)
#' @export
newlines_to_definition <- function(definition, line_length = 10) {
   definition = strsplit(definition, split = " ")[[1]]
   if (length(definition) > line_length) {
      remainder = length(definition) %% line_length
      n_new_lines = floor((length(definition)/line_length))
      new_line_index = seq(line_length,(n_new_lines*line_length),line_length)
      definition[new_line_index] = paste0("\n", definition[new_line_index])
   }
   definition = paste(definition,collapse = " ")
   return(definition)
}




#' Create adjacency matrix of HPO child-parent relationships
#'
#' This is needed for plots created using ggnetwork as it can coerce an
#' adjacency matrix into a directed graph object, and assigns
#' each datapoint an x,y coordinate. Maybe shouldn't use a for loop for this.
#' It also may be possible to use a hash table for ggnetwork, which may be more
#' efficent for the web app
#'
#' @param pheno_ids a character vector of HPO Ids
#' @param hpo ontology object (available in ontologyIndex package)
#'
#' @returns adjacency matrix with HPO Ids for col and row names.
#' If adjacency[i,j] == 1 then phenotype[i] is a parent of phenotype[j]
#'
#' @examples
#' adjacency_matrix(hpo_id_char_vector, hpo)
#' @export

 adjacency_matrix <- function(pheno_ids, hpo) {

   HPO_id = unique(pheno_ids)
   size = length(HPO_id)
   adjacency = data.frame(matrix(nrow = size, ncol = size))
   rownames(adjacency) = HPO_id
   colnames(adjacency) = HPO_id
   adjacency[is.na(adjacency)] = 0
   for (id in HPO_id) {
     children = hpo$children[id][[1]]
     adjacency[id, children] = 1
   }
   return(adjacency[HPO_id,HPO_id])
 }




#' Identify relative ontology level for HPO term in a subset of ontology
#'
#' When plotting subsets of the HPO data as a network plot, there is often more than
#' one connected component to be plotted (i.e. sections of the graph are not connected
#' by any edges). To map the ontology level to the size of nodes (such that high up terms
#' are bigger than low level terms). This was made to find the relative ontology level
#' with respect to other terms in its connected component so that each connected section
#' of the plot would have the same size for its root node. It would also be possible
#' to just use absolute ontology level, so that "Phenotypic abnormality" would always
#' be the largest datapoint etc.
#'
#' @param phenotype HPO term Id (string)
#' @param phenoAdj A adjacency matrix (produced by the adjacency_matrix function)
#' @param hpo The HPO ontology data object
#'
#' @returns A integer representing the relative ontology level of a term within
#' a connected component of a subset of the HPO.
#' @export

 find_parent <- function (phenotype,phenoAdj,hpo){
   pos_parents = hpo$parents[phenotype]
   phenotypes = rownames(phenoAdj)
   paths = list()
   for (p in phenotypes){
     if (phenoAdj[p,phenotype] == 1) {
       if (p %in% pos_parents) {
         paths[p] = 1 + find_parent(p,phenoAdj,hpo) # <- recursion
       }
     }
   }
   if (length(paths) == 0) {
     return (0)
   } else {
     parents = 0
     for (i in seq(length(paths))) {
       if (paths[[i]] > parents) {
         parents = paths[[i]]
       }
     }
   }
   return (parents)
 }


#' Get component relative ontology level of all terms within a subset of HPO
#'
#' This calls the \code{find_parent} function on all phenotypes in the subset of
#' the HPO to be plotted. The subest chosen when creating the phenoAdj from the main
#' adjacency matrix of all phenotypes. So, the phenotypes to be plotted can be
#' found in the row and column names of phenoAdj.
#'
#' @param phenoAdj A adjacency matrix of phenotypes where 1 represents i is parent of j
#' and 0 represents that i is not a parent of j. It is a subset of the main phenotype adjacency matrix
#' @param hpo The HPO ontology data object
#' @param reverse A boolean, if TRUE it will reverse the ontology level numbers so that
#' the parent terms are larger than the child terms.
#' @returns A named vector of relative ontology level, where names are HPO Ids and
#' value is relative ontology level.
#' @export
 get_heirarchy <- function (phenoAdj,hpo,reverse=TRUE) {
   heirarchy = c()
   for (p in rownames(phenoAdj)) {
     heirarchy[p] = find_parent(p,phenoAdj,hpo)
   }
   if (reverse) {
     heirarchy = max(heirarchy) - heirarchy
   }
   return (heirarchy)
 }


#' Subset RD EWCE Results
#'
#' This subsets  the Rare disease EWCE results by cell type, q threshold and fold change.
#'
#' @param phenotype_to_genes The list of HPO terms with their assocaited gene lists taken from HPO website
#' @param all_results_merged The dataframe of RD EWCE Results
#' @param hpo The HPO ontology data object
#' @param cell_type A string representing the cell type of interest.
#' @param q_threshold The q threshold. The subset of results will have a q lower than this
#' @param fold_threshold The fold change threshold. The subest of results will have a fold change greater than this.
#'
#' @returns A data frame of results taken from the main data frame of results
#' @export
subset_phenos = function(phenotype_to_genes, all_results_merged, hpo,cell_type = "Neurons", q_threshold =0.0005, fold_threshold = 1) {
   phenos = get_cell_ontology(cell_type,all_results_merged,q_threshold = q_threshold, fold_threshold = fold_threshold, phenotype_to_genes, hpo)
   phenos = phenos[!is.na(phenos$HPO_term_Id) & phenos$HPO_term_valid,]
   return (phenos)
}



#' Make network object
#'
#' This uses the network package to coerce the adjacency matrix into a
#' network object. It also adds the fold change, label, and relative ontology level
#' parameters to each node in the network.
#'
#' @param phenos The subset of the results to be plotted
#' @param adjacency The adjacency matrix of all HPO terms
#' @param hpo The HPO ontology data object
#'
#' @returns A ggnetowrk graph/ network object of a subset of the RD EWCE results.
#' @export
make_network_object = function(phenos, adjacency, hpo) {
   ValidTerms = phenos$HPO_term_Id
   phenoAdj = adjacency[ValidTerms,ValidTerms]
   # MAKE NETWORK OBJECT
   phenoNet = network(phenoAdj, directed = TRUE)
   # To add a another value to the nodes do this
   phenoNet %v% "fold" = as.numeric(phenos$fold_change)
   phenoNet %v% "label" = as.character(phenos$list)
   phenoNet %v% "heirarchy" = as.numeric(get_heirarchy(phenoAdj,hpo,reverse=TRUE)+1)
   # (%v% is for setting vertex attributes, %e% is for edge, %n% is for network attributes)
  return(phenoNet)
 }



# PLOT USING GGNETWORK
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
#' @export
ggnetwork_plot <- function(phenoNet,phenos,disease_descriptions) {
  term_definitions = hpo_term_definition_list(phenos$HPO_term_Id,disease_descriptions)
  pheno_ggnetwork = ggnetwork(phenoNet, arrow.gap=0)
  pheno_ggnetwork$hover = paste(pheno_ggnetwork$label,
                                "\nId:",phenos$HPO_term_Id,
                                "\nFold:",phenos$fold_change,
                                "\nq:",phenos$q,
                                "\nDefinition:",term_definitions)

  network_plot <-  ggplot(pheno_ggnetwork, aes(x=x,y=y,xend=xend,yend=yend,text=hover)) +
    geom_edges(color = "darkgray")+
    geom_point(aes(colour = fold, size = heirarchy)) +  #, text= hover)) +

    geom_text(aes(label = label), color = "black",alpha = 0.7) +
    scale_colour_gradient2(low = "white", mid = "yellow", high = "red") +
    scale_size(trans = "exp") +
    guides(size = FALSE)+
     labs(colour="Fold")+
     theme_blank()#, tooltip = "hover") put ggplotly back before ggplot above
  return(network_plot)
}


#' Create interactive network plot start to finish (including subset data)
#'
#' This puts all the functions together from gettig the subest of results to creating
#' the final interactive plot.
#'
#'  @param penotype_to_genes The phenotype gene lists taken from the HPO
#'  @param all_results_merged The RD EWCE Results
#'  @param hpo The HPO ontology data object
#'  @param disease_descriptions The dataframe of all disease descriptions in the HPO
#'  @param cell_type The cell type of interest to be plotted
#'  @param q_threshold The q value threshold for the subset of results to be plotted
#'  @param fold_threshold The minimum fold change in specific expression for the subest of results to be plotted
#'
#'  @return A interactive network plot of the selected subset of results from RD EWCE analysis
#'  @export
ggnetwork_plot_full = function(phenotype_to_genes, all_results_merged, hpo, disease_descriptions,cell_type = "Neurons", q_threshold =0.0005, fold_threshold = 1){
  phenos = subset_phenos(phenotype_to_genes, all_results_merged, hpo,cell_type =cell_type, q_threshold =q_threshold, fold_threshold = fold_threshold)
  adjacency = adjacency_matrix(unique(phenos$HPO_term_Id), hpo)
  phenoNet = make_network_object(phenos,adjacency,hpo)
  network_plot = ggnetwork_plot(phenoNet, phenos,disease_descriptions)
  return(network_plot)
}


#' Get absolute ontology level
#'
#' This gets the absolute ontology level of a term (without consideration for the
#' particular subset of the data you are looking at, as in the find_parent function).
#'
#' @param hpo The HPO ontology data object
#' @param term_id HPO term ID <string>
#' @example get_ont_level(hpo,"HP:0000003")
#' @return returns the ontology level <numeric>
#' @export
get_ont_level = function(hpo,term_id) {
   children = unique(setdiff(unlist(hpo$children[term_id]), term_id))
   if (length(children) == 0) {
      return(0)
   } else {
      return(1 + get_ont_level(hpo,children)) #<- recursion..
   }
}



# a function that allows you to choose between the ggnetwork interactive plot and the ontoplot (better for printing maybe)
# ggnetwork_or_ontoplot_full = function(onto_or_ggnet, phenotype_to_genes, all_results_merged,
#                                       hpo, disease_descriptions, cell_type,
#                                       q_threshold = 0.0005, fold_threshold =1) {
#    if (onto_or_ggnet == "ontoplot") {
#       return(one_cell_ontology_plot_heatmap(all_results_merged,cell=cell_type,
#                                              heatmapped_value="fold change",q_threshold=q_threshold,
#                                              fold_threshold=fold_threshold,
#                                              phenotype_to_genes=phenotype_to_genes,hpo))
#    } else if (onto_or_ggnet == "ggnetwork") {
#       return(ggnetwork_plot_full(phenotype_to_genes, all_results_merged,
#                           hpo,disease_descriptions, cell_type = cell_type,
#                           q_threshold = q_threshold,
#                           fold_threshold = fold_threshold))
#    }
#
#    }




# TEST
# ggnetwork_plot_full(phenotype_to_genes, all_results_merged, hpo,cell_type = "Neurons", q_threshold =0.0005, fold_threshold = 1)
