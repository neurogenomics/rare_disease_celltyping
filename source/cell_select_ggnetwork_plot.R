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

hpo_get_term_definition <- function(ontologyId, disease_descriptions) {
   definition = disease_descriptions[ontologyId,"description"]
   definition = newlines_to_definition(definition)
   return (definition)
}

# applying hpo_get_term_definition to vector of Id's
hpo_term_definition_list <- function(ontologyId_list, disease_descriptions) {
   term_details <- c()
   for (term in ontologyId_list) {
      term_details[term] <- hpo_get_term_definition(term, disease_descriptions)
   }
   return (term_details)
}



# add newlines to the definition so it fits in hover box
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




# CREATE ADJACENCY MATRIX

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


 # find hierarchy for phenos



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



subset_phenos = function(phenotype_to_genes, all_results_merged, hpo,cell_type = "Neurons", q_threshold =0.0005, fold_threshold = 1) {
   phenos = get_cell_ontology(cell_type,all_results_merged,q_threshold = q_threshold, fold_threshold = fold_threshold, phenotype_to_genes, hpo)
   phenos = phenos[!is.na(phenos$HPO_term_Id) & phenos$HPO_term_valid,]
   return (phenos)
}

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




ggnetwork_plot_full = function(phenotype_to_genes, all_results_merged, hpo, disease_descriptions,cell_type = "Neurons", q_threshold =0.0005, fold_threshold = 1){
  phenos = subset_phenos(phenotype_to_genes, all_results_merged, hpo,cell_type =cell_type, q_threshold =q_threshold, fold_threshold = fold_threshold)
  adjacency = adjacency_matrix(unique(phenos$HPO_term_Id), hpo)
  phenoNet = make_network_object(phenos,adjacency,hpo)
  network_plot = ggnetwork_plot(phenoNet, phenos,disease_descriptions)
  return(network_plot)
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
