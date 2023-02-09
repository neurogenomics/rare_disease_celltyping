tfidf <- function(clusts,
                   col_name="dataset",
                   cluster_name="seurat_clusters",
                   terms_per_cluster=1,
                   replace_regex="[.]|_|-"){
  library(tidytext)
  data(stop_words)
  # From here: https://www.tidytextmining.com/tfidf.html
  clusts$cluster <- clusts[[cluster_name]]
  clusts$var <- gsub(replace_regex," ",clusts[[col_name]])
  clust_words <-clusts %>%
    unnest_tokens(word, var, drop = F) %>%
    anti_join(stop_words, keep=T) %>%
    dplyr::count(cluster, word, sort = TRUE)
  total_words <- clust_words %>%
    dplyr::group_by(cluster, .drop = F) %>%
    dplyr::summarize(total = sum(n, na.rm = T))
  total_samples <- clusts %>%
    dplyr::group_by(cluster, .drop = F) %>%
    dplyr::count(name = "samples")
  clust_words <- left_join(clust_words, total_words) %>%
    left_join(total_samples)
  clust_words <- clust_words %>%
    bind_tf_idf(word, cluster, n) %>%
    dplyr::group_by(cluster, .drop = F) %>%
    dplyr::slice_max(order_by = tf_idf, n=terms_per_cluster, with_ties = F) %>%
    subset(!word %in% c("cell","cells", stop_words))
  return(clust_words)
}


seurat_tfidf <- function(seurat,
                         label_var="label",
                         cluster_var="seurat_clusters",
                         replace_regex = " ",
                         terms_per_cluster=3,
                         force_new=F){
  # label_var="label";cluster_var="seurat_clusters";replace_regex = " ";terms_per_cluster=3;force_new=F
  input_dat <- seurat@meta.data
  if(any(c("enriched_words","tf_idf") %in% colnames(input_dat))){
    if(force_new){
      input_dat <- input_dat %>% dplyr::select(-c(enriched_words,tf_idf))
    }else {
      message("Previous TF-IDF results detected. Use force_new=T to re-run.")
      return(seurat)
    }
  }

  tfidf_df <- tfidf(clusts = input_dat,
                    col_name =  label_var,
                    cluster_name = cluster_var,
                    replace_regex = replace_regex,
                    terms_per_cluster = terms_per_cluster)

  new_metadata <- merge(x=input_dat,
                        y=tfidf_df %>%
                          dplyr::group_by(cluster) %>%
                          dplyr::summarise(enriched_words=paste(unique(word),collapse="; "),
                                           tf_idf=paste(unique(tf_idf),collapse="; ")),
                        all.x = T,
                        by.x="seurat_clusters", by.y = "cluster")
  # Make sure rows are in the right order
  new_metadata <- new_metadata[match(input_dat[[label_var]], new_metadata[[label_var]]),]
  if(sum(new_metadata[[label_var]]!=input_dat[[label_var]], na.rm = T)>0){
    stop("Seurat sample names and tfidf samples names are not aligned!")
  }
  seurat@meta.data <- data.frame(new_metadata, row.names = row.names(input_dat))
  return(seurat)
}


umap_tfidf <- function(seurat,
                       label_var="label",
                       cluster_var="seurat_clusters",
                       size_var=1,
                       color_var="cluster",
                       point_alpha=.7,
                       point_palette=c(unname(pals::alphabet()),rev(unname(pals::alphabet2()) )),
                       density_palette="Purples",
                       density_adjust=.2,
                       label_fill=alpha(c("white"),0.7),
                       show_plot=T,
                       replace_regex=" ",
                       terms_per_cluster=3,
                       ...){
  # label_var="label";cluster_var="seurat_clusters";size_var="genes";point_alpha=.7;density_palette="purples";density_adjust=.2;
  # label_fill=alpha(c("white"),0.7);show_plot=T;replace_regex=" ";terms_per_cluster=3; show_plot=T;

  # Prepare input_dat
  input_dat = seurat@meta.data[,!startsWith(colnames(seurat@meta.data), "UMAP")]
  input_dat  <- cbind(input_dat, seurat@reductions$umap@cell.embeddings)
  input_dat$cluster <-  input_dat[[cluster_var]]
  if(any(c("enriched_words","tf_idf") %in% colnames(input_dat))){
    input_dat <- input_dat %>% dplyr::select(-c(enriched_words,tf_idf))
  }
  # Prepare tfidf_df
  tfidf_df <- tfidf(clusts = input_dat,
                    col_name =  label_var,
                    cluster_name = cluster_var,
                    replace_regex = replace_regex,
                    terms_per_cluster = terms_per_cluster)


  cluster_centers <- input_dat %>%
    merge(tfidf_df %>%
            # Make sure non-enriched clusters get annotated
            dplyr::mutate(word=ifelse(is.na(word),"N/A",word)),
          by="cluster") %>%
    dplyr::group_by(cluster, word, .drop=F) %>%
    dplyr::summarise(x.mean=mean(UMAP_1, na.rm=T),
                     y.mean=mean(UMAP_2, na.rm=T),
                     size=mean(tf_idf, na.rm = T)) %>%
    dplyr::rename(term=word)


  cluster_number_centers <- cluster_centers %>%
    group_by(cluster, .drop=F) %>%
    dplyr::slice_head(n = 1)

  library(ggplot2)
  umap_tfidf <- ggplot(data = input_dat,
                       aes_string(x="UMAP_1", y="UMAP_2", size=size_var)) +
    stat_density_2d(aes(fill = ..level..),
                    adjust = density_adjust, contour = T,
                    geom = "polygon",contour_var = "ndensity") +
    scale_fill_distiller(palette=density_palette, direction=-1) +
    geom_point(aes_string(color=color_var,
                   # shape=Type,
                   size=size_var,
                   ...),
               alpha=point_alpha, show.legend = F) +
    scale_color_manual(values = point_palette) +
    # geom_label(data = cluster_centers,
    #            aes(x=x.mean, y=y.mean, label=cluster, size=size, color=cluster),
    #               fill = alpha(c("white"),0.8), inherit.aes = F, show.legend = F) +
    ggrepel::geom_label_repel(data = cluster_number_centers,
                              aes(x=x.mean, y=y.mean, label=cluster, size=size, color=cluster),
                              fill = alpha(c("black"),0.5),
                              ## Must set this outside of aes() due to bug that causes
                              ## conflict with geom_point (even when inherit.aes=F)
                              size=5,
                              min.segment.length = 0.1, box.padding = 3,
                              inherit.aes = F, max.overlaps = 30, show.legend = F) +
    ggrepel::geom_label_repel(data = cluster_centers,
                              aes(x=x.mean, y=y.mean, label=term, size=size, color=cluster),
                              fill = label_fill,
                              ## Must set this outside of aes() due to bug that causes
                              ## conflict with geom_point (even when inherit.aes=F)
                              # Also, log and rescale size so that all labels are still legible
                              size=scales::rescale(log10(cluster_centers$size), c(2,6)),
                              min.segment.length = 0.1, box.padding = .1,
                              inherit.aes = F, max.overlaps = 30, show.legend = F) +
    theme_void() +
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          legend.text = element_text(color="white"),
          legend.title = element_text(color="white")
          )

  if(show_plot) print(umap_tfidf)

  return(umap_tfidf)
}
