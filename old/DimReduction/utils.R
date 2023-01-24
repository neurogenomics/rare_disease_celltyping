binarize_ctd <- function(ctd,
                         level=2,
                         replace_nonzeros=T,
                         top_quantiles=2,
                         top_genes=NULL,
                         as_sparseMatrix=T){
  spec <- ctd[[level]]$specificity
  mat_cells <- spec
  specQ <- ctd[[level]]$specificity_quantiles

  if(is.null(top_genes)){
    print(paste("Selected top",top_quantiles,"quantiles per celltype."))
    for(x in colnames(mat_cells)){
      # print(x)
      quants <- unique(sort(specQ[,x], decreasing = T))
      fill0 <- which(!specQ[,x] %in% quants[1:top_quantiles])
      mat_cells[fill0,x] <- 0
    }
  } else {
    print(paste("Selected top",top_genes,"genes per celltype."))
    for(x in colnames(mat_cells)){
      # print(x)
      selected_genes <- sort(setNames(spec[,x], row.names(spec)), decreasing = T)[1:top_genes]
      mat_cells[!row.names(mat_cells) %in% selected_genes,x] <- 0
    }
  }
  if(replace_nonzeros){
    print("Replacing non-zero values with 1.")
    mat_cells[mat_cells>0] <- 1
  }

  if(as_sparseMatrix){
    mat_cells <- as(as.matrix(mat_cells), "sparseMatrix")
  }
  return(mat_cells)
}
