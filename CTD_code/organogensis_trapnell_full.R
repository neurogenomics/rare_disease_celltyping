# build ctd
library(monocle)
library(MAGMA.Celltyping)

df <- readRDS("cds_cleaned.RDS")
cell_annotations <- read.csv("cell_annotate.csv")

cell_names <- data.frame(unique(cell_annotations$Main_Cluster))
cell_names$types <- unique(cell_annotations$Main_cell_type)
cell_names <- cell_names[order(cell_names$unique.cell_annotations.Main_Cluster.),]
cell_names <- cell_names[1:37,]
rownames(cell_names) <- cell_names$unique.cell_annotations.Main_Cluster.

cell_names[,2] <- as.character(cell_names[,2])
cell_names[2,2] <- "Chondrocytes and osteoblasts"
cell_type <- cell_names[,2]
cell_names <- paste(cell_names[,2],cell_names[,1],sep="_")
#exp_matrix <- df@assayData$exprs

cell_cluster_1 <- df@phenoData@data$Cluster
cell_cluster_2 <- df@phenoData@data$sub_Cluster
cell_cluster_2 <- strtoi(cell_cluster_2)

cell_type <- cell_type[cell_cluster_1]

cell_level_1 <- cell_names[cell_cluster_1]
cell_level_2 <- df@phenoData@data$sub_cluster_id
#replace the "-" in the annotations
cell_level_2 <- sapply(cell_level_2,test <- function(x){return(sub("-","_",x))})
cell_level_2 <- paste(cell_type,cell_level_2,sep="_")

gene_mgi <- df@featureData@data$gene_short_name
rownames(df@assayData$exprs) = gene_mgi
annotLevels = list(level1class=cell_level_1,level2class=cell_level_2)
ctd <- makeCTDdata(df@assayData$exprs,annotLevels,groupName="trapnel")

############################

x_1 <- data.frame(gene_mgi)
x_2 <- data.frame(gene_mgi)
n_cluster_1 <- length(unique(cell_cluster_1))
for (i in 1:n_cluster_1){
  #level 1 cluster ids
  cid_1 <- which(cell_cluster_1 == i)
  mean_exp <- Matrix::rowMeans(df@assayData$exprs[,cid_1])
  x_1$new_cell <- mean_exp
  colnames(x_1)[which(names(x_1) == "new_cell")] <- cell_level_1[cid_1[1]]
  n_cluster_2 <- length(unique(cell_cluster_2[cid_1]))
  n_0 <- 0
  for (j in 1:n_cluster_2){
    #level 2 cluster ids
    k <- j + n_0
    while(length(which(cell_cluster_2[cid_1] == k)) == 0){
      n_0 <- n_0 + 1
      k <- j + n_0
    }
    cid_2 <- which(cell_cluster_2[cid_1] == k)
    mean_exp_2 <- Matrix::rowMeans(df@assayData$exprs[,cid_1[cid_2]])
    x_2$new_cell <- mean_exp_2
    colnames(x_2)[which(names(x_2) == "new_cell")] <- cell_level_2[cid_1[cid_2[1]]]
  }
}

# remove genes with duplicate gene symbols
y_1 <- x_1[which(duplicated(x_1$gene_mgi)==0),]
row.names(y_1) <- y_1$gene_mgi
y_1$gene_mgi <- NULL

y_2 <- x_2[which(duplicated(x_2$gene_mgi)==0),]
row.names(y_2) <- y_2$gene_mgi
y_2$gene_mgi <- NULL

normalised_meanExp_1 = t(t(y_1)*(1/colSums(y_1)))
specificity_1 = normalised_meanExp_1/(apply(normalised_meanExp_1,1,sum)+0.000000000001)

normalised_meanExp_2 = t(t(y_2)*(1/colSums(y_2)))
specificity_2 = normalised_meanExp_2/(apply(normalised_meanExp_2,1,sum)+0.000000000001)

ctd <- vector("list", 2) 
ctd[[1]]$mean_exp <- y_1
ctd[[1]]$specificity <- specificity_1
ctd[[1]]$annot <- cell_level_1

ctd[[2]]$mean_exp <- y_2
ctd[[2]]$specificity <- specificity_2
ctd[[2]]$annot <- cell_level_2
saveRDS(ctd,"ctd_cds_cleaned.RDS")
# ctd <- readRDS("ctd_cds_cleaned.RDS")
#################################



