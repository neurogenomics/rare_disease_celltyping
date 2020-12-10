library(stringr)

annot <- read.csv('annotations_facs.csv')
f_names <- dir(path = "FACS")
levels(annot$cell_ontology_class) <-
  c(levels(annot$cell_ontology_class), "NA")
annot$cell_ontology_class[which(annot$cell_ontology_class == "")] <-
  "NA"
alpha = 0.05 #parameter: non-zero rate
for (i in 1:20) {
  tissue <- str_sub(f_names[i], end = -12)
  annot1 <- annot[which(annot$tissue == tissue), ]
  data1 <- read.csv(paste("FACS/", f_names[i], sep = ""))
  rownames(data1) <- data1$X
  data1 <- data1[, -1]
  data2 <- data1[, as.character(annot1$cell)] # count data
  data2 <-  t( t(data2) / colSums(data2) *1e+06) # cpm
  data1 <- log2(data2+1) # log expression
  
  cells <- unique(annot1$cell_ontology_class)
  x <- data.frame(data2[, 1:length(cells)])
  names(x) <- paste(cells, tissue, sep = "_")
  x_log <- x # mean log expression
  x_log_nz <- x # mean log with non-zeros
  x_cpm_nz <- x # mean cpm with non-zeros
  for (k in 1:length(cells)) {
    cid <- which(annot1$cell_ontology_class == cells[k])
    x[, k] <- rowMeans(data2[, cid]) #mean cpm
    x_log[,k] <- rowMeans(data1[, cid]) # mean log expr
    x_log_nz[,k] <- x_log[,k]
    x_cpm_nz[,k] <- x[,k]
    nz_rate <- rowMeans(data1[,cid] > 0)
    c_nz <- which(nz_rate > alpha) # genes with dropout rate < 0.95
    x_log_nz[c_nz,k] <- x_log_nz[c_nz,k] / nz_rate[c_nz]
    x_cpm_nz[c_nz,k] <- x_cpm_nz[c_nz,k] / nz_rate[c_nz]
  }
  if (i == 1) {
    y = x
    y_log = x_log
    y_nz = x_log_nz
    y_cpm_nz = x_cpm_nz
  } else{
    y = cbind(y, x)
    y_log = cbind(y_log, x_log)
    y_nz = cbind(y_nz, x_log_nz)
    y_cpm_nz = cbind(y_cpm_nz, x_cpm_nz)
  }
}
x <- y
normalised_meanExp = t(t(x) * (1 / colSums(x)))
specificity = normalised_meanExp / (apply(normalised_meanExp, 1, sum) +
                                      0.000000000001)

ctd <- vector("list", 1)
ctd[[1]]$log_mean_exp_nz <- y_nz
ctd[[1]]$log_mean_exp <- y_log
ctd[[1]]$mean_exp <- x
ctd[[1]]$mean_exp_nz <- y_cpm_nz
ctd[[1]]$specificity <- specificity
ctd[[1]]$annot <- annot$cell_ontology_class

save(ctd, file = "ctd_tm_all.rda")