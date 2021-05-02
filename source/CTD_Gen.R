#Full code to generate ctd file containing level 1 and 2 cell annotations for the Tabula Muris dataset.
#Adapted from initial code written by Zijing Liu.


gen_ctd = function (output_ctd_file_name = "data/ctd_l1l2_nz.rda",
                    alpha = 0.05 ,
                    FACS_URL = "https://ndownloader.figshare.com/files/10700143" ,
                    annot_URL = "https://ndownloader.figshare.com/files/13088129" ,
                    L1_classifications_xlsx_path = "data/tm_level1classifications_full.xlsx"   ) {



  library(stringr)
  library(xlsx)
  library(readxl)


  #CHANGE: Download FACS data (if not already present) (https://tabula-muris.ds.czbiohub.org/)
  if (!file.exists(file.path("data/FACS"))){ #getw(), <- ?
    download.file(FACS_URL, "data/FACS.zip", mode = "wb")
    facs = unzip("data/FACS.zip",  exdir=paste0(getwd(), "/data")); unlink("__MACOSX", recursive = TRUE)  #to delete an extra file in the tm link
    }
  if (!file.exists("data/annotations_facs.csv")){
    download.file(annot_URL, "data/annotations_facs.csv")}

  f_names <- dir(path = "data/FACS") #Requires filepath "FACS" containing all tm FACS data
  annot <- read.csv("data/annotations_facs.csv") #Found on the original upload for the Tabula Muris data


  levels(annot$cell_ontology_class) <- c(levels(annot$cell_ontology_class), "NA")
  annot$cell_ontology_class[which(annot$cell_ontology_class == "")] <- "NA" #Unlabeled cells -> NA

  ###Generating level 2 data - remains unchanged from original code (written by Zijing Liu) aside from some label changes.

  alpha = alpha
  print("Generating L2 Data")
  for(i in 1:length(f_names)){
    cat("\r",paste("Files remaining:", length(f_names)-i, "    ")) # COUNTDOWN
    tissue <- str_sub(f_names[i], end = -12)
    annot1 <- annot[which(annot$tissue == tissue), ]
    data1 <- read.csv(paste("data/FACS/", f_names[i], sep = ""))
    rownames(data1) <- data1$X
    data1 <- data1[, -1]
    data2 <- data1[, as.character(annot1$cell)]
    data2 <- t(t(data2) / colSums(data2) * 1e+06)
    data1 <- log2(data2+1)

    uniquel2 <- unique(annot1$cell_ontology_class)
    x <- data.frame(data2[, 1:length(uniquel2)])
    names(x) <- paste(uniquel2, tissue, sep = "_")
    x_log <- x
    x_log_nz <- x
    x_cpm_nz <- x

    for(k in 1:length(uniquel2)){
      cid <- which(annot1$cell_ontology_class == uniquel2[k])
      x[, k] <- rowMeans(data2[, cid])
      x_log[, k] <- rowMeans(data1[, cid])
      x_log_nz[, k] <- x_log[, k]
      x_cpm_nz[, k] <- x[, k]
      nz_rate <- rowMeans(data1[, cid] > 0)
      c_nz <- which(nz_rate > alpha)
      x_log_nz[c_nz, k] <- x_log_nz[c_nz, k] / nz_rate[c_nz]
      x_cpm_nz[c_nz, k] <- x_cpm_nz[c_nz, k] / nz_rate[c_nz]
    }
    if(i == 1){
      y = x
      y_log = x_log
      y_nz = x_log_nz
      y_cpm_nz = x_cpm_nz
    }
    else{
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

  ctd <- vector("list", 2)
  ctd[[2]]$log_mean_exp_nz <- y_nz
  ctd[[2]]$log_mean_exp <- y_log
  ctd[[2]]$mean_exp <- x
  ctd[[2]]$mean_exp_nz <- y_cpm_nz
  ctd[[2]]$specificity <- specificity
  ctd[[2]]$annot <- annot$cell_ontology_class



  ###Generating level 1 data - done by pulling mean row values from the level 2 data according to groupings of cells found in "tm_level1classifications_full" excel.

  Level1Data <- read_excel(L1_classifications_xlsx_path)

  #Generating annotation labels by making a new annot column representing each individual cell's level 1 class.
  annot$cell_ontology_class_with_tissue <- annot$cell_ontology_class
  for(i in 1:length(annot$cell_ontology_class_with_tissue)){
    annot$cell_ontology_class_with_tissue[i] <- paste(annot$cell_ontology_class_with_tissue[i], annot$tissue[i], sep = "_")
  }

  annot$cell_ontology_class_l1 <- annot$cell_ontology_class_with_tissue

  for(i in 1:length(Level1Data$Level1Classification)){
    currentCells <- Level1Data$L1Combinations[i]
    currentCells <- unlist(strsplit(currentCells, ", "))
    annot$cell_ontology_class_l1[which(annot$cell_ontology_class_l1 %in% currentCells)] <- Level1Data$Level1Classification[i]
  }

  totall1 <- unique(Level1Data$Level1Classification) #Full list of level 1 classifications

  #Dataframe of level 2 data formed on a matrix with columns = level 1 groups (totall1) and rows = gene names (rownames(x))

  a <- matrix(0, ncol = length(totall1), nrow = nrow(x)) #(a, b) is essentially a level 1 version of (x, y) (above)
  b <- data.frame(a)
  colnames(b) <- totall1
  rownames(b) <- rownames(x)
  b_log <- b
  b_log_nz <- b
  b_cpm_nz <- b

  for(k in 1:ncol(b)){
    currentCells <- unlist(strsplit(Level1Data$L1Combinations[which(Level1Data$Level1Classification == colnames(b)[k])], ", "))
    cid <- which(colnames(x) %in% currentCells)
    if(length(currentCells) < 2){
      b[, k] <- x[, cid]
      b_log[, k] <- y_log[, cid]
      b_log_nz[, k] <- y_nz[, cid]
      b_cpm_nz[, k] <- y_cpm_nz[, cid]
    }
    else{
      b[, k] <- rowMeans(x[, cid])
      b_log[, k] <- rowMeans(y_log[, cid])
      b_log_nz[, k] <- rowMeans(y_nz[, cid])
      b_cpm_nz[, k] <- rowMeans(y_cpm_nz[, cid])
    }
  }

  normalised_meanExp_l1 = t(t(b) * (1/colSums(b)))
  specificityl1 = normalised_meanExp_l1 / (apply(normalised_meanExp_l1, 1, sum) + 0.000000000001)

  #Level 1 data becomes ctd[[1]]
  ctd[[1]]$mean_exp <- b
  ctd[[1]]$mean_exp_nz <- b_cpm_nz
  ctd[[1]]$log_mean_exp <- b_log
  ctd[[1]]$log_mean_exp_nz <- b_log_nz
  ctd[[1]]$specificity <- specificityl1
  ctd[[1]]$annot <- annot$cell_ontology_class_l1

  save(ctd, file = output_ctd_file_name) #Full ctd containing levels 1 + 2 saved.

}
