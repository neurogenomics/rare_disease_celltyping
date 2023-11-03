f <- list.files("/rds/general/project/neurogenomics-lab/ephemeral/rare_disease/",
                # pattern = ".rds$",
                # recursive = TRUE,
                full.names = TRUE)
f2 <- sapply(f, function(x){list.files(x,pattern = ".rds", full.names = TRUE)})
length(unlist(f2))
dat <- (
  lapply(unlist(f2), readRDS) |> 
    data.table::rbindlist(use.names = TRUE, idcol = "batch")
)[,batch:=basename(batch)][,q:=stats::p.adjust(p,method = "bonf")]


unique(dat$batch)
missing_batches <- basename(f)[!basename(f) %in% unique(dat$batch)]
i <- sort(as.numeric(stringr::str_split(missing_batches,"_", simplify = TRUE)[,3]))
summary(i)
dat[q<0.05,]
