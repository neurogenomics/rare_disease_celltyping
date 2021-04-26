
gen_results = function (CTD_file = "ctd_l1l2_nz.rda",
                       phenotype_to_genes_txt_file = "phenotype_to_genes.txt",
                       output_path = "Output",
                       reps = 100000,
                       level = 1,
                       mc.cores = 24,
                       p_adj_method = "BH",
                       output_merged_rda_filename = "Results_merged.rda",
                       output_unmerged_rda_filename = "Results_unmerged.rda") {
  library(EWCE)
  library(limma)
  library(rvest)
  library(rlist)
  library(stringr)
  library(parallel)
  load(file= CTD_file) #Load CTD

  if (!file.exists(output_path)) {dir.create(output_path)}
  if (!file.exists(phenotype_to_genes_txt_file)) {
    download.file("http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
                  phenotype_to_genes_txt_file)}

  #Some gene symbols need correcting, this is achieved in EWCE with the following:
  if(!file.exists("MRK_List2.rpt")){
    download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", destfile="MRK_List2.rpt")
  }
  ctd[[1]]$specificity <- fix.bad.mgi.symbols(ctd[[1]]$specificity, mrk_file_path="MRK_List2.rpt")
  ctd[[1]]$mean_exp <- fix.bad.mgi.symbols(ctd[[1]]$mean_exp, mrk_file_path="MRK_List2.rpt")
  genedata = read.delim(phenotype_to_genes_txt_file, skip = 1) #skip=1 ? #HPO annotation

  # CHANGE: shortening gene data for speeding up test runs
  #genedata = genedata[1:12000,]

  colnames(genedata) = c("ID", "Phenotype", "EntrezID", "Gene", "Additional", "Source", "LinkID")
  desiredPhenotypes <- as.list(unique(genedata$Phenotype)) #Just the whole list of unique phenotypes in the "phenotype_to_genes.txt" file.
  data("mouse_to_human_homologs") #Part of EWCE
  m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol", "MGI.symbol")])
  PullGenes = function(phenotype, genedata){
    geneListLoad = subset(genedata, Phenotype == phenotype)$Gene
    return(geneListLoad)
  }
  RunEWCE = function(phenotype){
    ctd = ctd
    reps = reps
    level = level
    geneList = PullGenes(phenotype = phenotype, genedata = genedata)
    y=1
    mouse.hits = c("")
    for(i in geneList){
      if(i %in% m2h$HGNC.symbol){
        mouse.hits[y] = i
        y = y + 1
      }
    }
    mouse.bg = unique(m2h$MGI.symbol)
    mouse.hits = unique(m2h[m2h$HGNC.symbol %in% mouse.hits, "MGI.symbol"])

    #EWCE only works using lists of 3 or more genes.
    if(length(mouse.hits) > 3){
      full_results = bootstrap.enrichment.test(sct_data = ctd, hits = mouse.hits, bg = mouse.bg, reps = reps, annotLevel = level)
      full_results$results = cbind(full_results$results, list = phenotype)
      return(full_results$results)
    }
  }

  EWCEWrapper = function(phenotype){
    #Replaced all forward slashes in phenotype names with "-SLASH-"
    phenotype_nospace = gsub("/", "-SLASH-", phenotype, fixed = TRUE)
    phenotype_nospace = gsub(" ", "-", phenotype_nospace, fixed = TRUE)

    # CHANGE - may need to replace "." symbols? (error occured with phenotype "EEG > 3.5 hz")
    phenotype_nospace = gsub(".","-DOT-", phenotype_nospace, fixed = TRUE)

    filename = paste(phenotype_nospace, ".rda", sep = "")
    full_extension = paste0(output_path,"/", filename)
    if(!file.exists(full_extension)){
      result = RunEWCE(phenotype)
      assign(paste(phenotype_nospace, "tm_ewce", sep = "_"), result)
      resultname <- paste(phenotype_nospace, "tm_ewce", sep = "_")
      save(list = resultname, file = full_extension)
      return(result)
    }
    else{
      result = load(file = full_extension)
      return(result)
    }
  }

  all_results = mclapply(desiredPhenotypes, FUN = EWCEWrapper, mc.cores = mc.cores)
  save(all_results, file = output_unmerged_rda_filename)
  #Because EWCE only allows gene lists of 3 or more, some phenotypes were skipped over in the final results list and instead saved as just "1". This loop simply pulls out all the results for phenotypes that did work and stores them in a new list before saving, to get around this error.
  #"list" is the format for results that came out fine.
  all_results_two = c("")
  j = 1
  for(i in 1:length(all_results)){
    if(((typeof(all_results[[i]]) == "list")) == TRUE){
      all_results_two[j] = all_results[i]
      j = j + 1
    }
  }
  all_results_merged = data.table::rbindlist(all_results_two)
  all_results_merged$q = p.adjust(all_results_merged$p, method = p_adj_method) #Went with BH for now
  print("Debug: saving merged results")
  save(all_results_merged, file = output_merged_rda_filename)

}
