# remove terms not related to 'phenotypic abnormality'
hpo_clinical_mod <- hpo$children["HP:0012823"]
hpo_inheritance <- hpo$children["HP:0000005"]

descendants = list()
for (b in hpo_clinical_mod[[1]]) {
  descendants = ontologyIndex::get_descendants(hpo, b)
  for(d in descendants){
    all_results_merged2 <- all_results_merged2[!grepl(d, all_results_merged2$HPO_id),]
  }
}

descendants = list()
for (b in hpo_inheritance[[1]]) {
  descendants = ontologyIndex::get_descendants(hpo, b)
  for(d in descendants){
    all_results_merged2 <- all_results_merged2[!grepl(d, all_results_merged2$HPO_id),]
  }
}

# get children terms of 'phenotypic abnormality'
hpo_branches = hpo$children["HP:0000118"]
hpo_branches_names = c()
for (b in hpo_branches) {
  hpo_branches_names = c(hpo_branches_names, hpo$name[b])
}
main_branch_df = data.frame(hpo_id = hpo_branches, phenotype = hpo_branches_names)
descendants = list()
for (b in hpo_branches[[1]]) {
  descendants[[hpo$name[b]]] = ontologyIndex::get_descendants(hpo, 
                                                              b)
}

if (!"branch" %in% colnames(all_results_merged2)) {
  phenotypes <- unique(all_results_merged2$Phenotype)
  
  # remove phenotypes related to age of onset 
  phenotypes <- gsub("^Adult onset$", "", phenotypes)
  phenotypes <- gsub("^Antenatal onset$", "", phenotypes)
  phenotypes <- gsub("^Pediatric onset$", "", phenotypes)
  phenotypes <- gsub("^Congenital onset$", "", phenotypes)
  phenotypes <- gsub("^Neonatal onset$", "", phenotypes)
  phenotypes <- gsub("^Puerpural onset$", "", phenotypes)
  
  phenotypes <- phenotypes[phenotypes != ""]
  
  all_results_merged2$branch = NA
  pheno_branches = c()
  for (p in phenotypes) {
    
    hpo_id = unique(all_results_merged2$HPO_id[all_results_merged2$Phenotype==p])
    if (!is.na(hpo_id)) {
      id = hpo_id
      for (i in seq(1, length(descendants))) {
        if (id %in% descendants[[i]]) {
          pheno_branches[p] = names(descendants[i])
        }
      }
    }
  }
  
  
  
  for (p in unique(all_results_merged2$Phenotype)) {
    all_results_merged2$branch[all_results_merged2$Phenotype == 
                                p] = pheno_branches[p]
  }
  
  sum <- sum(is.na(all_results_merged2$branch))
  
  
  if(sum>0){
    missing_phenos <- phenotypes[phenotypes %in% unique(all_results_merged2$Phenotype[is.na(all_results_merged2$branch)])]
    missing_phenos <- sub("level", "concentration", missing_phenos)
    meta <- HPOExplorer::hpo_meta
     missing_branches_df <- data.frame(Phenotype=NA, branch=NA)
    for(p in missing_phenos){
      column <- grep(p, meta)
      
      if(length(column)>0 & length(column)<2){
        
        index <- grep(p, meta[[column]])
        alt_id <- meta$HPO_ID[index]

      for (i in seq(1, length(descendants))) {
        if (alt_id %in% descendants[[i]]) {
          missing_branches <- data.frame(Phenotype=p, branch=names(descendants[i]))
          missing_branches_df <- rbind(missing_branches_df, missing_branches)
          }
        }
      }
  }
  
missing_branches_df <- missing_branches_df[-1,]
unique_key1 <- paste(missing_branches_df$Phenotype)
unique_key2 <- paste(all_results_merged2$Phenotype)
inds <- is.na(all_results_merged2$branch)
all_results_merged2$branch[inds] <- missing_branches_df$branch[match(unique_key2[inds], unique_key1)]
    
}
}

# manual fix to branch names
test <- test[!grepl("Gonosomal inheritance", test$Phenotype),]
test <- test[!grepl("Heterogeneous", test$Phenotype),]
test$branch[test$Phenotype=="Abnormal circulating insulin level"] <- "Abnormality of the endocrine system"
test$branch[test$Phenotype=="Abnormality of calvarial morphology"] <- "Abnormality of head or neck"
test$branch[test$Phenotype=="Abnormality of the metacarpal bones"] <- "Abnormality of limbs"
test$branch[test$Phenotype=="Abnormality of urine glucose concentration"] <- "Abnormality of the genitourinary system"
test$branch[test$Phenotype=="Decreased circulating aldosterone level"] <- "Abnormality of the endocrine system"
test$branch[test$Phenotype=="Abnormal echocardiogram"] <- "Abnormality of the cardiovascular system"
test$branch[test$Phenotype=="Abnormality of malar bones"] <- "Abnormality of head or neck"
test$branch[test$Phenotype=="Abnormality of the periosteum"] <- "Abnormality of the musculoskeletal system"
test$branch[test$Phenotype=="Distal sensory loss of all modalities"] <- "Abnormality of the nervous system"
test$branch[test$Phenotype=="Bilateral external ear deformity"] <- "Abnormality of the ear"
test$branch[test$Phenotype=="Glucose intolerance"] <- "Abnormality of metabolism/homeostasis"
test$branch[test$Phenotype=="Iron accumulation in globus pallidus"] <- "Abnormality of the nervous system"
test$branch[test$Phenotype=="Marked muscular hypertrophy"] <- "Abnormality of the musculoskeletal system "
test$branch[test$Phenotype=="Marked muscular hypertrophy"] <- "Abnormality of the musculoskeletal system"
test$branch[test$Phenotype=="Decreased pulmonary function"] <- "Abnormality of the respiratory system"
test$branch[test$Phenotype=="Large beaked nose"] <- "Abnormality of the nose"
test$branch[test$Phenotype=="Dysphasia"] <- "Abnormality of the nervous system"
test$branch[test$Phenotype=="Loss of ability to walk"] <- "Abnormality of the nervous system"
test$branch[test$Phenotype=="Peripheral arterial stenosis"] <- "Abnormality of the cardiovascular system"
test$branch[test$Phenotype=="Sparse and thin eyebrow"] <- "Abnormality of head or neck"
