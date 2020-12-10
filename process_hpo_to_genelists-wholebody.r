# The HPO data was downloaded on
download_date = "30thSept2019"
full_url = "http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt"
short_url = "goo.gl/1qvFo7"

# Load the data
library("data.table")
allhpo = fread("http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt",fill=TRUE,sep="\t",skip=1)
colnames(allhpo) = c("HPOid","PhenotypeName","EntrezGeneID","HGNC.symbol")
unique_terms = unique(allhpo$PhenotypeName)

# Get all fields of interest
termstouse = unique(allhpo$PhenotypeName) #read.csv("/Users/natske/OneDrive - Imperial College London/Datasets/Human Phenotype Ontology/TermsToUse.csv",stringsAsFactors = FALSE)

for(i in 1:length(termstouse)){
    termName = termstouse[i]
    termNameAdj = gsub("/"," or ",termName)
    genes = unique(allhpo[grep(termName,allhpo$PhenotypeName),]$HGNC.symbol)
    if(length(genes)>20){
        label1 = sprintf("HPO - %s",termName)
        label2 = sprintf("Downloaded on %s from %s and processed in process_hpo_to_genelists.r",download_date,short_url)
        fileData = c(label1,label2,genes)
        filePath = sprintf("/Users/natske/OneDrive - Imperial College London/Datasets/Human Phenotype Ontology/GeneLists-WholeBody/HPO - %s.txt",termNameAdj)
        write.table(fileData,file=filePath,quote=FALSE,row.names=FALSE,col.names=FALSE)
        #print(length())
    }
}
# 
# library("ontologyIndex")
# data(hpo)
# get_term_property(ontology=hpo, property="ancestors", term="HP:0001873", as_names=TRUE)
# 
