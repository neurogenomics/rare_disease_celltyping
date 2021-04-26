# The HPO data was downloaded on
download_date = "26042021"
full_url = "https://ewce.page.link/gene2phenotype"
short_url = "https://ewce.page.link/gene2phenotype"
output_dir = "GeneLists-WholeBody_new"


# Load the data
library("data.table")
if (!file.exists("genes_to_phenotype.csv")) {
    download.file(full_url, "genes_to_phenotype.csv")
    allhpo = read.csv("genes_to_phenotype.csv")}
if (!file.exists(output_dir)) {dir.create(output_dir)}

allhpo = read.csv("genes_to_phenotype.csv")
unique_terms = unique(allhpo$PhenotypeName)

# Get all fields of interest
termstouse = unique(allhpo$PhenotypeName) #read.csv("/Users/natske/OneDrive - Imperial College London/Datasets/Human Phenotype Ontology/TermsToUse.csv",stringsAsFactors = FALSE)
print("Getting gene lists")
for(i in 1:length(termstouse)){
    cat("\rTerms Remaining:",length(termstouse)-i, "      ") #countdown
    termName = termstouse[i]
    termNameAdj = gsub("/"," or ",termName)
    genes = unique(allhpo[grep(termName,allhpo$PhenotypeName),]$HGNC.symbol)
    if(length(genes)>20){
        label1 = sprintf("HPO - %s",termName)
        label2 = sprintf("Downloaded on %s from %s and processed in process_hpo_to_genelists.r",download_date,short_url)
        fileData = c(label1,label2,genes)
        filePath = sprintf("%s/HPO - %s.txt",output_dir,termNameAdj)
        write.table(fileData,file=filePath,quote=FALSE,row.names=FALSE,col.names=FALSE)
        #print(length())
    }
}
#
# library("ontologyIndex")
# data(hpo)
# get_term_property(ontology=hpo, property="ancestors", term="HP:0001873", as_names=TRUE)
#
