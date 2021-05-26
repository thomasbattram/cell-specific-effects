annotation.package <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
library(annotation.package, character.only=T)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("missMethyl")
BiocManager::install("getMappedEntrezIDs")
library(missMethyl)
BioManager::install("org.Hs.eg.db")
library('org.Hs.eg.db')
## if the package is not installed then do this:
##   install.packages("BiocManager")
##   BiocManager::install(annotation.package)

data(list=annotation.package)
data(Locations)
annotation <- as.data.frame(Locations)
annotation <- annotation[which(annotation$chr %in% paste0("chr", 1:22)),]
## function for annotating a set of CpG site identifiers

annotation.sites <- function(sites) {
   annotation[match(sites, rownames(annotation)),]
}

# load in cpg lists
tca1 <- read.csv("celltype1top100significant_TCA.csv")[,1]
# get locations
tca1_anno <- annotation.sites(tca1)
# get gene ids 
gene_ids <- getMappedEntrezIDs(tca1, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
tca1_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')
write(tca1_gene_names, file="tca1_genenames.txt", sep=",")

tca2 <- read.csv("celltype2top100significant_TCA.csv")[,1]
tca2_anno <- annotation.sites(tca2)
# get locations
tca2_anno <- annotation.sites(tca2)
# get gene ids 
gene_ids <- getMappedEntrezIDs(tca2, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
tca2_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')


tca3 <- read.csv("celltype3top100significant_TCA.csv")[,1]
tca3_anno <- annotation.sites(tca3)
# get locations
tca3_anno <- annotation.sites(tca3)
# get gene ids 
gene_ids <- getMappedEntrezIDs(tca3, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
tca3_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')


tca4 <- read.csv("celltype4top100significant_TCA.csv")[,1]
tca4_anno <- annotation.sites(tca4)
# get locations
tca4_anno <- annotation.sites(tca4)
# get gene ids 
gene_ids <- getMappedEntrezIDs(tca4, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
tca4_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')


tca5 <- read.csv("celltype5top100significant_TCA.csv")[,1]
tca5_anno <- annotation.sites(tca5)
# get locations
tca5_anno <- annotation.sites(tca5)
# get gene ids 
gene_ids <- getMappedEntrezIDs(tca5, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
tca5_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')


tca6 <- read.csv("celltype6top100significant_TCA.csv")[,1]
tca16_anno <- annotation.sites(tca6)
# get locations
tca6_anno <- annotation.sites(tca6)
# get gene ids 
gene_ids <- getMappedEntrezIDs(tca6, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
tca6_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')



hire1 <- read.csv("celltype1top100significant_HIRE.csv")[,1]
# get locations
hire1_anno <- annotation.sites(hire1)
# get gene ids 
gene_ids <- getMappedEntrezIDs(hire1, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids <- c(gene_ids$sig.eg)
# create gene list
hire1_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')
write(hire1_gene_names, file="hire1_genenames.txt", sep=",")



hire2 <- read.csv("celltype2top100significant_HIRE.csv")[,1]
hire2_anno <- annotation.sites(hire2)
# get gene ids 
gene_ids <- getMappedEntrezIDs(hire2, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
hire2_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')
write(hire2_gene_names, file="hire2_genenames.txt", sep=",")


hire3 <- read.csv("celltype3top100significant_HIRE.csv")[,1]
hire3_anno <- annotation.sites(hire3)
# get gene ids 
gene_ids <- getMappedEntrezIDs(hire3, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
hire3_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')
write(hire3_gene_names, file="hire3_genenames.txt", sep=",")


hire4 <- read.csv("celltype4top100significant_HIRE.csv")[,1]
hire4_anno <- annotation.sites(hire4)
# get gene ids 
gene_ids <- getMappedEntrezIDs(hire4, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
hire4_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')
write(hire4_gene_names, file="hire4_genenames.txt", sep=",")


hire5 <- read.csv("celltype5top100significant_HIRE.csv")[,1]
hire5_anno <- annotation.sites(hire5)
# get gene ids 
gene_ids <- getMappedEntrezIDs(hire5, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
hire5_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')
write(hire5_gene_names, file="hire5_genenames.txt", sep=",")


hire6 <- read.csv("celltype6top100significant_HIRE.csv")[,1]
hire6_anno <- annotation.sites(hire6)
# get gene ids 
gene_ids <- getMappedEntrezIDs(hire6, all.cpg = NULL, anno = NULL, array.type="450K")
gene_ids = c(gene_ids$sig.eg)
# create gene list
hire6_gene_names <- mapIds(org.Hs.eg.db, gene_ids, 'SYMBOL', 'ENTREZID')
write(hire6_gene_names, file="hire6_genenames.txt", sep=",")


write(hire1_anno, file="hire1.bed", sep='\t')
write(hire2_anno, file="hire2.bed", sep='\t')
write(hire3_anno, file="hire3.bed", sep='\t')
write(hire4_anno, file="hire4.bed", sep='\t')
write(hire5_anno, file="hire5.bed", sep='\t')
write(hire6_anno, file="hire6.bed", sep='\t')
write(tca1_anno, file="tca1.bed", sep='\t')
write(tca2_anno, file="tca2.bed", sep='\t')
write(tca3_anno, file="tca3.bed", sep='\t')
write(tca4_anno, file="tca4.bed", sep='\t')
write(tca5_anno, file="tca5.bed", sep='\t')
write(tca6_anno, file="tca6.bed", sep='\t')




