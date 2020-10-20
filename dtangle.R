###dtangle
#install dtangle
# packageurl <-'https://cran.r-project.org/bin/macosx/contrib/4.0/dtangle_2.0.9.tgz'
# install.packages(packageurl, repos=NULL, type="source")
require(dtangle)
library(limma)
library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Biobase)
library(reshape2)
library(Matrix)

options(error=recover)

#read in data
dat <- readMM("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_count_matrix.mtx")
genenames <- readLines("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_column_metadata.txt")
indgroup <- read_excel("indgroup.xlsx")
bulk_data=read.table(file='~/ROSMAP_RNAseq_FPKM_gene.tsv', sep = '\t', header = TRUE)
bulkmeta=read.csv('/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/Clinical/ROSMAP_IDkey.csv')
data=dat
colnames(data)=filtered.colMetadata$TAG
rownames(data)=genenames
bulk=bulk_data
bulk=bulk[,-1]
rownames(bulk)=bulk[,1]
bulk=bulk[,-1]
#convert gene id to gene symbol
keytypes(org.Hs.eg.db)
geneid=rownames(bulk)
geneid_modified <- sapply(strsplit(geneid,"\\."), function(x) x[1])
gene_symbol <- bitr(geneid_modified, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
# make a dataframe contains gene expression, probeid and annotation
ind <- match(geneid_modified, gene_symbol$ENSEMBL)
bulk$Gene_Symbol <- gene_symbol$SYMBOL[ind]
# drop na in gene symbol
bulk=bulk[!is.na(bulk$Gene_Symbol), ]
#drop duplicated gene symbols in bulk data(60 duplicated)
bulk<-bulk[!duplicated(bulk$Gene_Symbol), ] 
rownames(bulk)=bulk$Gene_Symbol
bulk=bulk[,-641]

saveRDS(bulk,'bulk')
#DEG analysis
# bulk <- backgroundCorrect(bulk, method='normexp')

#find the genes that are common to both datasets and join the datasets using them
commongenes <- intersect (rownames(bulk), rownames(data))
bulk <- bulk[pmatch(commongenes, rownames(bulk)), ]
sce <- data[pmatch(commongenes, rownames(data)), ]
y <- cbind(sce, bulk)
#apply quantile normalization in order to ensure that they are indeed comparable
y <- normalizeBetweenArrays(y)
y <- t(y)

#create an object that records can the pure samples of each cell type
all_cell_type <- unique(filtered.colMetadata$broad.cell.type)
pure_samples <- lapply(1:length(all_cell_type), function(i) {
  which(filtered.colMetadata$broad.cell.type == all_cell_type[i])
})
names(pure_samples) = all_cell_type

#find marker genes
marker_list = find_markers(y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')
#choose top 10% of all marker genes for each type
q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
n_markers

#perform dtangle
marks = marker_list$L
dc <- dtangle(y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'microarray-gene', markers = marks)

#extract the proportion estimates for mixture dataset
final_est <- dc$estimates[(dim(sce)[2]+1):dim(y)[1],]
colnames(final_est) <-  all_cell_type
head(final_est)

#plot the proportion estimates
library(ggplot2)
library(reshape2)
plot_data <- melt(final_est)
colnames(plot_data) <- c("Sample", "Cell Type", "Proportion") 
plot_data$Proportion <- as.numeric(plot_data$Proportion)
dtanglere=ggplot(plot_data, aes(x = `Cell Type`, y=Proportion))+geom_violin(aes(fill = `Cell Type`)) + geom_jitter(height = 0, width = 0.1)



