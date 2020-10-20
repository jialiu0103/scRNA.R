library(MuSiC)
library(xbioc)
require(Biobase)
library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Biobase)
library(reshape2)
library(Matrix)
library(tidyverse)

options(error=recover)

#read in data
dat <- readMM("~/filtered_count_matrix.mtx")
genenames <- readLines("~/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("~/filtered_column_metadata.txt")
indgroup <- read_excel("indgroup.xlsx")
bulk_data=read.table(file='~/ROSMAP_RNAseq_FPKM_gene.tsv', sep = '\t', header = TRUE)
bulkmeta=read.csv('~/Clinical/ROSMAP_IDkey.csv')
pid=as.data.frame(filtered.colMetadata$projid)
colnames(pid)='projid'
ad=inner_join(pid,indgroup, by = 'projid')
data=dat
data=as.matrix(data)
colnames(data)=filtered.colMetadata$broad.cell.type
rownames(data)=genenames

#folowing dataset is sc datset contains only control cells
my_data=data[,which(ad$ad==0)]
colnames(my_data)=filtered.colMetadata$TAG[ad$ad==0]
sccelltype=filtered.colMetadata$broad.cell.type[ad$ad==0]
sccelltag=filtered.colMetadata$TAG[ad$ad==0]

#following dataset is dataset with all cells
my_data=data
colnames(my_data)=filtered.colMetadata$TAG
sccelltype=filtered.colMetadata$broad.cell.type
sccelltag=filtered.colMetadata$TAG

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
#process bulk metadata
# bulkmt=bulkmeta[!is.na(bulkmeta$rnaseq_id), ]

#convert dataframe to expressionset 
aa=data_frame(sccelltype=sccelltype,sccelltag=sccelltag)
aa=as.data.frame(aa)
colnames(my_data)=sccelltag
rownames(aa)=sccelltag
##
metadata <- data.frame(labelDescription= c("Cell Type Name","Sample ID" ), row.names=c("sccelltype", "sccelltag"))
SC.eset = ExpressionSet(assayData = data.matrix(my_data), phenoData =  new("AnnotatedDataFrame", data = aa, varMetadata = metadata) )
sc_data <- SC.eset
##
### ANOTHER WAY TO CREATE EXPRESSIONSET, COME ACROSS THE BUG
# aa=as.data.frame(filtered.colMetadata)
# pdata <-AnnotatedDataFrame(aa)
# sc_data <- new("ExpressionSet", exprs=as.matrix(data))
# sc_data@phenoData<-pdata
###
# sc_data@phenoData<-pdata
bulk_data<-new("ExpressionSet", exprs=as.matrix(bulk))

###Bulk Tissue Cell Type Estimation
# Estimate cell type proportions
Est.prop.bulk = music_prop(bulk.eset = bulk_data, sc.eset = sc_data, clusters ='sccelltype',
                           samples = 'sccelltag', verbose = F)
names(Est.prop.bulk)
saveRDS(Est.prop.bulk,'estprop_music')

# Jitter plot of estimated cell type proportions
jitter.fig = Jitter_Est(list(data.matrix(Est.prop.bulk$Est.prop.weighted),
                             data.matrix(Est.prop.bulk$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')


###Estimation of cell type proportions with pre-grouping of cell types
#To deal with collinearity, MuSiC employs a tree-guided procedure 
#that recursively zooms in on closely related cell types.

##Clustering single cell data
# Produce the first step information
sub.basis = music_basis(sc_data, clusters = 'sccelltype', samples ='sccelltag')
# Plot the dendrogram of design matrix and cross-subject mean of realtive abundance
par(mfrow = c(1, 2))
d <- dist(t(log(sub.basis$Disgn.mtx + 1e-6)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d <- dist(t(log(sub.basis$M.theta + 1e-8)), method = "euclidean")
# Hierarchical clustering using Complete Linkage
# hc2 <- hclust(d, method = "complete" )
hc2 <- hclust(d, method = "complete")
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1, main = 'Cluster log(Mean of RA)')
saveRDS(sub.basis,'sub.basis2_music')

##Bulk Tissue Cell Type Estimation with Pre-grouping of Cell Types
clusters.type = list(C1 = 'Mic', C2 = c('Ex', 'In', 'Oli','Opc', 'Ast'), C3 = 'End', C4 = 'Per')
cl.type = as.character(sccelltype)
for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
pData(sc_data)$clusterType = factor(cl.type, levels = names(clusters.type))


#select DEGs for cluster2
ind <- clusterType=='C2'
m <- exprs(sc_data)
data_c2=m[,ind]

c2 <- CreateSeuratObject(counts = data_c2, min.cells = 3, min.features = 200)
c2 <- NormalizeData(c2, normalization.method = "LogNormalize", scale.factor = 10000)
c2 <- FindVariableFeatures(c2, selection.method = "vst", nfeatures = 2000)
c2 <- ScaleData(c2, features = all.genes)
c2 <- RunPCA(c2, features = VariableFeatures(object = c2))
c2 <- FindNeighbors(c2, dims = 1:10)
c2 <- FindClusters(c2, resolution = 0.5)
c2 <- RunUMAP(c2, dims = 1:10)
c2.markers <- FindAllMarkers(c2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- c2.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
topdeg=unique(top20$gene)
DEG=list(C1='',C2=topdeg,C3='',C4='')

#estimate proportion
# Est.pregrouping.bulk = music_prop.cluster(bulk.eset = bulk_data, sc.eset = sc_data, group.markers = DEG, clusters = 'sccelltype', group =  'clusterType', samples ='sccelltag', clusters.type =clusters.type)
# sub.basis = music_basis(sc_data, clusters = 'clusterType', samples = 'sccelltag',select.ct = aa)#aa=names(clusters.type)

Est.pregrouping.bulk = music_prop.cluster(bulk.eset = bulk_data, sc.eset = sc_data, 
                                          clusters ='sccelltype',group.markers = DEG, 
                                          groups =  'clusterType',samples ='sccelltag',
                                          clusters.type =clusters.type)

jitter.fig = Jitter_Est(Est.pregrouping.bulk,method.name = c('MuSiC', 'NNLS'),title = 'Jitter plot of Est Pregrouping Proportions')
saveRDS(Est.pregrouping.bulk,'callcell_music')

saveRDS(Est.pregrouping.bulk,'control_music')
