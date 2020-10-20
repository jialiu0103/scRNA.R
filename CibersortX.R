library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Biobase)
library(reshape2)
library(Matrix)
library(tidyverse)
library(xbioc)
require(Biobase)
library(dplyr)

dat <- readMM("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_count_matrix.mtx")
genenames <- readLines("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_column_metadata.txt")
indgroup <- read_excel("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/mathys_pathologygroup.xlsx")

bulk_data=read.table(file='~/ROSMAP_RNAseq_FPKM_gene.tsv', sep = '\t', header = TRUE)
bulkmeta=read.csv('~/Clinical/ROSMAP_IDkey.csv')
pid=as.data.frame(filtered.colMetadata$projid)
colnames(pid)='projid'
ad=inner_join(pid,indgroup, by = 'projid')
data=dat
data=as.matrix(data)
colnames(data)=filtered.colMetadata$broad.cell.type
rownames(data)=genenames

#convert gene id in row bulk dta to gene symbel, prepare bulk data
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


data1=rbind(as.character(filtered.colMetadata$broad.cell.type),data)
gn=c('names',genenames)
data2=cbind(gn,data1)


#M1: random select 10% cells
# write.table(data2, file = "scdt.txt", sep = "\t")
# my_data <- read.delim("scdt.txt")
# newsc=sample(my_data, 7100)
# write.table(newsc, file = "newsc.txt", sep = "\t")

#M2: select 30% cells based on cell type proportion 
pid=as.data.frame(filtered.colMetadata$projid)
colnames(pid)='projid'
ad=inner_join(pid,indgroup, by = 'projid')
dat <- readMM("~/filtered_count_matrix.mtx")
data=dat
data=as.matrix(data)
colnames(data)=filtered.colMetadata$broad.cell.type
rownames(data)=genenames
my_data=data[,which(ad$ad==0)]
gsum=rowSums(my_data)
my_dt=my_data[rowSums(my_data)>20,]
new_genename=rownames(my_dt)

ast=my_dt[,colnames(my_dt)=='Ast']
ast=as.data.frame(ast)
ast_rm=sample(ast,416)

mic=my_dt[,colnames(my_dt)=='Mic']
mic=as.data.frame(mic)
mic_rm=sample(mic,257)

ex=my_dt[,colnames(my_dt)=='Ex']
ex=as.data.frame(ex)
ex_rm=sample(ex,4559)

In=my_dt[,colnames(my_dt)=='In']
In=as.data.frame(In)
In_rm=sample(In,1286)

end=my_dt[,colnames(my_dt)=='End']
end=as.data.frame(end)
end_rm=sample(end,16)

opc=my_dt[,colnames(my_dt)=='Opc']
opc=as.data.frame(opc)
opc_rm=sample(opc,356)

oli=my_dt[,colnames(my_dt)=='Oli']
oli=as.data.frame(oli)
oli_rm=sample(oli,2453)

per=my_dt[,colnames(my_dt)=='Per']
per=as.data.frame(per)
per_rm=sample(per,24)

ct=c(rep('ast',416),rep('mic',257),rep('In',1286),rep('ex',4559),rep('per',24),rep('oli',2453),rep('opc',356),rep('end',16))
afterdata=cbind(ast_rm,mic_rm,In_rm,ex_rm,per_rm,oli_rm,opc_rm,end_rm)

data1=rbind(as.character(ct),afterdata)
gn=c('names',new_genename)
data2=cbind(gn,data1)
data2=as.data.frame(data2)
write.table(data2, file = "controlsc_v1.txt", sep = "\t",row.names = FALSE,col.names = FALSE)
#my_data <- read.delim("controlsc_v1.txt")

#M3: filter out low expressed genes, keep all cells
# pid=as.data.frame(filtered.colMetadata$projid)
# colnames(pid)='projid'
# ad=inner_join(pid,indgroup, by = 'projid')
# dat <- readMM("~/filtered_count_matrix.mtx")
# data=dat
# data=as.matrix(data)
# colnames(data)=filtered.colMetadata$broad.cell.type
# rownames(data)=genenames
# my_data=data[,which(ad$ad==0)]
# #my_data contains only control samples
# my_dt=my_data[rowSums(my_data)>20,]
# new_genename=rownames(my_dt)
# bb=colSums(my_dt)
# v2_dt=my_dt[,bb>4000]
# v2_dt=my_dt[,bb<9000]
# v2_dt=v2_dt[rowSums(v2_dt)>20,]
# v2_dt=v2_dt[,colSums(v2_dt)>4000]
# new_genename=rownames(v2_dt)
# 
# #make up cell types which have low cell number 
# ct=c(rep('mic',25),rep('per',25),rep('end',25))
# mic=my_data[,colnames(my_data)=='Mic']
# mic=as.data.frame(mic)
# mic_rm=sample(mic,25)
# per=my_data[,colnames(my_data)=='Per']
# per=as.data.frame(per)
# per_rm=sample(per,25)
# end=my_data[,colnames(my_data)=='End']
# end=as.data.frame(end)
# end_rm=sample(end,25)
# afterdata=cbind(mic_rm,per_rm,end_rm)
# data1=rbind(as.character(ct),afterdata)
# dt1=data1[new_genename,]
# mydt=cbind(v2_dt,dt1)
# 
# data1=rbind(c(as.character(colnames(v2_dt)),ct),mydt)
# 
# gn=c('names',new_genename)
# data2=cbind(gn,data1)
# data2=as.data.frame(data2)
# write.table(data2, file = "ctsc_v2.txt", sep = "\t",row.names = FALSE,col.names = FALSE)
#my_data <- read.delim("controlsc_v2.txt")


#M4 sum cells into individual level
data=dat
colnames(data)=filtered.colMetadata$broad.cell.type
data=rbind(filtered.colMetadata$projid,data)
getctind=function(ctname){
  ct_data=data[,colnames(data)==ctname]
  colnames(ct_data)=ct_data[1,]
  ct_data=ct_data[-1,]
  #here we have 17926 genes in ct_data
  ct_ind=data.frame(projectid=rep(999,17926))
  #aggregating all cells within each individual 
  mypid=indgroup$projid
  for (i in mypid){
    if (i %in% as.numeric(colnames(ct_data))){
    ind2=ct_data[,colnames(ct_data)==i]
    ind2=as.data.frame(ind2)
    if (dim(ind2)[2]>1){
    dat_ind2=apply(ind2,1,sum)}
    else{dat_ind2=dat_ind2}
    dat_ind2=as.data.frame(dat_ind2)
    colnames(dat_ind2)=i
    ct_ind=cbind(ct_ind,dat_ind2)
    }
  }
  ct_ind=ct_ind[,-1]
  return(ct_ind)
}

ast_ind=getctind('Ast')
mic_ind=getctind('Mic')
opc_ind=getctind('Opc')
in_ind=getctind('In')
oli_ind=getctind('Oli')
per_ind=getctind('Per')
end_ind=getctind('End')
ex_ind=getctind('Ex')

ct_ind=cbind(ast_ind,mic_ind,opc_ind,in_ind,oli_ind,per_ind,end_ind,ex_ind)
ct=c(rep('ast',48),rep('mic',48),rep('opc',48),rep('In',48),rep('oli',48),rep('per',37),rep('end',37),rep('ex',48))

data1=rbind(as.character(ct),ct_ind)
gn=c('names',rownames(ct_ind))
data2=cbind(gn,data1)
data2=as.data.frame(data2)
write.table(data2, file = "indsc.txt", sep = "\t",row.names = FALSE,col.names = FALSE)

#visualization
bmode_randomv1 <- read_csv("Desktop/cibersortxresult/bmode_randomv1.csv")
smode_randomv1 <- read_csv("Desktop/cibersortxresult/smode_randomv1.csv")
v1_bmode <- read_csv("Desktop/control_ctx/v1_bmode.csv")
v1_smode <- read_csv("Desktop/control_ctx/v1_smode.csv")
v2_smode <- read_csv("Desktop/control_ctx/v2_smode.csv")
v2_bmode <- read_csv("Desktop/control_ctx/v2_bmode.csv")
indBmode <- read_csv("Desktop/cibersortxresult/indBmode.csv")
indSmode <- read_csv("Desktop/cibersortxresult/indSmode.csv")
#bmode=read_csv('bmode_randomv1.csv')
#bmode=bmode_randomv1
#bmode=smode_randomv1
#bmode=v1_bmode
#bmode=v1_smode
#bmode=v2_bmode
#bmode=indBmode
#bmode=indSmode
bmode=as.data.frame(bmode)
library(ggplot2)
library(reshape2)
plot_data <- melt(bmode[,1:9])
colnames(plot_data) <- c("Sample", "Cell Type", "Proportion") 
plot_data$Proportion <- as.numeric(plot_data$Proportion)
ctx_bmode_re=ggplot(plot_data, aes(x = `Cell Type`, y=Proportion))+geom_violin(aes(fill = `Cell Type`)) + geom_jitter(height = 0, width = 0.1)
