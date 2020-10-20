library(tidyverse)
library(tidymodels)
library(dplyr)
library(Seurat)
library(patchwork)
library(scales)
library(viridis)
library(Matrix)
library(BRETIGEA)
library(edgeR)
library(plyr)

#input dataset and metadata
sc_data <- read.delim(file='scRNA_logCounts.tsv')
sc_meta <- read.delim(file='scRNA_metadata.tsv')
sc_data=as.matrix(sc_data)
sddt=sc_data
genenames=sddt[,1]
rownames(sddt)=genenames
sddt=sddt[,-1]
scdt=sddt
mode(scdt) = "numeric"
celltag=colnames(scdt)
#raw expression matrix was composed of 33,694 genes and 14,876 cells
#filtered matrix consisted of 10,850 genes and 13,214 cells

#as the input data is filtered scRNA data, create seurat obj without filtering
sc_dt <- CreateSeuratObject(counts = scdt, project = "scdt")
sc_dt

###preprocessing
#The matrix was normalized with a scale factor of 10,000 before FindVariableGenes was used to 
#FindVariableGenes was used to define variable genes with the parameters x.low.cutoff = 0.0125, 
#x.high.cutoff = 3, and y.cutoff = 0.5.(but seurat has updated to v3.1 and doesn't contain this function now)
#So I used FindVariableFeatures() instead
#ScaleData was used to center the gene expression.
sc_dt <- NormalizeData(sc_dt, normalization.method = "LogNormalize", scale.factor = 10000)#normalization
sc_dt <- FindVariableFeatures(sc_dt, mean.cutoff = c(0.0125, 3),dispersion.cutoff = c(0.5, Inf))#define variable genes
all.genes <- rownames(sc_dt)#scale data
sc_dt <- ScaleData(sc_dt, features = all.genes)



###cell type identification with BRETIGEA (BRain cEll Type specIfic Gene Expression Analysis) and seurat
##obtain the marker gene sets from BRETIGEA for all six cell types 
#markers_df_human_brain
table(markers_df_human_brain$cell)
ast=markers_df_human_brain$markers[markers_df_human_brain$cell=='ast']
end=markers_df_human_brain$markers[markers_df_human_brain$cell=='end']
mic=markers_df_human_brain$markers[markers_df_human_brain$cell=='mic']
neu=markers_df_human_brain$markers[markers_df_human_brain$cell=='neu']
oli=markers_df_human_brain$markers[markers_df_human_brain$cell=='oli']
opc=markers_df_human_brain$markers[markers_df_human_brain$cell=='opc']
human_features <- list(ast=ast,end=end,mic=mic,neu=neu,oli=oli,opc=opc)

##calculated a module score for each cell type (Seurat’s AddModuleScore function)
#calculates the average expression levels of each cell type gene set subtracted 
#by the aggregated expression of a background gene set
sc_dt <- AddModuleScore(object = sc_dt,features = list(ast),name = 'ast_Features')
sc_dt <- AddModuleScore(object = sc_dt,features = list(end),name = 'end_Features')
sc_dt <- AddModuleScore(object = sc_dt,features = list(mic),name = 'mic_Features')
sc_dt <- AddModuleScore(object = sc_dt,features = list(neu),name = 'neu_Features')
sc_dt <- AddModuleScore(object = sc_dt,features = list(oli),name = 'oli_Features')
sc_dt <- AddModuleScore(object = sc_dt,features = list(opc),name = 'opc_Features')

##cell type identification step1: identify hybrid and 6 cell types
#1.each cell was assigned a cell type based on the highest cell type score across all six cell types
#2.defined a cell as a hybrid cell if the difference between the first and second highest cell type scores 
#were within 20% of the highest cell type score
opc_score=sc_dt$opc_Features1
mic_score=sc_dt$mic_Features1
end_score=sc_dt$end_Features1
ast_score=sc_dt$ast_Features1
oli_score=sc_dt$oli_Features1
neu_score=sc_dt$neu_Features1
all_score=rbind(opc_score,mic_score,ast_score,neu_score,oli_score,end_score)

#define hybrid cells
selctct=function(allscore){
  step1ct=999
  check=999
  for (cell in 1:dim(allscore)[2]) {
    cells=allscore[,cell]
    maxsc=max(cells)
    maxname=names(cells[which(cells==maxsc,arr.ind=TRUE)])
    cell2=cells[-which(cells==maxsc,arr.ind=TRUE)]
    max2sc=max(cell2)
    ct=ifelse((maxsc-max2sc)/maxsc>=0.2,maxname,'hybrid') 
    my_check=ifelse((maxsc-max2sc)/maxsc>=0.2,maxsc,NA) 
    step1ct=c(step1ct,ct)
    check=c(check,my_check)
  }
  return(list(check,step1ct))
}

stepone=selctct(all_score)
scoreone=unlist(stepone[1])[-1]
ctone=unlist(stepone[2])[-1]
table(ctone)
#ast_score end_score hybrid mic_score neu_score oli_score opc_score 
#2361       100       420       499       687      8005      1142
names(scoreone)=ctone
scoreone_tag=colnames(all_score)
hybridcell=scoreone_tag[names(scoreone)=='hybrid']
##cell type identification step2: identify undefined cells
#for each cell type, apply z-score transformation to the gene score distribution. 
#cells with low cell type score (5th percentile and below) were relabeled as unidentified cells.

#reproduce supplyment1 (visualization)
# step2score=na.omit(scoreone)#drop out hybrid cells' score
# Astrocyte_Scores=step2score[which(names(step2score)=='ast_score')]
# hist(Astrocyte_Scores)
# Oligodendrocytes_Scores=step2score[which(names(step2score)=='oli_score')]
# hist(Oligodendrocytes_Scores)
# OPC_Scores=step2score[which(names(step2score)=='opc_score')]
# hist(OPC_Scores)
# Neurons_Scores=step2score[which(names(step2score)=='neu_score')]
# hist(Neurons_Scores)
# Endothelial_Scores=step2score[which(names(step2score)=='end_score')]
# hist(Endothelial_Scores)
# Microglia_Scores=step2score[which(names(step2score)=='mic_score')]
# hist(Microglia_Scores)

#match score to tag
step2score=as.data.frame(scoreone)
step2score$celltag=scoreone_tag
step2score$type=names(scoreone)
#step2score=na.omit(step2score)#drop out hybrid cells' score

#calculate z score for each cell type
calcz=function(cellscore){
  ctscore=step2score[which(step2score$type==cellscore),]
  my_ctscore=ctscore$scoreone
  # m = mean(my_ctscore)
  # s = sd(my_ctscore)
  # zs = (my_ctscore - m)/s
  zs=my_ctscore
  zsmtx=as.data.frame(zs)
  zsmtx$celltag=ctscore$celltag
  zsmtx$flag=ifelse(zsmtx$zs<quantile(zsmtx$zs, probs = 0.05),'low','ok')
  zsmtx$celltype=rep(substr(cellscore,1,3),dim(zsmtx)[1])
  return(zsmtx)
}
ast_z=calcz('ast_score')
mic_z=calcz('mic_score')
end_z=calcz('end_score')
neu_z=calcz('neu_score')
oli_z=calcz('oli_score')
opc_z=calcz('opc_score')

all_z=rbind(ast_z,mic_z,end_z,neu_z,oli_z,opc_z)
all_z$celltype=ifelse(all_z$flag=='low','undefined',all_z$celltype)
table(all_z$celltype)
#My result from cell type identification
#ast       end       mic       neu       oli       opc   undefined   hybrid
#2,243    95       474       652     7,604   1,084       642          420

#Author's result from cell type identification
#ast       end       mic       neu       oli       opc   undefined   hybrid
#2,171    98        449       656    7,432   1,078      925          405

###Human single-nuclei UMAP and clustering analysis.
##create iden for seurat obj (based on scoring systom, import defined cell types into seurat obj)
iden=all_z[,c('celltag','celltype')]
hybrididen=data.frame(celltag=hybridcell,celltype=rep('hybrid',length(hybridcell)))
alliden=rbind(iden,hybrididen)
all_iden=alliden[match(colnames(sc_dt),alliden$celltag),]
sc_dt$ident=as.factor(all_iden$celltype)

##scdt reduce dimention
sc_dt <- RunPCA(sc_dt, features = VariableFeatures(object = sc_dt))
ElbowPlot(sc_dt)
sc_dt <- FindNeighbors(sc_dt, dims = 1:19)
sc_dt <- FindClusters(sc_dt, resolution = 0.8)
sc_dt <- RunUMAP(sc_dt, dims = 1:19)
DimPlot(sc_dt, reduction = "umap")
#saveRDS(sc_dt,file = 'grubmans_seuratobj.rds')

##visualization based on AD or control
adct=sc_meta$batchCond
Idents(sc_dt) <- adct
#DimPlot(sc_dt, reduction = "umap", pt.size = 0.8,cols = c('AD' = 'purple','ct' = 'green'))

##visualization based on individual
ptag=sc_meta$patient
Idents(sc_dt) <- ptag
DimPlot(sc_dt, reduction = "umap", pt.size = 0.8,
        cols = c('AD1' = heat.colors(12)[1],'AD2' = heat.colors(12)[2],'AD3' = heat.colors(12)[3],'AD4' = heat.colors(12)[4],
                 'AD5' = heat.colors(12)[5],'AD6' = heat.colors(12)[6],'Ct5' = topo.colors(12)[5],'Ct6' = topo.colors(12)[6],
                 'AD-un' = 'black','Ct1' = topo.colors(12)[1],'Ct2' = topo.colors(12)[2],'Ct3' = topo.colors(12)[3],
                 'Ct4' = topo.colors(12)[4],'Ct-un' = 'grey'))

##visualization based on cell type
#sc_dt$CellType <- sc_dt$ident
Idents(sc_dt) <- sc_dt$ident
DimPlot(sc_dt, reduction = "umap", label = TRUE, pt.size = 0.8,
        cols = c('ast' = 'pink', 'end' = 'brown', 'mic' = 'blue','hybrid'='black','neu'='red',
                 'oli'='orange','opc'='purple','undefined'='grey')) # + NoLegend()


#cell type specific genes
sc_dt.markers <- FindAllMarkers(sc_dt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sc_dt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- sc_dt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#DoHeatmap(sc_dt, features = top10$gene) + NoLegend()


##subset the data
asttag=all_iden[all_iden$celltype=='ast',1]
mictag=all_iden[all_iden$celltype=='mic',1]
opctag=all_iden[all_iden$celltype=='opc',1]
olitag=all_iden[all_iden$celltype=='oli',1]
neutag=all_iden[all_iden$celltype=='neu',1]
endtag=all_iden[all_iden$celltype=='end',1]
undftag=all_iden[all_iden$celltype=='undefined',1]
hybtag=all_iden[all_iden$celltype=='hybrid',1]

#run subcluster on specific cell type
runsubclust=function(ctname){
  ct_dt <- subset(sc_dt, idents = ctname)
  ct_dt <- NormalizeData(ct_dt, normalization.method = "LogNormalize", scale.factor = 10000)#normalization
  ct_dt <- FindVariableFeatures(ct_dt, mean.cutoff = c(0.0125, 3),dispersion.cutoff = c(0.5, Inf))#define variable genes
  all.genes <- rownames(ct_dt)#scale data
  ct_dt <- ScaleData(ct_dt, features = all.genes)
  ct_dt <- RunPCA(ct_dt, features = VariableFeatures(object = ct_dt))
  return(ct_dt)
}

#annotate subclusters
annosub=function(ct_dt,ctname){
  subct=ct_dt$seurat_clusters
  subcttag=names(ct_dt$nCount_RNA)
  subct=data.frame(subct=subct,cttag=subcttag)
  subct$ctcluster=paste(ctname,subct$subct,sep = "")
  return(subct[,c(2,3)])
}

#mic=15; ast=16; neu=15; oli=15; OPC=16; end=18; unidentified=7; hybrid=9

#ast subclusters: 8
ast_dt=runsubclust('ast')
#ElbowPlot(ast_dt)
ast_dt <- FindNeighbors(ast_dt, dims = 1:16)
ast_dt <- FindClusters(ast_dt, resolution = 0.8)#got 8 subclusters
ast_dt <- RunUMAP(ast_dt, dims = 1:16)
DimPlot(ast_dt, reduction = "umap")
subast=annosub(ast_dt,'ast')
#subast$group=substr(subast$cttag,22,23)

#proportion plot (AD vs. control)
new_perct=tapply(subast[,3],subast[,2],count)
a=c(7,3,9,225,0,145,145,0)#ad freq
b=c(690,506,243,2,212,3,0,53)
astname=c('ast0','ast1','ast2','ast3','ast4','ast5','ast6','ast7')
aa=data.frame(subtype=astname,value=a,group=rep('AD',8),pct=a/(a+b))
bb=data.frame(subtype=astname,value=b,group=rep('control',8),pct=b/(a+b))
astpc=rbind(aa,bb)
ce = ddply(astpc, "subtype", transform, percent_pt = pct * 100)
ggplot(ce, aes(x = subtype, y = percent_pt, fill = group)) +
  geom_bar(stat = "identity", colour = "black")+scale_fill_manual(values = c(rainbow(34)[5],rainbow(28)[17]))+theme_bw()+
  labs(x = "Ast subtype", y = "cell proportion", title = "AD VS. Control cell proportion in Ast subtype")

#ast subtypes specific genes
ast_dt.markers <- FindAllMarkers(ast_dt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ast_dt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- ast_dt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(ast_dt, features = top10$gene) + NoLegend()

#mic subclusters: 4
mic_dt=runsubclust('mic')
#ElbowPlot(mic_dt)
mic_dt <- FindNeighbors(mic_dt, dims = 1:15)
mic_dt <- FindClusters(mic_dt, resolution = 0.8)#got 4 subclusters
mic_dt <- RunUMAP(mic_dt, dims = 1:15)
DimPlot(mic_dt, reduction = "umap")
submic=annosub(mic_dt,'mic')

#oli subclusters: 7
oli_dt=runsubclust('oli')
#ElbowPlot(oli_dt)
oli_dt <- FindNeighbors(oli_dt, dims = 1:15)
oli_dt <- FindClusters(oli_dt, resolution = 0.8)#got 7 subclusters
oli_dt <- RunUMAP(oli_dt, dims = 1:15)
DimPlot(oli_dt, reduction = "umap")
suboli=annosub(oli_dt,'oli')
#suboli$group=substr(suboli$cttag,22,23)

#oli subtype proportion plot (AD vs. control)
new_perct=tapply(suboli[,3],suboli[,2],count)
a=c(2029,23,1932,21,618,62,94)#ad
b=c(29,1936,1,693,45,87,34)#ct
oliname=c('oli0','oli1','oli2','oli3','oli4','oli5','oli6')
aa=data.frame(subtype=oliname,value=a,group=rep('AD',7),pct=a/(a+b))
bb=data.frame(subtype=oliname,value=b,group=rep('control',7),pct=b/(a+b))
olipc=rbind(aa,bb)
ce = ddply(olipc, "subtype", transform, percent_pt = pct * 100)
ggplot(ce, aes(x = subtype, y = percent_pt, fill = group)) +
  geom_bar(stat = "identity", colour = "black")+scale_fill_manual(values = c(rainbow(34)[5],rainbow(28)[17]))+theme_bw()+
  labs(x = "Oli subtype", y = "cell proportion", title = "AD VS. Control cell proportion in Oli subtype")

#oli subtypes specific genes
oli_dt.markers <- FindAllMarkers(oli_dt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
oli_dt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- oli_dt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(oli_dt, features = top10$gene) + NoLegend()

#neu subclusters: 7
neu_dt=runsubclust('neu')
ElbowPlot(neu_dt)
neu_dt <- FindNeighbors(neu_dt, dims = 1:15)
neu_dt <- FindClusters(neu_dt, resolution = 0.8)#got 7 subclusters
neu_dt <- RunUMAP(neu_dt, dims = 1:15)
DimPlot(neu_dt, reduction = "umap")
subneu=annosub(neu_dt,'neu')
#subneu$group=substr(subneu$cttag,22,23)

#neu subtype proportion plot (AD vs. control)
new_perct=tapply(subneu[,3],subneu[,2],count)
a=c(44,41,88,5,36,11,20)
b=c(153,65,10,80,25,44,30)
neuname=c('neu0','neu1','neu2','neu3','neu4','neu5','neu6')
aa=data.frame(subtype=neuname,value=a,group=rep('AD',7),pct=a/(a+b))
bb=data.frame(subtype=neuname,value=b,group=rep('control',7),pct=b/(a+b))
neupc=rbind(aa,bb)
ce = ddply(neupc, "subtype", transform, percent_pt = pct * 100)
ggplot(ce, aes(x = subtype, y = percent_pt, fill = group)) +
  geom_bar(stat = "identity", colour = "black")+scale_fill_manual(values = c(rainbow(34)[5],rainbow(28)[17]))+theme_bw()+
  labs(x = "Neu subtype", y = "cell proportion", title = "AD VS. Control cell proportion in Neu subtype")

#neu subtypes specific genes
neu_dt.markers <- FindAllMarkers(neu_dt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
neu_dt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- neu_dt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(neu_dt, features = top10$gene) + NoLegend()

#end subclusters: 3
end_dt=runsubclust('end')
#ElbowPlot(end_dt)
end_dt <- FindNeighbors(end_dt, dims = 1:18)
end_dt <- FindClusters(end_dt, resolution = 0.8)#got 3 subclusters
end_dt <- RunUMAP(end_dt, dims = 1:18)
DimPlot(end_dt, reduction = "umap")
subend=annosub(end_dt,'end')

#opc subclusters: 5
opc_dt=runsubclust('opc')
#ElbowPlot(opc_dt)
opc_dt <- FindNeighbors(opc_dt, dims = 1:16)
opc_dt <- FindClusters(opc_dt, resolution = 0.8)#got 5 subclusters
opc_dt <- RunUMAP(opc_dt, dims = 1:16)
DimPlot(opc_dt, reduction = "umap")
subopc=annosub(opc_dt,'opc')

#undefined subclusters: 6
undefined_dt=runsubclust('undefined')
#ElbowPlot(undefined_dt)
undefined_dt <- FindNeighbors(undefined_dt, dims = 1:7)
undefined_dt <- FindClusters(undefined_dt, resolution = 0.8)#got 6 subclusters
undefined_dt <- RunUMAP(undefined_dt, dims = 1:7)
DimPlot(undefined_dt, reduction = "umap")
subundefined=annosub(undefined_dt,'undefined')

#hybrid subclusters: 6
hybrid_dt=runsubclust('hybrid')
ElbowPlot(hybrid_dt)
hybrid_dt <- FindNeighbors(hybrid_dt, dims = 1:9)
hybrid_dt <- FindClusters(hybrid_dt, resolution = 0.8)#got 6 subclusters
hybrid_dt <- RunUMAP(hybrid_dt, dims = 1:9)
DimPlot(hybrid_dt, reduction = "umap")
subhybrid=annosub(hybrid_dt,'hybrid')

#add all subtype tag to sc_dt ident
all_sub=rbind(subundefined,subhybrid,subopc,subend,subast,submic,suboli,subneu)
all_subiden=all_sub[match(colnames(sc_dt),all_sub$cttag),]
sc_dt$subiden=as.factor(all_subiden$ctcluster)
#table(all_sub$ctcluster)
# ast0       ast1       ast2       ast3       ast4       ast5       ast6       ast7       end0 
# 697        509        252        227        212        148        145         53         39 
# end1       end2    hybrid0    hybrid1    hybrid2    hybrid3    hybrid4    hybrid5       mic0 
# 31         25        100         97         87         73         35         28        168 
# mic1       mic2       mic3       neu0       neu1       neu2       neu3       neu4       neu5 
# 104        103         99        197        106         98         85         61         55 
# neu6       oli0       oli1       oli2       oli3       oli4       oli5       oli6       opc0 
# 50       2058       1959       1933        714        663        149        128        315 
# opc1       opc2       opc3       opc4 undefined0 undefined1 undefined2 undefined3 undefined4 
# 249        227        187        106        226        153        128         88         37 
# undefined5 
# 10 

##visualization based on subast
ch=as.character(sc_dt$subiden)

vissubast=ifelse(sc_dt$subiden %in% c('ast0','ast1','ast2','ast3','ast4','ast5','ast6','ast7'),ch,'others')
Idents(sc_dt) <- vissubast
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('ast0' = rainbow(10)[1], 'ast1' = rainbow(10)[2],'ast2' = rainbow(10)[3],'ast3' = rainbow(10)[4],
                 'ast4' = rainbow(10)[5], 'ast5' = rainbow(10)[6],'ast6' = rainbow(10)[7],'ast7' = rainbow(10)[8],
                 'others'='grey')) #+ NoLegend()

vissubend=ifelse(sc_dt$subiden %in% c('end0','end1','end2'),ch,'others')
Idents(sc_dt) <- vissubend
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('end0' = rainbow(10)[1], 'end1' = rainbow(10)[2],'end2' = rainbow(10)[3],'others'='grey')) 

vissubhybrid=ifelse(sc_dt$subiden %in% c('hybrid0','hybrid1','hybrid2','hybrid3','hybrid4','hybrid5'),ch,'others')
Idents(sc_dt) <- vissubhybrid
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('hybrid0' = rainbow(10)[1], 'hybrid1' = rainbow(10)[2],'hybrid2' = rainbow(10)[3],'hybrid3' = rainbow(10)[4],
                 'hybrid4' = rainbow(10)[6], 'hybrid5' = rainbow(10)[8],'others'='grey'))

vissubmic=ifelse(sc_dt$subiden %in% c('mic0','mic1','mic2','mic3'),ch,'others')
Idents(sc_dt) <- vissubmic
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('mic0' = rainbow(10)[1], 'mic1' = rainbow(10)[2],'mic2' = rainbow(10)[5],'mic3' = rainbow(10)[7],'others'='grey'))

vissuboli=ifelse(sc_dt$subiden %in% c('oli0','oli1','oli2','oli3','oli4','oli5','oli6'),ch,'others')
Idents(sc_dt) <- vissuboli
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('oli0' = rainbow(10)[1], 'oli1' = rainbow(10)[2],'oli2' = rainbow(10)[3],'oli3' = rainbow(10)[4],
                 'oli4' = rainbow(10)[6], 'oli5' = rainbow(10)[7],'oli6' = rainbow(10)[8],'others'='grey'))

vissubopc=ifelse(sc_dt$subiden %in% c('opc0','opc1','opc2','opc3','opc4'),ch,'others')
Idents(sc_dt) <- vissubopc
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('opc0' = rainbow(10)[1], 'opc1' = rainbow(10)[2],'opc2' = rainbow(10)[3],'opc3' = rainbow(10)[4],
                 'opc4' = rainbow(10)[6], 'others'='grey'))

vissubundefined=ifelse(sc_dt$subiden %in% c('undefined0','undefined1','undefined2','undefined3','undefined4','undefined5'),ch,'others')
Idents(sc_dt) <- vissubundefined
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('undefined0' = rainbow(10)[1], 'undefined1' = rainbow(10)[2],'undefined2' = rainbow(10)[3],'undefined3' = rainbow(10)[4],
                 'undefined4' = rainbow(10)[6], 'undefined5' = rainbow(10)[7],'others'='grey'))

vissubneu=ifelse(sc_dt$subiden %in% c('neu0','neu1','neu2','neu3','neu4','neu5','neu6'),ch,'others')
Idents(sc_dt) <- vissubneu
DimPlot(sc_dt, reduction = "umap", label = FALSE, pt.size = 0.8,
        cols = c('neu0' = rainbow(10)[1], 'neu1' = rainbow(10)[2],'neu2' = rainbow(10)[3],'neu3' = rainbow(10)[4],
                 'neu4' = rainbow(10)[6], 'neu5' = rainbow(10)[7],'neu6'=rainbow(10)[9],'others'='grey'))

##visualize patients' cell proportion in cell types
sc_dt$pid=sc_meta$patient
patientct=data.frame(patient=sc_dt$pid,celltype=sc_dt$ident)

new_perct=tapply(patientct[,1],patientct[,2],count)
newpatientct=data.frame(celltype=rep('ast',14),patient=new_perct$ast$x,pcount=new_perct$ast$freq)
newpatientct$perct=newpatientct$pcount/sum(newpatientct$pcount)

newmic=data.frame(celltype=rep('mic',14),patient=new_perct$mic$x,pcount=new_perct$mic$freq)
newmic$perct=newmic$pcount/sum(newmic$pcount)

newoli=data.frame(celltype=rep('oli',14),patient=new_perct$oli$x,pcount=new_perct$oli$freq)
newoli$perct=newoli$pcount/sum(newoli$pcount)

newopc=data.frame(celltype=rep('opc',14),patient=new_perct$opc$x,pcount=new_perct$opc$freq)
newopc$perct=newopc$pcount/sum(newopc$pcount)

newneu=data.frame(celltype=rep('neu',14),patient=new_perct$neu$x,pcount=new_perct$neu$freq)
newneu$perct=newneu$pcount/sum(newneu$pcount)

newend=data.frame(celltype=rep('end',14),patient=c(as.character(new_perct$end$x),'AD-un','Ct-un'),pcount=c(new_perct$end$freq,0,0))
newend$perct=newend$pcount/sum(newend$pcount)

newundefined=data.frame(celltype=rep('undefined',14),patient=c(as.character(new_perct$undefined$x),'Ct4'),pcount=c(new_perct$undefined$freq,0))
newundefined$perct=newundefined$pcount/sum(newundefined$pcount)

newhybrid=data.frame(celltype=rep('hybrid',14),patient=c(as.character(new_perct$hybrid$x),'AD6'),pcount=c(new_perct$hybrid$freq,0))
newhybrid$perct=newhybrid$pcount/sum(newhybrid$pcount)

newpatientct=rbind(newmic,newpatientct,newneu,newoli,newopc,newend,newundefined,newhybrid)
ce = ddply(newpatientct, "celltype", transform, percent_pt = perct * 100)
ggplot(ce, aes(x = celltype, y = percent_pt, fill = patient)) +
  geom_bar(stat = "identity", colour = "black")+scale_fill_manual(values = c(rainbow(34)[1:7],rainbow(28)[15:22]))+theme_bw()+
  labs(x = "cell types", y = "cell proportion", title = "patients' cell proportion in cell types")


###Identification of individual- and sex-specific genes
#scdt is the sc count mtx, sc_meta is metadata
celltype=sc_dt$ident

#group cell type and ad/ct
ad_ast=scdt[,sc_meta$batchCond=='AD' & celltype=='ast']
control_ast=scdt[,sc_meta$batchCond=='ct' & celltype=='ast']
#ad_ast_tag=sc_meta$sampleID[sc_meta$batchCond=='AD' & celltype=='ast']
ad_ast_sex=sc_meta$sex[sc_meta$batchCond=='AD' & celltype=='ast']
control_ast_sex=sc_meta$sex[sc_meta$batchCond=='ct' & celltype=='ast']
#individuals specific genes
ad_ast_ind=sc_meta$patient[sc_meta$batchCond=='AD' & celltype=='ast']
#ad_ast_ind=ifelse(ad_ast_ind=='AD1','AD1','Others')
ad_ast_ind=ifelse(ad_ast_ind=='AD2','AD2','Others')
control_ast_ind=sc_meta$patient[sc_meta$batchCond=='ct' & celltype=='ast']
# control_ast_ind=ifelse(control_ast_ind=='Ct1','Ct1','Others')
control_ast_ind=ifelse(control_ast_ind=='Ct2','Ct2','Others')

ad_mic=scdt[,sc_meta$batchCond=='AD' & celltype=='mic']
control_mic=scdt[,sc_meta$batchCond=='ct' & celltype=='mic']
ad_mic_sex=sc_meta$sex[sc_meta$batchCond=='AD' & celltype=='mic']
control_mic_sex=sc_meta$sex[sc_meta$batchCond=='ct' & celltype=='mic']
#individuals specific genes
ad_mic_ind=sc_meta$patient[sc_meta$batchCond=='AD' & celltype=='mic']
# ad_mic_ind=ifelse(ad_mic_ind=='AD1','AD1','Others')
ad_mic_ind=ifelse(ad_mic_ind=='AD2','AD2','Others')
control_mic_ind=sc_meta$patient[sc_meta$batchCond=='ct' & celltype=='mic']
# control_mic_ind=ifelse(control_mic_ind=='Ct1','Ct1','Others')
control_mic_ind=ifelse(control_mic_ind=='Ct2','Ct2','Others')

ad_neu=scdt[,sc_meta$batchCond=='AD' & celltype=='neu']
control_neu=scdt[,sc_meta$batchCond=='ct' & celltype=='neu']
ad_neu_sex=sc_meta$sex[sc_meta$batchCond=='AD' & celltype=='neu']
control_neu_sex=sc_meta$sex[sc_meta$batchCond=='ct' & celltype=='neu']
#individuals specific genes
ad_neu_ind=sc_meta$patient[sc_meta$batchCond=='AD' & celltype=='neu']
# ad_neu_ind=ifelse(ad_neu_ind=='AD1','AD1','Others')
ad_neu_ind=ifelse(ad_neu_ind=='AD2','AD2','Others')
control_neu_ind=sc_meta$patient[sc_meta$batchCond=='ct' & celltype=='neu']
# control_neu_ind=ifelse(control_neu_ind=='Ct1','Ct1','Others')
control_neu_ind=ifelse(control_neu_ind=='Ct2','Ct2','Others')

ad_oli=scdt[,sc_meta$batchCond=='AD' & celltype=='oli']
control_oli=scdt[,sc_meta$batchCond=='ct' & celltype=='oli']
ad_oli_sex=sc_meta$sex[sc_meta$batchCond=='AD' & celltype=='oli']
control_oli_sex=sc_meta$sex[sc_meta$batchCond=='ct' & celltype=='oli']
#individuals specific genes
ad_oli_ind=sc_meta$patient[sc_meta$batchCond=='AD' & celltype=='oli']
# ad_oli_ind=ifelse(ad_oli_ind=='AD1','AD1','Others')
ad_oli_ind=ifelse(ad_oli_ind=='AD2','AD2','Others')
control_oli_ind=sc_meta$patient[sc_meta$batchCond=='ct' & celltype=='oli']
# control_oli_ind=ifelse(control_oli_ind=='Ct1','Ct1','Others')
control_oli_ind=ifelse(control_oli_ind=='Ct2','Ct2','Others')

ad_opc=scdt[,sc_meta$batchCond=='AD' & celltype=='opc']
control_opc=scdt[,sc_meta$batchCond=='ct' & celltype=='opc']
ad_opc_sex=sc_meta$sex[sc_meta$batchCond=='AD' & celltype=='opc']
control_opc_sex=sc_meta$sex[sc_meta$batchCond=='ct' & celltype=='opc']
#individuals specific genes
ad_opc_ind=sc_meta$patient[sc_meta$batchCond=='AD' & celltype=='opc']
# ad_opc_ind=ifelse(ad_opc_ind=='AD1','AD1','Others')
ad_opc_ind=ifelse(ad_opc_ind=='AD2','AD2','Others')
control_opc_ind=sc_meta$patient[sc_meta$batchCond=='ct' & celltype=='opc']
# control_opc_ind=ifelse(control_opc_ind=='Ct1','Ct1','Others')
control_opc_ind=ifelse(control_opc_ind=='Ct2','Ct2','Others')

##sex specific genes
#create DEGList obj
group <- factor(ad_ast_sex)
y_sex <- DGEList(counts = ad_ast, group = group)
#filtering out 0 expression gene (nothing is filtered out)
# keep <- rowSums(cpm(y_sex)>1) >= 2
# y_sex <- y_sex[keep, , keep.lib.sizes=FALSE]
#normalization
y_sex <- calcNormFactors(y_sex)
#plotMDS(y_sex)
#estimate dispersion
y_sex <- estimateCommonDisp(y_sex, verbose=TRUE)
#Disp = 1e-04 , BCV = 0.01 
y_sex <- estimateTagwiseDisp(y_sex)
#plotBCV(y_sex)
#find out differential genes
et <- exactTest(y_sex)
top <- topTags(et)
#defind DEGs
summary(de <- decideTestsDGE(et,p.value = 0.01,lfc = 0.5))
#       Male-Female
# Down           875
# NotSig        9516
# Up             204
detags <- rownames(y_sex)[as.logical(de)]
#plot
plotSmear(et, de.tags=detags);
abline(h=c(-0.5, 0.5), col="blue");

getdeg=function(metadata,countdata){
  group <- factor(metadata)
  y <- DGEList(counts = countdata, group = group)
  y <- calcNormFactors(y)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y, verbose=TRUE)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y)
  top <- topTags(et)
  de <- decideTestsDGE(et,p.value = 0.01,lfc = 0.5)
  detags <- rownames(y)[as.logical(de)]
  return(list(de,detags))
}

#ast AD SEX degs
astdeg=getdeg(ad_ast_sex,ad_ast)
ad_ast_de_sex=unlist(astdeg[1])
ad_ast_detags_sex=unlist(astdeg[2])
#       Male-Female
# Down           875
# NotSig        9516
# Up             204

#ast control SEX degs
astdeg=getdeg(control_ast_sex,control_ast)
control_ast_de_sex=unlist(astdeg[1])
control_ast_detags_sex=unlist(astdeg[2])
#       Male-Female
# Down           7
# NotSig        10649
# Up             194

#mic AD SEX degs
micdeg=getdeg(ad_mic_sex,ad_mic)
ad_mic_de_sex=unlist(micdeg[1])
ad_mic_detags_sex=unlist(micdeg[2])
#       Male-Female
# Down           139
# NotSig        10680
# Up             31

#mic control SEX degs
micdeg=getdeg(control_mic_sex,control_mic)
control_mic_de_sex=unlist(micdeg[1])
control_mic_detags_sex=unlist(micdeg[2])
#       Male-Female
# Down           2
# NotSig        10844
# Up             4

#neu AD SEX degs
neudeg=getdeg(ad_neu_sex,ad_neu)
ad_neu_de_sex=unlist(neudeg[1])
ad_neu_detags_sex=unlist(neudeg[2])
#       Male-Female
# Down           547
# NotSig        10220
# Up             83

#neu control SEX degs
neudeg=getdeg(control_neu_sex,control_neu)
control_neu_de_sex=unlist(neudeg[1])
control_neu_detags_sex=unlist(neudeg[2])
#       Male-Female
# Down           5
# NotSig        10841
# Up             4

#oli AD SEX degs
olideg=getdeg(ad_oli_sex,ad_oli)
ad_oli_de_sex=unlist(olideg[1])
ad_oli_detags_sex=unlist(olideg[2])
#       Male-Female
# Down          1356
# NotSig        9273
# Up            221

#oli control SEX degs
olideg=getdeg(control_oli_sex,control_oli)
control_oli_de_sex=unlist(olideg[1])
control_oli_detags_sex=unlist(olideg[2])
#       Male-Female
# Down          4
# NotSig        10338
# Up            508

#opc AD SEX degs
opcdeg=getdeg(ad_opc_sex,ad_opc)
ad_opc_de_sex=unlist(opcdeg[1])
ad_opc_detags_sex=unlist(opcdeg[2])
#       Male-Female
# Down           46
# NotSig        10796
# Up             8

#opc control SEX degs
opcdeg=getdeg(control_opc_sex,control_opc)
control_opc_de_sex=unlist(opcdeg[1])
control_opc_detags_sex=unlist(opcdeg[2])
#       Male-Female
# Down          6
# NotSig        10839
# Up            5

##individual specific genes
#ast AD1 vs. others ind degs
astdeg=getdeg(ad_ast_ind,ad_ast)
ad_ast_de_ind=unlist(astdeg[1])
ad_ast_detags_ind=unlist(astdeg[2])
# -1     0     1 
# 6 10791    53 
# -1     0     1 
# 209 10500   141 
#mic AD1 vs. others ind degs
micdeg=getdeg(ad_mic_ind,ad_mic)
ad_mic_de_ind=unlist(micdeg[1])
ad_mic_detags_ind=unlist(micdeg[2])
#    0 
#10850 
# -1     0     1 
# 26 10809    15
#neu AD1 vs. others ind degs
neudeg=getdeg(ad_neu_ind,ad_neu)
ad_neu_de_ind=unlist(neudeg[1])
ad_neu_detags_ind=unlist(neudeg[2])
# -1     0     1 
# 6 10799    45 
# -1     0     1 
# 19 10815    16 
#opc AD1 vs. others ind degs
opcdeg=getdeg(ad_opc_ind,ad_opc)
ad_opc_de_ind=unlist(opcdeg[1])
ad_opc_detags_ind=unlist(opcdeg[2])
# -1     0 
# 1 10849
# -1     0 
# 4 10846 
#oli AD1 vs. others ind degs
olideg=getdeg(ad_oli_ind,ad_oli)
ad_oli_de_ind=unlist(olideg[1])
ad_oli_detags_ind=unlist(olideg[2])
# -1     0     1 
# 39 10063   748
# -1     0     1 
# 120 10407   323 

##individual specific genes
#ast ct1 vs. others ind degs
astdeg=getdeg(control_ast_ind,control_ast)
control_ast_de_ind=unlist(astdeg[1])
control_ast_detags_ind=unlist(astdeg[2])
# -1     0     1 
# 6 10821    23 
# -1     0     1 
# 1 10758    91
#mic ct1 vs. others ind degs
micdeg=getdeg(control_mic_ind,control_mic)
control_mic_de_ind=unlist(micdeg[1])
control_mic_detags_ind=unlist(micdeg[2])
# -1     0 
# 2 10848
# 0 
# 10850
#neu ct1 vs. others ind degs
neudeg=getdeg(control_neu_ind,control_neu)
control_neu_de_ind=unlist(neudeg[1])
control_neu_detags_ind=unlist(neudeg[2])
# -1     0     1 
# 3 10846     1
# 0 
# 10850 
#oli ct1 vs. others ind degs
olideg=getdeg(control_oli_ind,control_oli)
control_oli_de_ind=unlist(olideg[1])
control_oli_detags_ind=unlist(olideg[2])
# -1     0     1 
# 6 10479   365 
# -1     0     1 
# 8 10607   235 
#opc ct1 vs. others ind degs
opcdeg=getdeg(control_opc_ind,control_opc)
control_opc_de_ind=unlist(opcdeg[1])
control_opc_detags_ind=unlist(opcdeg[2])
# -1     0     1 
# 2 10846     2 
# -1     0     1 
# 6 10841     3

#visualization in ad2
ind_ad_bar=data.frame(numDEGs=c(141,-209,15,-26,16,-19,323,-120,0,-4),
                      CellType=c('Ast','Ast','Mic','Mic','Neu','Neu','Oli','Oli','Opc','Opc'))
ind_ad_bar$RegulateDirection=ifelse(ind_ad_bar$numDEGs>=0,'Up','Down')
ggplot(ind_ad_bar, aes(x=CellType,y=numDEGs, fill = RegulateDirection)) +theme_bw()+
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"))+
  labs(x = "cell types", y = "Numbers of DEGs", title = "DEGs in AD2-otherADs")
#visualization in ct2
ind_ct_bar=data.frame(numDEGs=c(91,-1,0,0,0,0,235,-8,3,-6),
                      CellType=c('Ast','Ast','Mic','Mic','Neu','Neu','Oli','Oli','Opc','Opc'))
ind_ct_bar$RegulateDirection=ifelse(ind_ct_bar$numDEGs>=0,'Up','Down')
ggplot(ind_ct_bar, aes(x=CellType,y=numDEGs, fill = RegulateDirection)) +theme_bw()+
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"))+
  labs(x = "cell types", y = "Numbers of DEGs", title = "DEGs in Control2-otherControls")

#visualization in ct1
ind_ct_bar=data.frame(numDEGs=c(23,-6,0,-2,1,-3,365,-6,2,-2),
                      CellType=c('Ast','Ast','Mic','Mic','Neu','Neu','Oli','Oli','Opc','Opc'))
ind_ct_bar$RegulateDirection=ifelse(ind_ct_bar$numDEGs>=0,'Up','Down')
ggplot(ind_ct_bar, aes(x=CellType,y=numDEGs, fill = RegulateDirection)) +theme_bw()+
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"))+
  labs(x = "cell types", y = "Numbers of DEGs", title = "DEGs in Control1-otherControls")

#visualization in ad1
ind_ad_bar=data.frame(numDEGs=c(53,-6,0,0,45,-6,748,-39,0,-1),
                      CellType=c('Ast','Ast','Mic','Mic','Neu','Neu','Oli','Oli','Opc','Opc'))
ind_ad_bar$RegulateDirection=ifelse(ind_ad_bar$numDEGs>=0,'Up','Down')
ggplot(ind_ad_bar, aes(x=CellType,y=numDEGs, fill = RegulateDirection)) +theme_bw()+
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"))+
  labs(x = "cell types", y = "Numbers of DEGs", title = "DEGs in AD1-otherADs")

#visualization
sex_ad_bar=data.frame(numDEGs=c(204,-857,31,-139,83,-547,221,-1356,8,-46),
                      CellType=c('Ast','Ast','Mic','Mic','Neu','Neu','Oli','Oli','Opc','Opc'))
sex_ad_bar$RegulateDirection=ifelse(sex_ad_bar$numDEGs>=0,'Up','Down')
ggplot(sex_ad_bar, aes(x=CellType,y=numDEGs, fill = RegulateDirection)) +theme_bw()+
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"))+
  labs(x = "cell types", y = "Numbers of DEGs", title = "DEGs in AD male&female")

#control
sex_ct_bar=data.frame(numDEGs=c(194,-7,4,-2,4,-5,508,-4,5,-6),
                      CellType=c('Ast','Ast','Mic','Mic','Neu','Neu','Oli','Oli','Opc','Opc'))
sex_ct_bar$RegulateDirection=ifelse(sex_ct_bar$numDEGs>=0,'Up','Down')
ggplot(sex_ct_bar, aes(x=CellType,y=numDEGs, fill = RegulateDirection)) +theme_bw()+
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"))+
  labs(x = "cell types", y = "Numbers of DEGs", title = "DEGs in Control male&female")






###Human single-nuclei differential expression and gene set enrichment analysis
# find markers for every cluster compared to all remaining cells, report only the positive ones

