#Bisque
library(Biobase)
library(BisqueRNA)
library(Seurat)
library(readxl)

#mathys bulk data
bulk.matrix=bulk[-1,]
bulk.mtx=data.matrix(bulk.matrix)
colnames(bulk.mtx)=as.character(bulk[1,])
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.mtx)

#mathys sc data
dat <- readMM("~/filtered_count_matrix.mtx")
genenames <- readLines("~/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("~/filtered_column_metadata.txt")
indgroup <- read_excel("~/mathys_pathologygroup.xlsx")
data=dat
pid=as.data.frame(filtered.colMetadata$projid)
colnames(pid)='projid'
ad=inner_join(pid,indgroup, by = 'projid')
data=as.matrix(data)
colnames(data)=filtered.colMetadata$TAG
rownames(data)=genenames

sample.ids <- filtered.colMetadata$TAG
individual.labels=filtered.colMetadata$projid
cell.type.labels=filtered.colMetadata$broad.cell.type
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=individual.labels,
                       cellType=cell.type.labels)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=data,
                                  phenoData=sc.pdata)
# markers could be set here

###Reference-based decomposition
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)
ref.based.estimates <- res$bulk.props
refProp=t(ref.based.estimates)
referenceProp=refProp[,c(1,2,3,4,5,6,7)]
write.csv(referenceProp,'Bisque_refProp.csv')
boxplot(refProp,main='Bisque')

###Marker-based decomposition
# cell.types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells")
# avg.props <- c(.5, .2, .2, .07, .03)
# sim.data <- SimulateData(n.ind=10, n.genes=100, n.cells=500, cell.types=cell.types, avg.props=avg.props)
# markers <- sim.data$markers

mathysobj <- CreateSeuratObject(counts = sc.eset, min.cells = 3, min.features = 200)
mathysobj <- NormalizeData(mathysobj, normalization.method = "LogNormalize", scale.factor = 10000)
mathysobj <- FindVariableFeatures(mathysobj, selection.method = "vst", nfeatures = 2000)
mathysobj <- ScaleData(mathysobj)
mathysobj <- RunPCA(mathysobj, features = VariableFeatures(object = mathysobj))
mathysobj <- FindNeighbors(mathysobj, dims = 1:10)
mathysobj <- FindClusters(mathysobj, resolution = 0.5)
mathysobj <- RunUMAP(mathysobj, dims = 1:10)
DimPlot(mathysobj, reduction = "umap", label = TRUE)
mathysobj.markers <- FindAllMarkers(mathysobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(mathysobj,'mathysobj.rds')
saveRDS(mathysobj.markers,'mathysobj_markers.rds')
submarker <- read_excel("~/DeconvolutionMethods/Mathys_subclusterMarkers.xlsx")
submarker=submarker[,1:4]
submarker=as.data.frame(submarker)
#ex -'NRGN', cluster2
#in -'GAD1', cluster6
#ast -'AQP4', cluster9
#oli -'MBP', cluster0
#mic -'CSF1R and CD74', cluster 13
#opc -'VCAN', cluster11
#end -'FLT1',cluster NA
#per -'AMBP', cluster NA
exmk1=mathysobj.markers[mathysobj.markers$cluster==2,]
exmk2=mathysobj.markers[mathysobj.markers$cluster==6,]
exmk=rbind(exmk1,exmk2)
astmk=mathysobj.markers[mathysobj.markers$cluster==9,]
olimk=mathysobj.markers[mathysobj.markers$cluster==0,]
micmk=mathysobj.markers[mathysobj.markers$cluster==13,]
opcmk=mathysobj.markers[mathysobj.markers$cluster==11,]
endmk=mathysobj.markers[mathysobj.markers$cluster==17,]
save(exmk, astmk,olimk,micmk,opcmk,endmk, file = "Mathys_ctmarkers.RData")

#subclster
mathysobj$ident=as.factor(filtered.colMetadata$Subcluster)
bb=as.character(filtered.colMetadata$Subcluster)
aa=bb[which(bb %in% c('Mic0','Mic1','Mic2','Mic3'))]
submic=mathysobj[,filtered.colMetadata$Subcluster %in% c('Mic0','Mic1','Mic2','Mic3')]
submic <- FindNeighbors(submic, dims = 1:10)
submic <- FindClusters(submic, resolution = 0.5)
submic <- RunUMAP(submic, dims = 1:10)
Idents(submic) <- as.factor(aa)
DimPlot(submic, reduction = "umap",label = submic$ident)
submic0.markers <- FindMarkers(submic, ident.1 = 'Mic0', min.pct = 0.25)
submic1.markers <- FindMarkers(submic, ident.1 = 'Mic1', min.pct = 0.25)
submic2.markers <- FindMarkers(submic, ident.1 = 'Mic2', min.pct = 0.25)
submic3.markers <- FindMarkers(submic, ident.1 = 'Mic3', min.pct = 0.25)

submic.markers <- FindAllMarkers(submic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
submic.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



exdeg=data.frame(gene=exmk$gene[1:100],cluster=rep('Ex',100),avg_logFC=exmk$avg_logFC[1:100])
indeg=data.frame(gene=inmk$gene[1:100],cluster=rep('In',100),avg_logFC=inmk$avg_logFC[1:100])
astdeg=data.frame(gene=astmk$gene[1:100],cluster=rep('Ast',100),avg_logFC=astmk$avg_logFC[1:100])
olideg=data.frame(gene=olimk$gene[1:100],cluster=rep('Oli',100),avg_logFC=olimk$avg_logFC[1:100])
micdeg=data.frame(gene=micmk$gene[1:100],cluster=rep('Mic',100),avg_logFC=micmk$avg_logFC[1:100])
opcdeg=data.frame(gene=opcmk$gene[1:100],cluster=rep('OPC',100),avg_logFC=opcmk$avg_logFC[1:100])
enddeg=data.frame(gene=endmk$gene[1:100],cluster=rep('End',100),avg_logFC=endmk$avg_logFC[1:100])

mathysdeg=rbind(exdeg,indeg)
mathysdeg=rbind(mathysdeg,astdeg)
mathysdeg=rbind(mathysdeg,olideg)
mathysdeg=rbind(mathysdeg,micdeg)
mathysdeg=rbind(mathysdeg,opcdeg)
mathysdeg=rbind(mathysdeg,enddeg)

res_marker <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, mathysdeg, weighted=F)
marker.based.estimates <- res_marker$bulk.props
markProp=t(marker.based.estimates)
write.csv(markProp,'Bisque_markProp.csv')
