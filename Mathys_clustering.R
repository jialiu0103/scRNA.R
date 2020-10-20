require(Seurat)
require(edgeR)
require(Matrix)

###read in data
dat = readMM("~/filtered_count_matrix.mtx")
rownames(dat) <- readLines("~/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("~/filtered_column_metadata.txt")
colnames(dat) <- filtered.colMetadata$TAG
datobj <- CreateSeuratObject(counts = dat, min.cells = 3, min.features = 200)

###normalization
#as the data has been filtered, here I just perform normalization
scdata=NormalizeData(datobj, normalization.method = "LogNormalize", scale.factor = 10000)

###feature selection
scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scdata), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scdata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

###Scaling 
all.genes <- rownames(scdata)
scdata <- ScaleData(scdata, features = all.genes)

###dimentional reduction
scdata <- RunPCA(scdata)
VizDimLoadings(scdata, dims = 1:2, reduction = "pca")
DimPlot(scdata, reduction = "pca")

###determine dimention
ElbowPlot(scdata)

###cluster
scdata <- FindNeighbors(scdata, dims = 1:15)
scdata <- FindClusters(scdata, resolution = 0.5)

###umap
scdata <- RunUMAP(scdata, dims = 1:15)
umap_scdata=DimPlot(scdata, reduction = "umap",label = TRUE)

###tsne
scdata <- RunTSNE(object = scdata)
tsne_scdata <- DimPlot(object = scdata, reduction = "tsne")


###find cluster biomarkers
# find all markers of cluster 1
cluster1.markers <- FindMarkers(scdata, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
scdata.markers <- FindAllMarkers(scdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2mg <- scdata.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#random walk, some marker genes mentioned in the paper
#Ex4 neurons were marked by LINGO1, RASGEF1B, and SLC26A3---3
#Oli0 cells were marked by CADM2, QDPR, NLGN1, and CRYAB---4
#Ast1 showed preferential expression of GLUL and of the AD risk factor CLU---2
Ex=FeaturePlot(scdata, features = c('SLC17A7','CAMK2A','NRGN'))
Oli=FeaturePlot(scdata, features = c("MBP", "MOBP", 'PLP1'))
Ast=FeaturePlot(scdata, features = c("AQP4", "GFAP"))
Mic=FeaturePlot(scdata, features = c("CD74", "CSF1R",'C3'))
OPC=FeaturePlot(scdata, features = c("PDGFRA", "VCAN",'CSPG4'))
In=FeaturePlot(scdata, features = c("GAD1", 'GAD2'))
Neuron=FeaturePlot(scdata, features = c("SYT1", 'SNAP25','GRIN1'))
End=FeaturePlot(scdata, features = c("FLT1", 'CLDN5'))
FeaturePlot(scdata, features = 'APOE')

#assign cell type to marker genes
new.cluster.ids <- c("Oli", "Oli", "Ex", "Ex", "Ex", "Ex",'Ast','In','Ex','Ex','In',
                     'OPC','Ex','Mix','In',"Ex","Ex",'In',"Ex","Ex",'End')
names(new.cluster.ids) <- levels(scdata)
scdata <- RenameIdents(scdata, new.cluster.ids)
DimPlot(scdata, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


###subcluster of mic
mic<-SubsetData(object = scdata, ident.use="13")

#normalize, scale
mic <- NormalizeData(object = mic)
mic <- ScaleData(object = mic)

#dimentional reduction
mic <- RunPCA(mic)
VizDimLoadings(mic, dims = 1:2, reduction = "pca")
DimPlot(mic, reduction = "pca")

###determine dimention
ElbowPlot(mic)

###cluster
mic <- FindNeighbors(mic, dims = 1:15)
mic <- FindClusters(mic, resolution = 0.5)

###umap
mic <- RunUMAP(mic, dims = 1:15)
umap_mic=DimPlot(mic, reduction = "umap",label = TRUE)

# find markers for every cluster, report only the positive ones
mic.markers <- FindAllMarkers(mic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mic_marker=mic.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#heatmap
top10 <- mic.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(mic, features = top10$gene) + NoLegend()
