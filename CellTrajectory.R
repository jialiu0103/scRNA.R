#mathys, cell trajectory for mic. package: slingshot
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

dat <- readMM("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_count_matrix.mtx")
genenames <- readLines("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_column_metadata.txt")
indgroup <- read_excel("indgroup.xlsx")
data=dat
colnames(data)=filtered.colMetadata$projid

#To prevent endothelial cells and per cells from being mistaken as very differentiated cell types derived from neural stem cells,
#only keep 6 cell types.
#ind=filtered.colMetadata$Subcluster == 'Mic0'
# cells_use=filtered.colMetadata$projid[ind]
# data=data[,cells_use]
# rownames(data)=genenames
# cells_porjid=filtered.colMetadata$projid[ind]
data=dat
colnames(data)=filtered.colMetadata$Subcluster
data=rbind(filtered.colMetadata$projid,data)
celltype_data=data[,colnames(data)=='Mic0']
cname=as.character(celltype_data[1,])
cname=as.numeric(cname)
#find special projid in subtype, try to modify
colnames(celltype_data)=cname
celltype_data=celltype_data[-1,]
rownames(celltype_data)=genenames
pid=pathologygroup$projid
cname=as.data.frame(cname)
cname$projid=cname[,1]
adstage=left_join(cname,pathologygroup)
#here we have 17926 genes in data
#celltype_ind=data.frame(projectid=rep(999,17926))
#aggregating all cells in a specific celltype within each individual 
# for (i in pid){
#   if (i %in% colnames(celltype_data)){
#     ind2=celltype_data[,colnames(celltype_data)==i]
#     ind2=as.data.frame(ind2)
#     if (length(ind2)>1){
#       dat_ind2<-apply(ind2,1,sum)
#     } else {
#       dat_ind2<-ind2}
#     dat_ind2=as.data.frame(dat_ind2)
#     colnames(dat_ind2)=i
#     celltype_ind=cbind(celltype_ind,dat_ind2)
#   }
# }
# celltype_ind=celltype_ind[,-1]

###preprocessing

##QC
seu=CreateSeuratObject(celltype_data) %>% 
  SCTransform()#normalize and scale
# Add cell type annotation to metadata
seu <- AddMetaData(seu, setNames(filtered.colMetadata$Subcluster[ind], cells_use), 
                   col.name = "cell_type")
VlnPlot(seu, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 1, group.by = "cell_type")

##Dimension reduction
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
ElbowPlot(seu, ndims = 70)

# Need to use DimPlot due to weird workflowr problem with PCAPlot that calls seu[[wflow.build]]
# and eats up memory. I suspect this is due to the sys.call() in 
# Seurat:::SpecificDimPlot. 
DimPlot(seu, reduction = "pca",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE) +
  scale_color_brewer(type = "qual", palette = "Set2")

#tsne
seu <- RunTSNE(seu, dims = 1:30, verbose = FALSE,perplexity = 50)
DimPlot(seu, reduction = "tsne",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE) +
  scale_color_brewer(type = "qual", palette = "Set2")

#umap
seu <- RunUMAP(seu, dims = 1:30, seed.use = 4867,k=8)
DimPlot(seu, reduction = "umap",
        group.by = "cell_type", pt.size = 0.5, label = TRUE, repel = TRUE) +
  scale_color_brewer(type = "qual", palette = "Set2")

##find subcluater
names(seu@meta.data)
seu <- FindNeighbors(seu, verbose = FALSE, dims = 1:30)
seu <- FindClusters(seu, random.seed = 256, resolution = 1)
DimPlot(seu, pt.size = 0.5, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(seu, pt.size = 0.5, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)
DimPlot(seu, pt.size = 0.5, reduction = "pca", group.by = "seurat_clusters", label = TRUE)
#trajectory inference
sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, stretch = 0)
#visualization
#' Assign a color to each cell based on some value
#' 
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# cell_colors <- cell_pal(seu$cell_type, brewer_pal("qual", "Set2"))
# cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())
# list=999
# for (i in adstage$pathology...5){if (i==0){list=c(list,2)} else if (i==1){list=c(list,1)} else if (i==5){list=c(list,3)}}
# mycolor=list[-1]
ad_colors <- cell_pal(as.factor(adstage$pathology...5), brewer_pal("div", "Set1"))
# plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
# lines(sds, lwd = 2, type = 'lineages', col = 'black')
ind_colors<-cell_pal(as.factor(adstage$projid), hue_pal())
plot(reducedDim(sds), col = ad_colors, pch = 16, cex = 1.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')

plot(reducedDim(sds), col = ad_colors, pch = 16, cex = 1.5)
lines(sds, lwd = 2, col = 'black')

plot(reducedDim(sds), col = ind_colors, pch = 16, cex = 1.5)
lines(sds, lwd = 2, col = 'black')
#############################random walk
# ls=999
# for (i in pathologygroup$pathology...5){if (i==0){ls<-c(ls,'yellow')} else if (i==1) {ls<-c(ls,'green')} else if (i==2) {ls<-c(ls,'red')}}
# ls=ls[-1]
# #use kmeans to cluster
# km <- kmeans(sds@reducedDim, center = 30)
# plot(km$centers,col=)
# plot(reducedDim(sds), col = ad_colors, pch = 20, cex = 1)
# lines(sds, lwd = 2, type = 'lineages', col = 'black')
# legend(10, 35, legend=c("ad", "control"),
#        col=c("red", "green"), lty=1:2, cex=0.8)
#############################3

#plot the pseudotime values for each lineage.
nc <- 2
pt <- slingPseudotime(sds)
mypt <- as.matrix(cells_ad$pathology...5)
colnames(mypt)[1]<-'curve1'
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}

#differential expression genes
#look at which genes are differentially expressed along one of those lineages
# Get top highly variable genes
top_hvg <- HVFInfo(seu) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(residual_variance)) %>% 
  top_n(100, residual_variance) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(seu, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,2], dat_use) # Do curve 2, so 2nd columnn
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

#The subset of data is randomly split into training and validation; 
#the model fitted on the training set will be evaluated on the validation set.
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)

#fit model
model <- rand_forest(mtry = 20, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)

#ranger
val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
  select(truth = pseudotime, estimate)
#The model is evaluated on the validation set with 3 metrics: room mean squared error (RMSE), 
#coefficient of determination using correlation (rsq, between 0 and 1), and mean absolute error (MAE).
metrics(data = val_results, truth, estimate)
summary(dat_use_df$pseudotime)

#plot some genes deemed the most important to predicting pseudotime
var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)[1:6]
par(mfrow = c(3, 2))
for (i in seq_along(top_genes)) {
  colors <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]
  plot(reducedDim(sds), col = colors, 
       pch = 16, cex = 0.5, main = top_genes[i])
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}
