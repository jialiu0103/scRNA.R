require(tidyverse)
require(Seurat)
filtered.colMetadata <- read.delim("~/filtered_column_metadata.txt")
mathysobj <- readRDS("~/DeconvolutionMethods/mathysobj.rds")
mathys.celltype=filtered.colMetadata$broad.cell.type
DimPlot(mathysobj, reduction = "umap")
#set ident as cell types provided by the metadata
Idents(object = mathysobj) <- mathys.celltype
#marks = list(Ex='NRGN',In='GAD1',Ast='AQP4',Oli='MBP',Mic=c('CD74',CSF1R'),End='FLT1',OPC='VCAN',Per='AMBP')
FeaturePlot(mathysobj, features = 'VCAN',label = TRUE)#use marker genes provided by the paper to check the cell type annotation


#find markers
# find all markers 
mic.markers <- FindMarkers(mathysobj, ident.1 = 'Mic', ident.2 = c('Ast','End','Ex','In','Oli','Opc','Per'), min.pct = 0.5)
ast.markers <- FindMarkers(mathysobj, ident.1 = 'Ast', ident.2 = c('Mic','End','Ex','In','Oli','Opc','Per'), min.pct = 0.5)
end.markers <- FindMarkers(mathysobj, ident.1 = 'End', ident.2 = c('Ast','Mic','Ex','In','Oli','Opc','Per'), min.pct = 0.5)
oli.markers <- FindMarkers(mathysobj, ident.1 = 'Oli', ident.2 = c('Ast','End','Ex','In','Mic','Opc','Per'), min.pct = 0.5)
ex.markers <- FindMarkers(mathysobj, ident.1 = 'Ex', ident.2 = c('Ast','End','Mic','In','Oli','Opc','Per'), min.pct = 0.5)
in.markers <- FindMarkers(mathysobj, ident.1 = 'In', ident.2 = c('Ast','End','Mic','Ex','Oli','Opc','Per'), min.pct = 0.5)
per.markers <- FindMarkers(mathysobj, ident.1 = 'Per', ident.2 = c('Ast','End','Mic','In','Oli','Opc','Ex'), min.pct = 0.5)
opc.markers <- FindMarkers(mathysobj, ident.1 = 'In', ident.2 = c('Ast','End','Mic','Ex','Oli','In','Per'), min.pct = 0.5)

#Define marker genes as genes with adjusted p-value of <0.01, and a log2 mean fold change >0.5. 
mic.markers=mic.markers[mic.markers$avg_logFC>0.5,]#num of markers for each celll type
mic.markers=mic.markers[mic.markers$p_val_adj<0.01,]#mic: 42
ast.markers=ast.markers[ast.markers$avg_logFC>0.5,]
ast.markers=ast.markers[ast.markers$p_val_adj<0.01,]#ast: 
end.markers=end.markers[end.markers$avg_logFC>0.5,]
end.markers=end.markers[end.markers$p_val_adj<0.01,]#end: 
oli.markers=oli.markers[oli.markers$avg_logFC>0.5,]
oli.markers=oli.markers[oli.markers$p_val_adj<0.01,]#oli: 
ex.markers=ex.markers[ex.markers$avg_logFC>0.5,]
ex.markers=ex.markers[ex.markers$p_val_adj<0.01,]#ex: 
in.markers=in.markers[in.markers$avg_logFC>0.5,]
in.markers=in.markers[in.markers$p_val_adj<0.01,]#in: 
opc.markers=opc.markers[opc.markers$avg_logFC>0.5,]
opc.markers=opc.markers[opc.markers$p_val_adj<0.01,]#opc: 
per.markers=per.markers[per.markers$avg_logFC>0.5,]
per.markers=per.markers[per.markers$p_val_adj<0.01,]#per: 

save(mic.markers,ast.markers,oli.markers,in.markers,ex.markers,end.markers,per.markers,opc.markers, file = "mathys_marker.Rdata")




