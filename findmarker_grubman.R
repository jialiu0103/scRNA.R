#grubman marker
grubmans_seuratobj=readRDS('~/GrubmanReproduce/grubmans_seuratobj.rds')
grubmans_annot <- readRDS("~/GrubmanReproduce/grubmans_annot.rds")

# find all markers 
mic.markers <- FindMarkers(grubmans_seuratobj, ident.1 = 'mic', ident.2 = c('ast','end','neu','oli'), min.pct = 0.5)
ast.markers <- FindMarkers(grubmans_seuratobj, ident.1 = 'ast', ident.2 = c('mic','end','neu','oli'), min.pct = 0.5)
end.markers <- FindMarkers(grubmans_seuratobj, ident.1 = 'end', ident.2 = c('mic','ast','neu','oli'), min.pct = 0.5)
neu.markers <- FindMarkers(grubmans_seuratobj, ident.1 = 'neu', ident.2 = c('mic','ast','end','oli'), min.pct = 0.5)
oli.markers <- FindMarkers(grubmans_seuratobj, ident.1 = 'oli', ident.2 = c('mic','ast','end','neu'), min.pct = 0.5)

#Define marker genes as genes with adjusted p-value of <0.01, and a log2 mean fold change >0.5. 
mic.markers=mic.markers[mic.markers$avg_logFC>0.5,]#num of markers for each celll type
mic.markers=mic.markers[mic.markers$p_val_adj<0.01,]#mic: 50 
ast.markers=ast.markers[ast.markers$avg_logFC>0.5,]
ast.markers=ast.markers[ast.markers$p_val_adj<0.01,]#ast: 115
neu.markers=neu.markers[neu.markers$avg_logFC>0.5,]
neu.markers=neu.markers[neu.markers$p_val_adj<0.01,]#neu: 69
end.markers=end.markers[end.markers$avg_logFC>0.5,]
end.markers=end.markers[end.markers$p_val_adj<0.01,]#end: 12
oli.markers=oli.markers[oli.markers$avg_logFC>0.5,]
oli.markers=oli.markers[oli.markers$p_val_adj<0.01,]#oli: 70

save(mic.markers,ast.markers,oli.markers,neu.markers,end.markers, file = "grubman_marker.Rdata")
