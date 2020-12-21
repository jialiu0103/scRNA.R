#x as cell types, y as log2value
colnames(mygenedt)=filtered.colMetadata$broad.cell.type
data=rbind(filtered.colMetadata$projid,mygenedt)

getctfunc=function(celltype,data){
  ct_data=data[,colnames(data)==celltype]
  colnames(ct_data)=ct_data[1,]
  ct_data=ct_data[-1,]
  ct_ind=data.frame(projectid=rep(999,20))
  for (i in pid){
    ind2=ct_data[,colnames(ct_data)==i]
    dat_ind2=apply(ind2,1,sum)
    dat_ind2=as.data.frame(dat_ind2)
    colnames(dat_ind2)=i
    ct_ind=cbind(ct_ind,dat_ind2)
  }
  ct_ind=ct_ind[,-1]
  alldt=melt(t(ct_ind))
  pidlist=as.character(alldt$Var1)
  alldt.group=ifelse(groupinfo[pidlist,2]==0,'Control','AD')
  alldt$group=alldt.group
  colnames(alldt)=c('projid','Gene','value','group')
  alldt$group <- factor(alldt$group,levels = c("Control", "AD"))
  alldt$CellType=rep(celltype)
  return(alldt)
}

astplot=getctfunc('Ast',data)
explot=getctfunc('Ex',data)
inplot=getctfunc('In',data)
micplot=getctfunc('Mic',data)
oliplot=getctfunc('Oli',data)
opcplot=getctfunc('Opc',data)
allctdt=rbind(astplot,micplot,oliplot,explot,inplot,opcplot)

magenes=c('CR1','BIN1','HLA-DRB1','TREM2','CD2AP','NYAP1',
          'EPHA1','PTK2B','CLU','SPI1','MS4A2','PICALM','SORL1','FERMT2',
          'SLC24A4','ABCA7','APOE','CASS4','ECHDC3','ACE')

onegene=allctdt[allctdt$Gene=='CR1',]
onegene=allctdt[allctdt$Gene=='BIN1',]
onegene=allctdt[allctdt$Gene=='HLA-DRB1',]
onegene=allctdt[allctdt$Gene=='TREM2',]
onegene=allctdt[allctdt$Gene=='CD2AP',]
onegene=allctdt[allctdt$Gene=='NYAP1',]
onegene=allctdt[allctdt$Gene=='EPHA1',]
onegene=allctdt[allctdt$Gene=='PTK2B',]
onegene=allctdt[allctdt$Gene=='CLU',]
onegene=allctdt[allctdt$Gene=='SPI1',]
onegene=allctdt[allctdt$Gene=='MS4A2',]
onegene=allctdt[allctdt$Gene=='PICALM',]
onegene=allctdt[allctdt$Gene=='SORL1',]
onegene=allctdt[allctdt$Gene=='FERMT2',]
onegene=allctdt[allctdt$Gene=='SLC24A4',]
onegene=allctdt[allctdt$Gene=='ABCA7',]
onegene=allctdt[allctdt$Gene=='APOE',]
onegene=allctdt[allctdt$Gene=='CASS4',]
onegene=allctdt[allctdt$Gene=='ECHDC3',]
onegene=allctdt[allctdt$Gene=='ACE',]
onegene$log2value=log2(onegene$value+1)
#onegene$log2value=ifelse(abs(onegene$value)>1,log2(onegene$value),0)
ggplot(onegene, aes(x=CellType, y=log2value, fill=group)) + 
  geom_boxplot()


###mathys
colnames(mygenedt)=filtered.colMetadata$broad.cell.type
data=rbind(filtered.colMetadata$projid,mygenedt)
ct_data=data[,colnames(data)=='Ex']
colnames(ct_data)=ct_data[1,]
ct_data=ct_data[-1,]
ct_ind=data.frame(projectid=rep(999,20))
for (i in pid){
  ind2=ct_data[,colnames(ct_data)==i]
  dat_ind2=apply(ind2,1,sum)
  dat_ind2=as.data.frame(dat_ind2)
  colnames(dat_ind2)=i
  ct_ind=cbind(ct_ind,dat_ind2)
}
ct_ind=ct_ind[,-1]
alldt=melt(t(ct_ind))
pidlist=as.character(alldt$Var1)
alldt.group=ifelse(groupinfo[pidlist,2]==0,'Control','AD')
alldt$group=alldt.group
colnames(alldt)=c('projid','Gene','value','group')
alldt$group <- factor(alldt$group,levels = c("Control", "AD"))
