#IHC data (benchmark data)
ast_bm <- read.delim("~/DeconvolutionMethods/benchmark_DATA/IHC.astro.txt")
mic_bm <- read.delim("~/DeconvolutionMethods/benchmark_DATA/IHC.microglia.txt")
neu_bm <- read.delim("~/DeconvolutionMethods/benchmark_DATA/IHC.neuro.txt")
oli_bm <- read.delim("~/DeconvolutionMethods/benchmark_DATA/IHC.oligo.txt")
end_bm <- read.delim("~/DeconvolutionMethods/benchmark_DATA/IHC.endo.txt")
ihc_meta <- read_excel("~/DeconvolutionMethods/newprop/IHC_metadata.xlsx")
ihc_meta=ihc_meta[-1,]
colnames(ihc_meta)=ihc_meta[1,]
ihcmeta=data.frame(ihc_meta)
ihcmeta=ihcmeta[-1,]
#CONTROL
ihcmeta=ihcmeta[ihcmeta$clinicalAD==0,]
controlid=ihcmeta$projID
THEid=as.numeric(controlid)
controlid=as.data.frame(THEid)

#AD
ihcmeta=data.frame(ihc_meta)
ihcmeta=ihcmeta[-1,]
ihcAD=ihcmeta[ihcmeta$clinicalAD==1,]
ADid=ihcAD$projID
THEid=as.numeric(ADid)
ADid=as.data.frame(THEid)

all_bm=merge(ast_bm,mic_bm,all = T)
all_bm=merge(all_bm,neu_bm,all = T)
all_bm=merge(all_bm,oli_bm,all = T)
all_bm=merge(all_bm,end_bm,all = T)
rownames(all_bm)=c('ast','mic','neu','oli','end')
benchmark=t(all_bm)
alltp=na.omit(benchmark)#49 samples (ad and control)
alltp=as.data.frame(alltp)
sampleid=substr(rownames(alltp),2,10)
THEid=as.numeric(sampleid)
alltp$THEid=THEid

controlIHC=inner_join(controlid,alltp) #31 control samples have proportions
adIHC=inner_join(ADid,alltp) #18 ad samples have proportions

write.csv(controlIHC,'controlIHC.csv')
write.csv(adIHC,'adIHC.csv')
finalallprop=apply(alltp, 2, mean) 
# ast       mic       neu       oli       end 
# 0.1375989 0.1010994 0.1879012 0.1814057 0.5083996 

fourtp=benchmark[,c('ast','mic','neu','oli')]
fourtp=na.omit(fourtp)
foursum=apply(fourtp, 1, sum)
foursum=as.matrix(foursum)
fourtp=cbind(fourtp,foursum)
ast_4=fourtp[,1]/fourtp[,5]
mic_4=fourtp[,2]/fourtp[,5]
neu_4=fourtp[,3]/fourtp[,5]
oli_4=fourtp[,4]/fourtp[,5]
fourprop=cbind(as.data.frame(ast_4),as.data.frame(mic_4))
fourprop=cbind(fourprop,as.data.frame(neu_4))
fourprop=cbind(fourprop,as.data.frame(oli_4))             
final4prop=apply(fourprop, 2, mean)              

# ast_4     mic_4     neu_4     oli_4 
# 0.2263323 0.1648211 0.3124922 0.2963544 
