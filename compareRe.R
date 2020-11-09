#import proportion results from excel
stares=function(csvfile,celltypename=c('mic','ast','oli')){
  resulttb <- read_csv(csvfile)
  resulttb=as.data.frame(resulttb)
  rownames(resulttb)=resulttb[,1]
  resulttb=resulttb[,-1]
  #controlname is defined in the note at the bottom
  ctres=resulttb[controlname,]
  ctres=ctres[,celltypename]  
  ctre=apply(ctres, 2, mean)
  sub3=sum(ctre)
  micmean=ctre[celltypename[1]]/sub3
  astmean=ctre[celltypename[2]]/sub3
  olimean=ctre[celltypename[3]]/sub3
  meanre=c(micmean,astmean,olimean)
  
  #truere
  tmic=0.0585
  tast=0.1875
  toli=0.751
  
  tre=c(tmic,tast,toli)
  
  #start to compare
  #1. chi-square test
  ca = matrix(c(micmean,astmean,olimean,tmic,tast,toli), nrow = 2,byrow = T)
  csq=chisq.test(ca)#X-squared = 0.1176, df = 2, p-value = 0.9429
  
  #2. t-test
  ctres$sum=apply(ctres,1,sum)
  olit=ctres$oli/ctres$sum
  oli_ttest=t.test(olit,mu=toli)#t = 1.4324, df = 200, p-value = 0.1536
  astt=ctres$ast/ctres$sum
  ast_ttest=t.test(astt,mu=tast)#t = -215.05, df = 200, p-value < 2.2e-16
  mict=ctres$mic/ctres$sum
  mic_ttest=t.test(mict,mu=tmic)#t = 16.301, df = 200, p-value < 2.2e-16
  
  return(list(csq,mic_ttest,ast_ttest,oli_ttest,summary(data.frame(mict,astt,olit)),meanre,mict,astt,olit))
}

#dtangle1(dtangle_grubman'sMK)
csvfile='~/DeconvolutionMethods/resulttables/dtangle1.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
dtangle1=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#dtangle2(dtangle_dtangle'sMK)
csvfile='~/DeconvolutionMethods/resulttables/dtangle2.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
dtangle2=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#dtangle3(dtangle_mathys'sMK)
csvfile='~/DeconvolutionMethods/resulttables/dtangle3.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
dtangle3=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#dtangle4(dtangle_identcellMK)
csvfile='~/DeconvolutionMethods/resulttables/dtangle4.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
dtangle4=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))

#music1(music_allgene)
csvfile='~/DeconvolutionMethods/resulttables/music1.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
MuSiC_basic1=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#music2(music_grubmanMK)
csvfile='~/DeconvolutionMethods/resulttables/music2.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
MuSiC_basic2=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#music3(music_identcell)
csvfile='~/DeconvolutionMethods/resulttables/music3.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
MuSiC_basic3=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#music4(music_all_MK_pregroup)
csvfile='~/DeconvolutionMethods/resulttables/music4.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal (pmic_t=0.0001037)
MuSiC_pregroup1=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#music5(music_grubmanMK_pregroup)
csvfile='~/DeconvolutionMethods/resulttables/music5.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast<0.05 -- not equal toli=0.1536-- equal
MuSiC_pregroup2=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))

#ctx1(ctx_mkgene_Bmode)
csvfile='~/DeconvolutionMethods/resulttables/ctx1.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Bmode1=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx2(ctx_mkgene_Smode)
csvfile='~/DeconvolutionMethods/resulttables/ctx2.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Smode1=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx3(ctx_mk_ind_Bmode)
csvfile='~/DeconvolutionMethods/resulttables/ctx3.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Bmode2=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx4(ctx_mk_ind_Smode)
csvfile='~/DeconvolutionMethods/resulttables/ctx4.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Smode2=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx5(ctx_allgene_Bmode)
csvfile='~/DeconvolutionMethods/resulttables/ctx5.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Bmode3=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx6(ctx_allgene_Smode)
csvfile='~/DeconvolutionMethods/resulttables/ctx6.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Smode3=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx7(ctx_allgene_ind_Bmode)
csvfile='~/DeconvolutionMethods/resulttables/ctx7.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Bmode4=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx8(ctx_allgene_ind_Smode)
csvfile='~/DeconvolutionMethods/resulttables/ctx8.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Smode4=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx9(ctx_allgene_ctind_Bmode)
csvfile='~/DeconvolutionMethods/resulttables/ctx9.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Bmode5=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))
#ctx10(ctx_allgene_ctind_Smode)
csvfile='~/DeconvolutionMethods/resulttables/ctx10.csv'
aa=stares(csvfile)#p(chisq)>0.05. tmic,tast,toli<0.05 -- not equal
CibersortX_Smode5=data.frame(mic=unlist(aa[7]),ast=unlist(aa[8]),oli=unlist(aa[9]))

#plot
celltypetb=function(n){
  aa=as.data.frame(c(dtangle1[,n],dtangle2[,n],dtangle3[,n],dtangle4[,n],MuSiC_basic1[,n],MuSiC_basic2[,n],
                     MuSiC_basic3[,n],MuSiC_pregroup1[,n],MuSiC_pregroup2[,n],CibersortX_Bmode1[,n],CibersortX_Bmode2[,n],
                     CibersortX_Bmode3[,n],CibersortX_Bmode4[,n],CibersortX_Bmode5[,n],CibersortX_Smode1[,n],CibersortX_Smode2[,n],
                     CibersortX_Smode3[,n],CibersortX_Smode4[,n],CibersortX_Smode5[,n]))
  names(aa)='value'
  aa$name=c(rep('dtange1',201),rep('dtange2',201),rep('dtange3',201),rep('dtange4',201),
            rep('MuSiC_basic1',201),rep('MuSiC_basic2',201),rep('MuSiC_basic3',201),rep('MuSiC_pregroup1',201),
            rep('MuSiC_pregroup2',201),rep('CibersortX_Bmode1',201),rep('CibersortX_Bmode2',201),
            rep('CibersortX_Bmode3',201),rep('CibersortX_Bmode4',201),rep('CibersortX_Bmode5',201),rep('CibersortX_Smode1',201),
            rep('CibersortX_Smode2',201),rep('CibersortX_Smode3',201),rep('CibersortX_Smode4',201),rep('CibersortX_Smode5',201))
  return(aa)
}
tbmic=celltypetb(1)
tbast=celltypetb(2)
tboli=celltypetb(3)
# tmic=0.0585
# tast=0.1875
# toli=0.751
tboli %>% 
  ggplot(., aes(name, value)) +
  geom_violin(aes(fill = name, colour = name)) +
  geom_boxplot(fill = NA) +
  geom_hline(aes(yintercept = 0.751),colour='red',linetype='dashed')+
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Estimated Oli Proportions",x='')+
  coord_flip()













#note for get control sample names
#process bulk sample
rowbulk=read.table(file='~/ROSMAP_RNAseq_FPKM_gene.tsv', sep = '\t', header = TRUE)
bulkmeta=read.csv('ROSMAP_IDkey.csv')
keyID=bulkmeta
sampleIDs <- names(rowbulk)[-(1:2)]
keyID <- keyID %>% dplyr::select(SubjectID = projid, SampleID = rnaseq_id) %>% dplyr::filter(SampleID != "")
keyID <- unique(keyID)
#Remove duplicated samples
duplicateSampleIDs <- c("492_120515_0","492_120515_6")
keyID <- dplyr::filter(keyID, !(SampleID %in% duplicateSampleIDs))#sample-id and projid
data <- rowbulk[,!colnames(rowbulk) %in% duplicateSampleIDs]
sampleIDs <- names(data)[-(1:2)]
#In the keyID file the latest sequencing sample IDs (from plate 7 - 8) lack the "_X" tag present in the expression dataset
removeSampleIDTag <- function(sampleID) {
  gsub("^([^_]+_[0-9]+)(_[0-9]+)?$", "\\1", sampleID)
}
keyID$SampleID <- sapply(keyID$SampleID, removeSampleIDTag)
sampleIDs <- sapply(sampleIDs, removeSampleIDTag)
length(intersect(keyID$SampleID, names(data)[-(1:2)]))# Check full correspondance
#Extract subjectIDs in the same order as in the data
subjectIDs <- sapply(sampleIDs, function(x) {keyID[keyID$SampleID == x,]$SubjectID})
#Create a matrix with the data ready for limma analysis
mData <- as.matrix(data[,-(1:2)])
colnames(mData) <- names(subjectIDs)
#2. Process covariate file
covFile <- read.csv('~/ROSMAP_Clinical_2019-05-03.csv')
cov <- covFile[,c("projid", "cts_mmse30_lv", "braaksc", "ceradsc", "cogdx")]
setnames(cov, c("projid", "cts_mmse30_lv", "braaksc", "ceradsc", "cogdx"), c("SubjectID", "MMSE", "Braak", "CERAD", "ClinicalDiagnosis"))

# ROSMAP ClinicalDiagnosis levels: 1 = NCI, 2 = MCI, 3 = MCI + other cause of CI, 
#4 = AD, 5 = AD + other cause of CI, 6 = Non-AD Dementia
#disease:2,3,4,5,6
#normal:1
cov_ct=cov[cov$ClinicalDiagnosis==1,]
cov_ct=inner_join(cov_ct,keyID,by='SubjectID')#201 control samples
#convert sampleid in bulk data
bulkid <- sapply(colnames(rowbulk)[-c(1,2)], removeSampleIDTag)
bulkid=substr(bulkid,2,20)
bulkname=colnames(rowbulk)[-c(1,2)]
#bulk is gene symble converted bulk data from music.r
colnames(bulk)=bulkid
bulk[1,]=bulkname
controlbulk=bulk[,cov_ct$SampleID]
controlname=controlbulk[1,]
controlname=as.character(controlname)
