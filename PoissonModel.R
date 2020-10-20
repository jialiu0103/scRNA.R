#The model with RUV. RUV covariates were added in the regression.
#indgroup is a table containing 0/1 as indicator of AD and projid indicates 48 individuals

library(DESeq2)
library(RUVSeq)
library(lme4)
library(Matrix)
library(readxl)
library(tidyverse)
library(reshape2)
require(sandwich)
require(msm)

dat <- readMM("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_count_matrix.mtx")
genenames <- readLines("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_gene_row_names.txt")
filtered.colMetadata <- read.delim("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/snRNA-seq/filtered_column_metadata.txt")
indgroup <- read_excel("/Users/liujia/Dropbox (Partners HealthCare)/CATS-OMICS/Projects/Deconvolution/Data/ROSMAP/mathys_pathologygroup.xlsx")

###collect MIC cells into each individual
data=dat
colnames(data)=filtered.colMetadata$broad.cell.type
data=rbind(filtered.colMetadata$projid,data)
ct_data=data[,colnames(data)=='Mic']
colnames(ct_data)=ct_data[1,]
ct_data=ct_data[-1,]
pid=indgroup$projid
rownames(ct_data)=genenames
#here we have 17926 genes in data
ct_ind=data.frame(projectid=rep(999,17926))
#aggregating all MIC cells within each individual 
for (i in pid){
  ind2=ct_data[,colnames(ct_data)==i]
  dat_ind2=apply(ind2,1,sum)
  dat_ind2=as.data.frame(dat_ind2)
  colnames(dat_ind2)=i
  ct_ind=cbind(ct_ind,dat_ind2)
}
ct_ind=ct_ind[,-1]
#extract gene APOE
apoe=ct_ind['APOE',]


#create a data frame for apoe. includes: projid, gene expression value, ad indicator 0/1
apoe_ct=melt(apoe)
apoe_ct$projid=as.character(apoe_ct$variable)
apoe_ct$projid=as.numeric(apoe_ct$projid)
apoe_ct=inner_join(apoe_ct, indgroup, by = 'projid')

###RUV
#mic_ind is an individual-level count matrix by aggregating all cells within mic cell type and each individual.

##filtering
d_e <- DGEList(ct_ind, genes=rownames(ct_ind))
keep <- rowSums(cpm(d_e)>1) >= 3
d_e <- d_e[keep, , keep.lib.sizes=FALSE]

##normalization
#Calculate normalization factors to scale the raw library sizes.
#'TMM': weighted trimmed mean of M-values. the weights are from the delta method on Binomial data.
#Can be dominated by a few highly expressed genes(don't use raw library size itself for normalization)
d_e <- calcNormFactors(d_e, method="TMM")
design=model.matrix(~apoe_ct$ad)

##compute the residuals from the GLM fit
#Estimates a common negative binomial dispersion parameter for a DGE dataset with a general experimental design.
d_e <- estimateGLMCommonDisp(d_e, design)
#Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each tag
d_e <- estimateGLMTagwiseDisp(d_e, design)
#Fit a negative binomial generalized log-linear model to the read counts for each gene
fit1 <- glmFit(d_e, design)
#residuals() extracts model residuals from objects returned by modeling functions.
res1 <- residuals(fit1, type="deviance")

##Remove Unwanted Variation Using Residuals
ruvn <- 10#number of factors of unwanted variation to be estimated from the data
#RUVr uses residuals, (e.g:from a first-pass GLM regression of the counts on the covariates of interest.)
#to estimate the factors of unwanted variation W
ruv_cov <- RUVr(round(d_e$counts), as.character(rownames(d_e$counts)), k=ruvn, res1)

#calculate lib for individuals
ct_lib=colSums(ct_ind)

###poisson mixed model
#add factors of unwanted variation W as covariates in glmer model
#in APOE, MIC. paper:z=3.288 p=0.001008
ct_ruv=ruv_cov$W#factors of unwanted variation: ruv_cov$W
re=glmer(formula = value ~ ad + ct_ruv[,1]+ct_ruv[,2]+ct_ruv[,3]+ct_ruv[,4]+ct_ruv[,5]+ct_ruv[,6]+ct_ruv[,7]+ct_ruv[,8]+ct_ruv[,9] +ct_ruv[,10] + offset(log(ct_lib))+ (1 | projid),data=apoe_ct, family = poisson, nAGQ=10)
summary(re)#my result:z=3.288 p=0.001008

###poisson model on individual level
ct_ruv=ruv_cov$W#factors of unwanted variation: ruv_cov$W
m1 <- glm(value ~ ad+ ct_ruv[,1]+ct_ruv[,2]+ct_ruv[,3]+ct_ruv[,4]+ct_ruv[,5]+ct_ruv[,6]+ct_ruv[,7]+ct_ruv[,8]+ct_ruv[,9] +ct_ruv[,10] + offset(log(ct_lib)),data=apoe_ct,family="poisson")
summary(m1)#my result:z=6.021 p=1.73e-09
