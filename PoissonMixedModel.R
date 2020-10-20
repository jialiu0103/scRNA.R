#The model with RUV. RUV covariates were added in the regression.
#indgroup is a 2-cols table containing 0/1 as indicator of AD and projid indicates 48 individuals

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
indgroup <- read_excel("indgroup.xlsx")

###collect ast cells into each individual
data=dat
colnames(data)=filtered.colMetadata$broad.cell.type
data=rbind(filtered.colMetadata$projid,data)
ast_data=data[,colnames(data)=='Ast']
colnames(ast_data)=ast_data[1,]
ast_data=ast_data[-1,]
pid=indgroup$projid
rownames(ast_data)=genenames
#here we have 17926 genes in ast_data
ast_ind=data.frame(projectid=rep(999,17926))
check=99999
#aggregating all ast cells within each individual 
for (i in pid){
  ind2=ast_data[,colnames(ast_data)==i]
  check=append(check,length(colnames(ind2)))
  dat_ind2=apply(ind2,1,sum)
  dat_ind2=as.data.frame(dat_ind2)
  colnames(dat_ind2)=i
  ast_ind=cbind(ast_ind,dat_ind2)
}
ast_ind=ast_ind[,-1]
#extract gene APOE
apoe=ast_ind['APOE',]

#check
# ind3=ast_dat[,colnames(ast_dat)==11409232]
# dat_ind3=apply(ind3,1,sum)

#create a data frame for apoe. includes: projid, gene expression value, ad indicator 0/1
apoe_ast=melt(apoe)
apoe_ast$projid=as.character(apoe_ast$variable)
apoe_ast$projid=as.numeric(apoe_ast$projid)
apoe_ast=inner_join(apoe_ast, indgroup, by = 'projid')

###RUV
#ast_ind is an individual-level count matrix by aggregating all cells within ast cell type and each individual.

##filtering
d_e <- DGEList(ast_ind, genes=rownames(ast_ind))
keep <- rowSums(cpm(d_e)>1) >= 3
d_e <- d_e[keep, , keep.lib.sizes=FALSE]

##normalization
#Calculate normalization factors to scale the raw library sizes.
#'TMM': weighted trimmed mean of M-values. the weights are from the delta method on Binomial data.
#Can be dominated by a few highly expressed genes(don't use raw library size itself for normalization)
d_e <- calcNormFactors(d_e, method="TMM")
design=model.matrix(~apoe_ast$ad)

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

#calculate lib for induviduals
ast_lib=colSums(ast_ind)

###poisson mixed model
#add factors of unwanted variation W as covariates in glmer model
#in APOE, AST. paper:z=-1.956695178, p=0.050383317
ast_ruv=ruv_cov$W[match(pid, ds$projid),1:10]#factors of unwanted variation: ruv_cov$W
re=glmer(formula = value ~ ad + ast_ruv[,1]+ast_ruv[,2]+ast_ruv[,3]+ast_ruv[,4]+ast_ruv[,5]+ast_ruv[,6]+ast_ruv[,7]+ast_ruv[,8]+ast_ruv[,9] +ast_ruv[,10] + offset(log(ast_lib))+ (1 | projid),data=apoe_ast, family = poisson, nAGQ=10)
summary(re)#my result: -1.955 0.050568

###poisson model on individual level
ast_ruv=ruv_cov$W[match(pid, ds$projid),1:10]#factors of unwanted variation: ruv_cov$W
m1 <- glm(value ~ ad+ ast_ruv[,1]+ast_ruv[,2]+ast_ruv[,3]+ast_ruv[,4]+ast_ruv[,5]+ast_ruv[,6]+ast_ruv[,7]+ast_ruv[,8]+ast_ruv[,9] +ast_ruv[,10] + offset(log(ast_lib)),data=apoe_ast,family="poisson")
summary(m1)#my result:z=-6.925 p=4.36e-12
