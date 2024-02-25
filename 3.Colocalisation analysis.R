#####colocalisation analysis######

################################Parameters for coloc#####################################
#beta1 - effect sizes for trait 1                                                       #
#beta2 - effect sizes for trait 2                                                       #
#p1 - p-values for trait 1                                                              #
#p2 - p-values for trait 2                                                              #
#MAF1 - minor allele frequencies for SNPs of trait 1                                    #
#MAF2 - minor allele frequencies for SNPs of trait 2                                    #
#N1 - sample size for trait 1                                                           #
#N2 - sample size for trait 2                                                           #
#s - for a case control dataset, the proportion of samples in dataset that are cases    #
#########################################################################################
library(coloc)
coloc.analysis <- function(snp,beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,s){
  #Convert the inputs in order to run in coloc function.
  #type, quant (quantitative) for pQTL study
  dataset1 <- list(snp=snp, beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1)
  #type, cc (case coontrol) for binary study
  dataset2 <- list(snp=snp, beta=beta2, varbeta=se2^2, MAF=MAF2,type="cc",s=s,N=N2)
  #Run the coloc analysis, setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12) as 1E-5.
  #p1 prior probability a SNP is associated with trait 1, default 1e-4
  #p2 prior probability a SNP is associated with trait 2, default 1e-4
  #p12 prior probability a SNP is associated with both traits, default 1e-5
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  #Format the data to save out.
  #List into data frame.
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  #Label the columns in the data frame.
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  return(df)
}

##example: colocalization analysis of protein level of BGAT on T2D##
setwd("/Users/zheng/Desktop/mr/")
library(readxl)
#read merged exposure data and outcome data#
df <- read_excel("input.xlsx")
df <- df[complete.cases(df),]
#get maf data#
df$MAF.exposure <- NULL
df$MAF.outcome <- NULL
for (j in 1:nrow(df)){
  if(df$maf_exposure[j]>0.5){df$MAF_exposure[j]=1-df$maf_exposure[j]} else {df$MAF_exposure[j]=df$maf_exposure[j]}  
  if(df$eaf_outcome[j]>0.5){df$MAF.outcome[j]=1-df$eaf_outcome[j]} else {df$MAF.outcome[j]=df$eaf_outcome[j]}
}
#sample size for two studies
N1<-35373
N2<-57698+308252          
#the proportion of samples in study 2 that are cases
s<- 57698/(57698+308252)
#remove invalid data
df <- df[!duplicated(df$rsids),]
df <- df[which(df$MAF_exposure>0),]
df <- df[which(df$MAF.outcome>0),]
#colocalisation analysis
result <- coloc.analysis(df$rsids,df$beta_exposure, df$beta_outcome, df$se_exposure, df$se_outcome, df$MAF_exposure,df$MAF.outcome, N1, N2, s) 
result
##########Approximate Bayes Factor colocalisation analyses###########
#H0: neither trait has a genetic association in the region          #
#H1:only trait 1 has a genetic association in the region            #
#H2: only trait 2 has a genetic associaiton in the region           #
#H3: both traits are associated, but with different causal variants #
#H4: both traits are associated, and share a single causal variant  #
