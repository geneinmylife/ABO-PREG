###### MR analysis: estimate the effect of protein level of BGAT on pregnancy outcomes using genetic instruments from deCODE #####

library(readxl)
library(writexl)
library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)

#read genetic instruments from deCODE#
exposure <-read_excel("deCODE.xlsx")
exposure <- format_data(exposure,
                        header = TRUE,
                        type = "exposure",
                        phenotype_col = "Exposure", 
                        snp_col = "SNP",
                        beta_col = "beta.exposure",
                        se_col = "se.exposure",
                        eaf_col = "eaf.exposure",
                        effect_allele_col = "effect_allele.exposure",
                        other_allele_col = "other_allele.exposure",
                        pval_col = "pval.exposure",
                        chr_col = "CHR_SNP",
                        pos_col = "POS_SNP",
                        log_pval = FALSE
)
#read outcome data#
outcome <- fread("summary_stats-finngen_R9_T2D.gz", header = T)
outcome <- format_data(outcome,
                       type = "outcome",
                       header = TRUE,
                       snp_col = "rsids",
                       beta_col = "beta",
                       se_col = "sebeta",
                       eaf_col = "af_alt",
                       effect_allele_col = "alt",
                       other_allele_col = "ref",
                       pval_col = "pval",
                       chr_col = "#chrom",
                       pos_col = "pos",
                       log_pval = FALSE
)
#harmonisation#
dat_harmonised <- harmonise_data(exposure_dat = exposure,
                                 outcome_dat = outcome)
#format harmonised data for MR#
dat_harmonised <- dat_to_MRInput(dat_harmonised, get_correlations = F, pop = "EUR")
#create LD matrix of instruments#
dat_harmonised$ABO.outcome@correlation <- matrix(c(1,-0.213,-0.213,1),nrow = 2,ncol = 2)
#MR#
MR_result <- MendelianRandomization:: mr_ivw(dat_harmonised$`ABO.outcome`,
                                           model= "default",
                                           robust = FALSE,
                                           penalized = FALSE,
                                           correl = T,
                                           weights ="simple",
                                           psi = 0,
                                           distribution = "normal",
                                           alpha = 0.05)
str(MR_result) 


