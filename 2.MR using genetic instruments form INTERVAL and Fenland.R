###### MR analysis: estimate the effect of protein level of BGAT on pregnancy outcomes using genetic instruments from INTERVAL and Fenland #####
###part 1:using genetic instruments from INTERVAL###
###part 2:using genetic instruments from Fenland###

library(readxl)
library(writexl)
library(data.table)
library(TwoSampleMR)
library(MendelianRandomization)

###part 1:using genetic instruments from INTERVAL###

#read exposure data from INTERVAL#
exposure_xlsx <-read_excel("INTERVAL.xlsx")
exposure <- format_data(exposure_xlsx,
                        phenotype_col = "Target",
                        type = "exposure",
                        header = TRUE,
                        snp_col = "Conditional variant",
                        beta_col = "Univariable beta",
                        se_col = "Univariable SE",
                        eaf_col = "Conditional variant EAF",
                        effect_allele_col = "Conditional variant EA",
                        other_allele_col = "Conditional variant OA",
                        pval_col = "Univariable p",
                        chr_col = "Conditional variant Chr",
                        pos_col = "Conditional variant Pos",
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
#find proxy snp for rs505922#
outcome<- outcome[which(outcome$SNP=="rs576123"),]
outcome[1,5]="rs505922"
#harmonisation#
dat_harmonised <- harmonise_data(exposure_dat = exposure,
                                 outcome_dat = outcome)
#MR#
MR_result <- mr(dat_harmonised)
str(MR_result)


###part 2:using genetic instruments from Fenland###

#read exposure data from Fenland#
exposure_xlsx <-read_excel("Fenland.xlsx")
exposure <- format_data(exposure_xlsx,
                        phenotype_col = "Target",
                        type = "exposure",
                        header = TRUE,
                        snp_col = "rsID",
                        beta_col = "Effect",
                        se_col = "SE",
                        eaf_col = "EAF",
                        effect_allele_col = "EA",
                        other_allele_col = "NEA",
                        pval_col = "P.value",
                        chr_col = "Chr",
                        pos_col = "Position",
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
#find proxy snp for rs576125#
outcome<- outcome[c(which(outcome$SNP=="rs576123"),which(outcome$SNP=="rs2013075")),]
outcome[1,3]="A"
outcome[1,4]="G"
outcome[1,5]="rs576125"
#harmonisation#
dat_harmonised <- harmonise_data(exposure_dat = exposure,
                                 outcome_dat = outcome)
#format harmonised data for MR#
dat_harmonised <- dat_to_MRInput(dat_harmonised, get_correlations = F, pop = "EUR")
#create LD matrix of instruments#
dat_harmonised$`'BGAT.outcome`@correlation <- matrix(c(1,0.398402,0.398402,1),nrow = 2,ncol = 2) 
#MR#
MR_result <- MendelianRandomization:: mr_ivw(dat_harmonised$`'BGAT.outcome`,
                                             model= "default",
                                             robust = FALSE,
                                             penalized = FALSE,
                                             correl = T,
                                             weights ="simple",
                                             psi = 0,
                                             distribution = "normal",
                                             alpha = 0.05)
str(MR_result) 
