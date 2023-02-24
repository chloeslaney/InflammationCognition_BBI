# Project: inflammation and cognition
# Script for two-sample MR - adapted from IEU MR short-course code
# Chloe Slaney 
# R version 4.1.1 & TwoSampleMR _0.5.6

#####################
# PRELIMINARY STEPS #
#####################
library(TwoSampleMR)
library(ieugwasr)
library(tidyverse)
library(readr)
library(dplyr)
library(LDlinkR)
library(utils)
library(httr)
library(rlist)
library(haven)
library(AER)
library("data.table")
library(MRPRESSO)
library(R.utils)
library(rlist)
library(MendelianRandomization)
library(gsmr)
library(survey)

# set working directory
setwd("SET PATH HERE!")

# Clear the work environment
rm(list = ls())

######################################################
# PART 1: Read in exposure and outcome data for 2SMR #
######################################################
EXP_LIGTHART_CIS <- read_exposure_data(filename = "EXP_LIGTHART_CIS_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_LIGTHART_GWAS_CORRECTED <- read_exposure_data(filename = "EXP_LIGTHART_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_HAN_CIS <- read_exposure_data(filename = "EXP_HAN_CIS_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_HAN_GWAS_CORRECTED <- read_exposure_data(filename = "EXP_HAN_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_AHLUWALIA_CIS <- read_exposure_data(filename = "EXP_AHLUWALIA_CIS_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_AHLUWALIA_GWAS_CORRECTED <- read_exposure_data(filename = "EXP_AHLUWALIA_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_BORGES_GWAS_CORRECTED <- read_exposure_data(filename = "EXP_BORGES_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_KETTUNEN_GWAS_CORRECTED <- read_exposure_data(filename = "EXP_KETTUNEN_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_ROSA_INSTRUMENT_CORRECTED <- read_exposure_data(filename = "EXP_ROSA_INSTRUMENT_CORRECTED_FOR2SMR_pval.txt", sep= "\t")
EXP_SWERDLOW_INSTRUMENT_CORRECTED <- read_exposure_data(filename = "EXP_SWERDLOW_INSTRUMENT_CORRECTED_FOR2SMR_pval.txt", sep= "\t")
EXP_SARWAR_INSTRUMENT_CORRECTED <- read_exposure_data(filename = "EXP_SARWAR_INSTRUMENT_CORRECTED_FOR2SMR_pval.txt", sep= "\t")
OUTCOME_LAM_GWAS <- read_outcome_data(filename = "OUTCOME_LAM_GWAS_FOR2SMR_pval.txt", sep= "\t")
OUTCOME_HATOUM_GWAS <- fread("EFForChelsie.txt.gz")

################################
#  PART 2: HARMONIZE DATASETS  #
################################
# Extract instrument SNPs from outcome GWAS (if not available, uses LD proxies instead).
# action 1 = assumes all alleles coded on forward strand, doesn't flip
# action 2 = tries infer positive strand alleles using EAF for palindromes (default, conservative)
# action 3 = correct strand for non-palindromic snps, and drops all palindromic snps (more conservative)

# ligthart genome-wide (78 snps) - 74 available + 3 proxies = 77
ligthart_gw_lam<- harmonise_data( 
  exposure_dat = EXP_LIGTHART_GWAS_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

#ligthart cis (6 snps)
Ligthart_cis_lam <- harmonise_data( 
  exposure_dat = EXP_LIGTHART_CIS,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# han cis (20 snps) - 11 available + 2 proxies = 13
han_cis_lam<- harmonise_data( 
  exposure_dat = EXP_HAN_CIS,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# han genome-wide (552 snps) - 444 available + 50 proxies = 494
han_gw_lam<- harmonise_data( 
  exposure_dat = EXP_HAN_GWAS_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# ahluwalia cis (2 snps)
ahluwalia_cis_lam<- harmonise_data( 
  exposure_dat = EXP_AHLUWALIA_CIS,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# ahluwalia genome-wide (3 snps)
ahluwalia_gw_lam<- harmonise_data( 
  exposure_dat = EXP_AHLUWALIA_GWAS_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# borges (87 snps) - 73 available + 9 proxies = 82
borges_lam<- harmonise_data( 
  exposure_dat = EXP_BORGES_GWAS_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# kettunen (10 snps) - 9 available + 1 proxy = 10
kettunen_lam<- harmonise_data( 
  exposure_dat = EXP_KETTUNEN_GWAS_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# rosa (34 snps) - 27 available + 0 proxies = 27
# removed 5 snps for being palindromic with intermediate allele frequencies
rosa_lam<- harmonise_data( 
  exposure_dat = EXP_ROSA_INSTRUMENT_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# swerdlow (3 snps) 
swerdlow_lam<- harmonise_data( 
  exposure_dat = EXP_SWERDLOW_INSTRUMENT_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

# sarwar (1 snp) 
sarwar_lam<- harmonise_data( 
  exposure_dat = EXP_SARWAR_INSTRUMENT_CORRECTED,
  outcome_dat = OUTCOME_LAM_GWAS,
  action = 2
)

########################################
## PART 3: Check for palindromic snps ##
#######################################
# ligthart - gw = 4
palindromic_at<-subset(ligthart_gw_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(ligthart_gw_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(ligthart_gw_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(ligthart_gw_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# han - gw = 41
palindromic_at<-subset(han_gw_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(han_gw_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(han_gw_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(han_gw_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# ligthart - cis = 0
palindromic_at<-subset(Ligthart_cis_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(Ligthart_cis_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(Ligthart_cis_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(Ligthart_cis_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# han - cis = 2
palindromic_at<-subset(han_cis_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(han_cis_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(han_cis_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(han_cis_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# ahluwalia - gw = 1
palindromic_at<-subset(ahluwalia_gw_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(ahluwalia_gw_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(ahluwalia_gw_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(ahluwalia_gw_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# ahluwalia - cis = 1
palindromic_at<-subset(ahluwalia_cis_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(ahluwalia_cis_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(ahluwalia_cis_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(ahluwalia_cis_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# borges = 12
palindromic_at<-subset(borges_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(borges_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(borges_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(borges_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# kettunen = 2
palindromic_at<-subset(kettunen_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(kettunen_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(kettunen_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(kettunen_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# rosa = 5
palindromic_at<-subset(rosa_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(rosa_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(rosa_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(rosa_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# swerdlow = 0
palindromic_at<-subset(swerdlow_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(swerdlow_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(swerdlow_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(swerdlow_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

# sarwar = 0
palindromic_at<-subset(sarwar_lam,effect_allele.exposure %in% "A"&other_allele.exposure %in% "T")
palindromic_ta<-subset(sarwar_lam,effect_allele.exposure %in% "T"&other_allele.exposure %in% "A")
palindromic_gc<-subset(sarwar_lam,effect_allele.exposure %in% "G"&other_allele.exposure %in% "C")
palindromic_cg<-subset(sarwar_lam,effect_allele.exposure %in% "C"&other_allele.exposure %in% "G")

#################################################################################
#  PART 4: ESTIMATE CAUSAL EFFECT OF INFLAMMATION ON GENERAL COGNITIVE ABILITY  #
#################################################################################

#Set seed - SE obtained using random bootstraps, seed will ensure results are the same.
#IVW regression, mr egger, weighted median and weighted mode
set.seed(1234)
mr_results_ligthart_gw_lam <- mr(ligthart_gw_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_ligthart_cis_lam <- mr(Ligthart_cis_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_han_gw_lam <- mr(han_gw_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_han_cis_lam <- mr(han_cis_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_ahluwalia_gw_lam <- mr(ahluwalia_gw_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_ahluwalia_cis_lam <- mr(ahluwalia_cis_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_borges_lam <- mr(borges_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_kettunen_lam <- mr(kettunen_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_rosa_lam <- mr(rosa_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_swerdlow_lam <- mr(swerdlow_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_sarwar_lam <- mr(sarwar_lam, method_list=c("mr_wald_ratio"))

#dataframe of results
#results_ligthart_gw_lam <- mr_results_ligthart_gw_lam[, c("rsid","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype")]
#results_ligthart_gw_lam <- rename(results_ligthart_gw_lam,SNP = rsid,beta = beta.corrected,se = se,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

# Get 95% CI
mr_results_ligthart_gw_lam$lower <- mr_results_ligthart_gw_lam$b-1.96*mr_results_ligthart_gw_lam$se
mr_results_ligthart_gw_lam$upper <- mr_results_ligthart_gw_lam$b+1.96*mr_results_ligthart_gw_lam$se
mr_results_ligthart_cis_lam$lower <- mr_results_ligthart_cis_lam$b-1.96*mr_results_ligthart_cis_lam$se
mr_results_ligthart_cis_lam$upper <- mr_results_ligthart_cis_lam$b+1.96*mr_results_ligthart_cis_lam$se
mr_results_han_gw_lam$lower <- mr_results_han_gw_lam$b-1.96*mr_results_han_gw_lam$se
mr_results_han_gw_lam$upper <- mr_results_han_gw_lam$b+1.96*mr_results_han_gw_lam$se
mr_results_han_cis_lam$lower <- mr_results_han_cis_lam$b-1.96*mr_results_han_cis_lam$se
mr_results_han_cis_lam$upper <- mr_results_han_cis_lam$b+1.96*mr_results_han_cis_lam$se
mr_results_ahluwalia_gw_lam$lower <- mr_results_ahluwalia_gw_lam$b-1.96*mr_results_ahluwalia_gw_lam$se
mr_results_ahluwalia_gw_lam$upper <- mr_results_ahluwalia_gw_lam$b+1.96*mr_results_ahluwalia_gw_lam$se
mr_results_ahluwalia_cis_lam$lower <- mr_results_ahluwalia_cis_lam$b-1.96*mr_results_ahluwalia_cis_lam$se
mr_results_ahluwalia_cis_lam$upper <- mr_results_ahluwalia_cis_lam$b+1.96*mr_results_ahluwalia_cis_lam$se
mr_results_borges_lam$lower <- mr_results_borges_lam$b-1.96*mr_results_borges_lam$se
mr_results_borges_lam$upper <- mr_results_borges_lam$b+1.96*mr_results_borges_lam$se
mr_results_kettunen_lam$lower <- mr_results_kettunen_lam$b-1.96*mr_results_kettunen_lam$se
mr_results_kettunen_lam$upper <- mr_results_kettunen_lam$b+1.96*mr_results_kettunen_lam$se
mr_results_rosa_lam$lower <- mr_results_rosa_lam$b-1.96*mr_results_rosa_lam$se
mr_results_rosa_lam$upper <- mr_results_rosa_lam$b+1.96*mr_results_rosa_lam$se
mr_results_swerdlow_lam$lower <- mr_results_swerdlow_lam$b-1.96*mr_results_swerdlow_lam$se
mr_results_swerdlow_lam$upper <- mr_results_swerdlow_lam$b+1.96*mr_results_swerdlow_lam$se
mr_results_sarwar_lam$lower <- mr_results_sarwar_lam$b-1.96*mr_results_sarwar_lam$se
mr_results_sarwar_lam$upper <- mr_results_sarwar_lam$b+1.96*mr_results_sarwar_lam$se

# # save results
# results_ligthart_gw_lam <-cbind.data.frame(mr_results_ligthart_gw_lam$outcome,mr_results_ligthart_gw_lam$nsnp,mr_results_ligthart_gw_lam$method,mr_results_ligthart_gw_lam$b,mr_results_ligthart_gw_lam$se,mr_results_ligthart_gw_lam$pval, mr_results_ligthart_gw_lam$lower, mr_results_ligthart_gw_lam$upper)
# write.csv(results_ligthart_gw_lam, "./ligthart_gw_lam_results.csv")
# 
# results_ligthart_cis_lam <-cbind.data.frame(mr_results_ligthart_cis_lam$outcome,mr_results_ligthart_cis_lam$nsnp,mr_results_ligthart_cis_lam$method,mr_results_ligthart_cis_lam$b,mr_results_ligthart_cis_lam$se,mr_results_ligthart_cis_lam$pval, mr_results_ligthart_cis_lam$lower, mr_results_ligthart_cis_lam$upper)
# write.csv(results_ligthart_cis_lam, "./ligthart_cis_lam_results.csv")
# 
# results_han_gw_lam <-cbind.data.frame(mr_results_han_gw_lam$outcome,mr_results_han_gw_lam$nsnp,mr_results_han_gw_lam$method,mr_results_han_gw_lam$b,mr_results_han_gw_lam$se,mr_results_han_gw_lam$pval, mr_results_han_gw_lam$lower, mr_results_han_gw_lam$upper)
# write.csv(results_han_gw_lam, "./han_gw_lam_results.csv")
# 
# results_han_cis_lam <-cbind.data.frame(mr_results_han_cis_lam$outcome,mr_results_han_cis_lam$nsnp,mr_results_han_cis_lam$method,mr_results_han_cis_lam$b,mr_results_han_cis_lam$se,mr_results_han_cis_lam$pval, mr_results_han_cis_lam$lower, mr_results_han_cis_lam$upper)
# write.csv(results_han_cis_lam, "./han_cis_lam_results.csv")
# 
# results_ahluwalia_gw_lam <-cbind.data.frame(mr_results_ahluwalia_gw_lam$outcome,mr_results_ahluwalia_gw_lam$nsnp,mr_results_ahluwalia_gw_lam$method,mr_results_ahluwalia_gw_lam$b,mr_results_ahluwalia_gw_lam$se,mr_results_ahluwalia_gw_lam$pval, mr_results_ahluwalia_gw_lam$lower, mr_results_ahluwalia_gw_lam$upper)
# write.csv(results_ahluwalia_gw_lam, "./ahluwalia_gw_lam_results.csv")
# 
# results_ahluwalia_cis_lam <-cbind.data.frame(mr_results_ahluwalia_cis_lam$outcome,mr_results_ahluwalia_cis_lam$nsnp,mr_results_ahluwalia_cis_lam$method,mr_results_ahluwalia_cis_lam$b,mr_results_ahluwalia_cis_lam$se,mr_results_ahluwalia_cis_lam$pval, mr_results_ahluwalia_cis_lam$lower, mr_results_ahluwalia_cis_lam$upper)
# write.csv(results_ahluwalia_cis_lam, "./ahluwalia_cis_lam_results.csv")
# 
# results_borges_lam <-cbind.data.frame(mr_results_borges_lam$outcome,mr_results_borges_lam$nsnp,mr_results_borges_lam$method,mr_results_borges_lam$b,mr_results_borges_lam$se,mr_results_borges_lam$pval, mr_results_borges_lam$lower, mr_results_borges_lam$upper)
# write.csv(results_borges_lam, "./borges_lam_results.csv")
# 
# results_kettunen_lam <-cbind.data.frame(mr_results_kettunen_lam$outcome,mr_results_kettunen_lam$nsnp,mr_results_kettunen_lam$method,mr_results_kettunen_lam$b,mr_results_kettunen_lam$se,mr_results_kettunen_lam$pval, mr_results_kettunen_lam$lower, mr_results_kettunen_lam$upper)
# write.csv(results_kettunen_lam, "./kettunen_lam_results.csv")
# 
# results_swerdlow_lam <-cbind.data.frame(mr_results_swerdlow_lam$outcome,mr_results_swerdlow_lam$nsnp,mr_results_swerdlow_lam$method,mr_results_swerdlow_lam$b,mr_results_swerdlow_lam$se,mr_results_swerdlow_lam$pval, mr_results_swerdlow_lam$lower, mr_results_swerdlow_lam$upper)
# write.csv(results_swerdlow_lam, "./swerdlow_lam_results.csv")
# 
# results_sarwar_lam <-cbind.data.frame(mr_results_sarwar_lam$outcome,mr_results_sarwar_lam$nsnp,mr_results_sarwar_lam$method,mr_results_sarwar_lam$b,mr_results_sarwar_lam$se,mr_results_sarwar_lam$pval, mr_results_sarwar_lam$lower, mr_results_sarwar_lam$upper)
# write.csv(results_sarwar_lam, "./sarwar_lam_results.csv")
#
# results_rosa_lam <-cbind.data.frame(mr_results_rosa_lam$outcome,mr_results_rosa_lam$nsnp,mr_results_rosa_lam$method,mr_results_rosa_lam$b,mr_results_rosa_lam$se,mr_results_rosa_lam$pval, mr_results_rosa_lam$lower, mr_results_rosa_lam$upper)
# write.csv(results_rosa_lam, "./rosa_lam_results.csv")

###############################
#  PART 5: SENSITIVITY TESTS  #
###############################

##################################################
## Heterogeneity between SNPs using Cochran's Q ##
##################################################
#Calculated as weighted sum of squared differences between snps and pooled across snps
#Low power if number snps small
mr_het_ligthart_cis_lam <- mr_heterogeneity(Ligthart_cis_lam)
mr_het_ligthart_gw_lam <- mr_heterogeneity(ligthart_gw_lam)
mr_het_han_gw_lam <- mr_heterogeneity(han_gw_lam)
mr_het_han_cis_lam <- mr_heterogeneity(han_cis_lam)
mr_het_ahluwalia_gw_lam <- mr_heterogeneity(ahluwalia_gw_lam) 
mr_het_ahluwalia_cis_lam <- mr_heterogeneity(ahluwalia_cis_lam) 
mr_het_borges_lam <- mr_heterogeneity(borges_lam)
mr_het_kettunen_lam <- mr_heterogeneity(kettunen_lam)
mr_hett_rosa_lam <- mr_heterogeneity(rosa_lam)
mr_het_swerdlow_lam <- mr_heterogeneity(swerdlow_lam)

#################################
## Check horizontal pleiotropy ##
#################################
#Intercept term of MR-Egger indicates possible horizontal pleiotropy for ligthart GW and han GW
mr_pleiotropy_ligthart_cis_lam <- mr_pleiotropy_test(Ligthart_cis_lam)
mr_pleiotropy_ligthart_gw_lam <- mr_pleiotropy_test(ligthart_gw_lam)
mr_pleiotropy_han_gw_lam <- mr_pleiotropy_test(han_gw_lam)
mr_pleiotropy_han_cis_lam <- mr_pleiotropy_test(han_cis_lam)
mr_pleiotropy_ahluwalia_gw_lam <- mr_pleiotropy_test(ahluwalia_gw_lam)
mr_pleiotropy_ahluwalia_cis_lam <- mr_pleiotropy_test(ahluwalia_cis_lam)
mr_pleiotropy_borges_lam <- mr_pleiotropy_test(borges_lam)
mr_pleiotropy_kettunen_lam <- mr_pleiotropy_test(kettunen_lam)
mr_pleiotropy_rosa_lam <- mr_pleiotropy_test(rosa_lam)
mr_pleiotropy_swerdlow_lam <- mr_pleiotropy_test(swerdlow_lam)

#######################
## Check single snps ##
#######################
res_single_ligthart_cis <- mr_singlesnp(Ligthart_cis_lam)
res_single_ligtahart_gw <- mr_singlesnp(ligthart_gw_lam)
res_single_han_gw <- mr_singlesnp(han_gw_lam)
res_single_han_cis <- mr_singlesnp(han_cis_lam)
res_single_ahluwalia_cis <- mr_singlesnp(ahluwalia_cis_lam)
res_single_ahluwalia_gw <- mr_singlesnp(ahluwalia_gw_lam)
res_single_borges <- mr_singlesnp(borges_lam)
res_single_kettunen <- mr_singlesnp(kettunen_lam)
res_single_rosa <- mr_singlesnp(rosa_lam)
res_single_swerdlow <- mr_singlesnp(swerdlow_lam)

##############################
#  PART 6: VISUALIZE EFFECT  #
##############################
##################################################
## 1. scatter plots comparing different methods ##
##################################################
#mr_report(dat)
png("./ligthart_gw_lam_scatter.png")
mr_scatter_plot(mr_results_ligthart_gw_lam, ligthart_gw_lam)
dev.off()

png("./ligthart_cis_lam_scatter.png")
mr_scatter_plot(mr_results_ligthart_cis_lam, Ligthart_cis_lam)
dev.off()

png("./han_cis_lam_scatter.png")
mr_scatter_plot(mr_results_han_cis_lam, han_cis_lam)
dev.off()

png("./han_gw_lam_scatter.png")
mr_scatter_plot(mr_results_han_gw_lam, han_gw_lam)
dev.off()

png("./ahluwalia_cis_lam_scatter.png")
mr_scatter_plot(mr_results_ahluwalia_cis_lam, ahluwalia_cis_lam)
dev.off()

png("./ahluwalia_gw_lam_scatter.png")
mr_scatter_plot(mr_results_ahluwalia_gw_lam, ahluwalia_gw_lam)
dev.off()

png("./borges_lam_scatter.png")
mr_scatter_plot(mr_results_borges_lam, borges_lam)
dev.off()

png("./kettunen_lam_scatter.png")
mr_scatter_plot(mr_results_kettunen_lam, kettunen_lam)
dev.off()

png("./rosa_lam_scatter.png")
mr_scatter_plot(mr_results_rosa_lam, rosa_lam)
dev.off()

png("./swerdlow_lam_scatter.png")
mr_scatter_plot(mr_results_swerdlow_lam, swerdlow_lam)
dev.off()

#########################################
## 2. Forest plots of each snp effects ##
#########################################
#forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./ligthart_gw_lam_forest.png")
res_single_ligthart_gw_lam <- mr_singlesnp(ligthart_gw_lam)
mr_forest_plot(res_single_ligthart_gw_lam)#was res_single
dev.off()

png("./ligthart_cis_lam_forest.png")
res_single_ligthart_cis_lam <- mr_singlesnp(Ligthart_cis_lam)
mr_forest_plot(res_single_ligthart_cis_lam)#was res_single
dev.off()

png("./han_gw_lam_forest.png")
res_single_han_gw_lam <- mr_singlesnp(han_gw_lam)
mr_forest_plot(res_single_han_gw_lam)#was res_single
dev.off()

png("./han_cis_lam_forest.png")
res_single_han_cis_lam <- mr_singlesnp(han_cis_lam)
mr_forest_plot(res_single_han_cis_lam)#was res_single
dev.off()

png("./ahluwalia_gw_lam_forest.png")
res_single_ahluwalia_gw_lam <- mr_singlesnp(ahluwalia_gw_lam)
mr_forest_plot(res_single_ahluwalia_gw_lam)#was res_single
dev.off()

png("./ahluwalia_cis_lam_forest.png")
res_single_ahluwalia_cis_lam <- mr_singlesnp(ahluwalia_cis_lam)
mr_forest_plot(res_single_ahluwalia_cis_lam)#was res_single
dev.off()

png("./borges_lam_forest.png")
res_single_borges_lam <- mr_singlesnp(borges_lam)
mr_forest_plot(res_single_borges_lam)#was res_single
dev.off()

png("./kettunen_lam_forest.png")
res_single_kettunen_lam <- mr_singlesnp(kettunen_lam)
mr_forest_plot(res_single_kettunen_lam)#was res_single
dev.off()

png("./rosa_lam_forest.png")
res_single_rosa_lam <- mr_singlesnp(rosa_lam)
mr_forest_plot(res_single_rosa_lam)#was res_single
dev.off()

png("./swerdlow_lam_forest.png")
res_single_swerdlow_lam <- mr_singlesnp(swerdlow_lam)
mr_forest_plot(res_single_swerdlow_lam)#was res_single
dev.off()

########################################
## 3. Funnel plots to check asymmetry ##
########################################
png("./ligthart_gw_lam_funnel.png")
mr_funnel_plot(res_single_ligthart_gw_lam)
dev.off()

png("./ligthart_cis_lam_funnel.png")
mr_funnel_plot(res_single_ligthart_cis_lam)
dev.off()

png("./han_gw_lam_funnel.png")
mr_funnel_plot(res_single_han_gw_lam)
dev.off()

png("./han_cis_lam_funnel.png")
mr_funnel_plot(res_single_han_cis_lam)
dev.off()

png("./ahluwalia_gw_lam_funnel.png")
mr_funnel_plot(res_single_ahluwalia_gw_lam)
dev.off()

png("./ahluwalia_cis_lam_funnel.png")
mr_funnel_plot(res_single_ahluwalia_cis_lam)
dev.off()

png("./borges_lam_funnel.png")
mr_funnel_plot(res_single_borges_lam)
dev.off()

png("./kettunen_lam_funnel.png")
mr_funnel_plot(res_single_kettunen_lam)
dev.off()

png("./rosa_lam_funnel.png")
mr_funnel_plot(res_single_rosa_lam)
dev.off()

png("./swerdlow_lam_funnel.png")
mr_funnel_plot(res_single_swerdlow_lam)
dev.off()

###############################
## 4. Leave-one-out analysis ##
###############################
#test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo_ligthart_gw_lam <- mr_leaveoneout(ligthart_gw_lam)
png("./ligthart_gw_lam_loo.png")
mr_leaveoneout_plot(res_loo_ligthart_gw_lam)
dev.off()

res_loo_ligthart_cis_lam <- mr_leaveoneout(ligthart_cis_lam)
png("./ligthart_cis_lam_loo.png")
mr_leaveoneout_plot(res_loo_ligthart_cis_lam)
dev.off()

res_loo_han_gw_lam <- mr_leaveoneout(han_gw_lam)
png("./han_gw_lam_loo.png")
mr_leaveoneout_plot(res_loo_han_gw_lam)
dev.off()

res_loo_han_cis_lam <- mr_leaveoneout(han_cis_lam)
png("./han_cis_lam_loo.png")
mr_leaveoneout_plot(res_loo_han_cis_lam)
dev.off()

res_loo_ahluwalia_gw_lam <- mr_leaveoneout(ahluwalia_gw_lam)
png("./ahluwalia_gw_lam_loo.png")
mr_leaveoneout_plot(res_loo_ahluwalia_gw_lam)
dev.off()

res_loo_ahluwalia_cis_lam <- mr_leaveoneout(ahluwalia_cis_lam)
png("./ahluwalia_cis_lam_loo.png")
mr_leaveoneout_plot(res_loo_ahluwalia_cis_lam)
dev.off()

res_loo_borges_lam <- mr_leaveoneout(borges_lam)
png("./borges_lam_loo.png")
mr_leaveoneout_plot(res_loo_borges_lam)
dev.off()

res_loo_kettunen_lam <- mr_leaveoneout(kettunen_lam)
png("./kettunen_lam_loo.png")
mr_leaveoneout_plot(res_loo_kettunen_lam)
dev.off()

res_loo_rosa_lam <- mr_leaveoneout(rosa_lam)
png("./rosa_lam_loo.png")
mr_leaveoneout_plot(res_loo_rosa_lam)
dev.off()

res_loo_swerdlow_lam <- mr_leaveoneout(swerdlow_lam)
png("./swerdlow_lam_loo.png")
mr_leaveoneout_plot(res_loo_swerdlow_lam)
dev.off()

#############################################
#  PART 6: ADDITIONAL SENSITIVITY ANALYSIS  #
#############################################

###################
#### MR-PRESSO ####
###################

#Sd stands for estimated standard error of causal estimate (https://github.com/rondolab/MR-PRESSO/issues/15)
#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se_outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = ligthart_gw_lam, NbDistribution = 1000,  SignifThreshold = 0.05)
#NbDistribution = number of bootstrap replications, default is 1000. Outlier sig threshold, default is 0.05
presso_ligthart_gw_lam_2 <- run_mr_presso(ligthart_gw_lam, NbDistribution = 10000, SignifThreshold = 0.05)
presso_ligthart_cis_lam_2 <- run_mr_presso(Ligthart_cis_lam, NbDistribution = 10000, SignifThreshold = 0.05)
presso_han_gw_lam_2 <- run_mr_presso(han_gw_lam, NbDistribution = 10000, SignifThreshold = 0.05)
presso_han_cis_lam_2 <- run_mr_presso(han_cis_lam, NbDistribution = 10000, SignifThreshold = 0.05)
presso_borges_lam_2 <- run_mr_presso(borges_lam, NbDistribution = 10000, SignifThreshold = 0.05)
presso_kettunen_lam_2 <- run_mr_presso(kettunen_lam, NbDistribution = 10000, SignifThreshold = 0.05)
presso_rosa_lam_2 <- run_mr_presso(rosa_lam, NbDistribution = 10000, SignifThreshold = 0.05)

#saves data
#list.save(presso_ligthart_gw_lam_2, 'mrpresso_ligthart_gw_lam.rdata')
#list.save(presso_ligthart_cis_lam_2, 'mrpresso_ligthart_cis_lam.rdata')
#list.save(presso_han_gw_lam_2, 'mrpresso_han_gw_lam.rdata')
#list.save(presso_han_cis_lam_2, 'mrpresso_han_cis_lam.rdata')
#list.save(presso_borges_lam_2, 'mrpresso_borges_lam.rdata')
#list.save(presso_kettunen_lam_2, 'mrpresso_kettunen_lam.rdata')
#list.save(presso_rosa_lam_2, 'mrpresso_rosa_lam.rdata')

# 95% CI to MR-PRESSO causal estimates
# must load in above results separately
rm(list = ls())
presso_estimate <- paste(x[[1]][["Main MR results"]][["Causal Estimate"]])
presso_se <- paste(x[[1]][["Main MR results"]][["Sd"]])
presso_results <- data.frame(presso_estimate,presso_se) 
presso_results[] <- lapply(presso_results, function(x) as.numeric(as.character(x)))
presso_results$lower <- presso_results$presso_estimate-1.96*presso_results$presso_se
presso_results$upper <- presso_results$presso_estimate+1.96*presso_results$presso_se

###########################
#### Steiger filtering ####
###########################
#Identifies if stronger bidirectional effects (e.g., A on B or B on A) to select valid IVs from very large GWAS that can be used in MR. 
#Assumes that a valid IV should explain more variance in the exposure than the outcome and removes those genetic variants that do not satisfy this criterion

#Step 1: run Steiger filter on each exposure
#Ligthart GW = 77 true
ligthart_gw_lam$samplesize.exposure <- 204402
ligthart_gw_lam$samplesize.outcome <- 373617
steiger_ligthart_gw_lam <- steiger_filtering(ligthart_gw_lam)
table(steiger_ligthart_gw_lam$steiger_dir)

#Ligthart cis = 6 true
Ligthart_cis_lam$samplesize.exposure <- 204402
Ligthart_cis_lam$samplesize.outcome <- 373617
steiger_ligthart_cis_lam <- steiger_filtering(Ligthart_cis_lam)
table(steiger_ligthart_cis_lam$steiger_dir)

#Han GW = 490 true, 4 false. (if run with SD - beta and se = 6 false)
#used p-val and sample size to calculate R2 (as exposure and outcome are not standardised...)
han_gw_lam$samplesize.exposure <- 418642
han_gw_lam$samplesize.outcome <- 373617
steiger_han_gw_lam <- steiger_filtering(han_gw_lam)
table(steiger_han_gw_lam$steiger_dir)

#Han cis = 13 true
han_cis_lam$samplesize.exposure <- 418642
han_cis_lam$samplesize.outcome <- 373617
steiger_han_cis_lam <- steiger_filtering(han_cis_lam)
table(steiger_han_cis_lam$steiger_dir)

#Ahluwalia cis = 2 true
ahluwalia_cis_lam$samplesize.exposure <- 52654
ahluwalia_cis_lam$samplesize.outcome <- 373617
steiger_ahluwalia_cis_lam <- steiger_filtering(ahluwalia_cis_lam)
table(steiger_ahluwalia_cis_lam$steiger_dir)

#Ahluwalia GW = 3 TRUE
ahluwalia_gw_lam$samplesize.exposure <- 52654
ahluwalia_gw_lam$samplesize.outcome <- 373617
steiger_ahluwalia_gw_lam <- steiger_filtering(ahluwalia_gw_lam)
table(steiger_ahluwalia_gw_lam$steiger_dir)

#Rosa = 27 true
rosa_lam$samplesize.exposure <- 3301
rosa_lam$samplesize.outcome <- 373617
steiger_rosa_lam <- steiger_filtering(rosa_lam)
table(steiger_rosa_lam$steiger_dir)

#Borges = 82 true
borges_lam$samplesize.exposure <- 115078
borges_lam$samplesize.outcome <- 373617
steiger_borges_lam <- steiger_filtering(borges_lam)
table(steiger_borges_lam$steiger_dir)

#Kettunen = 10 true
kettunen_lam$samplesize.exposure <- 19270
kettunen_lam$samplesize.outcome <- 373617
steiger_kettunen_lam <- steiger_filtering(kettunen_lam)
table(steiger_kettunen_lam$steiger_dir)

#Swerdlow = 3 true
swerdlow_lam$samplesize.exposure <- 4479
swerdlow_lam$samplesize.outcome <- 373617
steiger_swerdlow_lam <- steiger_filtering(swerdlow_lam)
table(steiger_swerdlow_lam$steiger_dir)

#Sarwar = 1 true
sarwar_lam$samplesize.exposure <- 27185
sarwar_lam$samplesize.outcome <- 373617
steiger_sarwar_lam <- steiger_filtering(sarwar_lam)
table(steiger_sarwar_lam$steiger_dir)

#Step 2: run MR with "false" SNPs removed
steiger_han_gw_lam <- subset(steiger_han_gw_lam, steiger_dir=="TRUE", select=SNP:mr_keep)
set.seed(1234)
steiger_mr_results_han_gw_lam <- mr(steiger_han_gw_lam, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
#write.csv(steiger_mr_results_han_gw_lam, "./steiger_han_gw_lam_results.csv") #save data
