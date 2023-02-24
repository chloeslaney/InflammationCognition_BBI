# Project: inflammation and cognition 
# Script for two-sample MR - adapted from MRC IEU MR short-course code
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
EXP_LAM_FOR_LIGT <- read_exposure_data(filename = "EXP_LAM_FOR_LIGTHART_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_LAM_FOR_HAN <- read_exposure_data(filename = "EXP_LAM_FOR_HAN_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_LAM_FOR_AHLUWALIA <- read_exposure_data(filename = "EXP_LAM_FOR_AHLUWALIA_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_LAM_FOR_BORGES <- read_exposure_data(filename = "EXP_LAM_FOR_BORGES_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
EXP_LAM_FOR_KETTUNEN <- read_exposure_data(filename = "EXP_LAM_FOR_KETTUNEN_FOR2SMR_PROXIESADDED_pval.txt", sep= "\t")
OUTCOME_LIGTHART_GWAS <- read_outcome_data(filename = "OUTCOME_LIGTHART_GWAS_FOR2SMR.txt", sep= "\t")
OUTCOME_HAN_GWAS <- read_outcome_data(filename = "OUTCOME_HAN_GWAS_FOR2SMR.txt", sep= "\t")
OUTCOME_AHLUWALIA_GWAS <- read_outcome_data(filename = "OUTCOME_AHLUWALIA_GWAS_FOR2SMR.txt", sep= "\t")
OUTCOME_BORGES_GWAS <- read_outcome_data(filename = "OUTCOME_BORGES_GWAS_FOR2SMR.txt", sep= "\t")
OUTCOME_KETTUNEN_GWAS <- read_outcome_data(filename = "OUTCOME_KETTUNEN_GWAS_FOR2SMR.txt", sep= "\t")

################################
#  PART 2: HARMONIZE DATASETS  #
################################
# 250 snps meet criteria for Lam GWAS.
# Extract instrument SNPs from outcome GWAS (if not available, uses LD proxies instead).
# action 1 = assumes all alleles coded on forward strand, doesn't flip
# action 2 = tries infer positive strand alleles using EAF for palindromes (default, conservative)
# action 3 = correct strand for non-palindromic snps, and drops all palindromic snps (more conservative)
# whilst the exposures not have "corrected" in names, they are corrected (i.e., beta/eaf flipped)

# lam - ligthart (219) - 129 + 90 proxies = 219
lam_ligthart <- harmonise_data( 
  exposure_dat = EXP_LAM_FOR_LIGT,
  outcome_dat = OUTCOME_LIGTHART_GWAS,
  action = 2
)

# lam - han (249) - no proxy available for one missing snp
lam_han <- harmonise_data( 
  exposure_dat = EXP_LAM_FOR_HAN,
  outcome_dat = OUTCOME_HAN_GWAS,
  action = 2
)

# lam - ahluwalia (222) - 137 + 85 proxies = 222
lam_ahluwalia <- harmonise_data( 
  exposure_dat = EXP_LAM_FOR_AHLUWALIA,
  outcome_dat = OUTCOME_AHLUWALIA_GWAS,
  action = 2
)

# lam - borges (250)
lam_borges <- harmonise_data( 
  exposure_dat = EXP_LAM_FOR_BORGES,
  outcome_dat = OUTCOME_BORGES_GWAS,
  action = 2
)

# lam - kettunen (250) - 248 + 2 proxies = 250 
lam_kettunen <- harmonise_data( 
  exposure_dat = EXP_LAM_FOR_KETTUNEN,
  outcome_dat = OUTCOME_KETTUNEN_GWAS,
  action = 2
)

#################################################################################
#  PART 3: ESTIMATE CAUSAL EFFECT OF INFLAMMATION ON GENERAL COGNITIVE ABILITY  #
#################################################################################
#Set seed - SE obtained using random bootstraps, seed will ensure results are the same.
#IVW regression, mr egger, weighted median and weighted mode
set.seed(1234)
mr_results_lam_ligthart <- mr(lam_ligthart, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_lam_han <- mr(lam_han, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_lam_ahluwalia <- mr(lam_ahluwalia, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_lam_borges <- mr(lam_borges, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
set.seed(1234)
mr_results_lam_kettunen <- mr(lam_kettunen, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

#Get 95% CI
mr_results_lam_ligthart$lower <- mr_results_lam_ligthart$b-1.96*mr_results_lam_ligthart$se
mr_results_lam_ligthart$upper <- mr_results_lam_ligthart$b+1.96*mr_results_lam_ligthart$se
mr_results_lam_han$lower <- mr_results_lam_han$b-1.96*mr_results_lam_han$se
mr_results_lam_han$upper <- mr_results_lam_han$b+1.96*mr_results_lam_han$se
mr_results_lam_ahluwalia$lower <- mr_results_lam_ahluwalia$b-1.96*mr_results_lam_ahluwalia$se
mr_results_lam_ahluwalia$upper <- mr_results_lam_ahluwalia$b+1.96*mr_results_lam_ahluwalia$se
mr_results_lam_borges$lower <- mr_results_lam_borges$b-1.96*mr_results_lam_borges$se
mr_results_lam_borges$upper <- mr_results_lam_borges$b+1.96*mr_results_lam_borges$se
mr_results_lam_kettunen$lower <- mr_results_lam_kettunen$b-1.96*mr_results_lam_kettunen$se
mr_results_lam_kettunen$upper <- mr_results_lam_kettunen$b+1.96*mr_results_lam_kettunen$se

#save results
results_lam_ligthart <-cbind.data.frame(mr_results_lam_ligthart$outcome,mr_results_lam_ligthart$nsnp,mr_results_lam_ligthart$method,mr_results_lam_ligthart$b,mr_results_lam_ligthart$se,mr_results_lam_ligthart$pval, mr_results_lam_ligthart$lower, mr_results_lam_ligthart$upper)
#write.csv(results_lam_ligthart, "./lam_ligthart_reversedir_results.csv")
results_lam_han <-cbind.data.frame(mr_results_lam_han$outcome,mr_results_lam_han$nsnp,mr_results_lam_han$method,mr_results_lam_han$b,mr_results_lam_han$se,mr_results_lam_han$pval, mr_results_lam_han$lower, mr_results_lam_han$upper)
#write.csv(results_lam_han, "./lam_han_reversedir_results.csv")
results_lam_ahluwalia <-cbind.data.frame(mr_results_lam_ahluwalia$outcome,mr_results_lam_ahluwalia$nsnp,mr_results_lam_ahluwalia$method,mr_results_lam_ahluwalia$b,mr_results_lam_ahluwalia$se,mr_results_lam_ahluwalia$pval, mr_results_lam_ahluwalia$lower, mr_results_lam_ahluwalia$upper)
#write.csv(results_lam_ahluwalia, "./lam_ahluwalia_reversedir_results.csv")
results_lam_borges <-cbind.data.frame(mr_results_lam_borges$outcome,mr_results_lam_borges$nsnp,mr_results_lam_borges$method,mr_results_lam_borges$b,mr_results_lam_borges$se,mr_results_lam_borges$pval, mr_results_lam_borges$lower, mr_results_lam_borges$upper)
#write.csv(results_lam_borges, "./lam_borges_reversedir_results.csv")
results_lam_kettunen <-cbind.data.frame(mr_results_lam_kettunen$outcome,mr_results_lam_kettunen$nsnp,mr_results_lam_kettunen$method,mr_results_lam_kettunen$b,mr_results_lam_kettunen$se,mr_results_lam_kettunen$pval, mr_results_lam_kettunen$lower, mr_results_lam_kettunen$upper)
#write.csv(results_lam_kettunen, "./lam_kettunen_reversedir_results.csv")

##################################################
## Heterogeneity between SNPs using Cochran's Q ##
##################################################
#Calculated as weighted sum of squared differences between snps and pooled across snps
#Low power if number snps small
mr_het_lam_ligthart <- mr_heterogeneity(lam_ligthart)
mr_het_lam_han <- mr_heterogeneity(lam_han)
mr_het_lam_ahluwalia <- mr_heterogeneity(lam_ahluwalia)
mr_het_lam_borges <- mr_heterogeneity(lam_borges)
mr_het_lam_kettunen <- mr_heterogeneity(lam_kettunen)

#################################
## Check horizontal pleiotropy ##
#################################
mr_pleiotropy_lam_ligthart <- mr_pleiotropy_test(lam_ligthart)
mr_pleiotropy_lam_han <- mr_pleiotropy_test(lam_han)
mr_pleiotropy_lam_ahluwalia <- mr_pleiotropy_test(lam_ahluwalia)
mr_pleiotropy_lam_borges <- mr_pleiotropy_test(lam_borges)
mr_pleiotropy_lam_kettunen <- mr_pleiotropy_test(lam_kettunen)

#######################
## Check single snps ##
#######################
res_single_lam_ligthart <- mr_singlesnp(lam_ligthart)
res_single_lam_han <- mr_singlesnp(lam_han)
res_single_lam_ahluwalia <- mr_singlesnp(lam_ahluwalia)
res_single_lam_borges <- mr_singlesnp(lam_borges)
res_single_lam_kettunen <- mr_singlesnp(lam_kettunen)

##############################
#  PART 4: VISUALIZE EFFECT  #
##############################
##################################################
## 1. scatter plots comparing different methods ##
##################################################
#mr_report(dat)
png("./lam_ligthart_reversedir_scatter.png")
mr_scatter_plot(mr_results_lam_ligthart, lam_ligthart)
dev.off()

png("./lam_han_reversedir_scatter.png")
mr_scatter_plot(mr_results_lam_han, lam_han)
dev.off()

png("./lam_ahluwalia_reversedir_scatter.png")
mr_scatter_plot(mr_results_lam_ahluwalia, lam_ahluwalia)
dev.off()

png("./lam_borges_reversedir_scatter.png")
mr_scatter_plot(mr_results_lam_borges, lam_borges)
dev.off()

png("./lam_kettunen_reversedir_scatter.png")
mr_scatter_plot(mr_results_lam_kettunen, lam_kettunen)
dev.off()

#########################################
## 2. Forest plots of each snp effects ##
#########################################
png("./lam_ligthart_reversedir_forest.png")
res_single_lam_ligthart <- mr_singlesnp(lam_ligthart)
mr_forest_plot(res_single_lam_ligthart)
dev.off()

png("./lam_han_reversedir_forest.png")
res_single_lam_han<- mr_singlesnp(lam_han)
mr_forest_plot(res_single_lam_han)
dev.off()

png("./lam_ahluwalia_reversedir_forest.png")
res_single_lam_ahluwalia <- mr_singlesnp(lam_ahluwalia)
mr_forest_plot(res_single_lam_ahluwalia)
dev.off()

png("./lam_borges_reversedir_forest.png")
res_single_lam_borges <- mr_singlesnp(lam_borges)
mr_forest_plot(res_single_lam_borges)
dev.off()

png("./lam_kettunen_reversedir_forest.png")
res_single_lam_kettunen <- mr_singlesnp(lam_kettunen)
mr_forest_plot(res_single_lam_kettunen)
dev.off()

########################################
## 3. Funnel plots to check asymmetry ##
########################################
png("./lam_ligthart_reversedir_funnel.png")
mr_funnel_plot(res_single_lam_ligthart)
dev.off()

png("./lam_han_reversedir_funnel.png")
mr_funnel_plot(res_single_lam_han)
dev.off()

png("./lam_ahluwalia_reversedir_funnel.png")
mr_funnel_plot(res_single_lam_ahluwalia)
dev.off()

png("./lam_borges_reversedir_funnel.png")
mr_funnel_plot(res_single_lam_borges)
dev.off()

png("./lam_kettunen_reversedir_funnel.png")
mr_funnel_plot(res_single_lam_kettunen)
dev.off()

###############################
## 4. Leave-one-out analysis ##
###############################
res_loo_lam_ligthart <- mr_leaveoneout(lam_ligthart)
png("./lam_ligthart_reversedir_loo.png")
mr_leaveoneout_plot(res_loo_lam_ligthart)
dev.off()

res_loo_lam_han <- mr_leaveoneout(lam_han)
png("./lam_han_reversedir_loo.png")
mr_leaveoneout_plot(res_loo_lam_han)
dev.off()

res_loo_lam_ahluwalia <- mr_leaveoneout(lam_ahluwalia)
png("./lam_ahluwalia_reversedir_loo.png")
mr_leaveoneout_plot(res_loo_lam_ahluwalia)
dev.off()

res_loo_lam_borges <- mr_leaveoneout(lam_borges)
png("./lam_borges_reversedir_loo.png")
mr_leaveoneout_plot(res_loo_lam_borges)
dev.off()

res_loo_lam_kettunen <- mr_leaveoneout(lam_kettunen)
png("./lam_kettunen_reversedir_loo.png")
mr_leaveoneout_plot(res_loo_lam_kettunen)
dev.off()

#############################################
#  PART 5: ADDITIONAL SENSITIVITY ANALYSIS  #
#############################################

###################
#### MR-PRESSO ####
###################
#Sd stands for estimated standard error of causal estimate (https://github.com/rondolab/MR-PRESSO/issues/15)
#mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se_outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = ligthart_gw_lam, NbDistribution = 1000,  SignifThreshold = 0.05)
#NbDistribution = number of bootstrap replications, default is 1000. Outlier sig threshold, default is 0.05
presso_lam_ligthart_reversedir <- run_mr_presso(lam_ligthart, NbDistribution = 10000, SignifThreshold = 0.05)
presso_lam_han_reversedir <- run_mr_presso(lam_han, NbDistribution = 10000, SignifThreshold = 0.05)
presso_lam_ahluwalia_reversedir <- run_mr_presso(lam_ahluwalia, NbDistribution = 10000, SignifThreshold = 0.05)
presso_lam_borges_reversedir <- run_mr_presso(lam_borges, NbDistribution = 10000, SignifThreshold = 0.05)
presso_lam_kettunen_reversedir <- run_mr_presso(lam_kettunen, NbDistribution = 10000, SignifThreshold = 0.05)

#saves data
list.save(presso_lam_ligthart_reversedir, 'mrpresso_lam_ligthart_reversedir.rdata')
list.save(presso_lam_han_reversedir, 'mrpresso_lam_han_reversedir.rdata')
list.save(presso_lam_ahluwalia_reversedir, 'mrpresso_lam_ahluwalia_reversedir.rdata')
list.save(presso_lam_borges_reversedir, 'mrpresso_lam_borges_reversedir.rdata')
list.save(presso_lam_kettunen_reversedir, 'mrpresso_lam_kettunen_reversedir.rdata')

# 95% CI to MR-PRESSO causal estimates
# must load in above results separately
rm(list = ls())
presso_estimate <- paste(x[[1]][["Main MR results"]][["Causal Estimate"]])
presso_se <- paste(x[[1]][["Main MR results"]][["Sd"]])
presso_pval <- paste(x[[1]][["Main MR results"]][["P-value"]])
presso_results <- data.frame(presso_estimate,presso_se,presso_pval) 
presso_results[] <- lapply(presso_results, function(x) as.numeric(as.character(x)))
presso_results$lower <- presso_results$presso_estimate-1.96*presso_results$presso_se
presso_results$upper <- presso_results$presso_estimate+1.96*presso_results$presso_se

###########################
#### Steiger filtering ####
###########################
#Identifies if stronger bidirectional effects (e.g., A on B or B on A) to select valid IVs from very large GWAS that can be used in MR. 
#Assumes that a valid IV should explain more variance in the exposure than the outcome and removes those genetic variants that do not satisfy this criterion
#Step 1: run Steiger filter on each exposure

#Lam - ligthart = 218 TRUE, 1 FALSE
lam_ligthart$samplesize.outcome <- 204402
lam_ligthart$samplesize.exposure <- 373617
steiger_lam_ligthart_reversedir <- steiger_filtering(lam_ligthart)
table(steiger_lam_ligthart_reversedir$steiger_dir)

#Lam - han = 246 TRUE, 3 FALSE
lam_han$samplesize.outcome <- 418642
lam_han$samplesize.exposure <- 373617
steiger_lam_han_reversedir <- steiger_filtering(lam_han)
table(steiger_lam_han_reversedir$steiger_dir)

#Lam - ahluwalia = 215 TRUE, 7 FALSE
lam_ahluwalia$samplesize.outcome <- 52654
lam_ahluwalia$samplesize.exposure <- 373617
steiger_lam_ahluwalia_reversedir <- steiger_filtering(lam_ahluwalia)
table(steiger_lam_ahluwalia_reversedir$steiger_dir)

#Lam - borges = 239 TRUE, 11 FALSE
lam_borges$samplesize.outcome <- 115078
lam_borges$samplesize.exposure <- 373617
steiger_lam_borges_reversedir <- steiger_filtering(lam_borges)
table(steiger_lam_borges_reversedir$steiger_dir)

#Lam - kettunen = 192 TRUE, 58 FALSE
#only 191 included in MR - unclear why
lam_kettunen$samplesize.outcome <- 19270
lam_kettunen$samplesize.exposure <- 373617
steiger_lam_kettunen_reversedir <- steiger_filtering(lam_kettunen)
table(steiger_lam_kettunen_reversedir$steiger_dir)

#Step 2: run MR with "false" SNPs removed
steiger_lam_ligthart_reversedir <- subset(steiger_lam_ligthart_reversedir, steiger_dir=="TRUE", select=SNP:mr_keep)
set.seed(1234)
steiger_mr_results_lam_ligthart <- mr(steiger_lam_ligthart_reversedir, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

steiger_lam_han_reversedir <- subset(steiger_lam_han_reversedir, steiger_dir=="TRUE", select=SNP:mr_keep)
set.seed(1234)
steiger_mr_results_lam_han <- mr(steiger_lam_han_reversedir, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

steiger_lam_ahluwalia_reversedir <- subset(steiger_lam_ahluwalia_reversedir, steiger_dir=="TRUE", select=SNP:mr_keep)
set.seed(1234)
steiger_mr_results_lam_ahluwalia <- mr(steiger_lam_ahluwalia_reversedir, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

steiger_lam_borges_reversedir <- subset(steiger_lam_borges_reversedir, steiger_dir=="TRUE", select=SNP:mr_keep)
set.seed(1234)
steiger_mr_results_lam_borges <- mr(steiger_lam_borges_reversedir, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

steiger_lam_kettunen_reversedir <- subset(steiger_lam_kettunen_reversedir, steiger_dir=="TRUE", select=SNP:mr_keep)
set.seed(1234)
steiger_mr_results_lam_kettunen <- mr(steiger_lam_kettunen_reversedir, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

#Get 95% CI
steiger_mr_results_lam_ligthart$lower <- steiger_mr_results_lam_ligthart$b-1.96*steiger_mr_results_lam_ligthart$se
steiger_mr_results_lam_ligthart$upper <- steiger_mr_results_lam_ligthart$b+1.96*steiger_mr_results_lam_ligthart$se
steiger_mr_results_lam_han$lower <- steiger_mr_results_lam_han$b-1.96*steiger_mr_results_lam_han$se
steiger_mr_results_lam_han$upper <- steiger_mr_results_lam_han$b+1.96*steiger_mr_results_lam_han$se
steiger_mr_results_lam_ahluwalia$lower <- steiger_mr_results_lam_ahluwalia$b-1.96*steiger_mr_results_lam_ahluwalia$se
steiger_mr_results_lam_ahluwalia$upper <- steiger_mr_results_lam_ahluwalia$b+1.96*steiger_mr_results_lam_ahluwalia$se
steiger_mr_results_lam_borges$lower <- steiger_mr_results_lam_borges$b-1.96*steiger_mr_results_lam_borges$se
steiger_mr_results_lam_borges$upper <- steiger_mr_results_lam_borges$b+1.96*steiger_mr_results_lam_borges$se
steiger_mr_results_lam_kettunen$lower <- steiger_mr_results_lam_kettunen$b-1.96*steiger_mr_results_lam_kettunen$se
steiger_mr_results_lam_kettunen$upper <- steiger_mr_results_lam_kettunen$b+1.96*steiger_mr_results_lam_kettunen$se

#saves data
#write.csv(steiger_mr_results_lam_ligthart, "./steiger_lam_ligthart_reversedir_results.csv")
#write.csv(steiger_mr_results_lam_han, "./steiger_lam_han_reversedir_results.csv")
#write.csv(steiger_mr_results_lam_ahluwalia, "./steiger_lam_ahluwalia_reversedir_results.csv")
#write.csv(steiger_mr_results_lam_borges, "./steiger_lam_borges_reversedir_results.csv")
#write.csv(steiger_mr_results_lam_kettunen, "./steiger_lam_kettunen_reversedir_results.csv")
