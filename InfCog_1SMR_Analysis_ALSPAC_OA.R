# Project: inflammation and cognition
# Script for one-sample Mendelian randomization (parts code adapted from IEU short course)
# Chloe Slaney

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

################################################
########## PART 1: Load and check data #########
################################################
# set working directory
setwd("INSERT FILE LOCATION")

# Clear the work environment
rm(list = ls())

##################################
# ALSPAC pheno data= 15,645 obs. #
##################################
# note: observational data includes related individuals, genetic data does not.
ALSPAC_pheno <- read_stata("INSERT ALSPAC DATA FILE")
ALSPAC_pheno <- ALSPAC_pheno[, c('aln', 'qlet', 'rCRP_age24', 'rGlycA_age24', 'dp_24new', 'rERT', 'rSSRT', 'rBMI_age24', 'rIQ_age8', 'rCRP_age9', 'rCRP_age15', 'rCRP_age17', 'rGp_age7', 'rGp_age15', 'rGp_age17', 'rIL6_age9', 'ralcohol_age24', 'rsex', 'rethnicity', 'rmated', 'rmatSEP', 'rsmoke_age24', 'rCRP_age24_nlog', 'rGp_age24_nlog', 'zscoreCRP_age24', 'zscoreGlycA_age24', 'zscoreERT', 'zscoredp_24new', 'zscoreSSRT', 'rCRP_age24_sensitivity', 'rGlycA_age24_sensitivity')]
ALSPAC_pheno$alnqlet <- paste(ALSPAC_pheno$aln,ALSPAC_pheno$qlet, sep="")
table(ALSPAC_pheno$rCRP_age24 == 0)#false

# 1. check observations not missing, mean, sd, min, max (to 2dp)
# all match stata file
#crp - age 24
sum(! is.na(ALSPAC_pheno$rCRP_age24)) #3015
mean(ALSPAC_pheno$rCRP_age24, na.rm=TRUE) #2.28
sd(ALSPAC_pheno$rCRP_age24, na.rm=TRUE) #6.52
min(ALSPAC_pheno$rCRP_age24, na.rm=TRUE) #0.1
max(ALSPAC_pheno$rCRP_age24, na.rm=TRUE) #224.72

#glyca - age 24
sum(! is.na(ALSPAC_pheno$rGlycA_age24)) #3,258
mean(ALSPAC_pheno$rGlycA_age24, na.rm=TRUE) #1.23
sd(ALSPAC_pheno$rGlycA_age24, na.rm=TRUE) #0.17
min(ALSPAC_pheno$rGlycA_age24, na.rm=TRUE) #0.84
max(ALSPAC_pheno$rGlycA_age24, na.rm=TRUE) #2.25

#wm
sum(! is.na(ALSPAC_pheno$dp_24new)) #3,478
mean(ALSPAC_pheno$dp_24new, na.rm=TRUE) #2.76
sd(ALSPAC_pheno$dp_24new, na.rm=TRUE) #0.80
min(ALSPAC_pheno$dp_24new, na.rm=TRUE) #0
max(ALSPAC_pheno$dp_24new, na.rm=TRUE) #3.78

#ert
sum(! is.na(ALSPAC_pheno$rERT)) #3,613
mean(ALSPAC_pheno$rERT, na.rm=TRUE) #66.36
sd(ALSPAC_pheno$rERT, na.rm=TRUE) #7.89
min(ALSPAC_pheno$rERT, na.rm=TRUE) #25
max(ALSPAC_pheno$rERT, na.rm=TRUE) #88

#ssrt
sum(! is.na(ALSPAC_pheno$rSSRT)) #3,430
mean(ALSPAC_pheno$rSSRT, na.rm=TRUE) #258.72
sd(ALSPAC_pheno$rSSRT, na.rm=TRUE) #53.11
min(ALSPAC_pheno$rSSRT, na.rm=TRUE) #67
max(ALSPAC_pheno$rSSRT, na.rm=TRUE)#508

#sex (0=male; 1=female)
table(ALSPAC_pheno$rsex) #7,690 male, 7,348 female
 
#ethnicity (0=white; 1=non-white)
table(ALSPAC_pheno$rethnicity) #11,523 white, 613 non-white

#bmi_24
sum(! is.na(ALSPAC_pheno$rBMI_age24)) #3,974
mean(ALSPAC_pheno$rBMI_age24, na.rm=TRUE) #24.92
sd(ALSPAC_pheno$rBMI_age24, na.rm=TRUE) #5.08
min(ALSPAC_pheno$rBMI_age24, na.rm=TRUE) #13.68
max(ALSPAC_pheno$rBMI_age24, na.rm=TRUE) #63.74

#maternal education (0=uni degree; 1=A-level; 2=O-level; 3=vocational; 4=CSE)
# uni degree (1,608), A-level (2,793), O-level (4,323), vocational (1,229), CSE (1,750)
table(ALSPAC_pheno$rmated) 

#maternal SEP (0=professional; 1=intermediate; 2=skilled-nonmanual; 3=skilled-manual; 4=part-skilled; 5=unskilled)
#prof (595); intermediate (3,180); skilled-nonman (4,322); skilled-manual (790); part-skilled (997); unskilled (222)
table(ALSPAC_pheno$rmatSEP)

#smoke 24 (0=never smoked whole cigarette; 1=not smoked last 30 days; 2=not daily smoker; 3=daily smoker)
# never smoked whole cigarette (1,435); not last 30 days (1,390); not daily (642); daily (486)
table(ALSPAC_pheno$rsmoke_age24)

#alcohol 24
sum(! is.na(ALSPAC_pheno$ralcohol_age24)) #3,928
mean(ALSPAC_pheno$ralcohol_age24, na.rm=TRUE) #5.16
sd(ALSPAC_pheno$ralcohol_age24, na.rm=TRUE) #2.51
min(ALSPAC_pheno$ralcohol_age24, na.rm=TRUE) #0
max(ALSPAC_pheno$ralcohol_age24, na.rm=TRUE) #12

#IQ 8
sum(! is.na(ALSPAC_pheno$rIQ_age8)) #7,346
mean(ALSPAC_pheno$rIQ_age8, na.rm=TRUE) #103.97
sd(ALSPAC_pheno$rIQ_age8, na.rm=TRUE) #16.54
min(ALSPAC_pheno$rIQ_age8, na.rm=TRUE) #45
max(ALSPAC_pheno$rIQ_age8, na.rm=TRUE) #151

#crp - age 9
sum(! is.na(ALSPAC_pheno$rCRP_age9)) #5,080
mean(ALSPAC_pheno$rCRP_age9, na.rm=TRUE) #0.80
sd(ALSPAC_pheno$rCRP_age9, na.rm=TRUE) #2.72
min(ALSPAC_pheno$rCRP_age9, na.rm=TRUE) #0.01
max(ALSPAC_pheno$rCRP_age9, na.rm=TRUE) #67.44

#crp - age 15
sum(! is.na(ALSPAC_pheno$rCRP_age15)) #3,488
mean(ALSPAC_pheno$rCRP_age15, na.rm=TRUE) #1.24
sd(ALSPAC_pheno$rCRP_age15, na.rm=TRUE) #3.79
min(ALSPAC_pheno$rCRP_age15, na.rm=TRUE) #0.07
max(ALSPAC_pheno$rCRP_age15, na.rm=TRUE) #72.55

#crp - age 17
sum(! is.na(ALSPAC_pheno$rCRP_age17)) #3,285
mean(ALSPAC_pheno$rCRP_age17, na.rm=TRUE) #1.61
sd(ALSPAC_pheno$rCRP_age17, na.rm=TRUE) #4.97
min(ALSPAC_pheno$rCRP_age17, na.rm=TRUE) #0.02
max(ALSPAC_pheno$rCRP_age17, na.rm=TRUE) #176.1

#glyca - age 7 
sum(! is.na(ALSPAC_pheno$rGp_age7)) #5,518
mean(ALSPAC_pheno$rGp_age7, na.rm=TRUE)  #1.23
sd(ALSPAC_pheno$rGp_age7, na.rm=TRUE) #0.14
min(ALSPAC_pheno$rGp_age7, na.rm=TRUE) #0.86
max(ALSPAC_pheno$rGp_age7, na.rm=TRUE) #2.31

#glyca - age 15
sum(! is.na(ALSPAC_pheno$rGp_age15)) #3,363
mean(ALSPAC_pheno$rGp_age15, na.rm=TRUE) #1.21
sd(ALSPAC_pheno$rGp_age15, na.rm=TRUE) #0.13
min(ALSPAC_pheno$rGp_age15, na.rm=TRUE) #0.87
max(ALSPAC_pheno$rGp_age15, na.rm=TRUE) #1.90

#glyca - age 17
sum(! is.na(ALSPAC_pheno$rGp_age17)) #3,173 
mean(ALSPAC_pheno$rGp_age17, na.rm=TRUE) #1.22
sd(ALSPAC_pheno$rGp_age17, na.rm=TRUE) ##0.14
min(ALSPAC_pheno$rGp_age17, na.rm=TRUE) #0.87
max(ALSPAC_pheno$rGp_age17, na.rm=TRUE) #1.93

#il-6 - age 9
sum(! is.na(ALSPAC_pheno$rIL6_age9)) #5,070
mean(ALSPAC_pheno$rIL6_age9, na.rm=TRUE) #1.29
sd(ALSPAC_pheno$rIL6_age9, na.rm=TRUE) #1.59
min(ALSPAC_pheno$rIL6_age9, na.rm=TRUE) #0.007
max(ALSPAC_pheno$rIL6_age9, na.rm=TRUE) #20.05

#########################################################
#      Correlations between CRP/GlycA over time         #
#########################################################
# spearman's for crp as not normal distrib - note: not restricted to those who have genetics data - to get this value re-run using ALSPAC_merged.
cor(ALSPAC_pheno$rCRP_age24, ALSPAC_pheno$rCRP_age17, method="spearman", use="complete.obs") #0.32
cor(ALSPAC_pheno$rCRP_age24, ALSPAC_pheno$rCRP_age15, method="spearman", use="complete.obs") #0.28
cor(ALSPAC_pheno$rCRP_age24, ALSPAC_pheno$rCRP_age9, method="spearman", use="complete.obs") #0.29
cor(ALSPAC_pheno$rGlycA_age24, ALSPAC_pheno$rGp_age17, method="pearson", use="complete.obs") #0.39
cor(ALSPAC_pheno$rGlycA_age24, ALSPAC_pheno$rGp_age15, method="pearson", use="complete.obs") #0.36
cor(ALSPAC_pheno$rGlycA_age24, ALSPAC_pheno$rGp_age7, method="pearson", use="complete.obs") #0.22

##########################################################
# ALSPAC PRS data - with proxies added where appropriate #
##########################################################
# ligthart - cis
# identical(ligthartcis_PRS$FID,ligthartcis_PRS$IID) #yes
ligthartcis_PRS <- read.table("ligthartcis_score_sum.profile", header=TRUE)
ligthartcis_PRS <- rename(ligthartcis_PRS, ligthartcis_PRS=SCORESUM, alnqlet=IID)
LIGTHARTCIS_PRS <- ligthartcis_PRS[, c('alnqlet', 'ligthartcis_PRS')]

# ligthart - genomewide (proxies added)
ligthart_GW_PRS <- read.table("ligthart_genomewide_proxiesadded_score_sum.profile", header=TRUE)
identical(ligthart_GW_PRS$FID, ligthart_GW_PRS$IID)
ligthart_GW_PRS <- rename(ligthart_GW_PRS, ligthart_GW_PRS=SCORESUM, alnqlet=IID)
LIGTHART_GW_PRS <- ligthart_GW_PRS[, c('alnqlet','ligthart_GW_PRS')]

# han - cis (proxies added)
hancis_PRS <- read.table("hancis_proxiesadded_score_sum.profile", header=TRUE)
identical(hancis_PRS$FID, hancis_PRS$IID)
hancis_PRS <- rename(hancis_PRS, hancis_PRS=SCORESUM, alnqlet=IID)
HANCIS_PRS <- hancis_PRS[, c('alnqlet','hancis_PRS')]

# han - genomewide (proxies added)
han_GW_PRS <- read.table("han_genomewide_proxiesadded_score_sum.profile", header=TRUE)
identical(han_GW_PRS$FID, han_GW_PRS$IID)
han_GW_PRS <- rename(han_GW_PRS, han_GW_PRS=SCORESUM, alnqlet=IID)
HAN_GW_PRS <- han_GW_PRS[, c('alnqlet','han_GW_PRS')]

# ahluwalia - cis
ahluwaliacis_PRS <- read.table("ahluwaliacis_score_sum.profile", header=TRUE)
identical(ahluwaliacis_PRS$FID, ahluwaliacis_PRS$IID)
ahluwaliacis_PRS <- rename(ahluwaliacis_PRS, ahluwaliacis_PRS=SCORESUM, alnqlet=IID)
AHLUWALIACIS_PRS <- ahluwaliacis_PRS[, c('alnqlet','ahluwaliacis_PRS')]

# ahluwalia - genomewide
ahluwalia_GW_PRS <- read.table("ahluwalia_genomewide_score_sum.profile", header=TRUE)
identical(ahluwalia_GW_PRS$FID, ahluwalia_GW_PRS$IID)
ahluwalia_GW_PRS <- rename(ahluwalia_GW_PRS, ahluwalia_GW_PRS=SCORESUM, alnqlet=IID)
AHLUWALIA_GW_PRS <- ahluwalia_GW_PRS[, c('alnqlet','ahluwalia_GW_PRS')]

# borges (proxies added)
borges_PRS <- read.table("borges_eafchecked_proxiesadded_score_sum.profile", header=TRUE)
identical(borges_PRS$FID, borges_PRS$IID)
borges_PRS <- rename(borges_PRS, borges_PRS=SCORESUM, alnqlet=IID)
BORGES_PRS <- borges_PRS[, c('alnqlet','borges_PRS')]

# kettunen
kettunen_PRS <- read.table("kettunen_score_sum.profile", header=TRUE)
identical(kettunen_PRS$FID, kettunen_PRS$IID)
kettunen_PRS <- rename(kettunen_PRS, kettunen_PRS=SCORESUM, alnqlet=IID)
KETTUNEN_PRS <- kettunen_PRS[, c('alnqlet','kettunen_PRS')]

# rosa
rosa_PRS <- read.table("rosa_score_sum.profile", header=TRUE)
identical(rosa_PRS$FID, rosa_PRS$IID)
rosa_PRS <- rename(rosa_PRS, rosa_PRS=SCORESUM, alnqlet=IID)
ROSA_PRS <- rosa_PRS[, c('alnqlet','rosa_PRS')]

# swerdlow
swerdlow_PRS <- read.table("swerdlow_score_sum.profile", header=TRUE)
identical(swerdlow_PRS$FID, swerdlow_PRS$IID)
swerdlow_PRS <- rename(swerdlow_PRS, swerdlow_PRS=SCORESUM, alnqlet=IID)
SWERDLOW_PRS <- swerdlow_PRS[, c('alnqlet','swerdlow_PRS')]

# sarwar
sarwar_PRS <- read.table("sarwar_score_sum.profile", header=TRUE)
identical(sarwar_PRS$FID, sarwar_PRS$IID)
sarwar_PRS <- rename(sarwar_PRS, sarwar_PRS=SCORESUM, alnqlet=IID)
SARWAR_PRS <- sarwar_PRS[, c('alnqlet','sarwar_PRS')]

# mahedy wm
mahedy_wm_PRS <- read.table("mahedy_wm_score_sum.profile", header=TRUE)
identical(mahedy_wm_PRS$FID, mahedy_wm_PRS$IID)
mahedy_wm_PRS <- rename(mahedy_wm_PRS, mahedy_wm_PRS=SCORESUM, alnqlet=IID)
MAHEDY_WM_PRS <- mahedy_wm_PRS[, c('alnqlet','mahedy_wm_PRS')]

# mahedy ert
mahedy_ert_PRS <- read.table("mahedy_ert_score_sum.profile", header=TRUE)
identical(mahedy_ert_PRS$FID, mahedy_ert_PRS$IID)
mahedy_ert_PRS <- rename(mahedy_ert_PRS, mahedy_ert_PRS=SCORESUM, alnqlet=IID)
MAHEDY_ERT_PRS <- mahedy_ert_PRS[, c('alnqlet','mahedy_ert_PRS')]

# mahedy inhibition
mahedy_inhib_PRS <- read.table("mahedy_inhibition_score_sum.profile", header=TRUE)
identical(mahedy_inhib_PRS$FID, mahedy_inhib_PRS$IID)
mahedy_inhib_PRS <- rename(mahedy_inhib_PRS, mahedy_inhib_PRS=SCORESUM, alnqlet=IID)
MAHEDY_INHIB_PRS <- mahedy_inhib_PRS[, c('alnqlet','mahedy_inhib_PRS')]

#principal components
PC <- read.table("children_eigenvec.txt", header=FALSE) #available from ALSPAC
PC_10 <- rename(PC, alnqlet=V1, PC1=V3, PC2=V4, PC3=V5, PC4=V6, PC5=V7, PC6=V8, PC7=V9, PC8=V10, PC9=V11, PC10=V12)
PC_10 <- PC_10[, c('alnqlet','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')]

####################################################
# merging phenotype and genetic data = 8,130 observations #
####################################################
# genetic (8,253 obs) once merged with phenotype data (8,130) = removes 33 (no pheno data) and rest are related.
ALSPAC_merged <- merge(ALSPAC_pheno, LIGTHARTCIS_PRS, by = "alnqlet") #add ligthartcis
ALSPAC_merged <- merge(ALSPAC_merged, LIGTHART_GW_PRS, by = "alnqlet") #add ligthart GW
ALSPAC_merged <- merge(ALSPAC_merged, HANCIS_PRS, by = "alnqlet") #add hancis
ALSPAC_merged <- merge(ALSPAC_merged, HAN_GW_PRS, by = "alnqlet") #add han GW
ALSPAC_merged <- merge(ALSPAC_merged, AHLUWALIACIS_PRS, by = "alnqlet") #add ahluwaliacis
ALSPAC_merged <- merge(ALSPAC_merged, AHLUWALIA_GW_PRS, by = "alnqlet") #add ahluwalia GW
ALSPAC_merged <- merge(ALSPAC_merged, BORGES_PRS, by = "alnqlet")  #add borges
ALSPAC_merged <- merge(ALSPAC_merged, KETTUNEN_PRS, by = "alnqlet") #add kettunen
ALSPAC_merged <- merge(ALSPAC_merged, ROSA_PRS, by = "alnqlet") #add rosa
ALSPAC_merged <- merge(ALSPAC_merged, SWERDLOW_PRS, by = "alnqlet") #add swerdlow
ALSPAC_merged <- merge(ALSPAC_merged, SARWAR_PRS, by = "alnqlet") #add sarwar
ALSPAC_merged <- merge(ALSPAC_merged, MAHEDY_WM_PRS, by = "alnqlet") # add mahedy wm
ALSPAC_merged <- merge(ALSPAC_merged, MAHEDY_ERT_PRS, by = "alnqlet") #add mahedy ert
ALSPAC_merged <- merge(ALSPAC_merged, MAHEDY_INHIB_PRS, by = "alnqlet") #add mahedy inhibition
ALSPAC_merged <- merge(ALSPAC_merged, PC_10, by = "alnqlet") #add top 10 principal components

###############################################
## Check distributions: exposure and outcome ##
###############################################
#inflammation - crp and il-6 are highly skewed, glyca fine.
hist(ALSPAC_merged$rCRP_age24) #crp highly skewed
hist(ALSPAC_merged$rCRP_age24_nlog) #log crp ok
hist(ALSPAC_merged$rGlycA_age24) #glyca ok
hist(ALSPAC_merged$rIL6_age9) #il6 highly skewed
ALSPAC_merged$rIL6_age9_nlog <- log(ALSPAC_merged[,c("rIL6_age9")]) #create log for IL-6
hist(ALSPAC_merged$rIL6_age9_nlog) #log il6 ok

# sensitivity analysis - removing individuals CRP > 10 at age 24
hist(ALSPAC_merged$rCRP_age24_sensitivity) #highly skewed
ALSPAC_merged$rCRP_age24_sensitivity_nlog <- log(ALSPAC_merged[,c("rCRP_age24_sensitivity")])
hist(ALSPAC_merged$rGlycA_age24_sensitivity) #ok

# Outcomes
hist(ALSPAC_merged$dp_24new) #working memory outcome skewed
ALSPAC_merged$dp_24new_nlog <- log(ALSPAC_merged[,c("dp_24new")])
hist(ALSPAC_merged$dp_24new_nlog) #log working memory outcome still skewed
hist(ALSPAC_merged$rERT) #emotion recognition outcome ok
hist(ALSPAC_merged$rSSRT) #response inhibition outcome ok

# PRS
# some skewed, likely due to small N snps to create PRS.
hist(ALSPAC_merged$ligthart_GW_PRS) #ok
hist(ALSPAC_merged$hancis_PRS) #ok
hist(ALSPAC_merged$han_GW_PRS) #ok
hist(ALSPAC_merged$rosa_PRS) #ok
hist(ALSPAC_merged$borges_PRS) #ok
hist(ALSPAC_merged$kettunen_PRS) #ok

## sqrt transform - not substantially improve so used original
#ligthart - cis
hist(ALSPAC_merged$ligthartcis_PRS)
ALSPAC_merged$ligthartcis_PRS_sqrtlog <- sqrt(ALSPAC_merged[,c("ligthartcis_PRS")])
hist(ALSPAC_merged$ligthartcis_PRS_sqrtlog)

# ahluwalia
hist(ALSPAC_merged$ahluwaliacis_PRS) 
hist(ALSPAC_merged$ahluwalia_GW_PRS) 
ALSPAC_merged$ahluwalia_GW_PRS_sqrtlog <- sqrt(ALSPAC_merged[,c("ahluwalia_GW_PRS")])
ALSPAC_merged$ahluwaliacis_PRS_sqrtlog <- sqrt(ALSPAC_merged[,c("ahluwaliacis_PRS")])

#swerdlow and sarwar
hist(ALSPAC_merged$swerdlow_PRS) 
hist(ALSPAC_merged$sarwar_PRS) 
ALSPAC_merged$swerdlow_PRS_sqrtlog <- sqrt(ALSPAC_merged[,c("swerdlow_PRS")])
ALSPAC_merged$sarwar_PRS_sqrtlog <- sqrt(ALSPAC_merged[,c("sarwar_PRS")])

###############################
### Create standardized PRS ###
###############################
# Only used for cognitive outcomes.
ALSPAC_merged$ligthart_GW_PRS_sd <- scale(ALSPAC_merged$ligthart_GW_PRS)
ALSPAC_merged$ligthartcis_PRS_sd <- scale(ALSPAC_merged$ligthartcis_PRS)
ALSPAC_merged$hancis_PRS_sd <- scale(ALSPAC_merged$hancis_PRS)
ALSPAC_merged$han_GW_PRS_sd <- scale(ALSPAC_merged$han_GW_PRS)
ALSPAC_merged$borges_PRS_sd <- scale(ALSPAC_merged$borges_PRS)
ALSPAC_merged$ahluwalia_GW_PRS_sd <- scale(ALSPAC_merged$ahluwalia_GW_PRS)
ALSPAC_merged$ahluwaliacis_PRS_sd <- scale(ALSPAC_merged$ahluwaliacis_PRS)
ALSPAC_merged$swerdlow_PRS_sd <- scale(ALSPAC_merged$swerdlow_PRS)
ALSPAC_merged$sarwar_PRS_sd <- scale(ALSPAC_merged$sarwar_PRS)
ALSPAC_merged$rosa_PRS_sd <- scale(ALSPAC_merged$rosa_PRS)
ALSPAC_merged$kettunen_PRS_sd <- scale(ALSPAC_merged$kettunen_PRS)
# cognitive outcomes
ALSPAC_merged$dp_24new_sd <- scale(ALSPAC_merged$dp_24new)
ALSPAC_merged$rERT_sd <- scale(ALSPAC_merged$rERT)
ALSPAC_merged$rSSRT_sd <- scale(ALSPAC_merged$rSSRT)

#################################################
# PART 2 - are PRS' associated with phenotypes? #
#################################################

# INFLAMMATION
# log CRP - age 24
summary(lm(rCRP_age24_nlog ~ ligthart_GW_PRS, ALSPAC_merged)) #f-stat(99.94, r2=0.043)
summary(lm(rCRP_age24_nlog ~ ligthartcis_PRS, ALSPAC_merged)) #f-stat (26.72, r2=0.012)
summary(lm(rCRP_age24_nlog ~ hancis_PRS, ALSPAC_merged)) #f-stat(15.75,r2=.007)
summary(lm(rCRP_age24_nlog ~ han_GW_PRS, ALSPAC_merged)) #f-stat(79.9, r2=0.035)

#log IL6 - age 9
summary(lm(rIL6_age9_nlog ~ ahluwaliacis_PRS, ALSPAC_merged)) #f-stat (87.57, r2=0.021)
summary(lm(rIL6_age9_nlog ~ ahluwalia_GW_PRS, ALSPAC_merged)) #f-stat (77.27, r2=0.018)
summary(lm(rIL6_age9_nlog ~ rosa_PRS, ALSPAC_merged)) #f-stat (76.11, r2=0.018)
summary(lm(rIL6_age9_nlog ~ swerdlow_PRS, ALSPAC_merged)) #f-stat (95.41, r2=0.022)
summary(lm(rIL6_age9_nlog ~ sarwar_PRS, ALSPAC_merged)) #f-stat (87.87, r2=0.021)

#glyca
summary(lm(rGlycA_age24 ~ borges_PRS, ALSPAC_merged)) #f-stat (65.69, r2=0.027)
summary(lm(rGlycA_age24 ~ kettunen_PRS, ALSPAC_merged)) #f-stat (51.75, r2=0.021)

##### sensitivity analysis - removing > 10 mg/L in CRP at age 24 #####
summary(lm(rCRP_age24_sensitivity_nlog ~ ligthart_GW_PRS, ALSPAC_merged)) #f-stat (93.13, r2=0.042)
summary(lm(rCRP_age24_sensitivity_nlog ~ ligthartcis_PRS, ALSPAC_merged)) #f-stat (24.36, r2=0.011)
summary(lm(rCRP_age24_sensitivity_nlog ~ hancis_PRS, ALSPAC_merged)) #f-stat (12.03, r2=0.006)
summary(lm(rCRP_age24_sensitivity_nlog ~ han_GW_PRS, ALSPAC_merged)) #f-stat (64.75, r2=0.029)
summary(lm(rGlycA_age24_sensitivity ~ borges_PRS, ALSPAC_merged)) #f-stat (66.72, r2=0.030)
summary(lm(rGlycA_age24_sensitivity ~ kettunen_PRS, ALSPAC_merged)) #f-stat (50.37, r2=0.023)

# COGNITION - note: ALSPAC was used as GWAS to get SNPs
summary(lm(dp_24new ~ mahedy_wm_PRS, ALSPAC_merged)) #f-stat (59.67, r2=0.023)
summary(lm(rERT ~ mahedy_ert_PRS, ALSPAC_merged)) #f-stat (112.1, r2=0.041)
summary(lm(rSSRT ~ mahedy_inhib_PRS, ALSPAC_merged)) #f-stat (112.1, r2=0.043)

#################################################
# PART 4 - are PRS associated with confounders? #
#################################################
# potential confounds: sex, ethnicity, BMI, maternal education, maternal SEP, smoking, alcohol, IQ
# summary: PRS not strong evidence of association with confounds (except for maternal education)

#ligthart - cis
# no evidence
summary(lm(ligthartcis_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthartcis_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#ligthart - genome-wide
# maternal education (p=0.075)
summary(lm(ligthart_GW_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ligthart_GW_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

# han - cis
# no evidence
summary(lm(hancis_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(hancis_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

# han - genome-wide
# BMI (p=0.064), maternal education (p=0.00006), alcohol use (p=0.015), IQ (p=0.004)
summary(lm(han_GW_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(han_GW_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

# ahluwalia - cis
# no evidence
summary(lm(ahluwaliacis_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwaliacis_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#ahluwalia - genome-wide
# no evidence
summary(lm(ahluwalia_GW_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(ahluwalia_GW_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#borges
#alcohol (p=0.074)
summary(lm(borges_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(borges_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#kettunen
#no evidence
summary(lm(kettunen_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(kettunen_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#rosa
#no evidence
summary(lm(rosa_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(rosa_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#swerdlow
#no evidence
summary(lm(swerdlow_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(swerdlow_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#sarwar
#no evidence
summary(lm(sarwar_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~rIQ_age8 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(sarwar_PRS ~PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#mahedy - wm
#no evidence
summary(lm(mahedy_wm_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_wm_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_wm_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_wm_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_wm_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_wm_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_wm_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#mahedy - ert
#BMI (p=0.061), maternal SEP (p=0.013)
summary(lm(mahedy_ert_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_ert_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_ert_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_ert_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_ert_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_ert_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_ert_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#mahedy - ssrt
#no evidence
summary(lm(mahedy_inhib_PRS ~rsex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_inhib_PRS ~rethnicity + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_inhib_PRS ~rBMI_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_inhib_PRS ~rmated + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_inhib_PRS ~rmatSEP + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_inhib_PRS ~rsmoke_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))
summary(lm(mahedy_inhib_PRS ~ralcohol_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ALSPAC_merged))

#################################################################
# PART 5 - Does inflammation have a causal effect on cognition? #
#################################################################
# no evidence to suggest inflammation has a causal effect on cognition.
# conclusions same when 2SLS regressions include top 10 PCs

#######################
### primary analysis ##
#######################
#crp - no evidence
#ligthart - cis
summary(tsls_ligthartcis_wm <- ivreg(dp_24new_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ligthartcis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ligthartcis_ert <- ivreg(rERT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ligthartcis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ligthartcis_ssrt <- ivreg(rSSRT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ligthartcis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#han - cis
summary(tsls_hancis_wm <- ivreg(dp_24new_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | hancis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_hancis_ert <- ivreg(rERT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | hancis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_hancis_ssrt <- ivreg(rSSRT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | hancis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#glyca - no evidence
summary(tsls_borges_wm <- ivreg(dp_24new_sd ~ rGlycA_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | borges_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_borges_ert <- ivreg(rERT_sd ~ rGlycA_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | borges_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_borges_ssrt <- ivreg(rSSRT_sd ~ rGlycA_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | borges_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#il6 - ahluwalia - no evidence
summary(tsls_ahluwaliacis_wm <- ivreg(dp_24new_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ahluwaliacis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ahluwaliacis_ert <- ivreg(rERT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ahluwaliacis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ahluwaliacis_ssrt <- ivreg(rSSRT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ahluwaliacis_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#rosa - no evidence
summary(tsls_rosa_wm <- ivreg(dp_24new_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | rosa_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_rosa_ert <- ivreg(rERT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | rosa_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_rosa_ssrt <- ivreg(rSSRT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | rosa_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

##########################
### secondary analysis ###
##########################
# no evidence that inflammation has a causal effect on cognition
#crp
#ligthart - genome-wide
summary(tsls_ligthart_GW_wm <- ivreg(dp_24new_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ligthart_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ligthart_GW_ert <- ivreg(rERT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ligthart_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ligthart_GW_ssrt <- ivreg(rSSRT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ligthart_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#han - genome-wide
summary(tsls_han_GW_wm <- ivreg(dp_24new_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | han_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_han_GW_ert <- ivreg(rERT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | han_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_han_ssrt <- ivreg(rSSRT_sd ~ rCRP_age24_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | han_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#glyca
summary(tsls_kett_wm <- ivreg(dp_24new_sd ~ rGlycA_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | kettunen_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_kett_ert <- ivreg(rERT_sd ~ rGlycA_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | kettunen_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_kett_ssrt <- ivreg(rSSRT_sd ~ rGlycA_age24 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | kettunen_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#il6
summary(tsls_ahluwalia_GW_wm <- ivreg(dp_24new_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ahluwalia_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ahluwalia_GW_ert <- ivreg(rERT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ahluwalia_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ahluwalia_GW_ssrt <- ivreg(rSSRT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | ahluwalia_GW_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#swerdlow
summary(tsls_swerdlow_wm <- ivreg(dp_24new_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10| swerdlow_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_swerdlow_ert <- ivreg(rERT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | swerdlow_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_swerdlow_ssrt <- ivreg(rSSRT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | swerdlow_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#sarwar
summary(tsls_sarwar_wm <- ivreg(dp_24new_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | sarwar_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_sarwar_ert <- ivreg(rERT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | sarwar_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_sarwar_ssrt <- ivreg(rSSRT_sd ~ rIL6_age9_nlog + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | sarwar_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#################################################################
# PART 6 - Does cognition have a causal effect on inflammation? #
#################################################################
#wm
summary(tsls_wm_crp <- ivreg(rCRP_age24_nlog_sd ~ dp_24new + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_wm_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_wm_glyca <- ivreg(rGlycA_age24_sd ~ dp_24new + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_wm_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_wm_il6 <- ivreg(rIL6_age9_nlog_sd ~ dp_24new + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_wm_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#ert
summary(tsls_ert_crp <- ivreg(rCRP_age24_nlog_sd ~ rERT + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_ert_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ert_glyca <- ivreg(rGlycA_age24_sd ~ rERT + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_ert_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_ert_il6 <- ivreg(rIL6_age9_nlog_sd ~ rERT + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_ert_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)

#inhibition
summary(tsls_inhib_crp <- ivreg(rCRP_age24_nlog_sd ~ rSSRT + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_inhib_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_inhib_glyca <- ivreg(rGlycA_age24_sd ~ rSSRT + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_inhib_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
summary(tsls_inhib_il6 <- ivreg(rIL6_age9_nlog_sd ~ rSSRT + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | mahedy_inhib_PRS + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=ALSPAC_merged), diagnostics=T)
