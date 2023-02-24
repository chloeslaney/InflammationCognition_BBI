# Project: inflammation and cognition
# Script extracts snps for creating PRS for cognitive domains
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

#####################################################
########### Load in available data ##################
#####################################################
# set working directory
setwd("INSERT FILE PATH HERE")

# Clear the work environment
rm(list = ls())

#######################
###### Cognition ###### 
#######################

# 1. Working memory = 3 SNPs
mahedy_wm <- read.table("2020.02.27_GWAS_Working_Memory_Data_Sheet.txt", sep = " ", header=TRUE)
mahedy_wm_pval <- mahedy_wm %>% filter (P < 5e-06)
mahedy_wm_pval <- rename(mahedy_wm_pval, pval = P, rsid = SNP, effect_allele = EFFECT_ALLELE, other_allele = OTHER_ALLELE, beta = BETA)
MAHEDY_WM_GWAS <- ld_clump(mahedy_wm_pval,clump_kb=1000, clump_r2=0.01)
#sum(MAHEDY_WM_GWAS$EAF < 0.01 | MAHEDY_WM_GWAS$EAF > 0.99) #0

# 2. Emotion recognition = 6 SNPs
mahedy_ERT <- read.table("2020.02.27_GWAS _Emotion_Recognition_Data_Sheet.txt", sep = " ", header=TRUE)
mahedy_ERT_pval <- mahedy_ERT %>% filter (P < 5e-06)
mahedy_ERT_pval <- rename(mahedy_ERT_pval, pval = P, rsid = SNP, effect_allele = EFFECT_ALLELE, other_allele = OTHER_ALLELE, beta = BETA)
MAHEDY_ERT_GWAS <- ld_clump(mahedy_ERT_pval,clump_kb=1000, clump_r2=0.01)
#sum(MAHEDY_ERT_GWAS$EAF < 0.01 | MAHEDY_ERT_GWAS$EAF > 0.99) #0

# 3. Response inhibition = 6 SNPs
mahedy_inhibition <- read.table("2020.02.27_GWAS_Response_Inhibition_Data_Sheet.txt", sep = " ", header=TRUE)
mahedy_inhibition_pval <- mahedy_inhibition %>% filter (P < 5e-06)
mahedy_inhibition_pval <- rename(mahedy_inhibition_pval, pval = P, rsid = SNP, effect_allele = EFFECT_ALLELE, other_allele = OTHER_ALLELE, beta = BETA)
MAHEDY_INHIBITION_GWAS <- ld_clump(mahedy_inhibition_pval,clump_kb=1000, clump_r2=0.01)
#sum(MAHEDY_INHIBITION_GWAS$EAF < 0.01 | MAHEDY_INHIBITION_GWAS$EAF > 0.99) #0

################################
###### Beta flip for PGRS ######
################################
#Adapted from Christina Dardani code - flips all beta values to be positive i.e. abs(i$beta) and alleles to agree with this change
#Does not flip EAF as this is not kept/saved
COGGWAS<- list(MAHEDY_WM_GWAS,MAHEDY_ERT_GWAS,MAHEDY_INHIBITION_GWAS) # list of each GWAS hits to apply beta flip to

#Applies beta flipping to each GWAS listed above
for(i in 1:length(COGGWAS)){
  
  COGGWAS[[i]]$effect_allele.corrected<- NA
  COGGWAS[[i]]$effect_allele.corrected<- ifelse(COGGWAS[[i]]$beta<0, COGGWAS[[i]]$other_allele, COGGWAS[[i]]$effect_allele)
  
  COGGWAS[[i]]$other_allele.corrected<- NA
  COGGWAS[[i]]$other_allele.corrected<- ifelse(COGGWAS[[i]]$effect_allele.corrected==COGGWAS[[i]]$effect_allele, COGGWAS[[i]]$other_allele, COGGWAS[[i]]$effect_allele)
  
  COGGWAS[[i]]$beta.corrected<- NA
  COGGWAS[[i]]$beta.corrected<- abs(COGGWAS[[i]]$beta)
}

#Saves output above in new data frame for each GWAS - ensure order is same as order in list above**
MAHEDY_WM_GWAS_CORRECTED <- COGGWAS[[1]]
MAHEDY_ERT_GWAS_CORRECTED <- COGGWAS[[2]]
MAHEDY_INHIBITION_GWAS_CORRECTED <- COGGWAS[[3]]

###########################################
###### Save GWAS as txt file for HPC ######
###########################################
ALL_COG_INSTRUMENT_SNP <- data.frame(rsid = c(MAHEDY_WM_GWAS_CORRECTED[,"rsid"], MAHEDY_ERT_GWAS_CORRECTED[,"rsid"], MAHEDY_INHIBITION_GWAS_CORRECTED[,"rsid"]))
ALL_COG_INSTRUMENT_SNP_DUPLICATESREMOVED <- distinct(ALL_COG_INSTRUMENT_SNP)
#write.table(ALL_COG_INSTRUMENT_SNP_DUPLICATESREMOVED, "mahedy_cognition_SNPs_duplicatesremoved.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#mahedy - working memory
MAHEDY_WM_OUTPUT <- MAHEDY_WM_GWAS_CORRECTED[, c('rsid', 'effect_allele.corrected', 'beta.corrected')]
#write.table(MAHEDY_WM_OUTPUT, "mahedy_wm_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#mahedy - emotion recognition task
MAHEDY_ERT_OUTPUT <- MAHEDY_ERT_GWAS_CORRECTED[, c('rsid', 'effect_allele.corrected', 'beta.corrected')]
#write.table(MAHEDY_ERT_OUTPUT, "mahedy_ert_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#mahedy - response inhibition
MAHEDY_INHIBITION_OUTPUT <- MAHEDY_INHIBITION_GWAS_CORRECTED[, c('rsid', 'effect_allele.corrected', 'beta.corrected')]
#write.table(MAHEDY_INHIBITION_OUTPUT, "mahedy_inhibition_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

###############################################################################
# For strand checking - list of all SNPs, effect allele and non-effect allele #
###############################################################################
mahedy_wm_output_checkstrand <- MAHEDY_WM_GWAS_CORRECTED[, c('rsid', 'effect_allele.corrected', 'other_allele.corrected')]
mahedy_ert_output_checkstrand <- MAHEDY_ERT_GWAS_CORRECTED[, c('rsid', 'effect_allele.corrected', 'other_allele.corrected')]
mahedy_inhibition_output_checkstrand <- MAHEDY_INHIBITION_GWAS_CORRECTED[, c('rsid', 'effect_allele.corrected', 'other_allele.corrected')]
CHECKSTRANDS_COGNITION <- rbind (mahedy_ert_output_checkstrand, mahedy_wm_output_checkstrand, mahedy_inhibition_output_checkstrand)
#write.table(CHECKSTRANDS_COGNITION, "CHECKSTRANDS_COGNITION.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#################################
##### remove ambiguous snps #####
#################################
ambiguous_snps <- read.table("alspac_ambiguous_snps_88distinct_proxiesadded.txt", sep = "\t", header=FALSE)
MAHEDY_WM_OUTPUT_AMBIGREMOVED <- MAHEDY_WM_OUTPUT[ ! MAHEDY_WM_OUTPUT$rsid %in% ambiguous_snps$V1, ]
MAHEDY_ERT_OUTPUT_AMBIGREMOVED <- MAHEDY_ERT_OUTPUT[ ! MAHEDY_ERT_OUTPUT$rsid %in% ambiguous_snps$V1, ]
MAHEDY_INHIBITION_OUTPUT_AMBIGREMOVED <- MAHEDY_INHIBITION_OUTPUT[ ! MAHEDY_INHIBITION_OUTPUT$rsid %in% ambiguous_snps$V1, ]

########################
###### Save files ######
########################
write.table(MAHEDY_WM_OUTPUT_AMBIGREMOVED, "mahedy_wm_ambigremoved_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(MAHEDY_ERT_OUTPUT_AMBIGREMOVED, "mahedy_ert_ambigremoved_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(MAHEDY_INHIBITION_OUTPUT_AMBIGREMOVED, "mahedy_inhibition_ambigremoved_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
