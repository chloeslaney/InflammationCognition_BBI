# Project: inflammation and cognition
# Script creates instruments for two-sample MR (cognition on inflammation)
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

################################################
########## PART 1: Load and check data #########
################################################
# set working directory
setwd("INSERT PATH HERE!")

# Clear the work environment
rm(list = ls())

#######################################################################
#  PART 2: GET SNP-EXPOSURE (cognition) AND SNP-OUTCOME SUMMARY DATA  #
#######################################################################
# criteria applied (p value< 5x10-8, clump criteria(R2=0.01, kb=1000), MAF>0.01.
# beta-flipped (EA related to increase in cognition)

#Exposure: Lam - General cognitive ability GWAS  - 250 SNPS
lam <- fread("Cognition_MTAG_allsnps_mtag_meta.txt")
EXPOSURE_LAM_GWAS <- lam %>% filter(mtag_pval < 5e-08)
EXPOSURE_LAM_GWAS <- rename(EXPOSURE_LAM_GWAS, rsid = SNP, effect_allele = A1, other_allele = A2, pval = mtag_pval)
EXPOSURE_LAM_GWAS <- ld_clump(EXPOSURE_LAM_GWAS,clump_kb=1000, clump_r2=0.01)
sum(EXPOSURE_LAM_GWAS$meta_freq < 0.01 | EXPOSURE_LAM_GWAS$meta_freq > 0.99) #0

#Outcome GWAS
LIGTHART_GWAS <- fread("GWAS_ligthart_crp.txt", header=TRUE) 
HAN_GWAS <- fread("GWAS_han_crp.txt", sep = ",", header=TRUE)
AHLUWALIA_GWAS <- fread("Ahluwalia_IL6_summary_results_discovery_GWAS_MA_CHARGE.txt", sep = "\t", header=TRUE)
BORGES_GWAS <- fread("GWAS_borges_glyca.txt", header=TRUE) 
KETTUNEN_GWAS <- fread("GWAS_kettunen_glyca.txt", header=TRUE)

#Rename columns
LIGTHART_GWAS <- rename(LIGTHART_GWAS, SNP = variant_id, se = standard_error, eaf = effect_allele_frequency, pval = p_value)
HAN_GWAS <- rename(HAN_GWAS, beta = BETA, se = SE, effect_allele = ALLELE1, other_allele = ALLELE0, eaf = A1FREQ, pval = P)
AHLUWALIA_GWAS <- rename(AHLUWALIA_GWAS, effect_allele = EA, other_allele = OA, eaf = EAF, pval = p)
BORGES_GWAS <- rename(BORGES_GWAS, SNP = variant_id, se = standard_error, eaf = effect_allele_frequency, pval = p_value)
KETTUNEN_GWAS <- rename(KETTUNEN_GWAS, SNP = variant_id, se = standard_error, eaf = effect_allele_frequency, pval = p_value)
EXPOSURE_LAM_GWAS <- rename(EXPOSURE_LAM_GWAS, SNP = rsid)

###################################################
# PART 3: IDENTIFY SNPS AVAILABLE IN OUTCOME GWAS #
###################################################
#GWAS - Borges has no missing SNPs
available_snps_ligthart <- merge(LIGTHART_GWAS,EXPOSURE_LAM_GWAS, by="SNP") #129/250
available_snps_han <- merge(HAN_GWAS,EXPOSURE_LAM_GWAS, by="SNP") #249/250
available_snps_borges <- merge(BORGES_GWAS, EXPOSURE_LAM_GWAS, by="SNP") #250/250
available_snps_kettunen <- merge(KETTUNEN_GWAS,EXPOSURE_LAM_GWAS, by="SNP") #248/250
available_snps_ahluwalia <- merge(AHLUWALIA_GWAS,EXPOSURE_LAM_GWAS, by="SNP") #137/250

#missing
ligthart_missing <- data.frame(SNP = c(setdiff(EXPOSURE_LAM_GWAS$SNP, available_snps_ligthart$SNP))) #121
han_missing <- data.frame(SNP = c(setdiff(EXPOSURE_LAM_GWAS$SNP, available_snps_han$SNP))) #1
kettunen_missing <- data.frame(SNP = c(setdiff(EXPOSURE_LAM_GWAS$SNP, available_snps_kettunen$SNP))) #2
ahluwalia_missing <- data.frame(SNP = c(setdiff(EXPOSURE_LAM_GWAS$SNP, available_snps_ahluwalia$SNP))) #113

#combined list of missing snps - 237 total - 121 distinct
total_missing_snps <- data.frame(SNP=c(ligthart_missing$SNP,han_missing$SNP,kettunen_missing$SNP,ahluwalia_missing$SNP))
total_missing_snps <- distinct(total_missing_snps)

###################################################################
# PART 4: FIND PROXY SNPS FOR THOSE NOT AVAILABLE IN OUTCOME GWAS #
###################################################################
#LDproxy_batch(total_missing_snps, pop="EUR", r2d = "r2", token = "INSERT TOKEN HERE!!", append=TRUE) #saves txt file as combined_query_snp_list
proxySNPs <- read.table("combined_query_snp_list.txt", header=TRUE, row.names=NULL)
proxySNPs <- filter(proxySNPs, RS_Number != ".") #removes snps with no rsid
proxySNPs <- proxySNPs %>% filter(R2 > 0.8) #keeps snps which have R2 > 0.8
proxySNPs_available_ligt <- merge(proxySNPs, LIGTHART_GWAS, by.x="RS_Number", by.y="SNP")#only keep proxies that are also available in outcome GWAS
proxySNPs_available_han <- merge(proxySNPs, HAN_GWAS, by.x="RS_Number", by.y="SNP")
proxySNPs_available_kett <- merge(proxySNPs, KETTUNEN_GWAS, by.x="RS_Number", by.y="SNP")
proxySNPs_available_ahluw <- merge(proxySNPs, AHLUWALIA_GWAS, by.x="RS_Number", by.y="SNP")

#################################################################
# PART 5: FIND PROXY SNPS FOR MISISNG SNPs IN EACH OUTCOME GWAS #
#################################################################
##############
## LIGTHART ##
##############
# 90 proxies available for 121 mising snps
# A. identifies proxy snps available in lam then keeps those with highest R2
proxySNPs_meetscriteria_ligthart <- merge(proxySNPs_available_ligt, lam, by.x="RS_Number", by.y="SNP")
proxySNPs_meetscriteria_ligthart <- proxySNPs_meetscriteria_ligthart %>% 
  group_by(query_snp) %>% #group by missing snps
  slice(c(which.max(R2))) #keep snp with highest R2

# B. extract list of proxy snps (of missing, how many found proxies for) = 121 snps not in ligthart - 90 proxies
ligthart_replace <- merge(ligthart_missing, proxySNPs_meetscriteria_ligthart, by.x="SNP", by.y="query_snp") #list of missing snps for which proxies available
ligthart_proxies <- data.frame(ligthart_replace$RS_Number) #list of proxy snps

# C. delete missing snps in lam gwas = 250 (lam GWAS) - 121 (snps not in ligt) = 129 snps
EXP_LAM_GWAS_LIGTHART_DELETEDMISS <- EXPOSURE_LAM_GWAS[ ! EXPOSURE_LAM_GWAS$SNP %in% ligthart_missing$SNP, ] #remove ligthart missing snps
EXP_LAM_GWAS_LIGTHART_DELETEDMISS$id <- NULL 

# D. find proxies in lam full sum stats = 90 proxies available
LAM_GWAS_PROXIES_LIGTHART <- merge(ligthart_proxies, lam, by.x="ligthart_replace.RS_Number", by.y="SNP")
LAM_GWAS_PROXIES_LIGTHART <- rename(LAM_GWAS_PROXIES_LIGTHART, SNP = ligthart_replace.RS_Number, effect_allele = A1, other_allele = A2, pval = mtag_pval)

# E. add proxies to Lam GWAS (with snps missing in outcome deleted) = 129+90=219
EXPOSURE_LAM_GWAS_LIGTHART_PROXIESADDED <- rbind(EXP_LAM_GWAS_LIGTHART_DELETEDMISS, LAM_GWAS_PROXIES_LIGTHART)
EXPOSURE_LAM_GWAS_LIGTHART_PROXIESADDED$id.exposure <- NA
EXPOSURE_LAM_GWAS_LIGTHART_PROXIESADDED <- rename(EXPOSURE_LAM_GWAS_LIGTHART_PROXIESADDED, EAF = meta_freq, beta = mtag_beta, se = mtag_se)

#########
## HAN ##
#########
# no proxy available for one missing snp
# A. identifies proxy snps available in lam then keeps those with highest R2
proxySNPs_meetscriteria_han <- merge(proxySNPs_available_han, lam, by.x="RS_Number", by.y="SNP")
proxySNPs_meetscriteria_han <- proxySNPs_meetscriteria_han %>% 
  group_by(query_snp) %>% #group by missing snps
  slice(c(which.max(R2))) #keep snp with highest R2

# B. extract list of proxy snps (of missing, how many found proxies for) = 1 snp not in han - 0 proxies available!
han_replace <- merge(han_missing, proxySNPs_meetscriteria_han, by.x="SNP", by.y="query_snp") #list of missing snps for which proxies available
han_proxies <- data.frame(han_replace$RS_Number) #list of proxy snps

# C. delete missing snps in lam gwas = 250 (lam GWAS) - 1 (snps not in han) = 249 snps
EXP_LAM_GWAS_HAN_DELETEDMISS <- EXPOSURE_LAM_GWAS[ ! EXPOSURE_LAM_GWAS$SNP %in% han_missing$SNP, ]
EXP_LAM_GWAS_HAN_DELETEDMISS$id <- NULL 
EXPOSURE_LAM_GWAS_HAN_PROXIESADDED <- EXP_LAM_GWAS_HAN_DELETEDMISS
EXPOSURE_LAM_GWAS_HAN_PROXIESADDED <- rename(EXPOSURE_LAM_GWAS_HAN_PROXIESADDED, EAF = meta_freq, beta = mtag_beta, se = mtag_se)

###############
## AHLUWALIA ##
###############
# 85 proxies available for 113 missing snps
# A. identifies proxy snps available in lam then keeps those with highest R2
proxySNPs_meetscriteria_ahluwalia <- merge(proxySNPs_available_ahluw, lam, by.x="RS_Number", by.y="SNP")
proxySNPs_meetscriteria_ahluwalia <- proxySNPs_meetscriteria_ahluwalia %>% 
  group_by(query_snp) %>% #group by missing snps
  slice(c(which.max(R2))) #keep snp with highest R2

# B. extract list of proxy snps (of missing, how many found proxies for) = 113 snps not in ahluwalia - 85 proxies
ahluwalia_replace <- merge(ahluwalia_missing, proxySNPs_meetscriteria_ahluwalia, by.x="SNP", by.y="query_snp")
ahluwalia_proxies <- data.frame(ahluwalia_replace$RS_Number) #list of proxy snps

# C. delete missing snps in lam gwas = 250 (lam GWAS) - 113 (snps not in ahluwalia) = 137 snps
EXP_LAM_GWAS_AHLUWALIA_DELETEDMISS <- EXPOSURE_LAM_GWAS[ ! EXPOSURE_LAM_GWAS$SNP %in% ahluwalia_missing$SNP, ]
EXP_LAM_GWAS_AHLUWALIA_DELETEDMISS$id <- NULL 

# D. find proxies in lam full sum stats = 85 proxies available
LAM_GWAS_PROXIES_AHLUWALIA <- merge(ahluwalia_proxies, lam, by.x="ahluwalia_replace.RS_Number", by.y="SNP")
LAM_GWAS_PROXIES_AHLUWALIA <- rename(LAM_GWAS_PROXIES_AHLUWALIA, SNP = ahluwalia_replace.RS_Number, effect_allele = A1, other_allele = A2, pval = mtag_pval)

# E. add proxies to Lam GWAS (with snps missing in outcome deleted) = 137+85=222
EXPOSURE_LAM_GWAS_AHLUWALIA_PROXIESADDED <- rbind(EXP_LAM_GWAS_AHLUWALIA_DELETEDMISS, LAM_GWAS_PROXIES_AHLUWALIA)
EXPOSURE_LAM_GWAS_AHLUWALIA_PROXIESADDED$id.exposure <- NA
EXPOSURE_LAM_GWAS_AHLUWALIA_PROXIESADDED <- rename(EXPOSURE_LAM_GWAS_AHLUWALIA_PROXIESADDED, EAF = meta_freq, beta = mtag_beta, se = mtag_se)

############
## BORGES ##
############
# no snps missing 
EXPOSURE_LAM_GWAS_BORGES_PROXIESADDED <- EXPOSURE_LAM_GWAS
EXPOSURE_LAM_GWAS_BORGES_PROXIESADDED <- rename(EXPOSURE_LAM_GWAS_BORGES_PROXIESADDED, EAF = meta_freq, beta = mtag_beta, se = mtag_se)

##############
## KETTUNEN ##
##############
# 2 proxies available for 2 missing snps
# A. identifies proxy snps available in lam then keeps those with highest R2
proxySNPs_meetscriteria_kettunen <- merge(proxySNPs_available_kett, lam, by.x="RS_Number", by.y="SNP")
proxySNPs_meetscriteria_kettunen <- proxySNPs_meetscriteria_kettunen %>% 
  group_by(query_snp) %>% #group by missing snps
  slice(c(which.max(R2))) #keep snp with highest R2

# B. extract list of proxy snps (of missing, how many found proxies for) = 2 snps not in kettunen - 2 proxies
kettunen_replace <- merge(kettunen_missing, proxySNPs_meetscriteria_kettunen, by.x="SNP", by.y="query_snp")
kettunen_proxies <- data.frame(kettunen_replace$RS_Number) #list of proxy snps

# C. delete missing snps in lam gwas = 250 (lam GWAS) - 2 (snps not in kettunen) = 248 snps
EXP_LAM_GWAS_KETTUNEN_DELETEDMISS <- EXPOSURE_LAM_GWAS[ ! EXPOSURE_LAM_GWAS$SNP %in% kettunen_missing$SNP, ]
EXP_LAM_GWAS_KETTUNEN_DELETEDMISS$id <- NULL 

# D. find proxies in lam full sum stats = 2 proxies available
LAM_GWAS_PROXIES_KETTUNEN <- merge(kettunen_proxies, lam, by.x="kettunen_replace.RS_Number", by.y="SNP")
LAM_GWAS_PROXIES_KETTUNEN <- rename(LAM_GWAS_PROXIES_KETTUNEN, SNP = kettunen_replace.RS_Number, effect_allele = A1, other_allele = A2, pval = mtag_pval)

# E. add proxies to Lam GWAS (with snps missing in outcome deleted) =248+2=250
EXPOSURE_LAM_GWAS_KETTUNEN_PROXIESADDED <- rbind(EXP_LAM_GWAS_KETTUNEN_DELETEDMISS, LAM_GWAS_PROXIES_KETTUNEN)
EXPOSURE_LAM_GWAS_KETTUNEN_PROXIESADDED$id.exposure <- NA
EXPOSURE_LAM_GWAS_KETTUNEN_PROXIESADDED <- rename(EXPOSURE_LAM_GWAS_KETTUNEN_PROXIESADDED, EAF = meta_freq, beta = mtag_beta, se = mtag_se)

######################################################################################
# PART 6: Beta flip - all beta values positive and alleles/eaf to agree with change  #
######################################################################################
#Code adapted from Christina Dardani's code - flips all beta values to be positive i.e. abs(i$beta) and alleles and eaf to agree with this change
#Positive value = better cognitive functioning
#QC:checked beta/eaf correct in top 10 of each instrument 
GWAS <- list(EXPOSURE_LAM_GWAS_LIGTHART_PROXIESADDED, EXPOSURE_LAM_GWAS_HAN_PROXIESADDED, EXPOSURE_LAM_GWAS_AHLUWALIA_PROXIESADDED, EXPOSURE_LAM_GWAS_BORGES_PROXIESADDED, EXPOSURE_LAM_GWAS_KETTUNEN_PROXIESADDED) # list of each GWAS hits

# applies beta flipping to each GWAS hits listed above
for(i in 1:length(GWAS)){
  
  GWAS[[i]]$effect_allele.corrected<- NA
  GWAS[[i]]$effect_allele.corrected<- ifelse(GWAS[[i]]$beta<0, GWAS[[i]]$other_allele, GWAS[[i]]$effect_allele)
  
  GWAS[[i]]$other_allele.corrected<- NA
  GWAS[[i]]$other_allele.corrected<- ifelse(GWAS[[i]]$effect_allele.corrected==GWAS[[i]]$effect_allele, GWAS[[i]]$other_allele, GWAS[[i]]$effect_allele)
  
  GWAS[[i]]$beta.corrected<- NA
  GWAS[[i]]$beta.corrected<- abs(GWAS[[i]]$beta)
  
  GWAS[[i]]$eaf.corrected<- NA
  GWAS[[i]]$eaf.corrected<- ifelse(GWAS[[i]]$effect_allele.corrected==GWAS[[i]]$effect_allele, GWAS[[i]]$EAF, (1 - GWAS[[i]]$EAF))
}

#Saves output above in new data frame for each GWAS
#Important: ensure order is same as order in list above!
EXP_LAM_LIGT_CORRECTED <- GWAS[[1]]
EXP_LAM_HAN_CORRECTED <- GWAS[[2]]
EXP_LAM_AHLUWALIA_CORRECTED <- GWAS[[3]]
EXP_LAM_BORGES_CORRECTED <- GWAS[[4]]
EXP_LAM_KETTUNEN_CORRECTED <- GWAS[[5]]

####################################
# PART 7: Re-format data for 2SMR #
####################################
#Minimum required info: SNP, beta, se, effect_allele (other useful: other allele, eaf, phenotype)
#Exposure (Lam) - for outcome (ligthart)
#QC: ensured CORRECTED beta/eaf used.
EXP_LAM_LIGT_CORRECTED$Phenotype <- "General_Cognition"
EXP_LAM_LIGT_CORRECTED <- EXP_LAM_LIGT_CORRECTED[, c("SNP","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_LAM_LIGT_CORRECTED <- rename(EXP_LAM_LIGT_CORRECTED,beta = beta.corrected,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Exposure (Lam) - for outcome (han)
EXP_LAM_HAN_CORRECTED$Phenotype <- "General_Cognition"
EXP_LAM_HAN_CORRECTED <- EXP_LAM_HAN_CORRECTED[, c("SNP","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_LAM_HAN_CORRECTED <- rename(EXP_LAM_HAN_CORRECTED,beta = beta.corrected,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Exposure (Lam) - for outcome (ahluwalia)
EXP_LAM_AHLUWALIA_CORRECTED$Phenotype <- "General_Cognition"
EXP_LAM_AHLUWALIA_CORRECTED <- EXP_LAM_AHLUWALIA_CORRECTED[, c("SNP","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_LAM_AHLUWALIA_CORRECTED <- rename(EXP_LAM_AHLUWALIA_CORRECTED,beta = beta.corrected,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Exposure (Lam) - for outcome (borges)
EXP_LAM_BORGES_CORRECTED$Phenotype <- "General_Cognition"
EXP_LAM_BORGES_CORRECTED <- EXP_LAM_BORGES_CORRECTED[, c("SNP","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_LAM_BORGES_CORRECTED <- rename(EXP_LAM_BORGES_CORRECTED,beta = beta.corrected,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Exposure (Lam) - for outcome (kettunen)
EXP_LAM_KETTUNEN_CORRECTED$Phenotype <- "General_Cognition"
EXP_LAM_KETTUNEN_CORRECTED <- EXP_LAM_KETTUNEN_CORRECTED[, c("SNP","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_LAM_KETTUNEN_CORRECTED <- rename(EXP_LAM_KETTUNEN_CORRECTED,beta = beta.corrected,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Outcomes
LIGTHART_GWAS$Phenotype <- "CRP (Ligthart)"
LIGTHART_GWAS <- LIGTHART_GWAS[, c("SNP","beta","se","effect_allele","other_allele","eaf","Phenotype","pval")]
HAN_GWAS$Phenotype <- "CRP (Han)"
HAN_GWAS <- HAN_GWAS[, c("SNP","beta","se","effect_allele","other_allele","eaf","Phenotype","pval")]
AHLUWALIA_GWAS$Phenotype <- "IL-6 (Ahluwalia)"
AHLUWALIA_GWAS <- AHLUWALIA_GWAS[, c("SNP","beta","se","effect_allele","other_allele","eaf","Phenotype","pval")]
BORGES_GWAS$Phenotype <- "GlycA (Borges)"
BORGES_GWAS <- BORGES_GWAS[, c("SNP","beta","se","effect_allele","other_allele","eaf","Phenotype","pval")]
KETTUNEN_GWAS$Phenotype <- "GlycA (Kettunen)"
KETTUNEN_GWAS <- KETTUNEN_GWAS[, c("SNP","beta","se","effect_allele","other_allele","eaf","Phenotype","pval")]

#####################################################
# Save outputs to be used for 2SMR - pvals included #
#####################################################
write.table(EXP_LAM_LIGT_CORRECTED, "EXP_LAM_FOR_LIGTHART_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_LAM_HAN_CORRECTED, "EXP_LAM_FOR_HAN_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_LAM_AHLUWALIA_CORRECTED, "EXP_LAM_FOR_AHLUWALIA_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_LAM_BORGES_CORRECTED, "EXP_LAM_FOR_BORGES_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_LAM_KETTUNEN_CORRECTED, "EXP_LAM_FOR_KETTUNEN_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(LIGTHART_GWAS, "OUTCOME_LIGTHART_GWAS_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(HAN_GWAS, "OUTCOME_HAN_GWAS_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(AHLUWALIA_GWAS, "OUTCOME_AHLUWALIA_GWAS_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(BORGES_GWAS, "OUTCOME_BORGES_GWAS_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(KETTUNEN_GWAS, "OUTCOME_KETTUNEN_GWAS_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
