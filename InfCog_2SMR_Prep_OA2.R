# Project: inflammation and cognition
# Script creates instruments for two-sample MR 
# Packages: R version 4.1.1 & TwoSampleMR _0.5.6

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
setwd("INSERT PATH LOCATION HERE")

# Clear the work environment
rm(list = ls())

#################################################################################
#  PART 2: GET SNP-EXPOSURE (Inflammatory markers) AND SNP-OUTCOME SUMMARY DATA #
#################################################################################
# criteria applied (p value < 5x10-8, clump criteria(R2=0.01, kb=1000), MAF>0.01.
# beta-flipped (EA related to increase in inflammation)
# divided into cis variants (within 1mB of protein-coding gene) where possible
# section taken and adapted from script used to create SNP list for PRS (Inflammation_cognition)

###############################
# C-reactive protein = 2 GWAS #
###############################
# 1.Ligthart CRP GWAS - full summary statistics = 78 SNPs
#data source = VCF file from IEU database converted to txt using bcftools (separate script for this)
#ligthart_IEU <- tophits(id="ieu-b-35", pval=5e-08, r2=0.01, kb=1000) 
ligthart <- read.table("GWAS_ligthart_crp.txt", header=TRUE) #reads data
ligthart_pval <- ligthart %>% filter(p_value < 5e-08) #filter to include SNPs p < 5*10-8
ligthart_pval <- rename(ligthart_pval, pval = p_value, rsid = variant_id, CHR=chromosome, BP=base_pair_location, EAF = effect_allele_frequency) #rename for ld_clump
LIGTHART_GWAS <- ld_clump(ligthart_pval,clump_kb=1000, clump_r2=0.01) #clump SNPs
sum(LIGTHART_GWAS$EAF < 0.01 | LIGTHART_GWAS$EAF > 0.99) #0
#setdiff(ligthart_IEU$rsid, LIGTHART_GWAS$rsid)

# 2. Han CRP GWAS - full summary statistics (EA=ALLELE1) = 552 SNPs
#data source = requested from authors and sent link (used curl on linux to download to txt file)
#due to concerns that curl does not download recursively, re-downloaded using wget = same # snps
han <- read.table("GWAS_han_crp.txt", sep = ",", header=TRUE)
#han_wget <- read.table("CRP_all_SNPs.txt", sep = ",", header=TRUE) #used wget instead of curl - same # snps
han_pval <- han %>% filter (P < 5e-08)
han_pval <- rename(han_pval, pval = P, rsid = SNP, effect_allele = ALLELE1, other_allele = ALLELE0, beta = BETA, EAF = A1FREQ)
HAN_GWAS <- ld_clump(han_pval,clump_kb=1000, clump_r2=0.01)
sum(HAN_GWAS$EAF < 0.01 | HAN_GWAS$EAF > 0.99) #0

##########################
# Interleukin-6 = 1 GWAS #
##########################
# 1. Ahluwalia IL-6 GWAS - full summary statistics = 3 SNPs
#data source = requested from authors and sent as txt file
ahluwalia <- read.table("Ahluwalia_IL6_summary_results_discovery_GWAS_MA_CHARGE.txt", sep = "\t", header=TRUE)
ahluwalia_pval <- ahluwalia %>% filter(p < 5e-08)
ahluwalia_pval <- rename(ahluwalia_pval, pval = p, rsid = SNP, effect_allele = EA, other_allele = OA)
AHLUWALIA_GWAS <- ld_clump(ahluwalia_pval,clump_kb=1000, clump_r2=0.01)
sum(AHLUWALIA_GWAS$EAF < 0.01 | AHLUWALIA_GWAS$EAF > 0.99) #0

#################################
# Glycoprotein Acetyls = 2 GWAS #
#################################
# 1. Borges - full summary statistics = 88 SNPs (or 87 if removed SNP EAF < 0.01)
#data source = VCF file from IEU database converted to txt using bcftools (separate script for this)
#borges_IEU <- tophits(id="met-d-GlycA", r2=0.01, kb=1000)
borges <- read.table("GWAS_borges_glyca.txt", header=TRUE) 
borges_pval <- borges %>% filter(p_value < 5e-08) 
borges_pval <- rename(borges_pval, pval = p_value, rsid = variant_id, CHR = chromosome, BP = base_pair_location, EAF = effect_allele_frequency)
BORGES_GWAS <- ld_clump(borges_pval,clump_kb=1000, clump_r2=0.01)
sum(BORGES_GWAS$EAF < 0.01 | BORGES_GWAS$EAF > 0.99) #1 SNP! 
BORGES_GWAS <- BORGES_GWAS %>% filter(EAF > 0.01) # only keep SNPs with EAF > 0.01

# 2. Kettunen - full summary statistics = 10 SNPs
#data source = VCF file from IEU database converted to txt using bcftools (separate script for this)
#kettunen_IEU <- tophits(id="met-c-863", r2 = 0.01, kb = 1000)
kettunen <- read.table("GWAS_kettunen_glyca.txt", header=TRUE)
kettunen_pval <- kettunen %>% filter(p_value < 5e-08)
kettunen_pval <- rename(kettunen_pval, pval = p_value, rsid = variant_id, CHR = chromosome, BP = base_pair_location, EAF = effect_allele_frequency)
KETTUNEN_GWAS <- ld_clump(kettunen_pval,clump_kb=1000, clump_r2=0.01)
sum(KETTUNEN_GWAS$EAF < 0.01 | KETTUNEN_GWAS$EAF > 0.99) #0

#####################
# OTHER INSTRUMENTS #
#####################
# 1. Rosa = 34 SNPs (sIL6R)
#data source = taken from supplementary table in paper
#EAF not available
ROSA_INSTRUMENT <- read.table("Rosa_Instrument_34SNPs.txt", sep = "\t", header=TRUE)
ROSA_INSTRUMENT <- rename(ROSA_INSTRUMENT, rsid = SNP_ID, CHR=Chr, effect_allele=EA, other_allele=NEA, beta=Beta, BP = Position)
ROSA_INSTRUMENT$EAF <- NA

# 2. Swerdlow - 3 SNPs (IL6)
#data source = taken from paper and available on Nils Rek OSF
SWERDLOW_INSTRUMENT <- read.table("Swerdlow_Instrument_3SNPs_pvalsincluded.txt", sep = "\t", header=TRUE, colClasses = c("character", "character", "character", "numeric", "numeric", "numeric", "character", "numeric"))
SWERDLOW_INSTRUMENT <- rename(SWERDLOW_INSTRUMENT, rsid=SNP, CHR = chr, EAF = eaf)
sum(SWERDLOW_INSTRUMENT$EAF < 0.01 | SWERDLOW_INSTRUMENT$EAF > 0.99) #0

# 3. Sarwar - 1 SNP (IL6)
#data source = taken from Nils Rek OSF
#EAF not available
SARWAR_INSTRUMENT <- read.table("Sarwar_Instrument_1SNP_pvalsincluded.txt", sep = "\t", header=TRUE)
SARWAR_INSTRUMENT <- rename(SARWAR_INSTRUMENT, rsid=SNP, CHR = chr)
SARWAR_INSTRUMENT$EAF <- NA

# 4. Lam - General cognition ability GWAS 
#snps = 8,990,900
OUTCOME_LAM_GWAS <- fread("Cognition_MTAG_allsnps_mtag_meta.txt")
OUTCOME_LAM_GWAS <- rename(OUTCOME_LAM_GWAS, rsid=SNP)

###################################################
# PART 3: IDENTIFY SNPS AVAILABLE IN OUTCOME GWAS #
###################################################
# ahluwalia, swerdlow, sarwar not have missing SNPs
available_snps_ligthart <- merge(LIGTHART_GWAS,OUTCOME_LAM_GWAS, by="rsid") #74/78 snps
available_snps_han <- merge(HAN_GWAS,OUTCOME_LAM_GWAS, by="rsid") #444/552 snps
available_snps_borges <- merge(BORGES_GWAS, OUTCOME_LAM_GWAS, by="rsid") #73/87 snps
available_snps_kettunen <- merge(KETTUNEN_GWAS,OUTCOME_LAM_GWAS, by="rsid") #9/10 snps
available_snps_rosa <- merge(ROSA_INSTRUMENT,OUTCOME_LAM_GWAS, by="rsid") #27/34 snps
available_snps_ahluwalia <- merge(AHLUWALIA_GWAS,OUTCOME_LAM_GWAS, by="rsid") #3/3 snps
available_snps_swerdlow <- merge(SWERDLOW_INSTRUMENT,OUTCOME_LAM_GWAS, by="rsid") #3/3 snps
available_snps_sarwar <- merge(SARWAR_INSTRUMENT, OUTCOME_LAM_GWAS, by="rsid") #1/1 snp

# make list of all available snps = 634 SNPs (611 are distinct)
available_snps_combined <- data.frame(SNP = c(available_snps_ligthart[,"rsid"], available_snps_han[,"rsid"], available_snps_ahluwalia[,"rsid"], available_snps_borges[,"rsid"], available_snps_kettunen[,"rsid"], available_snps_rosa[,"rsid"], available_snps_swerdlow[,"rsid"], available_snps_sarwar[,"rsid"]))
available_snps_combined <- distinct(available_snps_combined)

# all inflammation snps requested - 768 SNPs (745 are distinct)
all_inflammation_snps <- data.frame(SNP = c(LIGTHART_GWAS[,"rsid"], HAN_GWAS[,"rsid"], AHLUWALIA_GWAS[,"rsid"], BORGES_GWAS[,"rsid"], KETTUNEN_GWAS[,"rsid"], ROSA_INSTRUMENT[,"rsid"], SWERDLOW_INSTRUMENT[,"rsid"], SARWAR_INSTRUMENT[,"rsid"]))
all_inflammation_snps <- distinct(all_inflammation_snps)

####################################################################
#  PART 4: FIND PROXY SNPS FOR THOSE NOT AVAILABLE IN OUTCOME GWAS #
####################################################################
# 134 missing snps
# code provides access to 'LDlink' API using R console. Enables perform batch queries in 1000 Genomes Project using 'LDlink'
missingSNPs <- data.frame(setdiff(all_inflammation_snps$SNP, available_snps_combined$SNP)) #snps requested but not available in Lam GWAS = 134 snps
#LDproxy_batch(missingSNPs, pop="EUR", r2d = "r2", token = "INSERT TOKEN HERE!!", append=TRUE) #saves into txt file in folder (combined_query_snp_list)
proxySNPs <- read.table("combined_query_snp_list_2smr_inflammtocog.txt", header=TRUE, row.names=NULL)#note: renamed output file from LDlink 
proxySNPs_availableOutcome <- merge(proxySNPs, OUTCOME_LAM_GWAS, by.x="RS_Number", by.y="rsid")#only keep proxies available in Lam GWAS
proxySNPs_removemissing <- filter(proxySNPs_availableOutcome, RS_Number != ".") # removes snps with no rsid=none
proxySNPs_LD <- proxySNPs_removemissing %>% filter(R2 > 0.8) #keeps snps R2 > 0.8 = 1980

# list of missing snps
SNPs_missing <- rename(missingSNPs, SNP = setdiff.all_inflammation_snps.SNP..available_snps_combined.SNP.)

#############################################################
#  PART 5: FIND PROXY SNPS IN ORIGINAL GWAS AND ADD TO LIST #
#############################################################
# Finds proxies for GWAS sources, cannot get proxy for Rosa as only instrument.

###########################################################
## LIGTHART CRP GWAS - replace missing snps with proxies ##
###########################################################
# A. identifies proxy snps also available in ligthart then keep those with highest R2
proxySNPs_meetscriteria_ligthart <- merge(proxySNPs_LD, ligthart, by.x="RS_Number", by.y="variant_id")
proxySNPs_meetscriteria_ligthart <- proxySNPs_meetscriteria_ligthart %>% 
  group_by(query_snp) %>%
  slice(c(which.max(R2)))

# B. extract list of missing snps from ligthart = 4 snps not in Lam - 3 proxies
ligthart_missing <- merge(SNPs_missing, LIGTHART_GWAS, by.x="SNP", by.y="rsid") #merges list of missing snps with ligthart gwas snps = obtain ligthart GWAS missing snps
ligthart_replace <- merge(ligthart_missing, proxySNPs_meetscriteria_ligthart, by.x="SNP", by.y="query_snp") #list of missing snps for which proxies available
ligthart_proxies <- data.frame(ligthart_replace$RS_Number) #list of proxy snps

# C.delete missing snps in ligthart gwas = 78 (ligthart GWAS) - 4 (snps not in Lam) = 74 snps
LIGTHART_GWAS_DELETEDMISS <- LIGTHART_GWAS[ ! LIGTHART_GWAS$rsid %in% ligthart_missing$SNP, ] #remove ligthart missing snps
LIGTHART_GWAS_DELETEDMISS$id <- NULL 

# D.find proxies in full sum stats = 3 proxies available
LIGTHART_GWAS_PROXIES <- merge(ligthart_proxies, ligthart, by.x="ligthart_replace.RS_Number", by.y="variant_id")
LIGTHART_GWAS_PROXIES <- rename(LIGTHART_GWAS_PROXIES, rsid = ligthart_replace.RS_Number, EAF = effect_allele_frequency, BP = base_pair_location, CHR = chromosome, pval = p_value)

# E. add proxies to ligthart (74+3=77)
LIGTHART_GWAS_PROXIESADDED <- rbind(LIGTHART_GWAS_DELETEDMISS, LIGTHART_GWAS_PROXIES)
LIGTHART_GWAS_PROXIESADDED$id.exposure <- NA

######################################################
## HAN CRP GWAS - replace missing snps with proxies ##
######################################################
# A. identifies proxy snps also available in han then keep those with highest R2
proxySNPs_meetscriteria_han <- merge(proxySNPs_LD, han, by.x="RS_Number", by.y="SNP")
proxySNPs_meetscriteria_han <- proxySNPs_meetscriteria_han %>% 
  group_by(query_snp) %>%
  slice(c(which.max(R2)))

# B. extract list of missing snps from han = 108 missing - 50 proxies
han_missing <- merge(SNPs_missing, HAN_GWAS, by.x="SNP", by.y="rsid") 
han_replace <- merge(han_missing, proxySNPs_meetscriteria_han, by.x="SNP", by.y="query_snp") 
han_proxies <- data.frame(han_replace$RS_Number) #list of proxy snps

# C. delete missing snps in han gwas = 552 - 108 = 444 snps
HAN_GWAS_DELETEDMISS <- HAN_GWAS[ ! HAN_GWAS$rsid %in% han_missing$SNP, ]
HAN_GWAS_DELETEDMISS$id <- NULL 

# D. find proxies in full sum stats = 50 proxies available
HAN_GWAS_PROXIES <- merge(han_proxies, han, by.x="han_replace.RS_Number", by.y="SNP")
HAN_GWAS_PROXIES <- rename(HAN_GWAS_PROXIES, rsid = han_replace.RS_Number, pval = P, effect_allele = ALLELE1, other_allele = ALLELE0, beta = BETA, EAF = A1FREQ)

# E. add proxies to han (444+50=494)
HAN_GWAS_PROXIESADDED <- rbind(HAN_GWAS_DELETEDMISS, HAN_GWAS_PROXIES)
HAN_GWAS_PROXIESADDED$id.exposure <- NA

###########################################################
## BORGES GlyCA GWAS - replace missing snps with proxies ##
###########################################################
# A. identifies proxy snps also available in borges then keep those with highest R2
proxySNPs_meetscriteria_borges <- merge(proxySNPs_LD, borges, by.x="RS_Number", by.y="variant_id")
proxySNPs_meetscriteria_borges <- proxySNPs_meetscriteria_borges %>% 
  group_by(query_snp) %>%
  slice(c(which.max(R2)))

# B. extract list of missing snps from borges = 14 missing - 9 proxies
borges_missing <- merge(SNPs_missing, BORGES_GWAS, by.x="SNP", by.y="rsid") 
borges_replace <- merge(borges_missing, proxySNPs_meetscriteria_borges, by.x="SNP", by.y="query_snp") 
borges_proxies <- data.frame(borges_replace$RS_Number) #list of proxy snps

# C. delete missing snps in borges gwas = 87 - 14 = 73
BORGES_GWAS_DELETEDMISS <- BORGES_GWAS[ ! BORGES_GWAS$rsid %in% borges_missing$SNP, ]
BORGES_GWAS_DELETEDMISS$id <- NULL 

# D. find proxies in full sum stats = 9 proxies
BORGES_GWAS_PROXIES <- merge(borges_proxies, borges, by.x="borges_replace.RS_Number", by.y="variant_id")
BORGES_GWAS_PROXIES <- rename(BORGES_GWAS_PROXIES, pval = p_value, rsid = borges_replace.RS_Number, CHR = chromosome, BP = base_pair_location, EAF = effect_allele_frequency)

# E. add proxies to borges (73+9=82)
BORGES_GWAS_PROXIESADDED <- rbind(BORGES_GWAS_DELETEDMISS, BORGES_GWAS_PROXIES)
BORGES_GWAS_PROXIESADDED$id.exposure <- NA

#############################################################
## KETTUNEN GlyCA GWAS - replace missing snps with proxies ##
#############################################################
# A. identifies proxy snps also available in kettunen then keep those with highest R2
proxySNPs_meetscriteria_kettunen <- merge(proxySNPs_LD, kettunen, by.x="RS_Number", by.y="variant_id")
proxySNPs_meetscriteria_kettunen <- proxySNPs_meetscriteria_kettunen %>% 
  group_by(query_snp) %>%
  slice(c(which.max(R2)))

# B. extract list of missing snps from kettunen = 1 missing - 1 proxy
kettunen_missing <- merge(SNPs_missing, KETTUNEN_GWAS, by.x="SNP", by.y="rsid") 
kettunen_replace <- merge(kettunen_missing, proxySNPs_meetscriteria_kettunen, by.x="SNP", by.y="query_snp") 
kettunen_proxies <- data.frame(kettunen_replace$RS_Number) #list of proxy snps

# C. delete missing snps in kettunen gwas = 10 - 1 = 9
KETTUNEN_GWAS_DELETEDMISS <- KETTUNEN_GWAS[ ! KETTUNEN_GWAS$rsid %in% kettunen_missing$SNP, ]
KETTUNEN_GWAS_DELETEDMISS$id <- NULL 

# D. find proxies in full sum stats = 1 proxy
KETTUNEN_GWAS_PROXIES <- merge(kettunen_proxies, kettunen, by.x="kettunen_replace.RS_Number", by.y="variant_id")
KETTUNEN_GWAS_PROXIES <- rename(KETTUNEN_GWAS_PROXIES, pval = p_value, rsid = kettunen_replace.RS_Number, CHR = chromosome, BP = base_pair_location, EAF = effect_allele_frequency)

# E. add proxies to kettunen (9+1-10)
KETTUNEN_GWAS_PROXIESADDED <- rbind(KETTUNEN_GWAS_DELETEDMISS, KETTUNEN_GWAS_PROXIES)
KETTUNEN_GWAS_PROXIESADDED$id.exposure <- NA

########################################################################################
#   PART 6: Beta flip - all beta values positive and alleles/eaf to agree with change  #
########################################################################################
# flips all beta values to be positive i.e. abs(i$beta) and alleles and eaf to agree with this change
GWAS <- list(LIGTHART_GWAS_PROXIESADDED,HAN_GWAS_PROXIESADDED,AHLUWALIA_GWAS,BORGES_GWAS_PROXIESADDED,KETTUNEN_GWAS_PROXIESADDED,ROSA_INSTRUMENT,SWERDLOW_INSTRUMENT,SARWAR_INSTRUMENT) # list of each GWAS hits

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

# saves output above in new data frame for each GWAS
# important: ensure order is same as order in list above**
EXP_LIGTHART_GWAS_CORRECTED <- GWAS[[1]]
EXP_HAN_GWAS_CORRECTED <- GWAS[[2]]
EXP_AHLUWALIA_GWAS_CORRECTED <- GWAS[[3]]
EXP_BORGES_GWAS_CORRECTED <- GWAS[[4]]
EXP_KETTUNEN_GWAS_CORRECTED <- GWAS[[5]]
EXP_ROSA_INSTRUMENT_CORRECTED <- GWAS[[6]]
EXP_SWERDLOW_INSTRUMENT_CORRECTED <- GWAS[[7]]
EXP_SARWAR_INSTRUMENT_CORRECTED <- GWAS[[8]]

# Remove uncorrected GWAS
rm(GWAS,LIGTHART_GWAS, HAN_GWAS, AHLUWALIA_GWAS, BORGES_GWAS, KETTUNEN_GWAS, ROSA_INSTRUMENT, SWERDLOW_INSTRUMENT, SARWAR_INSTRUMENT)

########################################################
########### PART 7: Extract cis variants ###############
########################################################

######################
#   CRP: CRP Gene    #
######################
# examine whether the variants are cis acting for gene: CRP
# Based on GRCh37 assembly as this assembly was used by both ligthart and han GWAS.
# 1. based on Gene Cards (GRCh37), CRP is located: chr1:159,682,079-159,684,379
# 2. two columns with the limits of the cis region i.e +-1,000,000 bp
CRP.LOWER <- 159682079-1000000
CRP.UPPER <- 159684379+1000000

# 3. identify which ones are cis
EXP_HAN_GWAS_CORRECTED$CHRPOS<- ifelse(EXP_HAN_GWAS_CORRECTED$CHR==1 & (EXP_HAN_GWAS_CORRECTED$BP>=CRP.LOWER & EXP_HAN_GWAS_CORRECTED$BP<=CRP.UPPER), "Cis", "Trans")
EXP_HAN_CIS <- filter (EXP_HAN_GWAS_CORRECTED, CHRPOS == "Cis")

EXP_LIGTHART_GWAS_CORRECTED$CHRPOS<- ifelse(EXP_LIGTHART_GWAS_CORRECTED$CHR==1 & (EXP_LIGTHART_GWAS_CORRECTED$BP>=CRP.LOWER & EXP_LIGTHART_GWAS_CORRECTED$BP<=CRP.UPPER), "Cis", "Trans")
EXP_LIGTHART_CIS <- filter (EXP_LIGTHART_GWAS_CORRECTED, CHRPOS == "Cis")

#####################
#    IL6:IL6R gene  #
#####################
# Ahluwalia GWAS based on GRCh36 assembly, and so new BP were found for variants based on latest GRCh38 assembly.
# examine whether the variants are cis acting for gene: IL6R
# 1. based on Gene Cards (GRCh38), IL6R is located:chr1:154,405,193-154,469,450
IL6R.LOWER<- 154405193-1000000
IL6R.UPPER<- 154469450+1000000

# 2. create new BP for SNPs (using GRCh38) - extracted manually
ahluwalia_grch38 <- read.table("ahluwalia_3SNPs_GRCH38location.txt", sep = "\t", header=TRUE)
EXP_AHLUWALIA_GWAS_CORRECTED <- merge(EXP_AHLUWALIA_GWAS_CORRECTED, ahluwalia_grch38, by.x="rsid", by.y="rsid")

# 3. identify which are cis
EXP_AHLUWALIA_GWAS_CORRECTED$CHRPOS<- ifelse(EXP_AHLUWALIA_GWAS_CORRECTED$CHR==1 & (EXP_AHLUWALIA_GWAS_CORRECTED$GRCH38>=IL6R.LOWER & EXP_AHLUWALIA_GWAS_CORRECTED$GRCH38<=IL6R.UPPER), "Cis", "Trans")
EXP_AHLUWALIA_CIS <- filter (EXP_AHLUWALIA_GWAS_CORRECTED, CHRPOS == "Cis")

#####################################
## PART 8: Re-format data for 2SMR ##
#####################################
#Minimum required info: SNP, beta, se, effect_allele (other useful: other allele, eaf, phenotype)
#Ligthart - cis
EXP_LIGTHART_CIS$Phenotype <- "CRP (Ligthart_cis)"
EXP_LIGTHART_CIS <- EXP_LIGTHART_CIS[, c("rsid","beta.corrected","standard_error","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_LIGTHART_CIS <- rename(EXP_LIGTHART_CIS,SNP = rsid,beta = beta.corrected,se = standard_error,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Ligthart - GW
EXP_LIGTHART_GWAS_CORRECTED$Phenotype <- "CRP (Ligthart_gw)"
EXP_LIGTHART_GWAS_CORRECTED <- EXP_LIGTHART_GWAS_CORRECTED[, c("rsid","beta.corrected","standard_error","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_LIGTHART_GWAS_CORRECTED <- rename(EXP_LIGTHART_GWAS_CORRECTED,SNP = rsid,beta = beta.corrected,se = standard_error,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Han - cis
EXP_HAN_CIS$Phenotype <- "CRP (Han_cis)"
EXP_HAN_CIS <- EXP_HAN_CIS[, c("rsid","beta.corrected","SE","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_HAN_CIS <- rename(EXP_HAN_CIS,SNP = rsid,beta = beta.corrected,se = SE,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Han - GW
EXP_HAN_GWAS_CORRECTED$Phenotype <- "CRP (Han_gw)"
EXP_HAN_GWAS_CORRECTED <- EXP_HAN_GWAS_CORRECTED[, c("rsid","beta.corrected","SE","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_HAN_GWAS_CORRECTED <- rename(EXP_HAN_GWAS_CORRECTED,SNP = rsid,beta = beta.corrected,se = SE,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Ahluwalia - cis
EXP_AHLUWALIA_CIS$Phenotype <- "IL-6 (Ahluwalia_cis)"
EXP_AHLUWALIA_CIS <- EXP_AHLUWALIA_CIS[, c("rsid","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_AHLUWALIA_CIS <- rename(EXP_AHLUWALIA_CIS,SNP = rsid,beta = beta.corrected,se = se,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Ahluwalia - GW
EXP_AHLUWALIA_GWAS_CORRECTED$Phenotype <- "IL-6 (Ahluwalia_gw)"
EXP_AHLUWALIA_GWAS_CORRECTED <- EXP_AHLUWALIA_GWAS_CORRECTED[, c("rsid","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_AHLUWALIA_GWAS_CORRECTED <- rename(EXP_AHLUWALIA_GWAS_CORRECTED,SNP = rsid,beta = beta.corrected,se = se,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Borges
EXP_BORGES_GWAS_CORRECTED$Phenotype <- "GLYCA (Borges)"
EXP_BORGES_GWAS_CORRECTED <- EXP_BORGES_GWAS_CORRECTED[, c("rsid","beta.corrected","standard_error","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_BORGES_GWAS_CORRECTED <- rename(EXP_BORGES_GWAS_CORRECTED,SNP = rsid,beta = beta.corrected,se = standard_error,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Kettunen
EXP_KETTUNEN_GWAS_CORRECTED$Phenotype <- "GLYCA (kettunen)"
EXP_KETTUNEN_GWAS_CORRECTED <- EXP_KETTUNEN_GWAS_CORRECTED[, c("rsid","beta.corrected","standard_error","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_KETTUNEN_GWAS_CORRECTED <- rename(EXP_KETTUNEN_GWAS_CORRECTED,SNP = rsid,beta = beta.corrected,se = standard_error,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Rosa
EXP_ROSA_INSTRUMENT_CORRECTED$Phenotype <- "sIL6 (Rosa)"
EXP_ROSA_INSTRUMENT_CORRECTED <- EXP_ROSA_INSTRUMENT_CORRECTED[, c("rsid","beta.corrected","SE","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","p.value")]
EXP_ROSA_INSTRUMENT_CORRECTED <- rename(EXP_ROSA_INSTRUMENT_CORRECTED,SNP = rsid,beta = beta.corrected,se = SE,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected,pval = p.value)

#Swerdlow
EXP_SWERDLOW_INSTRUMENT_CORRECTED$Phenotype <- "IL-6 (Swerdlow)"
EXP_SWERDLOW_INSTRUMENT_CORRECTED <- EXP_SWERDLOW_INSTRUMENT_CORRECTED[, c("rsid","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_SWERDLOW_INSTRUMENT_CORRECTED <- rename(EXP_SWERDLOW_INSTRUMENT_CORRECTED,SNP = rsid,beta = beta.corrected,se = se,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Sarwar
EXP_SARWAR_INSTRUMENT_CORRECTED$Phenotype <- "IL-6 (Sarwar)"
EXP_SARWAR_INSTRUMENT_CORRECTED <- EXP_SARWAR_INSTRUMENT_CORRECTED[, c("rsid","beta.corrected","se","effect_allele.corrected","other_allele.corrected","eaf.corrected","Phenotype","pval")]
EXP_SARWAR_INSTRUMENT_CORRECTED <- rename(EXP_SARWAR_INSTRUMENT_CORRECTED,SNP = rsid,beta = beta.corrected,se = se,effect_allele = effect_allele.corrected,other_allele = other_allele.corrected,eaf = eaf.corrected)

#Outcome - Lam
OUTCOME_LAM_GWAS$Phenotype <- "GeneralCognitiveAbility (Lam)"
OUTCOME_LAM_GWAS <- OUTCOME_LAM_GWAS[, c("rsid","mtag_beta","mtag_se","A1","A2","meta_freq","Phenotype","mtag_pval")]
OUTCOME_LAM_GWAS <- rename(OUTCOME_LAM_GWAS,SNP = rsid,beta = mtag_beta,se = mtag_se,effect_allele = A1,other_allele = A2,eaf = meta_freq,pval = mtag_pval)

######################################
## Save outputs to be used for 2SMR ##
######################################
# write.table(EXP_LIGTHART_CIS, "EXP_LIGTHART_CIS_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_LIGTHART_GWAS_CORRECTED, "EXP_LIGTHART_GWAS_CORRECTED_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_HAN_CIS, "EXP_HAN_CIS_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_HAN_GWAS_CORRECTED, "EXP_HAN_GWAS_CORRECTED_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_AHLUWALIA_CIS, "EXP_AHLUWALIA_CIS_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_AHLUWALIA_GWAS_CORRECTED, "EXP_AHLUWALIA_GWAS_CORRECTED_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_BORGES_GWAS_CORRECTED, "EXP_BORGES_GWAS_CORRECTED_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_KETTUNEN_GWAS_CORRECTED, "EXP_KETTUNEN_GWAS_CORRECTED_FOR2SMR_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_ROSA_INSTRUMENT_CORRECTED, "EXP_ROSA_INSTRUMENT_CORRECTED_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_SWERDLOW_INSTRUMENT_CORRECTED, "EXP_SWERDLOW_INSTRUMENT_CORRECTED_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(EXP_SARWAR_INSTRUMENT_CORRECTED, "EXP_SARWAR_INSTRUMENT_CORRECTED_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
# write.table(OUTCOME_LAM_GWAS, "OUTCOME_LAM_GWAS_FOR2SMR.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

#####################################################
# Save outputs to be used for 2SMR - pvals included #
#####################################################
write.table(EXP_LIGTHART_CIS, "EXP_LIGTHART_CIS_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_LIGTHART_GWAS_CORRECTED, "EXP_LIGTHART_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_HAN_CIS, "EXP_HAN_CIS_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_HAN_GWAS_CORRECTED, "EXP_HAN_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_AHLUWALIA_CIS, "EXP_AHLUWALIA_CIS_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_AHLUWALIA_GWAS_CORRECTED, "EXP_AHLUWALIA_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_BORGES_GWAS_CORRECTED, "EXP_BORGES_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_KETTUNEN_GWAS_CORRECTED, "EXP_KETTUNEN_GWAS_CORRECTED_FOR2SMR_PROXIESADDED_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_ROSA_INSTRUMENT_CORRECTED, "EXP_ROSA_INSTRUMENT_CORRECTED_FOR2SMR_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_SWERDLOW_INSTRUMENT_CORRECTED, "EXP_SWERDLOW_INSTRUMENT_CORRECTED_FOR2SMR_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(EXP_SARWAR_INSTRUMENT_CORRECTED, "EXP_SARWAR_INSTRUMENT_CORRECTED_FOR2SMR_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(OUTCOME_LAM_GWAS, "OUTCOME_LAM_GWAS_FOR2SMR_pval.txt", sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
