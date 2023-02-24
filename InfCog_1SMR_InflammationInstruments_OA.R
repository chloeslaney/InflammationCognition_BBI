# Project: inflammation and cognition
# Script extracts snps to create PRS for different inflammatory markers
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

###############################
# C-reactive protein = 2 GWAS #
###############################
#Code extracts snps for PRS (uses C+T approach) and checks all SNPs EAF > 0.01

# 1. Ligthart CRP GWAS - full summary statistics = 78 SNPs
#Data source = VCF file from IEU database converted to txt using bcftools (separate script for this)
#ligthart_IEU <- tophits(id="ieu-b-35", pval=5e-08, r2=0.01, kb=1000) #78
ligthart <- read.table("GWAS_ligthart_crp.txt", header=TRUE) #reads data
ligthart_pval <- ligthart %>% filter(p_value < 5e-08) #filter to include SNPs p < 5*10-8
ligthart_pval <- rename(ligthart_pval, pval = p_value, rsid = variant_id, CHR=chromosome, BP=base_pair_location, EAF = effect_allele_frequency) #rename for ld_clump
LIGTHART_GWAS <- ld_clump(ligthart_pval,clump_kb=1000, clump_r2=0.01) #clump SNPs
#Quality checks
#sum(LIGTHART_GWAS$EAF < 0.01 | LIGTHART_GWAS$EAF > 0.99) #0
#sum(duplicated(LIGTHART_GWAS$rsid)) #0 duplicate rsids

# 2. Han CRP GWAS - full summary statistics (EA=ALLELE1) = 552 SNPs
#data source = requested from authors and sent link (used curl on linux to download to txt file - checked wget gives same snp #)
han <- read.table("GWAS_han_crp.txt", sep = ",", header=TRUE)
han_pval <- han %>% filter (P < 5e-08)
han_pval <- rename(han_pval, pval = P, rsid = SNP, effect_allele = ALLELE1, other_allele = ALLELE0, beta = BETA, EAF = A1FREQ)
HAN_GWAS <- ld_clump(han_pval,clump_kb=1000, clump_r2=0.01)
#Quality checks
#sum(HAN_GWAS$EAF < 0.01 | HAN_GWAS$EAF > 0.99) #0
#sum(duplicated(HAN_GWAS$rsid)) #0 duplicate rsids

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
#Quality check
#sum(duplicated(AHLUWALIA_GWAS$rsid)) #0 duplicate rsids

#################################
# Glycoprotein Acetyls = 2 GWAS #
#################################
# 1. Borges - full summary statistics = 88 SNPs (87 if remove SNP EAF < 0.01)
#data source = VCF file from IEU database converted to txt using bcftools (separate script for this)
#borges_IEU <- tophits(id="met-d-GlycA", r2=0.01, kb=1000) #88
borges <- read.table("GWAS_borges_glyca.txt", header=TRUE) 
borges_pval <- borges %>% filter(p_value < 5e-08) 
borges_pval <- rename(borges_pval, pval = p_value, rsid = variant_id, CHR = chromosome, BP = base_pair_location, EAF = effect_allele_frequency)
BORGES_GWAS <- ld_clump(borges_pval,clump_kb=1000, clump_r2=0.01)
#Quality checks
sum(BORGES_GWAS$EAF < 0.01 | BORGES_GWAS$EAF > 0.99) #1 SNP! 
BORGES_GWAS <- BORGES_GWAS %>% filter(EAF > 0.01) # only keep SNPs with EAF > 0.01
#sum(duplicated(HAN_GWAS$rsid)) #0 duplicate rsids

# 2. Kettunen - full summary statistics = 10 SNPs
#data source = VCF file from IEU database converted to txt using bcftools (separate script for this)
#kettunen_IEU <- tophits(id="met-c-863", r2 = 0.01, kb = 1000) #10
kettunen <- read.table("GWAS_kettunen_glyca.txt", header=TRUE)
kettunen_pval <- kettunen %>% filter(p_value < 5e-08)
kettunen_pval <- rename(kettunen_pval, pval = p_value, rsid = variant_id, CHR = chromosome, BP = base_pair_location, EAF = effect_allele_frequency)
KETTUNEN_GWAS <- ld_clump(kettunen_pval,clump_kb=1000, clump_r2=0.01)
#Quality checks
sum(KETTUNEN_GWAS$EAF < 0.01 | KETTUNEN_GWAS$EAF > 0.99) #0
#sum(duplicated(KETTUNEN_GWAS$rsid)) #0 duplicate rsids

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
#data source = taken from paper and checked against Nils Rek OSF
SWERDLOW_INSTRUMENT <- read.table("Swerdlow_Instrument_3SNPs.txt", sep = "\t", header=TRUE, colClasses = c("character", "character", "character", "numeric", "numeric", "numeric", "character", "numeric"))
SWERDLOW_INSTRUMENT <- rename(SWERDLOW_INSTRUMENT, rsid=SNP, CHR = chr, EAF = eaf)
sum(SWERDLOW_INSTRUMENT$EAF < 0.01 | SWERDLOW_INSTRUMENT$EAF > 0.99) #0

# 3. Sarwar - 1 SNP (IL6)
#data source = taken from Nils Rek OSF
#EAF not available
SARWAR_INSTRUMENT <- read.table("Sarwar_Instrument_1SNP.txt", sep = "\t", header=TRUE)
SARWAR_INSTRUMENT <- rename(SARWAR_INSTRUMENT, rsid=SNP, CHR = chr)
SARWAR_INSTRUMENT$EAF <- NA

###########################################################
#              Finding proxy SNPs                         #
###########################################################
#Code provides access to 'LDlink' API using R console. Enables perform batch queries in 1000 Genomes Project using 'LDlink'
#Load available SNPs in ALSPAC (must first check which SNPs from those extracted below are available in ALSPAC)
#availableSNPs <- read.table("available_snps.txt", header=FALSE) 

#Finds SNP in highest LD with missing SNP - identified proxies for 32/54 missing SNPs
#missingSNPs <- data.frame(setdiff(ALL_GWAS_SNP_DUPLICATESREMOVED$rsid, availableSNPs$V1)) #extracts snps that requested from all GWAS but not available in ALSPAC
#LDproxy_batch(missingSNPs, pop="EUR", r2d = "r2", token = "INSERT TOKEN HERE", append=TRUE) #using missing snps extracts proxies and saves into txt file (combined_query_snp_list)
proxySNPs <- read.table("combined_query_snp_list_1smr.txt", header=TRUE, row.names=NULL)
proxySNPs_highestld <- subset(proxySNPs, row.names == "2" | row.names == "3" | row.names == "4") #keep top three snps in highest ld
proxySNPs_removemissing <- filter (proxySNPs_highestld, RS_Number != ".") # removes snps with no rsid
proxySNPs_LD <- proxySNPs_removemissing %>% # groups by rsid and keeps ld snp with lowest row.names (i.e., higher R2)
    group_by(query_snp) %>%
    slice(c(which.min(row.names)))
proxySNPs_meetscriteria <- proxySNPs_LD %>% filter(R2 > 0.8) #keep only ld snps with r2 > 0.8

#New SNP list = add proxy snps to available snps list
#availablesnps_andproxies <- data.frame(list.append(availableSNPs$V1,proxySNPs_meetscriteria$RS_Number))

###########################################################
## LIGTHART CRP GWAS - replace missing snps with proxies ##
###########################################################
# A.extract list of missing snps from ligthart = 2 snps not in alspac 
ligthart_missing <- merge(proxySNPs_LD, LIGTHART_GWAS, by.x="query_snp", by.y="rsid") #merges list of missing snps with ligthart gwas snps = obtain ligthart GWAS missing snps in ALSPAC
ligthart_replace <- merge(proxySNPs_meetscriteria, LIGTHART_GWAS, by.x="query_snp", by.y="rsid") #list of missing snps for which proxies available
ligthart_proxies <- data.frame(ligthart_replace$RS_Number) #list of proxy snps

# B.delete missing snps in ligthart gwas = 78 (ligthart GWAS) - 2 (snps not in alspac) = 76 snps
LIGTHART_GWAS_DELETEDMISS <- LIGTHART_GWAS[ ! LIGTHART_GWAS$rsid %in% ligthart_missing$query_snp, ] #remove ligthart missing snps
LIGTHART_GWAS_DELETEDMISS$id <- NULL

# C.find proxies in full sum stats (only 1 of 2 snps available) 
LIGTHART_GWAS_PROXIES <- merge(ligthart_proxies, ligthart, by.x="ligthart_replace.RS_Number", by.y="variant_id")
LIGTHART_GWAS_PROXIES <- rename(LIGTHART_GWAS_PROXIES, pval = p_value, rsid = ligthart_replace.RS_Number, CHR=chromosome, BP=base_pair_location, EAF=effect_allele_frequency)

# D.add proxies to ligthart (76+1=77)
LIGTHART_GWAS_PROXIESADDED <- rbind(LIGTHART_GWAS_DELETEDMISS, LIGTHART_GWAS_PROXIES)

######################################################
## HAN CRP GWAS - replace missing snps with proxies ##
######################################################
# Note: must run proxy section below before running this code.
# A.extract list of missing snps from han = 43 snps not in alspac, only 24 with proxies available (r2 >0.8)
han_missing <- merge(proxySNPs_LD, HAN_GWAS, by.x="query_snp", by.y="rsid")
han_replace <- merge(proxySNPs_meetscriteria, HAN_GWAS, by.x="query_snp", by.y="rsid")
han_proxies <- data.frame(han_replace$RS_Number)

# B.delete missing snps in han gwas = 552 (han GWAS) - 43 (snps not in alspac) = 509 snps
HAN_GWAS_DELETEDMISS <- HAN_GWAS[ ! HAN_GWAS$rsid %in% han_missing$query_snp, ]
HAN_GWAS_DELETEDMISS$id <- NULL

# C.find proxies in full sum stats (only 20 of 24 snps available) 
HAN_GWAS_PROXIES <- merge(han_proxies, han, by.x="han_replace.RS_Number", by.y="SNP")
HAN_GWAS_PROXIES <- rename(HAN_GWAS_PROXIES, pval = P, rsid = han_replace.RS_Number, effect_allele = ALLELE1, other_allele = ALLELE0, beta = BETA, EAF = A1FREQ)

# D.add proxies to han (509+20=529)
HAN_GWAS_PROXIESADDED <- rbind(HAN_GWAS_DELETEDMISS, HAN_GWAS_PROXIES)

###########################################################
## BORGES GlycA GWAS - replace missing snps with proxies ##
###########################################################
# Note: SNP with EAF < 0.01 is excluded here
# A.extract list of missing snps from borges = 9 snps not in alspac, only 6 with proxies available (r2 >0.8)
borges_missing <- merge(proxySNPs_LD, BORGES_GWAS, by.x="query_snp", by.y="rsid")
borges_replace <- merge(proxySNPs_meetscriteria, BORGES_GWAS, by.x="query_snp", by.y="rsid")
borges_proxies <- data.frame(borges_replace$RS_Number)

# B.delete missing snps in borges gwas = 87 (borges GWAS) - 9 (snps not in alspac) = 78 snps
BORGES_GWAS_DELETEDMISS <- BORGES_GWAS[ ! BORGES_GWAS$rsid %in% borges_missing$query_snp, ]
BORGES_GWAS_DELETEDMISS$id <- NULL

# C.find proxies in full sum stats (6 of 6 snps available) 
BORGES_GWAS_PROXIES <- merge(borges_proxies, borges, by.x="borges_replace.RS_Number", by.y="variant_id")
BORGES_GWAS_PROXIES <- rename(BORGES_GWAS_PROXIES, pval = p_value, rsid = borges_replace.RS_Number, BP = base_pair_location, CHR = chromosome, EAF = effect_allele_frequency)

# D.add proxies to borges (78+6=84)
BORGES_GWAS_PROXIESADDED <- rbind(BORGES_GWAS_DELETEDMISS, BORGES_GWAS_PROXIES)

#########################################################
#               Beta flip for PRS                       #
#########################################################
#Code adapted from Christina Dardani - flips all beta values to be positive i.e. abs(i$beta) and alleles to agree with this change
#Not flip EAF as this is not kept/used for 1SMR
GWAS <- list(LIGTHART_GWAS,HAN_GWAS,AHLUWALIA_GWAS,BORGES_GWAS,KETTUNEN_GWAS,ROSA_INSTRUMENT,SWERDLOW_INSTRUMENT,SARWAR_INSTRUMENT,HAN_GWAS_PROXIESADDED,LIGTHART_GWAS_PROXIESADDED,BORGES_GWAS_PROXIESADDED) # list of each GWAS hits

# applies beta flipping to each GWAS hits listed above
for(i in 1:length(GWAS)){
    
    GWAS[[i]]$effect_allele.corrected<- NA
    GWAS[[i]]$effect_allele.corrected<- ifelse(GWAS[[i]]$beta<0, GWAS[[i]]$other_allele, GWAS[[i]]$effect_allele)
    
    GWAS[[i]]$other_allele.corrected<- NA
    GWAS[[i]]$other_allele.corrected<- ifelse(GWAS[[i]]$effect_allele.corrected==GWAS[[i]]$effect_allele, GWAS[[i]]$other_allele, GWAS[[i]]$effect_allele)
    
    GWAS[[i]]$beta.corrected<- NA
    GWAS[[i]]$beta.corrected<- abs(GWAS[[i]]$beta)
}

#Saves output above in new data frame for each GWAS
#Important: ensure order is same as order in list above!
LIGTHART_GWAS_CORRECTED <- GWAS[[1]]
HAN_GWAS_CORRECTED <- GWAS[[2]]
AHLUWALIA_GWAS_CORRECTED <- GWAS[[3]]
BORGES_GWAS_CORRECTED <- GWAS[[4]]
KETTUNEN_GWAS_CORRECTED <- GWAS[[5]]
ROSA_INSTRUMENT_CORRECTED <- GWAS[[6]]
SWERDLOW_INSTRUMENT_CORRECTED <- GWAS[[7]]
SARWAR_INSTRUMENT_CORRECTED <- GWAS[[8]]
HAN_GWAS_PROXIESADDED_CORRECTED <- GWAS[[9]]
LIGTHART_GWAS_PROXIESADDED_CORRECTED <- GWAS[[10]]
BORGES_GWAS_PROXIESADDED_CORRECTED <- GWAS[[11]]

# Remove uncorrected GWAS #
rm(GWAS,LIGTHART_GWAS, HAN_GWAS, AHLUWALIA_GWAS, BORGES_GWAS, KETTUNEN_GWAS, ROSA_INSTRUMENT, SWERDLOW_INSTRUMENT, SARWAR_INSTRUMENT, HAN_GWAS_PROXIESADDED, LIGTHART_GWAS_PROXIESADDED)

#####################################################
########### Extract cis/trans variants ##############
#####################################################

#####################
#   CRP: CRP Gene   #
#####################
#Examine whether the variants are cis acting for gene: CRP
#Based on GRCh37 assembly as this assembly was used by both ligthart and han GWAS.
# 1. Based on Gene Cards (GRCh37), CRP is located: chr1:159,682,079-159,684,379
# 2. Two columns with the limits of the cis region i.e +-1,000,000 bp
CRP.LOWER <- 159682079-1000000
CRP.UPPER <- 159684379+1000000

# 3. Identify which ones are cis/trans
HAN_GWAS_CORRECTED$CHRPOS<- ifelse(HAN_GWAS_CORRECTED$CHR==1 & (HAN_GWAS_CORRECTED$BP>=CRP.LOWER & HAN_GWAS_CORRECTED$BP<=CRP.UPPER), "Cis", "Trans")
HAN_CIS <- filter (HAN_GWAS_CORRECTED, CHRPOS == "Cis")
HAN_TRANS <- filter (HAN_GWAS_CORRECTED, CHRPOS == "Trans")

LIGTHART_GWAS_CORRECTED$CHRPOS<- ifelse(LIGTHART_GWAS_CORRECTED$CHR==1 & (LIGTHART_GWAS_CORRECTED$BP>=CRP.LOWER & LIGTHART_GWAS_CORRECTED$BP<=CRP.UPPER), "Cis", "Trans")
LIGTHART_CIS <- filter (LIGTHART_GWAS_CORRECTED, CHRPOS == "Cis")
LIGTHART_TRANS <- filter (LIGTHART_GWAS_CORRECTED, CHRPOS == "Trans")

HAN_GWAS_PROXIESADDED_CORRECTED$CHRPOS<- ifelse(HAN_GWAS_PROXIESADDED_CORRECTED$CHR==1 & (HAN_GWAS_PROXIESADDED_CORRECTED$BP>=CRP.LOWER & HAN_GWAS_PROXIESADDED_CORRECTED$BP<=CRP.UPPER), "Cis", "Trans")
HAN_PROXIESADDED_CIS <- filter (HAN_GWAS_PROXIESADDED_CORRECTED, CHRPOS == "Cis")
HAN_PROXIESADDED_TRANS <- filter (HAN_GWAS_PROXIESADDED_CORRECTED, CHRPOS == "Trans")

LIGTHART_GWAS_PROXIESADDED_CORRECTED$CHRPOS<- ifelse(LIGTHART_GWAS_PROXIESADDED_CORRECTED$CHR==1 & (LIGTHART_GWAS_PROXIESADDED_CORRECTED$BP>=CRP.LOWER & LIGTHART_GWAS_PROXIESADDED_CORRECTED$BP<=CRP.UPPER), "Cis", "Trans")
LIGTHART_PROXIESADDED_CIS <- filter (LIGTHART_GWAS_PROXIESADDED_CORRECTED, CHRPOS == "Cis")
LIGTHART_PROXIESADDED_TRANS <- filter (LIGTHART_GWAS_PROXIESADDED_CORRECTED, CHRPOS == "Trans")

#####################
#   IL6:IL6R gene   #
#####################
#Ahluwalia GWAS based on GRCh36 assembly, and so new BP were found for variants based on latest GRCh38 assembly.
#Examine whether the variants are cis acting for gene: IL6R
# 1. Based on Gene Cards (GRCh38), IL6R is located:chr1:154,405,193-154,469,450
IL6R.LOWER<- 154405193-1000000
IL6R.UPPER<- 154469450+1000000

# 2. Get new BP for SNPs (using GRCh38) - extracted manually
ahluwalia_grch38 <- read.table("ahluwalia_3SNPs_GRCH38location.txt", sep = "\t", header=TRUE)
AHLUWALIA_GWAS_CORRECTED <- merge(AHLUWALIA_GWAS_CORRECTED, ahluwalia_grch38, by.x="rsid", by.y="rsid")

# 3. Identify which are cis/trans variants
AHLUWALIA_GWAS_CORRECTED$CHRPOS<- ifelse(AHLUWALIA_GWAS_CORRECTED$CHR==1 & (AHLUWALIA_GWAS_CORRECTED$GRCH38>=IL6R.LOWER & AHLUWALIA_GWAS_CORRECTED$GRCH38<=IL6R.UPPER), "Cis", "Trans")
AHLUWALIA_CIS <- filter (AHLUWALIA_GWAS_CORRECTED, CHRPOS == "Cis")
AHLUWALIA_TRANS <- filter (AHLUWALIA_GWAS_CORRECTED, CHRPOS == "Trans")

##########################################################
#              Save GWAS as txt files for HPC            #
##########################################################
#To create polygenic risk scores, need: (1) list all "sig" variants from all GWAS (2) variant, effect allele, beta (each individual GWAS)
############################
# A. Full list of all SNPs #
############################
#Inflammation and cognition SNPs (with proxies added) [USED]
COG_SNPs <- read.table("mahedy_cognition_SNPs_duplicatesremoved.txt", sep = "\t", header=FALSE)
INFLAMMATION_COGNITION_SNPs_PROXIESADDED <- data.frame(rsid = c(COG_SNPs[,"V1"], ALL_SNP_PROXIESADDED_DUPLICATESREMOVED[,"rsid"]))
INFLAMMATION_COGNITION_SNPs_PROXIESADDED_DUPLICATESREMOVED <- distinct(INFLAMMATION_COGNITION_SNPs_PROXIESADDED)
write.table(INFLAMMATION_COGNITION_SNPs_PROXIESADDED_DUPLICATESREMOVED, "Inflammation_Cognition_SNPs_duplicatesremoved_eafchecked_proxiesadded.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#######################################################
# B. List of SNPs for each GWAS/instrument separately #
#######################################################
#han with proxies added - cis and genome wide
HAN_PROXIESADDED_CIS_OUTPUT <- HAN_PROXIESADDED_CIS[, c('rsid', 'effect_allele.corrected', 'beta.corrected')]
write.table(HAN_PROXIESADDED_CIS_OUTPUT, "hancis_proxiesadded_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

HAN_PROXIESADDED_OUTPUT <- HAN_GWAS_PROXIESADDED_CORRECTED[, c('rsid', 'effect_allele.corrected', 'beta.corrected')]
write.table(HAN_PROXIESADDED_OUTPUT, "han_genomewide_proxiesadded_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#ligthart with proxy added
LIGTHART_PROXIESADDED_OUTPUT <- LIGTHART_GWAS_PROXIESADDED_CORRECTED[, c('rsid', 'effect_allele.corrected', 'beta.corrected')]
write.table(LIGTHART_PROXIESADDED_OUTPUT, "ligthart_genomewide_proxiesadded_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#borges with proxies added
BORGES_PROXIESADDED_OUTPUT <- BORGES_GWAS_PROXIESADDED_CORRECTED[, c('rsid', 'effect_allele.corrected', 'beta.corrected')]
write.table(BORGES_PROXIESADDED_OUTPUT, "borges_eafchecked_proxiesadded_for_plink.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

##############################################################################################
###### For strand checking (different script) - list SNPs from all GWAS with SNP, EA, OA #####
##############################################################################################
#No proxies
AHLUWALIA_SNPs <- AHLUWALIA_GWAS_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]
KETTUNEN_SNPs <- KETTUNEN_GWAS_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]
ROSA_SNPs <- ROSA_INSTRUMENT_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]
SWERDLOW_SNPs <- SWERDLOW_INSTRUMENT_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]
SARWAR_SNPs <- SARWAR_INSTRUMENT_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]

#Proxies added
LIGTHART_SNPs_PROXIESADDED <- LIGTHART_GWAS_PROXIESADDED_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]
HAN_SNPs_PROXIESADDED <- HAN_GWAS_PROXIESADDED_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]
BORGES_SNPs_PROXIESADDED <- BORGES_GWAS_PROXIESADDED_CORRECTED[,c("rsid", "effect_allele.corrected", "other_allele.corrected")]
CHECKSTRANDS_INFLAMMATION_PROXIESADDED <- rbind(LIGTHART_SNPs_PROXIESADDED,HAN_SNPs_PROXIESADDED,AHLUWALIA_SNPs,BORGES_SNPs_PROXIESADDED,KETTUNEN_SNPs,ROSA_SNPs,SWERDLOW_SNPs,SARWAR_SNPs)
write.table(CHECKSTRANDS_INFLAMMATION_PROXIESADDED, "CHECKSTRANDS_INFLAMMATION_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

#Merges list of inflammation snps (proxies added) & cognition snps
CHECKSTRANDS_COGNITION <- read.table("CHECKSTRANDS_COGNITION.txt", sep = "\t", header=FALSE)
CHECKSTRANDS_COGNITION <- rename(CHECKSTRANDS_COGNITION, rsid=V1, effect_allele.corrected=V2, other_allele.corrected=V3)
CHECKSTRANDS_INFLAMMATION_COGNITION_PROXIESADDED_SNPs <- rbind(CHECKSTRANDS_INFLAMMATION_PROXIESADDED, CHECKSTRANDS_COGNITION)
write.table(CHECKSTRANDS_INFLAMMATION_COGNITION_PROXIESADDED_SNPs, "CHECKSTRANDS_INFLAMMATION_COGNITION_SNPs_PROXIESADDED.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)#used in checkstrands R file

###################################
##### Ambiguous snps removed #####
##################################
ambiguous_snps <- read.table("alspac_ambiguous_snps_88distinct_proxiesadded.txt", sep = "\t", header=FALSE)

#Remove ambiguous snps from full snp list
INFLAM_COG_SNPs_PROXIESADD_DUPLICATESANDAMBIG_REMOVED <- data.frame(INFLAMMATION_COGNITION_SNPs_PROXIESADDED_DUPLICATESREMOVED[ ! INFLAMMATION_COGNITION_SNPs_PROXIESADDED_DUPLICATESREMOVED$rsid %in% ambiguous_snps$V1, ]) 
#write.table(INFLAM_COG_SNPs_PROXIESADD_DUPLICATESANDAMBIG_REMOVED, "INFLAMMATION_COGNITION_SNPs_PROXIESADDED_AMBIGUOUSREMOVED.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

# ligthart - cis and genomewide
LIGTHART_PROXIESADDED_OUTPUT_AMBIG_REMOVED <- LIGTHART_PROXIESADDED_OUTPUT[ ! LIGTHART_PROXIESADDED_OUTPUT$rsid %in% ambiguous_snps$V1, ]
LIGTHART_CIS_OUTPUT_AMBIG_REMOVED <- LIGTHART_CIS_OUTPUT[ ! LIGTHART_CIS_OUTPUT$rsid %in% ambiguous_snps$V1, ]

# han - cis and genomewide
HAN_PROXIESADDED_OUTPUT_AMBIG_REMOVED <- HAN_PROXIESADDED_OUTPUT[ ! HAN_PROXIESADDED_OUTPUT$rsid %in% ambiguous_snps$V1, ]
HAN_PROXIESADDED_CIS_OUTPUT_AMBIG_REMOVED <- HAN_PROXIESADDED_CIS_OUTPUT[ ! HAN_PROXIESADDED_CIS_OUTPUT$rsid %in% ambiguous_snps$V1, ]

# ahluwalia - cis and genomewide
AHLUWALIA_OUTPUT_AMBIG_REMOVED <- AHLUWALIA_OUTPUT[ ! AHLUWALIA_OUTPUT$rsid %in% ambiguous_snps$V1, ]
AHLUWALIA_CIS_OUTPUT_AMBIG_REMOVED <- AHLUWALIA_CIS_OUTPUT[ ! AHLUWALIA_CIS_OUTPUT$rsid %in% ambiguous_snps$V1, ]

# borges
BORGES_PROXIESADDED_OUTPUT_AMBIG_REMOVED <- BORGES_PROXIESADDED_OUTPUT[ ! BORGES_PROXIESADDED_OUTPUT$rsid %in% ambiguous_snps$V1, ]

# kettunen
KETTUNEN_OUTPUT_AMBIG_REMOVED <- KETTUNEN_OUTPUT[ ! KETTUNEN_OUTPUT$rsid %in% ambiguous_snps$V1, ]

# rosa
ROSA_OUTPUT_AMBIG_REMOVED <- ROSA_OUTPUT[ ! ROSA_OUTPUT$rsid %in% ambiguous_snps$V1, ]

# swerdlow
SWERDLOW_OUTPUT_AMBIG_REMOVED <- SWERDLOW_OUTPUT[ ! SWERDLOW_OUTPUT$rsid %in% ambiguous_snps$V1, ]

# sarwar
SARWAR_OUTPUT_AMBIG_REMOVED <- SARWAR_OUTPUT[ ! SARWAR_OUTPUT$rsid %in% ambiguous_snps$V1, ]

#Save above files
write.table(LIGTHART_PROXIESADDED_OUTPUT_AMBIG_REMOVED, "ligthart_genomewide_proxiesadded_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(LIGTHART_CIS_OUTPUT_AMBIG_REMOVED, "ligthartcis_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(HAN_PROXIESADDED_OUTPUT_AMBIG_REMOVED, "han_genomewide_proxiesadded_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(HAN_PROXIESADDED_CIS_OUTPUT_AMBIG_REMOVED, "hancis_proxiesadded_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(AHLUWALIA_OUTPUT_AMBIG_REMOVED, "ahluwalia_genomewide_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(AHLUWALIA_CIS_OUTPUT_AMBIG_REMOVED, "ahluwaliacis_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(BORGES_PROXIESADDED_OUTPUT_AMBIG_REMOVED, "borges_eafchecked_proxiesadded_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(KETTUNEN_OUTPUT_AMBIG_REMOVED, "kettunen_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(ROSA_OUTPUT_AMBIG_REMOVED, "rosa_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(SWERDLOW_OUTPUT_AMBIG_REMOVED, "swerdlow_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(SARWAR_OUTPUT_AMBIG_REMOVED, "sarwar_ambigremoved_for_plink.txt", sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

###############################################################################
##### Quality Check: Check GWAS have no duplicate rsids within instrument #####
###############################################################################
#Ligthart = 78 - no duplicates
LIGTHART_GWAS_DUPLICATECHECK <- data.frame(rsid=c(LIGTHART_GWAS_CORRECTED[,"rsid"]))
LIGTHART_GWAS_DUPLICATECHECK <-distinct(LIGTHART_GWAS_DUPLICATECHECK)
#Han = 552 - no duplicates
HAN_GWAS_DUPLICATECHECK <- data.frame(rsid=c(HAN_GWAS_CORRECTED[,"rsid"]))
HAN_GWAS_DUPLICATECHECK <-distinct(HAN_GWAS_DUPLICATECHECK)
#Ahluwalia = 3 - no duplicates
AHLUWALIA_GWAS_DUPLICATECHECK <- data.frame(rsid=c(AHLUWALIA_GWAS_CORRECTED[,"rsid"]))
AHLUWALIA_GWAS_DUPLICATECHECK <-distinct(AHLUWALIA_GWAS_DUPLICATECHECK)
#Borges = 87 - no duplicates
BORGES_GWAS_DUPLICATECHECK <- data.frame(rsid=c(BORGES_GWAS_CORRECTED[,"rsid"]))
BORGES_GWAS_DUPLICATECHECK <-distinct(BORGES_GWAS_DUPLICATECHECK)
#Kettunen = 10 - no duplicates
KETTUNEN_GWAS_DUPLICATECHECK <- data.frame(rsid=c(KETTUNEN_GWAS_CORRECTED[,"rsid"]))
KETTUNEN_GWAS_DUPLICATECHECK <-distinct(KETTUNEN_GWAS_DUPLICATECHECK)