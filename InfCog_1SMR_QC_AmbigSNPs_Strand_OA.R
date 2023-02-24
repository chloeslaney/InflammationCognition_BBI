#Project: inflammation and cognition
#Script used for strand checking and palindromic snps between base GWAS and target data (ALSPAC)
#Chloe Slaney

#####################
# PRELIMINARY STEPS #
#####################
library(tidyverse)
library(readr)
library(dplyr)
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

#########################################################
###### Load in base and target SNPs (both alleles) ######
#########################################################
# read in target data - proxies added
ALSPAC_snps <- read.table("all_child.bim") #snps requested and available in alspac = 721 distinct (all inflammation and cognition snps with proxies added)
ALSPAC_snps <- ALSPAC_snps[,c("V2", "V5", "V6")] #keep these columns
ALSPAC_snps <- rename(ALSPAC_snps, rsid=V2, allele1_target=V5, allele2_target=V6)
# read in base data - proxies added
basedata_snps_proxiesadded <- read.table("CHECKSTRANDS_INFLAMMATION_COGNITION_SNPs_PROXIESADDED.txt", sep = "\t", header=FALSE) #all snps = 756
basedata_snps_proxiesadded <- rename(basedata_snps_proxiesadded, rsid=V1, allele1_base=V2, allele2_base=V3)
basedata_rsids <- data.frame(rsid = c(basedata_snps_proxiesadded[,"rsid"]))#for checking
basedata_rsids_distinct <- distinct(basedata_rsids)#for checking = 733 distinct

###############################
# merge base and target data #
##############################
# some duplicates as some snps in multiple GWAS but may have different alleles reported so must check both
merged_file <- merge(ALSPAC_snps, basedata_snps_proxiesadded, by="rsid")

#########################
# no identical alleles #
########################
# checks base nor target data have identical snps (e.g., CC or TT)
sum(merged_file$allele1_base==merged_file$allele2_base) #0
sum(merged_file$allele1_target==merged_file$allele2_target) #0

###########################
# check strand mismatches #
###########################
# checks alleles reported in base and target match - if not, may be strand differences. No mismatches here.
merged_file$match_check <- NA
merged_file$match_check <- ifelse((merged_file$allele1_base==merged_file$allele1_target & merged_file$allele2_base==merged_file$allele2_target) | (merged_file$allele1_base==merged_file$allele2_target & merged_file$allele2_base==merged_file$allele1_target), "match", "mismatch")
sum(merged_file$match_check=="mismatch")

#########################################
# check for ambiguous snps in base data #
#########################################
# ambiguous snps are snps which are A/T or C/G (as same combination on forward and backward strand). 
merged_file$ambiguity <- NA
merged_file$ambiguity <- ifelse(merged_file$allele1_base=="A" & merged_file$allele2_base=="T" | (merged_file$allele1_base=="T" & merged_file$allele2_base=="A") | (merged_file$allele1_base=="C" & merged_file$allele2_base=="G") | (merged_file$allele1_base=="G" & merged_file$allele2_base=="C"), "ambiguous", "notambiguous") # identifies palindromic snps (A/T or C/G)
sum(merged_file$ambiguity=="ambiguous") #90
sum(merged_file$ambiguity=="notambiguous") #654
list_ambiguous_snps <- subset(merged_file, ambiguity == "ambiguous") #list of ambiguous snps
list_ambiguous_snps_distinct <- data.frame(rsid = c(list_ambiguous_snps[,"rsid"]))
list_ambiguous_snps_distinct <- distinct(list_ambiguous_snps_distinct) #88 distinct snps

###############################
# save list of ambiguous snps #
###############################
write.table(list_ambiguous_snps_distinct, "alspac_ambiguous_snps_88distinct_proxiesadded.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

