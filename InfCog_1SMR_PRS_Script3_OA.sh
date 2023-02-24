#!/bin/bash 
#SBATCH --job-name=EBI
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1 --tasks-per-node=1 --time=12:00:00
#SBATCH -o PRS_ALSPAC_script3.out
#---------------------------------------------

cd "INSERT WORKING DIRECTORY PATH HERE"

#---------------------------------------------
## module load
module load apps/plink/1.90

## NEW SCRIPT RUN - run next line for each score
## Generate sum scores = file inc var id, effect allele, beta
plink --bfile all_child  --score-no-mean-imputation --score ligthart_genomewide_proxiesadded_for_plink.txt sum --out ligthart_genomewide_proxiesadded_score_sum
plink --bfile all_child  --score-no-mean-imputation --score ligthartcis_for_plink.txt sum --out ligthartcis_score_sum
plink --bfile all_child  --score-no-mean-imputation --score han_genomewide_proxiesadded_for_plink.txt sum --out han_genomewide_proxiesadded_score_sum
plink --bfile all_child  --score-no-mean-imputation --score hancis_proxiesadded_for_plink.txt sum --out hancis_proxiesadded_score_sum
plink --bfile all_child  --score-no-mean-imputation --score ahluwalia_genomewide_for_plink.txt sum --out ahluwalia_genomewide_score_sum
plink --bfile all_child  --score-no-mean-imputation --score ahluwaliacis_for_plink.txt sum --out ahluwaliacis_score_sum
plink --bfile all_child  --score-no-mean-imputation --score borges_eafchecked_proxiesadded_for_plink.txt sum --out borges_eafchecked_proxiesadded_score_sum
plink --bfile all_child  --score-no-mean-imputation --score kettunen_for_plink.txt sum --out kettunen_score_sum
plink --bfile all_child  --score-no-mean-imputation --score rosa_for_plink.txt sum --out rosa_score_sum
plink --bfile all_child  --score-no-mean-imputation --score swerdlow_for_plink.txt sum --out swerdlow_score_sum
plink --bfile all_child  --score-no-mean-imputation --score sarwar_for_plink.txt sum --out sarwar_score_sum
plink --bfile all_child  --score-no-mean-imputation --score mahedy_wm_for_plink.txt sum --out mahedy_wm_score_sum
plink --bfile all_child  --score-no-mean-imputation --score mahedy_ert_for_plink.txt sum --out mahedy_ert_score_sum
plink --bfile all_child  --score-no-mean-imputation --score mahedy_inhibition_for_plink.txt sum --out mahedy_inhibition_score_sum