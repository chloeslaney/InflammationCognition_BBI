#!/bin/bash 
#SBATCH --job-name=EBI
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1 --tasks-per-node=1 --time=12:00:00
#SBATCH -o PRS_ALSPAC_script1.out
#---------------------------------------------

cd "INSERT WORKING DIRECTORY HERE"

#---------------------------------------------
## all PRS code adapted from Hannah Sallis code
## module load
module load apps/plink/1.90

## location of genetic data - 1000 genomes
BESTGUESS= "INSERT LOCATION OF GENETIC DATA" ## Path to the genotype file in plink binary format to be used in score construction - leave out the plink file extensions (*.bed, *.bim, *.fam).

## unrelated IDs - independent individuals
UNREL="INSERT LOCATION OF LIST OF UNRELATED IDs"

## extract snps
for i in {01..22}
do
plink --bfile ${BESTGUESS}/data_chr${i} --extract INFLAMMATION_COGNITION_SNPs_PROXIESADDED.txt --keep ${UNREL}/children_unrelated.txt --remove withdrawn_consent_20211111_for_plink.txt --out bg_chr${i}_child --make-bed
done
