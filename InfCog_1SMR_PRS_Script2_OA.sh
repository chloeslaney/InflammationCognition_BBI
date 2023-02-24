#!/bin/bash 
#SBATCH --job-name=EBI
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1 --tasks-per-node=1 --time=12:00:00
#SBATCH -o PRS_ALSPAC_script2.out

#---------------------------------------------

cd "INSERT WORKING DIRECTORY PATH"

#---------------------------------------------
## module load
module load apps/plink/1.90

## NEW SCRIPT - RUN UNTIL HERE - only include chr taken from script1
## combine to single dataset
plink --bfile bg_chr01_child --merge-list allfiles_child.txt --make-bed --out all_child
