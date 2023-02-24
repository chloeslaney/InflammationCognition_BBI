*Project: inflammation and cognition
*What this script does: (1) runs multiple imputation and (2) runs regression models (inflammation and cognition) on imputed data
*Chloe Slaney

************************************************
************** Load in ALSPAC data *************
************************************************
clear
use "INSERT ALSPAC DATA FILE LOCATION - INCLUDING ONLY INDIVIDUALS WHO HAVE DATA ON ALL THREE COGNITIVE MEASURES AT AGE 24"
log using "MultipleImputation_100impute.log", replace

*********************
**** Prior to MI ****
*********************
* Step 1. Examine number and proportion missing values among vars interest.
ssc install mdesc 
mdesc zscoredp_24new zscoreERT zscoreSSRT rsex zscoreCRP_age24 zscoreGlycA_age24 zscoreIL6_age24 rethnicity rBMI_age24 rsmoke_age24 ralcohol_age24 rIQ_age8 rmated rmatSEP rCRP_age9 rCRP_age15 rCRP_age17 rGp_age7 rGp_age15 rGp_age17 rWM_age10 ralcohol_age17 rmat_financediff rmatage rBMI_age17 rBMI_age7 rBMI_age13 rBMI_age15 rEPDS_18wk rMonocytes_age24 rLymphocytes_age24 rpaternalSEP rIL6_age9

* Step 2. Examine missing data patterns among vars interest and declare dataset as "mi" deatset and style store data
*mi set wide
*mi misstable summ dp_24new rERT rSSRT rCRP_age24 rGlycA_age24 rsex rethnicity rBMI_age24 rmated rmatSEP rsmoke_age24 ralcohol_age24 rIQ_age8 
*mi misstable patterns dp_24new rERT rSSRT rCRP_age24 rGlycA_age24 rsex rethnicity rBMI_age24 rmated rmatSEP rsmoke_age24 ralcohol_age24 rIQ_age8 

* Step 3. Identify auxiliary variables
*pwcorr zscoreCRP_age24 zscoreGlycA_age24 rethnicity rBMI_age24 rsmoke_age24 ralcohol_age24 rIQ_age8 rmated rmatSEP rCRP_age9 rCRP_age15 rCRP_age17 rGp_age7 rGp_age15 rGp_age17 rWM_age10 ralcohol_age17 rmat_financediff rmatage rBMI_age17 rBMI_age7 rBMI_age13 rBMI_age15 rEPDS_18wk rMonocytes_age24 rLymphocytes_age24 rpaternalSEP rIL6_age9

*********************************
****** Multiple Imputation ******
*********************************
* prepare datset for imputation
mi set wide

* variables included that do not need to be imputed (i.e., outcomes)
mi register regular zscoredp_24new zscoreERT zscoreSSRT rsex

* register variables that will need to be imputed
mi register imputed zscoreCRP_age24 zscoreGlycA_age24 zscoreIL6_age24 rethnicity rBMI_age24 rsmoke_age24 ralcohol_age24 rIQ_age8 rmated rmatSEP rCRP_age9 rCRP_age15 rCRP_age17 rGp_age7 rGp_age15 rGp_age17 rWM_age10 ralcohol_age17 rmat_financediff rmatage rBMI_age17 rBMI_age7 rBMI_age13 rBMI_age15 rEPDS_18wk rMonocytes_age24 rLymphocytes_age24 rpaternalSEP rIL6_age9

* chained imputations - 100 imputed models
mi impute chained (regress) zscoreCRP_age24 zscoreGlycA_age24 zscoreIL6_age24 rBMI_age24 ralcohol_age24 rIQ_age8 rCRP_age9 rCRP_age15 rCRP_age17 rGp_age7 rGp_age15 rGp_age17  rWM_age10 ralcohol_age17 rmat_financediff rmatage rBMI_age17 rBMI_age7 rBMI_age13 rBMI_age15 rEPDS_18wk rLymphocytes_age24 rMonocytes_age24 rIL6_age9(logit) rethnicity (ologit) rmated rmatSEP rpaternalSEP rsmoke_age24 = zscoredp_24new zscoreERT zscoreSSRT i.rsex, add(100) rseed(1234)

*************************************
******** Sensitivity analysis *******
*************************************
* predictive mean matching - 10 nearest neighbours - 100 imputed models (also changed this to 50 to check if this impacted SE)
*mi impute chained (pmm, knn(10)) zscoreCRP_age24 zscoreGlycA_age24 zscoreIL6_age24 rBMI_age24 ralcohol_age24 rIQ_age8 rCRP_age9 rCRP_age15 rCRP_age17 rGp_age7 rGp_age15 rGp_age17  rWM_age10 ralcohol_age17 rmat_financediff rmatage rBMI_age17 rBMI_age7 rBMI_age13 rBMI_age15 rEPDS_18wk rLymphocytes_age24 rMonocytes_age24 rIL6_age9(logit) rethnicity (ologit) rmated rmatSEP rpaternalSEP rsmoke_age24 = zscoredp_24new zscoreERT zscoreSSRT i.rsex, add(100) rseed(1234)

**********************************
******** Imputed analysis ********
**********************************

***********
****CRP****
***********
* model 1 - unadjusted 
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreCRP_age24

* model 2 - including sex, ethnicity and BMI
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24

* model 3 - additional include mated and maternal SEP
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

* model 4 - additional include smoking and alcohol use 
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

* model 5 - additional IQ
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8  

*************
****GlycA****
*************
* model 1 - unadjusted 
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreGlycA_age24

* model 2 - including sex, ethnicity and BMI
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24

* model 3 - additional include mated and maternal SEP
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

* model 4 - additional include smoking and alcohol use 
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

* model 5 - additional IQ
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8

***********
****IL6****
***********
* model 1 - unadjusted 
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreIL6_age24

* model 2 - including sex, ethnicity and BMI
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24

* model 3 - additional include mated and maternal SEP
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

* model 4 - additional include smoking and alcohol use 
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

* model 5 - additional IQ
mi estimate: mvreg zscoredp_24new zscoreERT zscoreSSRT = c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8

log close