*Project: inflammation and cognition in ALSPAC
*What this script does:
*(1) cleans ALSPAC variables
*(2) runs regression models to predict performance on three cognitive tasks (working memory, emotion recognition, response inhibition) from inflammation (CRP, GlycA, IL-6) at age 24 in ALSPAC.
*(3) creates a dataset for multiple imputation (including people who have data on all three cognitive tasks at age 24)
*(4) checks what variables predict missingness in cognitive data at age 24 
*Chloe Slaney

************************************************
************* Cleaning ALSPAC data *************
************************************************
clear	
use "INSERT ALSPAC DATA FILE LOCATION", clear
version 16
log using "Output_ObsALSPAC_InfCog.log",replace

******************
**** Outcomes ****
******************
* 1. Working memory
*code from Liam Mahedy
*For each variable, missing variable codes checked, graphs and min/max of each is as expected using code below for each var.
*tab FKEP2010, nolab 
*graph box FKEP2010 if FKEP2010 > -1
*summ FKEP2010 if FKEP2010 > -1

* Creates new variables replacing all minus values with .
recode FKEP2010 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP2010)
recode FKEP2020 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP2020)
recode FKEP2030 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP2030)
recode FKEP2040 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP2040)
recode FKEP2050 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP2050)

* Creates variables for number incorrect and number hits
* wrong = false alarms (out of 48 trials), missed (out of 8 trials), and no response (out of 48 trials)
* hits = correct (out of 48 trials) and hits (out of 8 - there were 8 targets)
gen numwrong = rFKEP2020 + rFKEP2040 + rFKEP2050
gen numhits = rFKEP2010 + rFKEP2030

* Creates new variables
* hit miss = hits (out of 8) + miss (out of 8)
* facorneg = false alarms (out of 48) + correct (out of 48)
gen hitmiss = rFKEP2030 + rFKEP2040
gen facorneg = rFKEP2020 + rFKEP2010
gen hitr = rFKEP2030/hitmiss
gen far = rFKEP2020/ facorneg
recode hitr (1=0.9375) (0=0.0625)
recode far (1=0.9875) (0=0.0125) 

* Creates d prime for age 24
gen dp_24new = invnorm(hitr) - invnorm(far)

* Creating new d prime  
*dp_24new1 = includes values above or equal to zero and individuals who responded on at least half trials 
gen dp_24new1 = dp_24new>=0 & FKEP2050<24
gen c_2bnew = -(invnorm(hitr) + invnorm(far))/2

* Replace d prime individuals as missing if no response on more than 24 trials OR d prime score less than 0
replace dp_24new =. if FKEP2050>24 | dp_24new<0
replace c_2bnew =. if dp_24new ==.

* Trying to identify people who hit the same button all of the time
*tab numwrong if rFKEP2030==8

* 2. Emotion Recognition 
* Checked all variables(FKEP1010-FKEP1075): missing labels, graphs and numbers as expect.
* Replace missing values as . 
recode FKEP1010 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1010)
recode FKEP1015 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1015)
recode FKEP1020 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1020)
recode FKEP1025 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1025)
recode FKEP1030 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1030)
recode FKEP1035 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1035)
recode FKEP1040 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1040)
recode FKEP1045 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1045)
recode FKEP1050 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1050)
recode FKEP1055 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1055)
recode FKEP1060 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1060)
recode FKEP1065 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1065)
recode FKEP1070 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rERT)
recode FKEP1075 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP1075)

* 3. Response Inhibition 
* Checked all variables as expect them to be (FKEP3010 to FKEP3060)
* Replace missing values as . (FKEP3050 not included - many neg vals and not necessary for analysis)
recode FKEP3010 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP3010)
recode FKEP3020 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP3020)
recode FKEP3030 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP3030)
recode FKEP3040 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rFKEP3040)
recode FKEP3060 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rSSRT)

*******************
**** Exposures ****
*******************
* 1. CRP at age 24 
* Checked missing value codes, graphs and values as expect and replaced missing as .
recode CRP_F24 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rCRP_age24)

* Create CRP var excluding vals > 10 (i.e., sign of current infection - sensitivity analysis). The forward slash denotes a range of values (includes beginning and end of range).
recode CRP_F24 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.) (10/max=.), gen(rCRP_age24_sensitivity)
egen zscoreCRP_age24_sensitivity = std(rCRP_age24_sensitivity)

* 2. GlycA at age 24
* Checked missing value codes, graphs and values as expect and replaced missing as .
recode Gp_F24 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rGlycA_age24)
gen rGlycA_age24_sensitivity = rGlycA_age24
replace rGlycA_age24_sensitivity = . if rCRP_age24 >= 10
egen zscoreGlycA_age24_sensitivity = std(rGlycA_age24_sensitivity)

*3. IL6 at age 24 (Olink NPX log2 Scale)
*missing value (-1 and -9999)
recode IL6_F24 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.) (-9999=.), gen(rIL6_age24)

* IL6 CRP > 10 removed but CRP missing kept
gen rIL6_age24_sensitivity_v2 = rIL6_age24
replace rIL6_age24_sensitivity_v2 = . if (rCRP_age24 >= 10 & rCRP_age24 < 300) 
egen zscoreIL6_age24_sensitivity_v2 = std(rIL6_age24_sensitivity_v2)

*****************************
**** Potential Confounds ****
*****************************

* 1. Participant sex - male to 0 and female to 1
recode kz021 (-1=.) (1=0) (2=1), gen(rsex)
label define rsex_lb 0 "Male" 1 "Female"
label values rsex rsex_lb

* 2. Ethnicity (c804)
recode c804 (-1=.) (1=0) (2=1), gen(rethnicity)
label define rethnicity_lb 0 "White" 1 "Non-white"
label values rethnicity rethnicity_lb

* 3. BMI at age 24 (FKMS1040)
recode FKMS1040 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rBMI_age24)

* 4. Maternal education (c645)
recode c645 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rmated)
recode rmated (1=4) (2=3) (3=2) (4=1) (5=0)
label define rmated_lb 4 "CSE" 3 "Vocational" 2 "O level" 1 "A level" 0 "University Degree"
label values rmated rmated_lb

* Create mated only three categories (note: higher value = lower education)
recode rmated (0=0) (1=0) (2=1) (3=2) (4=2), gen(rmated_collap)
label define rmated_collap_lb 0 "> O level" 1 "O level" 2 "< O level"
label values rmated_collap rmated_collap_lb

* 5. Maternal SEP
* Note: armed forces = .
recode c755 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rmatSEP)
recode rmatSEP (1=0) (2=1) (3=2) (4=3) (5=4) (6=5) (65=.)
label define mat_SEP_lb 0 "Professional" 1 "Intermediate" 2 "Skilled (non-manual)" 3 "Skilled (manual)" 4 "Partly skilled" 5 "Unskilled"
label values rmatSEP mat_SEP_lb
*Collapsed at skilled
recode rmatSEP (1=0) (2=1) (3=2) (4=2) (5=3) (6=4) (65=.), gen(rmatSEP_collap)
label define mat_SEP_lb_collap 0 "Professional" 1 "Intermediate" 2 "Skilled" 3 "Partly skilled" 4 "Unskilled"
label values rmatSEP_collap mat_SEP_lb_collap

* 6. Smoking at age 24 (Fagerstrom Test)
* Replace missing values as . and recoding variables creating four groups 
recode FKSM1150 (-99=.) (-10=.) (-9=.) (-7=.) (-5=2) (-4=1) (-3=0) (-2=.) (-1=.) (1=3) (2=3) (3=3) (4=3) (5=3) (6=3) (7=3) (8=3) (9=3), gen(rsmoke_age24)
label define smoking_lb 0 "Never smoked whole cigarette" 1 "Not smoked last 30 days" 2 "Not daily smoker" 3 "Daily smoker"
label values rsmoke_age24 smoking_lb

* 7. Alcohol at age 24 (AUDIT-C)
recode FKAL1500 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(ralcohol_age24)

* 8. IQ at age 8 (f8ws112)
recode f8ws112 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-3=.) (-2=.) (-1=.), gen(rIQ_age8)

*****************************
**** Auxiliary Variables ****
*****************************
* 1. CRP variables
*Checked all missing value codes (-1 only) and replaced missing values as .
recode CRP_f9 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rCRP_age9)
recode crp_TF3 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rCRP_age15)
recode CRP_TF4 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rCRP_age17)

* 2. GlycA variables
*Checked all missing value codes (-1 only) and replaced missing values as .
recode Gp_F7 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rGp_age7)
recode Gp_TF3 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rGp_age15)
recode Gp_TF4 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rGp_age17)

* 3. Working memory at 10
* Checked missing value codes (-9 and -2 only) and replaced missing values as .
recode fdcm110 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rWM_age10)

* 4. Asthma/eczema (kv1070)
* Checked missing value codes (-10 and -1 only) and replaced missing values as .
recode kv1070 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.) (4=0) (1=1) (2=1) (3=1), gen(rasthmaeczema_age10)
label define asthmaeczema_lb 0 "Not present" 1 "Yes, either or both"
label values rasthmaeczema_age10 asthmaeczema_lb

* 5. Inhibition at age 10
* Checked missing value codes (-9 and -2 only)
recode fdcm215 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rInhib_age10)

* 6. Smoking at 17 (number of cigarettes smoke each day)
* a. ever smoked
* Checked missing value codes (-10 and -4 only)
recode FJSM050 (-10=.) (-4=.) (2=0), gen(reversmoke_age17)
label define smoke_17_lb 0 "No" 1 "Yes"
label values reversmoke_age17 smoke_17_lb

* b. smokes daily
* Checked missing value codes (-10, -4 and -1 only)
recode FJSM350 (-10=.) (-4=.) (-1=.) (2=0), gen(rsmokesdaily_age17)
label define smokedaily_17_lb 0 "No" 1 "Yes"
label values rsmokesdaily_age17 smokedaily_17_lb

* 7. Alcohol at 17 (AUDIT)
* Checked missing value codes (-10 and -1 only)
recode FJAL4000 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(ralcohol_age17)

* 8. Financial difficulties - pregnancy
*Checked missing value codes (-7 and -1 only)
recode c525 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rmat_financediff)

* 9. Age of mother (delivery) 
* Checked missing values (-10, -4, -2) incorporated
recode mz028b (-10=.) (-9=.) (-8=.) (-7=.) (-4=.) (-99=.) (-2=.) (-1=.), gen(rmatage)

* 10. BMI at 7 
* missing (-9 and -1 only)
recode f7ms026a (-9=.) (-1=.), gen(rBMI_age7)

* 11. BMI at 17 
* missing (-10, -4 and -1 only)
recode FJMR022a (-10=.) (-4=.) (-1=.), gen(rBMI_age17)

*12. BMI at 9 (f9ms026a)
* missing (-9 and -1 only)
recode f9ms026a (-9=.) (-1=.), gen(rBMI_age9)

* 13. Depression at 24 
* Checked missing value codes (-99,-10,-9,-7,-2,-1 for all) and created binary variable 
recode FKDQ1000 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rmilddep_age24)
label define depression_age24_lb 0 "no" 1 "yes"
label values rmilddep_age24 depression_age24_lb
recode FKDQ1010 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rmoddep_age24)
label values rmoddep_age24 depression_age24_lb
recode FKDQ1020 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rsevdep_age24)
label values rsevdep_age24 depression_age24_lb

* 14. Mother told child has autism/aspergers/autistic spectrum
*missing (-10 and -1 only). Yes is 1.
recode ku360 (-10=.) (-1=.) (1=0) (2=1), gen(rASD)
label define diagnosisASD_lb 0 "No" 1 "Yes"
label values rASD diagnosisASD_lb

*15. Social difficulties score (ku710b - prorated)
* missing (-10, -6, -5 only)
recode ku710b (-10=.) (-6=.) (-5=.), gen(rSDQ)

*16. Any infection in 1st 3 months
* missing (-7, -1 only)
recode b065 (-7=.) (-1=.) (2=0), gen(rInfection_first3months)
label define infect_earlypreg_lb 0 "No" 1 "Yes"
label values rInfection_first3months infect_earlypreg_lb

*17. Birthweight
* missing (-10, -1 only)
recode kz030 (-10=.) (-1=.), gen(rbwt)

*18. Maternal - any infection
*missing (-7 and -1 only). Yes =1.
recode c064 (-7=.) (-1=.) (2=0), gen(rMatInfect)
label define infectpreg_lb 0 "No" 1 "Yes"
label values rMatInfect infectpreg_lb

*19.  Maternal depression (18wk and 32wk)
*missing (-7 and -1 only)
recode b370 (-7=.) (-1=.), gen(rEPDS_18wk)
* missing (-1 only)
recode e390 (-7=.) (-1=.), gen(rEPDS_32wk)

*20. Paternal Social Class
recode c765 (-10=.) (-9=.) (-8=.) (-7=.) (-99=.) (-2=.) (-1=.), gen(rpaternalSEP)
recode rpaternalSEP (1=0) (2=1) (3=2) (4=3) (5=4) (6=5) (65=.)
label define paternal_SEP_lb 0 "Professional" 1 "Intermediate" 2 "Skilled (non-manual)" 3 "Skilled (manual)" 4 "Partly skilled" 5 "Unskilled"
label values rpaternalSEP paternal_SEP_lb

*21. BMI at 13.5 years
* missing (-10, -6, -1)
recode fg3139 (-10=.) (-6=.) (-1=.), gen(rBMI_age13)

*22. BMI at 15.5 years
* missing (-10, -6, -1)
recode fh3019 (-10=.) (-6=.) (-1=.), gen(rBMI_age15)

*23. Maternal number cigarettes smoke (21m post delivery)
* Checked missing values (-1 only)
recode g820 (-1=.) (5=1) (10=2) (15=2) (20=3) (25=3) (30=4),gen(rmatsmoke_21m)
label define matsmoke_21m_lb 0 "Does not smoke" 1 "Smokes 1-9 cigs per day" 2 "Smokes 10-19 cigs per day" 3 "Smokes 20-29 cigs per day" 4 "Smokes > 30 cigs per day"
label values rmatsmoke_21m matsmoke_21m_lb

*24. Maternal current smoke (18 years)
* Checked missing values (-10 and -1)
recode t5520 (-10=.) (-1=.) (2=0), gen(rmatsmoke_18y)
label define matsmoke_age18_lb 0 "No" 1 "Yes"
label values rmatsmoke_18y matsmoke_age18_lb

*25. ImmunoglobinE at 7 years
* missing values (-1 only)
recode IGE_F7 (-1=.), gen(rIgE_age7)

*26. Lymphocytes at 24 years
*missing values (-99 and -1)
recode Lymphocytes_F24 (-99=.) (-1=.), gen(rLymphocytes_age24)

*27. Monocytes at 24 years
*missing values (-99 and -1)
recode Monocytes_F24 (-99=.) (-1=.), gen(rMonocytes_age24)

*28. IL6 at age 9
*missing value (-1 only)
recode IL6_f9 (-1=.), gen(rIL6_age9)

******************************
**** Create Log Variables ****
******************************
* Outcomes - working memory, emotion recognition and response inhibition
gen dp_24new_nlog=ln(dp_24new)
gen rERT_nlog=ln(rERT)
gen rSSRT_nlog=ln(rSSRT)

* Exposures
gen rCRP_age24_nlog=ln(rCRP_age24)
gen rGp_age24_nlog=ln(rGlycA_age24)
gen rIL6_age9_nlog=ln(rIL6_age9)
gen rIL6_age24_nlog=ln(rIL6_age24)

*****************************************
**** Creating standardized variables ****
*****************************************
egen zscoreCRP_age24 = std(rCRP_age24)
egen zscoreGlycA_age24 = std(rGlycA_age24)
egen zscoreIL6_age24 = std(rIL6_age24)
egen zscoreERT = std(rERT)
egen zscoredp_24new = std(dp_24new)
egen zscoreSSRT = std(rSSRT)

************************************
**** Descriptives/Visualisation ****
************************************
*Working Memory
summarize(dp_24new)
histogram(dp_24new)
histogram(dp_24new_nlog)

*Emotion Recognition Task
summarize(rERT)
histogram(rERT)

*Response Inhibition
summarize(rSSRT)
histogram(rSSRT)

*IL-6 age 9 
summarize(rIL6_age9)
histogram(rIL6_age9)
histogram(rIL6_age9_nlog)

*CRP age 24
summarize(rCRP_age24)
histogram(rCRP_age24)
histogram(rCRP_age24_nlog)

*GlycA age 24
summarize(rGlycA_age24)
histogram(rGlycA_age24)

****************************************************
************* Regression Models ********************
****************************************************
* Stata Handbook release 17: operator i. (specify indicators), c. (continuous), o. (omit or indicator variable), # (interactions), ## (full-factorial interactions).
* Note: regressions NOT restricted to N of final model (i.e., individuals who have data on all exposures, outcomes and covariates)
* Using natural log did not result in normal distribution so run on raw data.
* Includes standardised exposure and outcome

****************************
**** 1. CRP as exposure ****
****************************
*model 1 - unadjusted
regress zscoredp_24new c.zscoreCRP_age24
regress zscoreERT c.zscoreCRP_age24
regress zscoreSSRT c.zscoreCRP_age24

*model 2 - including sex, ethnicity and BMI
regress zscoredp_24new c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24
regress zscoreERT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24
regress zscoreSSRT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24

*model 3 - including socioeconomic factors 
regress zscoredp_24new c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreERT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreSSRT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

*model 4 - including smoking and alcohol use
regress zscoredp_24new c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreERT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreSSRT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

*model 5 - including IQ at age 8
regress zscoredp_24new c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreERT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreSSRT c.zscoreCRP_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8

******************************
**** 2. GlycA as exposure ****
******************************
*model 1 - unadjusted 
regress zscoredp_24new c.zscoreGlycA_age24
regress zscoreERT c.zscoreGlycA_age24
regress zscoreSSRT c.zscoreGlycA_age24

*model 2 - including sex, ethnicity and BMI
regress zscoredp_24new c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24
regress zscoreERT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24
regress zscoreSSRT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24

*model 3 - including socioeconomic factors 
regress zscoredp_24new c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreERT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreSSRT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

*model 4 - including smoking and alcohol use
regress zscoredp_24new c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreERT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreSSRT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

*model 5 - including IQ at age 8
regress zscoredp_24new c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreERT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreSSRT c.zscoreGlycA_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8

*****************************
**** 3. IL-6 as exposure ****
*****************************
*model 1 - unadjusted 
regress zscoredp_24new c.zscoreIL6_age24
regress zscoreERT c.zscoreIL6_age24
regress zscoreSSRT c.zscoreIL6_age24

*model 2 - including sex, ethnicity and BMI
regress zscoredp_24new c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24
regress zscoreERT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24
regress zscoreSSRT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24

*model 3 - including socioeconomic factors 
regress zscoredp_24new c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreERT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreSSRT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

*model 4 - including smoking and alcohol use
regress zscoredp_24new c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreERT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreSSRT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

*model 5 - including IQ at age 8
regress zscoredp_24new c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreERT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreSSRT c.zscoreIL6_age24 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8

**************************************************************************
******** Sensitivity Analysis - removing cases rCRP_F24 > 10 mg/L ********
**************************************************************************
*only CRP >10mg/L removed (CRP missing kept)
****************************
**** 1. CRP as exposure ****
****************************
*model 1 - unadjusted
regress zscoredp_24new c.zscoreCRP_age24_sensitivity
regress zscoreERT c.zscoreCRP_age24_sensitivity
regress zscoreSSRT c.zscoreCRP_age24_sensitivity

*model 2 - including sex, ethnicity and BMI
regress zscoredp_24new c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24
regress zscoreERT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24
regress zscoreSSRT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24

*model 3 - including socioeconomic factors 
regress zscoredp_24new c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreERT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreSSRT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

*model 4 - including smoking and alcohol use
regress zscoredp_24new c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreERT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreSSRT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

*model 5 - including IQ at age 8
regress zscoredp_24new c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreERT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreSSRT c.zscoreCRP_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8

******************************
**** 2. GlycA as exposure ****
****************************** 
*model 1 - unadjusted
regress zscoredp_24new c.zscoreGlycA_age24_sensitivity
regress zscoreERT c.zscoreGlycA_age24_sensitivity
regress zscoreSSRT c.zscoreGlycA_age24_sensitivity

*model 2 - including sex, ethnicity and BMI
regress zscoredp_24new c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24
regress zscoreERT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24
regress zscoreSSRT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24

*model 3 - including socioeconomic factors 
regress zscoredp_24new c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreERT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreSSRT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

*model 4 - including smoking and alcohol use
regress zscoredp_24new c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreERT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreSSRT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

*model 5 - including IQ at age 8
regress zscoredp_24new c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreERT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreSSRT c.zscoreGlycA_age24_sensitivity i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8


**************************
**** IL-6 sensitivity ****
**************************
*model 1 - unadjusted
regress zscoredp_24new c.zscoreIL6_age24_sensitivity_v2
regress zscoreERT c.zscoreIL6_age24_sensitivity_v2
regress zscoreSSRT c.zscoreIL6_age24_sensitivity_v2

*model 2 - including sex, ethnicity and BMI
regress zscoredp_24new c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24
regress zscoreERT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24
regress zscoreSSRT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24

*model 3 - including socioeconomic factors 
regress zscoredp_24new c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreERT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP
regress zscoreSSRT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP

*model 4 - including smoking and alcohol use
regress zscoredp_24new c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreERT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24
regress zscoreSSRT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24

*model 5 - including IQ at age 8
regress zscoredp_24new c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreERT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8
regress zscoreSSRT c.zscoreIL6_age24_sensitivity_v2 i.rsex i.rethnicity c.rBMI_age24 i.rmated i.rmatSEP i.rsmoke_age24 c.ralcohol_age24 c.rIQ_age8

******************************************************************************************
********************* Creating New Dataset for Multiple Imputation ***********************
******************************************************************************************
*Note: always ensure imputed model consistent analytic model  (include all vars in analytic model esp DV, recodes, transforms)
*keeping participants who have data on all three cognitive task
*gen allcog = 1 if dp_24new !=. & rERT !=. & rSSRT !=.
*keep if allcog == 1

*keep if !missing(dp_24new)
*keep if !missing(rSSRT)
*keep if !missing(rERT)

*keeping variables for new imputation dataset
*keep dp_24new rERT rSSRT rCRP_age24 rGlycA_age24 rIL6_age24 zscoredp_24new zscoreERT zscoreSSRT zscoreCRP_age24 zscore GlycA_age24 zscoreIL6_age24 rsex rethnicity rBMI_age24 rmated rmatSEP rsmoke_age24 ralcohol_age24 rIQ_age8 rCRP_age9 rCRP_age15 rCRP_age17 rGp_age7 rGp_age15 rGp_age17 rWM_age10 rasthmaeczema_age10 rInhib_age10 reversmoke_age17 rsmokesdaily_age17 ralcohol_age17 rmat_financediff rmatage rBMI_age17 rBMI_age7 rmilddep_age24 rmoddep_age24 rsevdep_age24 rInfection_first3months rbwt rMatInfect rEPDS_18wk rEPDS_32wk rpaternalSEP rBMI_age13 rBMI_age15 rmatsmoke_21m rmatsmoke_18y rIgE_age7 rLymphocytes_age24 rMonocytes_age24 rIL6_age9

**************************************************************************************
****************** Predict missingness for each outcome ******************************
**************************************************************************************
gen outcome_wm_missing = .
replace outcome_wm_missing = 1 if missing(dp_24new)
replace outcome_wm_missing = 0 if !missing(dp_24new)

gen outcome_ert_missing = .
replace outcome_ert_missing = 1 if missing(rERT)
replace outcome_ert_missing = 0 if !missing(rERT)

gen outcome_ssrt_missing = .
replace outcome_ert_missing = 1 if missing(rERT)
replace outcome_ert_missing = 0 if !missing(rERT)

gen outcome_allcog_missing = 1
replace outcome_allcog_missing = 0 if !missing(rERT) & !missing(dp_24new) & !missing(rSSRT)

*predict missingness in cognitive outcome data
logit outcome_allcog_missing c.rIQ_age8
logit outcome_allcog_missing i.rsex
logit outcome_allcog_missing i.rethnicity
logit outcome_allcog_missing c.rBMI_age24
logit outcome_allcog_missing i.rmated
logit outcome_allcog_missing i.rmatSEP
logit outcome_allcog_missing i.rsmoke_age24
logit outcome_allcog_missing c.ralcohol_age24
logit outcome_allcog_missing c.rCRP_age24
logit outcome_allcog_missing c.rGlycA_age24
logit outcome_allcog_missing c.rIL6_age24

*continued...early life inflammation
logit outcome_allcog_missing c.rCRP_age9
logit outcome_allcog_missing c.rCRP_age15
logit outcome_allcog_missing c.rCRP_age17
logit outcome_allcog_missing c.rGp_age7
logit outcome_allcog_missing c.rGp_age15
logit outcome_allcog_missing c.rGp_age17
logit outcome_allcog_missing c.rIL6_age9

*possible control
logit outcome_allcog_missing c.rmatage

log close