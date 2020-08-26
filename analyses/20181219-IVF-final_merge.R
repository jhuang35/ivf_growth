library(tidyverse)
library(ggplot2)
library(haven)
library(summarytools)
library(readxl)
library(broom)

options(tibble.print_max = 40, tibble.print_min = 30)

# 1. Set paths
plot_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/"
base_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

# 2. Merge
## DEMO AND ANTHRO
clean_data <- read_xlsx(paste0(base_path, "FormA428&429_20181001.xlsx"), na = c("N/A", "NA"," ", ""))
zbmi_annotate <- read_dta(paste0(base_path,"20181001-zBMI.dta")) %>% 
  select(SubjectID, visit, days, height, weight, bmi, zlen, zwei, zwfl, cbmi, zbmi)

merged <- clean_data %>% select(SubjectID, sex, mother_ethnicity, household_income, mother_highest_education, gdm_who_1999, gdm_treatment,
                                Dur_full_BF = Duration_full_BF_regrp_no_solids_new_4CAT, Dur_any_BF = Duration_anyBF_cat_new, ogtt_fasting_pw26, ogtt_2hour_pw26,
                                ppWeight, ppBMI, ppBMI_WHOclass, m_weight_pw26, m_height_pw26, mother_age_delivery, parity, marital_status, accommodation, 
                                GA, mode_of_delivery, HTN_clean, m_midarm_pw26, m_triceps_pw26, m_biceps_pw26, m_subscapular_pw26, m_suprailiac_pw26,
                                m_age_pw26, last_antenatal_weight, last_antenatal_weight_GA, weight1, weight2, weight3, weight4, weight5, weight6, weight7,
                                weight8, weight9, contains("gestationalWeek"), contains("Plasma"), contains("METPA"), mode_of_delivery) %>%
  left_join(zbmi_annotate, ., by = "SubjectID") %>% 
  mutate(IVF = substring(SubjectID,1,3) %in% c("019", "029"))
## PATERNAL ANTHROPOMETRICS
pat_anthro <- read_xlsx(paste0(base_path, "Paternal_anthropometry_20180815.xlsx"), na = c("N/A", "NA"," ", ""))
merged <- left_join(merged, pat_anthro)
## SMOKING VARIABLES
smoke <- read_excel(paste0(base_path, "Smoking variables_datateam0406.xls"), na = c("N/A", "NA"," ", ""))
merged <- smoke %>% mutate(SubjectID = subjectid, 
                           smk_prev = if_else(q11_3_smoking_before_pregnancy == "1_yes", 1, 
                                              if_else(q11_3_smoking_before_pregnancy == "0_no", 0, NA_real_)),
                           smk_home = if_else(q12_3_anyone_living_in_home_smok == "1_yes", 1, 
                                              if_else(q12_3_anyone_living_in_home_smok == "0_no", 0, NA_real_)),
                           smk_curr = if_else(q12_1_currently_smoking == "1_yes", 1, 
                                              if_else(q12_1_currently_smoking == "0_no", 0, NA_real_))) %>% 
  select(SubjectID, smk_prev, smk_curr, smk_home) %>% left_join(merged, .)
## QMR 
qmr_5 <- read_xlsx(paste0(base_path, "QMR_ext_year5_20181116_main variables.xlsx"), na = c("N/A", "NA"," ", ""))
qmr_6 <- read_xlsx(paste0(base_path, "QMR_ext_year6_20181116_main variables.xlsx"), na = c("N/A", "NA"," ", ""))
merged <- qmr_5 %>% select(SubjectID, fat_5_qmr = Fat, lean_5_qmr = Lean, wt_5_qmr = c_weight_yr5, age_5_qr = child_age_at_visit) %>% left_join(merged, .)
merged <- qmr_6 %>% select(SubjectID, fat_6_qmr = Fat_kg, lean_6_qmr = Lean_kg, wt_6_qmr = c_weight_yr6, age_6_qr = child_age_at_visit) %>% left_join(merged, .)
# MATERNAL PRENATAL BP
bp <- read_csv(paste0(base_path, "Delivery_CRF_5_Pregnancy_Data.csv"), na = c("N/A", "NA"," ", ""))
bp <- bp %>% rename(sbp1 = qn_1_antenatal_systolic_BP1, dbp1 = qn_1_antenatal_diastolic_BP1,
                    sbp2 = qn_1_antenatal_systolic_BP2, dbp2 = qn_1_antenatal_diastolic_BP2,
                    sbp3 = qn_1_antenatal_systolic_BP3, dbp3 = qn_1_antenatal_diastolic_BP3) %>% 
  mutate(pe = if_else(sbp1 >= 160 | sbp2 >= 160 | sbp3 >= 160 | dbp1 >= 110 | dbp2 >= 110 | dbp3 >= 110, 1, 0),
         hi_bp = if_else(sbp1 >= 140 | sbp2 >= 140 | sbp3 >= 140 | dbp1 >= 90 | dbp2 >= 90 | dbp3 >= 90, 1, 0))
merged <- bp %>% select(SubjectID = PSCID, contains("sbp"), contains("dbp"), pe, hi_bp) %>% left_join(merged, .)
# FETAL SIZE FROM U/S
fetal <- read_dta(paste0(base_path,"20181211-us_efw.dta"))
fetal <- fetal %>% select(-EDD, -DATE) %>% gather(var, val, -SubjectID, -visit) %>% mutate(var = paste0(var,"_",visit)) %>% select(-visit) %>% spread(var, val)
merged <- left_join(merged, fetal)
# FASTING GLUCOSE 6 YEARS 
glucose <- read_xlsx(paste0(base_path,"Child OGTT values as of 30 Sep 2017.xlsx")) %>% rename(SubjectID = SUBJID)
merged <- left_join(merged, glucose)
# MERGE SKINFOLDS (WIDE)
skin <- clean_data %>% select(SubjectID, contains("c_tricep"), contains("c_subscap"), 
                              contains("c_bicep"), contains("c_suprailiac"), -contains("day1"), -contains("wk1"), c_triceps_day1, c_subscapular_day1)
temp <- skin %>% rename(c_triceps_0 = c_triceps_day1, c_subscapular_0 = c_subscapular_day1,
                        c_triceps_18 = c_triceps_m18, c_triceps_24 = c_triceps_m24, c_triceps_36 = c_triceps_m36,
                        c_triceps_48 = c_triceps_yr4, c_triceps_54 = c_triceps_yr4.5, c_triceps_60 = c_triceps_yr5,
                        c_triceps_66 = c_triceps_yr5.5, c_triceps_72 = c_triceps_yr6, c_triceps_78 = c_triceps_yr6.5) %>%
  rename(c_subscapular_18 = c_subscapular_m18, c_subscapular_24 = c_subscapular_m24, c_subscapular_36 = c_subscapular_m36,
         c_subscapular_48 = c_subscapular_yr4, c_subscapular_54 = c_subscapular_yr4.5, c_subscapular_60 = c_subscapular_yr5,
         c_subscapular_66 = c_subscapular_yr5.5, c_subscapular_72 = c_subscapular_yr6, c_subscapular_78 = c_subscapular_yr6.5) %>%
  rename(c_biceps_18 = c_biceps_m18,c_biceps_24 = c_biceps_m24,c_biceps_36 = c_biceps_m36,c_biceps_48 = c_biceps_yr4,
         c_biceps_54 = c_biceps_yr4.5,c_biceps_60 = c_biceps_yr5,c_biceps_66 = c_biceps_yr5.5,c_biceps_72 = c_biceps_yr6,
         c_biceps_78 = c_biceps_yr6.5) %>%
  rename(c_suprailiac_48 = c_suprailiac_yr4,c_suprailiac_54 = c_suprailiac_yr4.5,c_suprailiac_60 = c_suprailiac_yr5,
         c_suprailiac_66 = c_suprailiac_yr5.5,c_suprailiac_72 = c_suprailiac_yr6,c_suprailiac_78 = c_suprailiac_yr6.5) %>%
  gather(measure, value, -SubjectID) %>% separate(measure, c("lead", "param", "visit")) %>% select(-lead) %>% spread(param, value)
temp <- temp %>% mutate(visit = as.double(visit))
merged <- left_join(merged, temp)
# ADD MRI
MRI_0 %>% select(TotalAbdVol_D_S_18Sept)
MRI_0 <- read_sav(paste0(base_path,"Neonatal MRI_2014-09-18_total abd fat_updated.sav")) %>% 
  select(SubjectID = PSCID, Age_MRI = Age_MRIday_D7, 
         SSAT = sSAT_18Sept, DSAT = dSAT_18Sept, IAT = IAT_18Sept, 
         AbVol = TotalAbdVol_D_S_18Sept) %>% add_column(visit = 0) %>% 
  mutate(SSAT_PCT = SSAT / AbVol, DSAT_PCT = DSAT/AbVol, IAT_PCT = IAT/AbVol)
merged <- left_join(merged, MRI_0)

MRI_54_1 <- read_sav(paste0(base_path,"Data_Transfer_MRS_4.5yr_9_june_2016.sav"))
MRI_54_2 <- read_sav(paste0(base_path,"MRI_54mth.sav")) %>% rename(SubjectID = subjectID)
MRI_54 <- full_join(MRI_54_1, MRI_54_2) %>% add_column(visit = 54)
merged <- left_join(merged, MRI_54)
# CHILD BLOOD PRESSURES
BP <- read_xlsx(paste0(base_path,"Child_Blood_Pressure_20180828.xlsx")) %>% 
  select(SubjectID, contains("sbp"), contains("dbp"), -contains("sitting"), -contains("yr7"))
BP <- BP %>% rename(SBP_36 = c_sbp_m36,
                    SBP_48 = c_sbp_yr4,
                    SBP_60 = c_sbp_yr5,
                    SBP_72 = c_sbp_yr6,
                    DBP_36 = c_dbp_m36,
                    DBP_48 = c_dbp_yr4,
                    DBP_60 = c_dbp_yr5,
                    DBP_72 = c_dbp_yr6) 
BP <- BP %>% gather(measure, value, -SubjectID) %>% separate(measure, c("param", "visit")) %>% spread(param, value) %>% mutate(visit = as.double(visit))
merged <- left_join(merged, BP)
# Generate Z-scores
# cbp_export <- merged %>% filter(visit %in% c(36,48,60,72)) %>%
#   select(id = SubjectID, sex, age = visit, agemon = days, height, sysbp = SBP, diabp = DBP) %>%
#   mutate(age = age/12, agemon = floor(agemon/365.25*12), sex = if_else(sex == "Male", 1, 2))
# write_dta(cbp_export, "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Blood pressure/20190104-CBP.dta")
# Data processing separately in: 20190104-cbp_pct
cbp_pct <- read_dta("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Blood pressure/results.dta") %>% 
  transmute(SubjectID = id, visit = age * 12, ZSBP = sys_per, ZDBP = dia_per, BP_out = outlier)
merged <- left_join(merged, cbp_pct, by = c("SubjectID", "visit"))

# 3. Save STATA file
write_dta(merged, paste0(base_path,"20181219-IVF-final.dta"))

# 4. Load STATA file
# merged <- read_dta(paste0(base_path,"20181219-IVF-final.dta"))


############################
############################
## DEPRECATED
############################
############################
# data merge verification
# clean_data <- read_xlsx(paste0(base_path, "FormA428&429_20180905.xlsx"), na = c("N/A", "NA"," ", ""))
# tmp <- left_join(merged3, clean_data, by = "SubjectID")
# tmp %>% select(contains(".y")) %>% names()
# 
# tmp %>% filter(visit == 1) %>% 
#   group_by(NOTMATCH = (mother_ethnicity.x != mother_ethnicity.y)) %>% 
#   count(mother_ethnicity.x, mother_ethnicity.y) #missing values present in clean_data
# 
# tmp %>% filter(visit == 1) %>% 
#   group_by(NOTMATCH = (mother_age_delivery.x != mother_age_delivery.y)) %>% 
#   count(NOTMATCH)
# 
# tmp %>% filter(visit == 1) %>% 
#   group_by(NOTMATCH = (GA.x != GA.x)) %>% 
#   count(NOTMATCH)
# 
# merged <- read_dta(paste0(base_path,"20181028-growth_visit.dta"))
# merged2 <- read_dta(paste0(base_path,"20181128-growth_visit_add.dta"))
# 
# merged %>% filter(visit == 0) %>% count(!is.na(GA), !is.na(mother_age_delivery), !is.na(mother_ethnicity))
# merged2 %>% filter(visit == 0) %>% count(!is.na(GA), !is.na(mother_age_delivery), !is.na(mother_ethnicity))
# merged3 %>% filter(visit == 0) %>% count(!is.na(GA), !is.na(mother_age_delivery), !is.na(mother_ethnicity))