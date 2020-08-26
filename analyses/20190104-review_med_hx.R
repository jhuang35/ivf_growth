library(tidyverse)
library(ggplot2)
library(haven)
library(summarytools)
library(readxl)
library(broom)

options(tibble.print_max = 80, tibble.print_min = 30)

plot_path <- "C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/"
base_path <- "C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

#plot_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/"
#base_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

merged <- read_dta(paste0(base_path,"20181219-IVF-final.dta"))

# Additional medical history data
add_data <- read_xlsx("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Clean dataset/forma428&429_20181224.xlsx", na = c("N/A", "NA"," ", ""))

add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% summarize_all(funs(sum(!is.na(.)))) %>% filter_all(funs(. != 0))
  
add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% filter(!is.na(medications_before_preg_pw11)) %>% 
  count(medications_before_preg_pw11) %>% arrange(-n)

# list of drugs for IVF
add_data %>% filter(participant_status != "ineligible") %>% 
  mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  filter(IVF == T) %>% 
  count(medications_before_preg_pw11, medication_other_pw11) %>% arrange(-n)

# list of drugs for non-IVF
add_data %>% #filter(participant_status != "ineligible") %>% 
  mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  filter(IVF == F) %>% 
  count(medications_before_preg_pw11, medication_other_pw11) %>% arrange(-n)

subfert <- add_data %>% filter(medication_other_pw11 %in% c("thyroxine", "thyroid, unspecified", "duromine", "progesterone", "aspirin",
                                                            "aspirin; hydroxychloroquine; prednisolone", "cabergoline", "carbimazole", "dexamethasone",
                                                            "clomifene", "cyproterone/ethinyl estradiol", "desogestrel/ethinyl estradiol; duromine", 
                                                            "dydrogesterone", "hormone pill, unspecified", "prednisolone; hydroxychloroquine", 
                                                            "propylthiouracil", "thiamazole")) %>% pull(SubjectID)

MISC <- add_data %>% 
  mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  filter(IVF == F) %>% 
  mutate(num_miss = 7 - (pw11_pregnancy_outcome_1 != "miscarriage" | is.na(pw11_pregnancy_outcome_1)) -
           (pw11_pregnancy_outcome_2 != "miscarriage" | is.na(pw11_pregnancy_outcome_2)) -
           (pw11_pregnancy_outcome_3 != "miscarriage" | is.na(pw11_pregnancy_outcome_3)) -
           (pw11_pregnancy_outcome_4 != "miscarriage" | is.na(pw11_pregnancy_outcome_4)) -
           (pw11_pregnancy_outcome_5 != "miscarriage" | is.na(pw11_pregnancy_outcome_5)) -
           (pw11_pregnancy_outcome_6 != "miscarriage" | is.na(pw11_pregnancy_outcome_6)) -
           (pw11_pregnancy_outcome_7 != "miscarriage" | is.na(pw11_pregnancy_outcome_7))) %>% 
  filter(num_miss >= 2) %>% pull(SubjectID)

mom_age <- merged %>% filter(visit == 0) %>% select(SubjectID, IVF, mother_age_delivery)

PCOS <- add_data %>% left_join(., mom_age) %>% filter((pcos_y4 == "1_yes") & (pcos_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
ENDO <- add_data %>% left_join(., mom_age) %>% filter((endometriosis_y4 == "1_yes") & (endometriosis_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
CYST <- add_data %>% left_join(., mom_age) %>% filter((ovariancyst_y4 == "1_yes") & (ovariancyst_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
FIBR <- add_data %>% left_join(., mom_age) %>% filter((uterinefibroid_y4 == "1_yes") & (uterinefibroid_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
THYR <- add_data %>% left_join(., mom_age) %>% filter( ((other_illness1_y4 == "Thyroid, Hyper" |
                                                           other_illness1_y4 == "Thyroid, Hypo" | 
                                                           other_illness1_y4 == "Thyroid, Unspecified") & 
                                                          (other_illness1_age_onset_y4 < mother_age_delivery)) |
    ((other_illness2_y4 == "Thyroid, Hypo") & (other_illness2_age_onset_y4 < mother_age_delivery)) ) %>% 
  filter(IVF == F) %>% pull(SubjectID)

# medications during pregnancy -- will not include; stick to pre-pregnancy conditions and medications
# add_data %>% filter(!(substring(SubjectID, 1,3) %in% c("019","029"))) %>% count(medication1_pw26)
# add_data %>% filter(!(substring(SubjectID, 1,3) %in% c("019","029"))) %>% count(medication2_pw26)
# add_data %>% filter(!(substring(SubjectID, 1,3) %in% c("019","029"))) %>% count(medication3_pw26)
# add_data %>% filter(!(substring(SubjectID, 1,3) %in% c("019","029"))) %>% count(medication4_pw26)



# look at match on characteristics
add_data %>% 
  mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  filter(IVF == T | SubjectID %in% subfert | 
           SubjectID %in% CYST | 
           SubjectID %in% ENDO | 
           SubjectID %in% FIBR |
           SubjectID %in% PCOS |
           SubjectID %in% MISC) %>% 
  group_by(IVF) %>% count(alcohol_consumption_preg)

# save the subsetted data
# not worthing including one additional individual with "Thyroid, Unspecified" dx per other_illness1_y4, but no medication record (010-21914)
subset <- merged %>% filter(IVF == T | 
                              SubjectID %in% subfert | 
                              SubjectID %in% CYST | 
                              SubjectID %in% ENDO | 
                              SubjectID %in% FIBR |
                              SubjectID %in% PCOS |
                              SubjectID %in% MISC)
subset_ID_2 <- subset %>% filter(visit == 0) %>% pull(SubjectID)
write_dta(subset, paste0(base_path,"20190107-IVF-final-SUB.dta"))

MI <- read_dta(paste0(base_path, "20181206-IVF_MI_Drop.dta"))
subset_MI <- MI %>% filter(IVF == T | 
                             SubjectID %in% subfert | 
                             SubjectID %in% CYST | 
                             SubjectID %in% ENDO | 
                             SubjectID %in% FIBR |
                             SubjectID %in% PCOS |
                             SubjectID %in% MISC)
write_dta(subset_MI, paste0(base_path, "20190107-IVF_MI_Drop-SUB.dta"))





## PLOTS
subset %>% filter(visit == 0, sex != "") %>% 
  select(contains("EFW"), IVF) %>% 
  gather(outcome, kg, -IVF) %>% 
  separate(outcome, c("outcome", "visit"), remove = T) %>% 
  ggplot() + 
  geom_jitter(aes(x = visit, y = kg, color = as.factor(IVF)), position = position_jitterdodge()) + 
  geom_boxplot(aes(x = visit, y = kg, color = as.factor(IVF)), alpha = 0.6) + 
  theme(legend.position = "bottom")

subset %>% ggplot(aes(x = days, y = height, color = as.factor(IVF))) + geom_point() + geom_smooth() + theme(legend.position = "bottom")
subset %>% ggplot(aes(x = days, y = zlen, color = as.factor(IVF))) + geom_point() + geom_smooth() + theme(legend.position = "bottom")
subset %>% ggplot(aes(x = days, y = weight, color = as.factor(IVF))) + geom_point() + geom_smooth() + theme(legend.position = "bottom")
subset %>% ggplot(aes(x = days, y = zwei, color = as.factor(IVF))) + geom_point() + geom_smooth() + theme(legend.position = "bottom")
subset %>% ggplot(aes(x = days, y = zwfl, color = as.factor(IVF))) + geom_point() + geom_smooth() + theme(legend.position = "bottom")
subset %>% ggplot(aes(x = days, y = zbmi, color = as.factor(IVF))) + geom_point() + geom_smooth() + theme(legend.position = "bottom")


# ANTHRO REGRESSIONS
#reg_out_sub <- read_excel("C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/ 8 Jan 2019-IVF-size_MI_SUB.xls", sheet = "adjusted")
reg_out_sub <- read_excel("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/ 8 Jan 2019-IVF-size_MI_SUB.xls", sheet = "adjusted")
out_all_visit <- reg_out_sub %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei"))
all_visit = rep(c(0, 1, 3, 6, 9, 12, 15, 18, 24, 36, 48, 54, 60, 66, 72, 78), 6)

out_all_visit %>% add_column(all_visit) %>% select(outcome, all_visit, N) %>% spread(outcome, N)
out_all_visit %>% group_by(outcome) %>% summarize(min = min(N), max = max(N))

out_all_visit %>% add_column(all_visit) %>% 
  filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  mutate(outcome = if_else(outcome == "zlen", "D. Height-for-age (SDS)",
                           if_else(outcome == "log_weight", "B. Weight (log(kg))",
                                   if_else(outcome == "height", "A. Height (cm)",
                                           if_else(outcome == "log_bmi", "C. BMI (log(kg/m^2))",
                                                   if_else(outcome == "zbmi", "F. BMI Z-score (SD)", 
                                                           if_else(outcome == "zwei", "E. Weight-for-age (SDS)", ""))))))) %>% 
  ggplot() + 
  geom_point(aes(x = all_visit, y = estimate)) +
  geom_errorbar(aes(x = all_visit, ymin = min95, ymax= max95)) + 
  geom_ribbon(aes(x = all_visit, ymin = min95, ymax= max95), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Associations between IVF status and anthropometry (83 IVF vs. 93 potentially subfertile subcohort).",
       x = "Age in Months",
       y = "") +
  facet_wrap(~outcome, scale = "free", labeller = label_bquote(.(outcome)))

# SKINFOLD REGRESSIONS
skinfold <- reg_out_sub %>% filter(outcome %in% c("log_sub", "log_tri", "log_bi", "log_supra"))
skin_visit <- c(0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                18, 24, 36, 48, 54, 60, 66, 72, 78, 
                0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                48, 54, 60, 66, 72, 78)
skinfold <- skinfold %>% add_column(skin_visit)
skinfold %>% group_by(outcome) %>% count(skin_visit, N) # Confirm against N for each measure / visit
skinfold %>% group_by(outcome) %>% summarize(min = min(N, na.rm = T), max = max(N, na.rm = T))

skinfold %>% mutate_at(c("estimate","min95","max95"),funs((exp(.)-1)*100)) %>% 
  mutate(outcome = if_else(outcome == "log_tri", "A. Triceps (N = 117 to 172)",
                           if_else(outcome == "log_sub", "B. Subscapular (N = 116 to 172)",
                                   if_else(outcome == "log_bi", "C. Biceps (N = 112 to 145)",
                                           if_else(outcome == "log_supra", "D. Suprailiac (N = 126 to 140)", ""))))) %>% 
  ggplot() + 
  geom_point(aes(x = skin_visit, y = estimate)) +
  geom_errorbar(aes(x = skin_visit, ymin = min95, ymax= max95)) + 
  geom_ribbon(aes(x = skin_visit, ymin = min95, ymax= max95), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Associations between IVF status and skinfold thickness, by location and visit (subcohort).",
       x = "Age in Months",
       y = "% difference in skinfold thickness") +
  facet_wrap(~outcome, labeller = label_bquote(.(outcome)))

# BLOOD PRESSURE REGRESSIONS
BP_visit = rep(c(36, 48, 60, 72), 4)
subset %>% group_by(visit) %>% summarize(N = sum(!is.na(SBP))) 
reg_out_sub %>% filter(outcome %in% c("SBP", "ZSBP", "DBP", "ZDBP")) %>% add_column(BP_visit) %>% 
  mutate(outcome = if_else(outcome == "SBP", "A. Systolic, absolute (mmHg)", 
                           if_else(outcome == "DBP", "B. Diastolic, absolute (mmHg)",
                                   if_else(outcome == "ZSBP", "C. Systolic, standardized (%ile rank)",
                                           if_else(outcome == "ZDBP", "D. Diastolic, standardized (%ile rank)", outcome))))) %>% 
  ggplot() + 
  geom_point(aes(x = BP_visit, y = estimate)) +
  geom_errorbar(aes(x = BP_visit, ymin = min95, ymax= max95), width = 1) + 
  geom_ribbon(aes(x = BP_visit, ymin = min95, ymax= max95), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~outcome) +
  labs(title = "Associations between IVF status and blood pressure (subcohort; N = 118 to 129).",
       x = "Age in Months",
       y = "Difference in absolute (mmHg) or standardized (%ile rank) blood pressure")

# EFW
EFW_visit = c(19, 26, 32)
reg_out_sub %>% filter(outcome %in% c("EFW_19", "EFW_26", "EFW_32")) %>% add_column(EFW_visit) %>% 
  ggplot() + 
  geom_point(aes(x = EFW_visit, y = estimate)) +
  geom_errorbar(aes(x = EFW_visit, ymin = min95, ymax= max95), width = 1) + 
  geom_ribbon(aes(x = EFW_visit, ymin = min95, ymax= max95), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Gestational Age (Weeks)",
       y = "Difference in weight (kg)",
       title = "Estimated Fetal Weight, By Ultrasound Visit Date.")


# alcohol usage
add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% select(IVF, contains("wine"), contains("beer")) %>% names()

add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% select(IVF, contains("wine")) %>% summarize_all(funs(sum(., na.rm = T)/n())) %>% 
  gather(var, val, -IVF)

add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% select(IVF, contains("beer")) %>% summarize_all(funs(sum(., na.rm = T)/n())) %>% 
  gather(var, val, -IVF)

add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% count(beer_times_per_month_pp, beer_times_per_month_preg)

add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% count(trad_wine_times_per_month_pp, trad_wine_times_per_month_preg)

add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% count(wine_times_per_month_pp, wine_times_per_month_preg)

add_data %>% mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029")) %>% 
  group_by(IVF) %>% count(wine_times_per_month_preg, trad_wine_times_per_month_preg, beer_times_per_month_preg)

# child health
health_03 <- read_csv("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Raw LORIS data/GUSTO child health/Child_Questionnaire_Month3.csv")
