library(ggbeeswarm)
library(vtable)

base_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"
new_path <- "D:/SICS/Data and Instruments/Archive/"
merge_cpg <- read_dta(paste0(base_path,"20190125-merge_CpG.dta"))
merged <- read_dta(paste0(base_path,"20181219-IVF-final.dta"))
merged_MI <- read_dta(paste0(new_path,"20190125-IVF_Gen_MI_Drop.dta"))
merged_MI <- merged_MI %>% filter(`_mi_m` == 1)

### DAD DATA
dad_hx <- read_excel("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/99 - SUPPORT/Data QC/Paternal SES/Father_demographics_20180913_wo_DOB.xlsx")
dad_hx %>% count(father_occupr)
dad_hx_clean <- dad_hx %>% select(-c(father_education, father_other_education, father_education_flags,
                                     father_occupation, father_other_occupation, father_occupation_flags,
                                     father_income_flags)) %>% 
  mutate(father_education_corrected = as.factor(substring(father_education_corrected,1,1)),
         father_occupation_corrected = as.factor(substring(father_occupation_corrected,1,1)),
         father_monthly_income = as.factor(substring(father_monthly_income,1,1)),
         father_household_income = as.factor(substring(father_household_income,1,1)),
         father_HBP = as.double(substring(father_HBP,1,1)),
         father_HBP_age = as.double(father_HBP_age),
         father_HBP_year = as.double(father_HBP_year),
         father_diabetes = as.double(substring(father_diabetes,1,1)),
         father_diabetes_age = as.double(father_diabetes_age),
         father_diabetes_year = as.double(father_diabetes_year),
         father_current_smoker = as.double(substring(father_current_smoker,1,1)),
         father_cigarette_sticks = as.double(father_cigarette_sticks),
         father_visit = as.double(substring(father_visit,2,3)))
vtable(dad_hx_clean, missing = T)
write_dta(dad_hx_clean, "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/99 - SUPPORT/Data QC/Paternal SES/dad_hx.dta")

dad_hx_clean %>% mutate(IVF = case_when(substring(SubjectID,1,3) %in% c("019", "029") ~ "IVF", TRUE ~ "SC")) %>% 
  group_by(IVF) %>% summarize(n = sum(!is.na(father_age_delivery)), 
                                      mu_age = mean(father_age_delivery, na.rm = T))

merged_dad <- left_join(merged, dad_hx_clean)
merge_cpg_dad <- left_join(merge_cpg, dad_hx_clean)
write_dta(merge_cpg_dad, paste0(base_path,"20190723-merge_CpG_dad.dta"))
merged_MI_dad <- left_join(merged_MI, dad_hx_clean)

merged_dad %>% ggplot() + 
  geom_density(aes(x = father_age_delivery, fill = as.factor(IVF)), alpha = 0.4)

sub_org  <- merged_dad %>% filter(SubjectID %in% sub_id | IVF == 1) %>% mutate(subset = 0)
sub_alt  <- merged_dad %>% filter(SubjectID %in% alt_sub | IVF == 1) %>% mutate(subset = 1)
sub_one  <- merged_dad %>% filter(SubjectID %in% one_cond | IVF == 1) %>% mutate(subset = 2)
sub_inv  <- merged_dad %>% filter(SubjectID %in% rand_sub | IVF == 1) %>% mutate(subset = 3)

all_sub <- bind_rows(sub_org, sub_alt, sub_one, sub_inv)
  
# check dad age by subgroup
all_sub %>% filter(visit == 0) %>% 
  mutate(IVF = factor(IVF, levels = c(1,0), labels = c("IVF","SC")),
         subset = factor(subset, levels = c(0,1,2,3), 
                         labels = c("A. Original", 
                                    "B. Original, <2 miscarriage",
                                    "C. 1 condition only",
                                    "D. No conditions"))) %>% 
  ggplot(aes(x = mother_age_delivery, 
             y = father_age_delivery, 
             color = as.factor(IVF))) + 
  geom_point() + 
  geom_smooth(method = glm) + 
  facet_wrap(~subset)

# height x age -- despite older age, no different trend by IVF status
merged_MI_dad %>% filter(visit == 72) %>% 
  mutate(IVF = factor(IVF, levels = c(1,0), labels = c("IVF","SC"))) %>% 
  ggplot(aes(x = father_age_delivery, y = f_height_m24, color = IVF)) + geom_point() + geom_smooth(method = glm)

merged_MI_dad %>% filter(visit == 72) %>% 
  group_by(IVF) %>% 
  summarize(age = mean(father_age_delivery, na.rm = T),
            edu = mean(as.double(father_education_corrected), na.rm = T),
            job = mean(as.double(father_occupation_corrected), na.rm = T),
            inc = mean(as.double(father_monthly_income), na.rm = T),
            hh_inc = mean(as.double(father_household_income), na.rm = T),
            HBP_pct = mean(father_HBP, na.rm = T),
            DM_pct = mean(father_diabetes, na.rm = T),
            smk_pct = mean(father_current_smoker, na.rm = T),
            cigs = mean(father_cigarette_sticks, na.rm = T))
  
merged_MI_dad %>% filter(visit == 72) %>% 
  lm(formula = zlen ~ IVF + mother_age_delivery + f_height_m24 + father_age_delivery + 
       father_HBP + father_diabetes) %>% 
  summary()

#### dad HTN
merged_MI_dad %>% filter(!is.na(father_HBP)) %>% 
  mutate(father_HBP = as.logical(father_HBP)) %>% 
  ggplot(aes(x = days, y = zlen, color = father_HBP)) + 
  geom_point() + geom_smooth() +
  labs(color = "Paternal HTN history") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_blank())
  
merged_MI_dad %>% filter(!is.na(father_HBP)) %>% 
  mutate(father_HBP = as.logical(father_HBP)) %>% 
  mutate(IVF = as.logical(IVF)) %>% 
  ggplot(aes(x = days, y = zlen, color = IVF)) + 
  geom_point() + geom_smooth() +
  labs(color = "IVF status") +
  facet_wrap(~father_HBP) +
  theme(legend.position = c(0.8, 0.2), legend.background = element_blank())


#### dad DM
merged_MI_dad %>% filter(!is.na(father_diabetes)) %>% 
  mutate(father_diabetes = as.logical(father_diabetes)) %>% 
  ggplot(aes(x = days, y = zlen, color = father_diabetes)) + 
  geom_point() + geom_smooth() +
  labs(color = "Paternal DM history") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_blank())


merged_MI_dad %>% filter(!is.na(father_diabetes)) %>% 
  mutate(father_diabetes = as.logical(father_diabetes)) %>% 
  mutate(IVF = as.logical(IVF)) %>% 
  ggplot(aes(x = days, y = zlen, color = IVF)) + 
  geom_point() + geom_smooth() +
  labs(color = "IVF status") +
  facet_wrap(~father_diabetes) +
  theme(legend.position = c(0.8, 0.2), legend.background = element_blank())

### DAD WEIGHT
merged_dad %>% filter(visit == 0) %>% count(!is.na(f_weight_m24), !is.na(f_weight_m36))
merged_MI_dad %>% mutate(DAD_BMI = f_weight_m24/(f_height_m24/100)^2) %>% 
  #filter(SubjectID %in% c(dad_id, sub_id) | IVF == 1) %>% 
  filter(SubjectID %in% inv_sub | IVF == 1) %>% 
  group_by(IVF) %>% summarize(mean(f_weight_m24/(f_height_m24/100)^2, na.rm = T))
  
### pull med hx (age > 35, HTN, DM)
dad_id <- merged_MI_dad %>% filter(visit == 0, IVF == 0) %>% 
  filter(father_age_delivery >= 35 | father_HBP == 1 | father_diabetes == 1) %>% pull(SubjectID)

## PLOT ZLEN INCLUDING DADS WITH HTN/DM Hx
merged_MI_dad %>% filter(SubjectID %in% c(dad_id, sub_id) | IVF == 1) %>% 
  mutate(DAD_BMI = f_weight_m24/(f_height_m24/100)^2) %>% 
  filter(DAD_BMI >= 25) %>% 
  mutate(IVF = as.logical(IVF)) %>% group_by(IVF) %>% 
  ggplot(aes(x = days, y = zlen, color = IVF)) + geom_point() + geom_smooth()

merged_MI_dad %>% filter(SubjectID %in% sub_id | IVF == 1) %>% 
  mutate(IVF = as.logical(IVF)) %>%
  filter(visit == 0) %>% 
  #group_by(IVF) %>% summarize(mean(f_weight_m24/(f_height_m24/100)^2, na.rm = T))
  ggplot(aes(x = f_weight_m24/(f_height_m24/100)^2, fill = IVF)) + 
  geom_density(alpha = 0.2)

merged_MI_dad %>% filter(SubjectID %in% sub_id | IVF == 1) %>% 
  filter(visit == 0) %>% 
  group_by(IVF) %>% 
  summarize(age = mean(father_age_delivery, na.rm = T),
            edu = mean(as.double(father_education_corrected), na.rm = T),
            job = mean(as.double(father_occupation_corrected), na.rm = T),
            inc = mean(as.double(father_monthly_income), na.rm = T),
            hh_inc = mean(as.double(father_household_income), na.rm = T),
            HBP_pct = mean(father_HBP, na.rm = T),
            DM_pct = mean(father_diabetes, na.rm = T),
            smk_pct = mean(father_current_smoker, na.rm = T),
            cigs = mean(father_cigarette_sticks, na.rm = T)) 

merged_MI_dad %>% #filter(SubjectID %in% sub_id | IVF == 1) %>% 
  filter(visit == 72) %>% 
  mutate(DAD_BMI = f_weight_m24/(f_height_m24/100)^2,
         BIAS_DOWN = case_when(IVF == 1 ~ (f_weight_m24-25)/(f_height_m24/100)^2,
                               TRUE ~ (f_weight_m24)/(f_height_m24/100)^2)) %>% 
  #group_by(IVF) %>% summarize(mean(BIAS_DOWN, na.rm = T))
  #ggplot() + geom_density(aes(x = BIAS_DOWN, fill = as.factor(IVF)), alpha = 0.2)
  # lm(formula = zlen ~ IVF + mother_age_delivery + father_age_delivery + 
  #      ppBMI + ogtt_fasting_pw26 + ogtt_2hour_pw26 + mother_ethnicity + parity + 
  #      hi_bp + f_weight_m24 + f_height_m24 + SCORE) %>% summary()
  lm(formula = zwei ~ IVF + mother_age_delivery + father_age_delivery + 
       ppBMI + ogtt_fasting_pw26 + ogtt_2hour_pw26 + mother_ethnicity + parity + 
       hi_bp + f_height_m24 + BIAS_DOWN + SCORE) %>% summary()

############ 
### INCLUDE DAD DATA REGRESSION RESULTS
############
### ANTHRO
reg_out_dad <- read_excel(paste0(out_path,"23 Jul 2019-IVF-dad_MI.xls"), sheet = "min")
out_dad <- reg_out_dad %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei"))
all_visit = rep(c(0, 1, 3, 6, 9, 12, 15, 18, 24, 36, 48, 54, 60, 66, 72, 78), 6)
out_dad <- out_dad %>% add_column(all_visit)
out_dad %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  group_by(outcome) %>% summarize(min = min(N), max = max(N))
out_dad %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  mutate(outcome = case_when(outcome == "height" ~ "A.~Height~(cm)",
                             outcome == "log_weight" ~ "B.~Weight~(log(kg))",
                             outcome == "log_bmi" ~ "C.~BMI~(log(kg/m^2))",
                             outcome == "zlen" ~ "D.~'Height-for-age'~'Z-score'~(SDS)",
                             outcome == "zwei" ~ "E.~'Weight-for-age'~'Z-score'~(SDS)",
                             outcome == "zbmi" ~ "F.~BMI~'Z-score'~(SDS)")) %>% 
  ggplot() + 
  geom_point(aes(x = all_visit, y = estimate)) +
  geom_errorbar(aes(x = all_visit, ymin = min95, ymax= max95)) + 
  geom_ribbon(aes(x = all_visit, ymin = min95, ymax= max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Age (Months)", y = "") +
  theme(legend.position = "none") +
  facet_wrap(~outcome, scale = "free", labeller = label_parsed)

### BP AND SKINFOLDS
BP_DAD <- reg_out_dad %>% filter(outcome %in% c("SBP", "ZSBP", "DBP", "ZDBP")) %>% 
  add_column(BP_visit = rep(c(36, 48, 60, 72), 4)) %>% 
  mutate(outcome = case_when(outcome == "SBP" ~ "A. Systolic, absolute (mmHg)", 
                             outcome == "DBP" ~ "B. Diastolic, absolute (mmHg)",
                             outcome == "ZSBP" ~ "C. Systolic, standardized (%ile rank)",
                             outcome == "ZDBP" ~ "D. Diastolic, standardized (%ile rank)")) %>% 
  ggplot() + geom_point(aes(x = BP_visit, y = estimate)) +
  geom_errorbar(aes(x = BP_visit, ymin = min95, ymax= max95), width = 1) + 
  geom_ribbon(aes(x = BP_visit, ymin = min95, ymax= max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(limits = c(-20, 15)) +
  scale_x_continuous(breaks = seq(36,72, 12)) +
  facet_wrap(~ outcome) +
  theme(legend.position = "none") +
  labs(x = "Age (Months)", y = "Absolute (mmHg) or standardized (%ile rank) 
       difference, (IVF - SC)")

SF_DAD <- reg_out_dad %>% filter(outcome %in% c("log_sub", "log_tri", "log_bi", "log_supra")) %>% 
  add_column(skin_visit = c(0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            48, 54, 60, 66, 72, 78)) %>% 
  mutate_at(c("estimate","min95","max95"),funs((exp(.)-1)*100)) %>% 
  mutate(outcome = case_when(outcome == "log_tri" ~ "A. Triceps (N = 800 to 1119)",
                             outcome == "log_sub" ~ "B. Subscapular (N = 765 to 1118)",
                             outcome == "log_bi" ~ "C. Biceps (N = 765 to 892)",
                             outcome == "log_supra" ~ "D. Suprailiac (N = 785 to 862)")) %>% 
  ggplot() + 
  geom_point(aes(x = skin_visit, y = estimate)) +
  geom_errorbar(aes(x = skin_visit, ymin = min95, ymax= max95)) + 
  geom_ribbon(aes(x = skin_visit, ymin = min95, ymax= max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-30,10)) +
  scale_x_continuous(breaks = c(0, 18, 24, 36, 48, 54, 60, 66, 72, 78))+
  theme(legend.position = "none") +
  labs(x = "", y = "% difference in skinfold thickness
       (IVF - SC)") +
  facet_wrap(~outcome)

ggarrange(SF_DAD, BP_DAD, labels = c("A", "B"), ncol = 1, nrow = 2)
