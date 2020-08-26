###################################################
# Load maternal medical history data
# Form original and alternative subcohorts
# Edited: 10 July 2019
###################################################

library(tidyverse)
library(ggplot2)
library(haven)
library(summarytools)
library(readxl)
library(broom)

options(tibble.print_max = 80, tibble.print_min = 30)

# plot_path <- "C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/"
# base_path <- "C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

plot_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/"
base_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

merged <- read_dta(paste0(base_path,"20181219-IVF-final.dta"))

# Additional medical history data
add_data <- read_xlsx("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Clean dataset/forma428&429_20181224.xlsx", na = c("N/A", "NA"," ", "")) %>% 
  mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029"))

add_data %>% group_by(IVF) %>% summarize_all(funs(sum(!is.na(.)))) %>% filter_all(funs(. != 0))

add_data %>% group_by(IVF) %>% filter(!is.na(medications_before_preg_pw11)) %>% 
  count(medications_before_preg_pw11) %>% arrange(-n)

# list of drugs for IVF
add_data %>% filter(participant_status != "ineligible") %>% 
  filter(IVF == T) %>% 
  count(medications_before_preg_pw11, medication_other_pw11) %>% arrange(-n)

# list of drugs for non-IVF
add_data %>% #filter(participant_status != "ineligible") %>% 
  filter(IVF == F) %>% 
  count(medications_before_preg_pw11, medication_other_pw11) %>% arrange(-n)

# n = 33 SC (35 total)
subfert <- add_data %>%  
  #filter(IVF == F) %>% 
  filter(medication_other_pw11 %in% c("thyroxine", "thyroid, unspecified", "duromine", "progesterone", "aspirin",
                                      "aspirin; hydroxychloroquine; prednisolone", "cabergoline", "carbimazole", "dexamethasone",
                                      "clomifene", "cyproterone/ethinyl estradiol", "desogestrel/ethinyl estradiol; duromine", 
                                      "dydrogesterone", "hormone pill, unspecified", "prednisolone; hydroxychloroquine", 
                                      "propylthiouracil", "thiamazole")) %>% pull(SubjectID)

# n = 14 SC (14 TOTAL)
add_data <- add_data %>%
  mutate(num_miss = 7 - 
           (pw11_pregnancy_outcome_1 != "miscarriage" | is.na(pw11_pregnancy_outcome_1)) -
           (pw11_pregnancy_outcome_2 != "miscarriage" | is.na(pw11_pregnancy_outcome_2)) -
           (pw11_pregnancy_outcome_3 != "miscarriage" | is.na(pw11_pregnancy_outcome_3)) -
           (pw11_pregnancy_outcome_4 != "miscarriage" | is.na(pw11_pregnancy_outcome_4)) -
           (pw11_pregnancy_outcome_5 != "miscarriage" | is.na(pw11_pregnancy_outcome_5)) -
           (pw11_pregnancy_outcome_6 != "miscarriage" | is.na(pw11_pregnancy_outcome_6)) -
           (pw11_pregnancy_outcome_7 != "miscarriage" | is.na(pw11_pregnancy_outcome_7)))
MISC <- add_data %>% 
  #filter(IVF == F) %>% 
  filter(num_miss >= 2) %>% pull(SubjectID)

mom_age <- merged %>% filter(visit == 0) %>% dplyr::select(SubjectID, mother_age_delivery)

# n = 16 PCOS (24 total)
#PCOS <- add_data %>% left_join(., mom_age) %>% filter((pcos_y4 == "1_yes") & (pcos_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
PCOS <- add_data %>% left_join(., mom_age) %>% filter((pcos_y4 == "1_yes") & (pcos_age_y4 < mother_age_delivery)) %>% pull(SubjectID)

# n = 10 ENDO (12 total)
#ENDO <- add_data %>% left_join(., mom_age) %>% filter((endometriosis_y4 == "1_yes") & (endometriosis_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
ENDO <- add_data %>% left_join(., mom_age) %>% filter((endometriosis_y4 == "1_yes") & (endometriosis_age_y4 < mother_age_delivery)) %>% pull(SubjectID)

# n = 19 CYST (24 total)
#CYST <- add_data %>% left_join(., mom_age) %>% filter((ovariancyst_y4 == "1_yes") & (ovariancyst_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
CYST <- add_data %>% left_join(., mom_age) %>% filter((ovariancyst_y4 == "1_yes") & (ovariancyst_age_y4 < mother_age_delivery)) %>% pull(SubjectID)

# n = 23 FIBROIDS (27 total)
#FIBR <- add_data %>% left_join(., mom_age) %>% filter((uterinefibroid_y4 == "1_yes") & (uterinefibroid_age_y4 < mother_age_delivery) & (IVF == F)) %>% pull(SubjectID)
FIBR <- add_data %>% left_join(., mom_age) %>% filter((uterinefibroid_y4 == "1_yes") & (uterinefibroid_age_y4 < mother_age_delivery)) %>% pull(SubjectID)

# n = 12 THYROID (13 total)
THYR <- add_data %>% left_join(., mom_age) %>% 
  filter( ((other_illness1_y4 == "Thyroid, Hyper" |
              other_illness1_y4 == "Thyroid, Hypo" | 
              other_illness1_y4 == "Thyroid, Unspecified") & 
             (other_illness1_age_onset_y4 < mother_age_delivery)) |
            ((other_illness2_y4 == "Thyroid, Hypo") & (other_illness2_age_onset_y4 < mother_age_delivery)) ) %>% 
  #filter(IVF == F) %>% 
  pull(SubjectID)

# n = 104 total
n_distinct(c(subfert, MISC, PCOS, ENDO, CYST, FIBR, THYR))

# sub_id_new <- merged %>% filter(visit == 0, !is.na(sex)) %>% 
#   dplyr::filter(SubjectID %in% c(subfert, MISC, PCOS, ENDO, CYST, FIBR, THYR) | IVF == 1) %>% 
#   group_by(IVF) %>% pull(SubjectID)

add_data <- add_data %>% 
  mutate(n_cond = (SubjectID %in% subfert) +
           (SubjectID %in% MISC) +
           (SubjectID %in% PCOS) +
           (SubjectID %in% ENDO) +
           (SubjectID %in% CYST) +
           (SubjectID %in% FIBR) +
           (SubjectID %in% THYR))

add_data %>% count(IVF, n_cond)
add_data %>% group_by(IVF) %>% 
  filter(SubjectID %in% sub_id) %>% 
  ggplot() + geom_density(aes(x = num_miss, fill = IVF), alpha = 0.4)
#ggplot() + geom_density(aes(x = n_cond, fill = IVF), alpha = 0.4)  
#summarize(mean(n_cond)) %>% 
#count(n_cond)
merge_cond <- add_data %>% select(SubjectID, n_cond, num_miss) %>% left_join(merged, .)
merge_cond %>% #filter(n_cond == 1) %>% 
  ggplot(aes(x = days, y = zlen, color = as.factor(IVF))) + 
  geom_point() + geom_smooth() +
  facet_wrap(~num_miss)
write_dta(merge_cond, paste0(base_path, "20190722-mom_hx_data.dta"))

one_cond <- add_data %>% filter(n_cond == 1 | IVF == TRUE) %>% pull(SubjectID)
alt_sub <- add_data %>% filter((num_miss == 1 | IVF == TRUE | SubjectID %in% sub_id) &
                                 num_miss < 2) %>% pull(SubjectID)
add_data %>% filter(SubjectID %in% alt_sub) %>% ggplot() + geom_density(aes(x = n_cond, fill = IVF), alpha = 0.4) 

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
# not worth including one additional individual with "Thyroid, Unspecified" dx per other_illness1_y4, but no medication record (010-21914)
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

### ALTERNATE SUBGROUPINGS
merged3 %>% filter(visit == 0, SubjectID %in% sub_id, !is.na(zlen)) %>% count(IVF)
merged3 %>% filter(visit == 0, SubjectID %in% alt_sub, !is.na(zlen)) %>% count(IVF)
merged3 %>% filter(visit == 0, SubjectID %in% one_cond, !is.na(zlen)) %>% count(IVF)
merged3 %>% filter(visit == 0, IVF == 1 | !(SubjectID %in% sub_id), !is.na(zlen)) %>% count(IVF)

inv_sub <- merged3 %>% filter(visit == 0, !(SubjectID %in% sub_id)) %>% pull(SubjectID)

set.seed(42782)
rand_sub <- tibble(sub_id = inv_sub, choose = rbinom(length(inv_sub), 1, 0.2)) %>% filter(choose == 1) %>% pull(sub_id)
merged3 %>% filter(visit == 0, !is.na(zlen), IVF == 1 | SubjectID %in% rand_sub) %>% count(IVF)

model_sub <- merged3 %>% select(days, zlen, IVF, SubjectID) %>% filter(SubjectID %in% sub_id) %>% mutate(model = "Original (N SC = 93)")
model_alt <- merged3 %>% select(days, zlen, IVF, SubjectID) %>% filter(SubjectID %in% alt_sub) %>% mutate(model = "One Miscarriage (N SC = 159)")
model_one <- merged3 %>% select(days, zlen, IVF, SubjectID) %>% filter(SubjectID %in% one_cond) %>% mutate(model = "One Indication Only (N SC = 76)")
model_inv <- merged3 %>% select(days, zlen, IVF, SubjectID) %>% filter(IVF == 1 | !(SubjectID %in% sub_id)) %>% mutate(model = "No Indications (N SC = 998")
model_rand_inv <- merged3 %>% select(days, zlen, IVF, SubjectID) %>% filter(IVF == 1 | (SubjectID %in% rand_sub)) %>% mutate(model = "No Indications (20% random sample; N SC = 204)")

all_models <- bind_rows(model_sub, model_alt, model_one, model_rand_inv)

all_models %>% mutate(IVF = factor(IVF, levels = c(1,0), labels = c("IVF", "SC"))) %>% 
  ggplot(aes(x = days, y = zlen, color = IVF), alpha = 0.4) +
  geom_point() + geom_smooth() + facet_wrap(~model) +
  labs(y = "length-for-age Z-score (SDS)", x = "age (days)", color = "IVF status") +
  theme(legend.position = "bottom")

### VERIFY INFERENCES DO NOT DIFFER
#
# All_MI <- read_dta(paste0(base_path,"20190125-IVF_Gen_MI_Drop.dta")) 
# All_MI_1 <- All_MI %>% filter(`_mi_m` == 1)
# 
# All_MI_1 %>% filter(IVF == 1 | !(SubjectID %in% rand_sub)) %>% 
#   lmer(zlen~IVF*as.factor(visit)+mother_age_delivery+as.factor(MOM_EDU)+as.factor(ETHNIC)+
#          as.factor(HH_INC)+m_height_pw26+ppBMI+f_height_m24+f_weight_m24+MALE+SCORE+parity+smk_home+days+
#          (1|SubjectID), data = .)  %>% summary()
