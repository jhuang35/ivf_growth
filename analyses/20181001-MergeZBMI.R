library(tidyverse)
library(ggplot2)
library(haven)
library(summarytools)
library(readxl)
# library(rms)
# library(mgcv)

options(tibble.print_max = 40, tibble.print_min = 30)

plot_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/PLOTS/"
base_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

# plot_path <- "C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/PLOTS/Sample Size/"
# base_path <- "C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

joint_bmi <- read_csv(paste0(base_path,"20180920-GDM_growth_visit.csv"))

zbmi_00 <- read_dta(paste0(base_path,"Zscores/WHOout_birth_z_rc.dta"))
zbmi_00 <- zbmi_00  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(0)) %>% filter(!is.na(cbmi))

zbmi_01 <- read_dta(paste0(base_path,"Zscores/WHOout_wk3_z_rc.dta"))
zbmi_01 <- zbmi_01  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(1)) %>% filter(!is.na(cbmi))

zbmi_03 <- read_dta(paste0(base_path,"Zscores/WHOout_m3_z_rc.dta"))
zbmi_03 <- zbmi_03  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(3)) %>% filter(!is.na(cbmi))

zbmi_06 <- read_dta(paste0(base_path,"Zscores/WHOout_m6_z_rc.dta"))
zbmi_06 <- zbmi_06  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(6)) %>% filter(!is.na(cbmi))

zbmi_09 <- read_dta(paste0(base_path,"Zscores/WHOout_m9_z_rc.dta"))
zbmi_09 <- zbmi_09  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(9)) %>% filter(!is.na(cbmi))

zbmi_12 <- read_dta(paste0(base_path,"Zscores/WHOout_m12_z_rc.dta"))
zbmi_12 <- zbmi_12  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(12)) %>% filter(!is.na(cbmi))

zbmi_15 <- read_dta(paste0(base_path,"Zscores/WHOout_m15_z_rc.dta"))
zbmi_15 <- zbmi_15  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(15)) %>% filter(!is.na(cbmi))

zbmi_18 <- read_dta(paste0(base_path,"Zscores/WHOout_m18_z_rc.dta"))
zbmi_18 <- zbmi_18  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(18)) %>% filter(!is.na(cbmi))

zbmi_24 <- read_dta(paste0(base_path,"Zscores/WHOout_m24_z_rc.dta"))
zbmi_24 <- zbmi_24  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(24)) %>% filter(!is.na(cbmi))

zbmi_36 <- read_dta(paste0(base_path,"Zscores/WHOout_m36_z_rc.dta"))
zbmi_36 <- zbmi_36  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(36)) %>% filter(!is.na(cbmi))

zbmi_48 <- read_dta(paste0(base_path,"Zscores/WHOout_yr4_z_rc.dta"))
zbmi_48 <- zbmi_48  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(48)) %>% filter(!is.na(cbmi))

zbmi_54 <- read_dta(paste0(base_path,"Zscores/WHOout_yr4_5_z_rc.dta"))
zbmi_54 <- zbmi_54  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi") %>% mutate(visit = as.integer(54)) %>% filter(!is.na(cbmi))

zbmi_60 <- read_dta(paste0(base_path,"Zscores/WHOout_yr5_z_rc.dta"))
zbmi_60 <- zbmi_60  %>% select(SubjectID, zlen = "_zlen", zwei = "_zwei", zwfl = "_zwfl", cbmi = "_cbmi", zbmi = "_zbmi")
zbmi_60_2 <- read_dta(paste0(base_path,"Zscores/WHOout_yr5_z.dta"))
zbmi_60_2 <- zbmi_60_2 %>% select(SubjectID, zlen_new = "_zhfa", zwei_new = "_zwfa", zbmi_new = "_zbfa")
zbmi_60 <- left_join(zbmi_60, zbmi_60_2, by = "SubjectID") %>% 
  mutate(zbmi = if_else(is.na(zbmi), zbmi_new, zbmi),
         zlen = if_else(is.na(zlen), zlen_new, zlen),
         zwei = if_else(is.na(zwei), zwei_new, zwei)) %>% 
  mutate(visit = as.integer(60)) %>% 
  select(-contains("_new")) %>% filter(!is.na(cbmi))

zbmi_66 <- read_dta(paste0(base_path,"Zscores/WHOout_yr5_5_z.dta"))
zbmi_66 <- zbmi_66  %>% select(SubjectID, zlen = "_zhfa", zwei = "_zwfa", cbmi = "_cbmi", zbmi = "_zbfa") %>% mutate(visit = as.integer(66)) %>% filter(!is.na(cbmi))

zbmi_72 <- read_dta(paste0(base_path,"Zscores/WHOout_yr6_z.dta"))
zbmi_72 <- zbmi_72  %>% select(SubjectID, zlen = "_zhfa", zwei = "_zwfa", cbmi = "_cbmi", zbmi = "_zbfa") %>% mutate(visit = as.integer(72)) %>% filter(!is.na(cbmi))

zbmi_78 <- read_dta(paste0(base_path,"Zscores/WHOout_yr6_5_z.dta"))
zbmi_78 <- zbmi_78  %>% select(SubjectID, zlen = "_zhfa", zwei = "_zwfa", cbmi = "_cbmi", zbmi = "_zbfa") %>% mutate(visit = as.integer(78)) %>% filter(!is.na(cbmi))

all_zbmi <- full_join(zbmi_00,
            full_join(zbmi_01,
            full_join(zbmi_03,
            full_join(zbmi_06,
            full_join(zbmi_09,
            full_join(zbmi_12,
            full_join(zbmi_15,
            full_join(zbmi_18,
            full_join(zbmi_24,
            full_join(zbmi_36,
            full_join(zbmi_48,
            full_join(zbmi_54,
            full_join(zbmi_60,
            full_join(zbmi_66,
            full_join(zbmi_72, zbmi_78)))))))))))))))

joint_zbmi <- left_join(joint_bmi, all_zbmi)

write_dta(joint_zbmi, paste0(base_path,"20181001-zBMI.dta"))

###################################
# Plot ZBMI by growth classifications
###################################
zbmi_annotate <- read_dta(paste0(base_path,"20181001-zBMI.dta"))

zbmi_annotate %>% group_by(SubjectID) %>% 
  mutate(early_norm = sum(zbmi < 0.5 & visit <= 24, na.rm = T) >= 5, #max = 9 points
         late_grow = sum(zbmi > 1 & visit > 24, na.rm = T) >= 3, #max = 7 point
         acc = early_norm & late_grow,
         all_high = sum(zbmi > 1, na.rm = T) >= 8,
         over_at_78 = sum(zbmi > 1 & visit == 78, na.rm = T) == 1,
         class = (!all_high & !acc & over_at_78) + 
           2*(!all_high & acc) + 
           3*(all_high & acc & over_at_78) +
           4*(all_high & !acc)) %>%
  ungroup() %>% filter(class != 0) %>% 
  mutate(class = factor(class, levels = c(1,2,3,4), 
                        labels = c("~6 years (n = 46)", "~3 years (n = 79)", "~2 years (n = 9)", "~6 months (n = 68)"))) %>% 
  #filter(visit == 0) %>% count(all_high, acc, over_at_78, class)
  ggplot(aes(x = days, y = zbmi)) + geom_point(aes(color = class)) + #geom_line() +
  scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)) +
  scale_x_continuous(breaks = c(0, 183, 365, 731, 1096, 1461, 1826, 2192, 2557),
                     labels = c("Birth", "6 months", "1 year", "2 years", "3 years", 
                                "4 years", "5 years", "6 years", "7 years")) +
  labs(x = "Age", y = "BMI Z-Score", color = "Approximate age @ > 1 SD: ",
       title = "General growth patterns for children who are overweight in pre-school ages (N = 202)") +
  geom_smooth(aes(x = days, y = zbmi, color = class), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = -1, linetype = "dotted") +
  theme(legend.position = "bottom")

# BY ETHNICITY
zbmi_annotate %>% group_by(SubjectID) %>% 
  mutate(early_norm = sum(zbmi < 0.5 & visit <= 24, na.rm = T) >= 5, #max = 9 points
         late_grow = sum(zbmi > 1 & visit > 24, na.rm = T) >= 3, #max = 7 point
         acc = early_norm & late_grow,
         all_high = sum(zbmi > 1, na.rm = T) >= 8,
         over_at_78 = sum(zbmi > 1 & visit == 78, na.rm = T) == 1,
         class = (!all_high & !acc & over_at_78) + 
           2*(!all_high & acc) + 
           3*(all_high & acc & over_at_78) +
           4*(all_high & !acc)) %>%
  ungroup() %>% filter(class != 0) %>% 
  mutate(ethnicity = factor(mother_ethnicity, levels = c("chinese", "indian","malay"),
                            labels = c("Chinese", "Indian", "Malay"))) %>% 
  mutate(class = factor(class, levels = c(1,2,3,4), 
                        labels = c("Gradual", "Late acceleration", "Early acceleration", "Consistently high"))) %>% 
  #filter(visit == 0) %>% count(all_high, acc, over_at_78, class)
  ggplot(aes(x = days, y = zbmi)) + geom_point(aes(color = class)) + #geom_line() +
  scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)) +
  scale_x_continuous(breaks = c(0, 365, 731, 1096, 1461, 1826, 2192, 2557),
                     labels = c("Birth", "1 year", "2 years", "3 years", 
                                "4 years", "5 years", "6 years", "7 years")) +
  labs(x = "Age", y = "BMI Z-Score", color = "Growth pattern: ",
       title = "General growth patterns for children who are overweight in pre-school ages, by ethnicity (N = 202)") +
  geom_smooth(aes(x = days, y = zbmi, color = class), alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = -1, linetype = "dotted") +
  theme(legend.position = "NONE") + facet_wrap(~ ethnicity)

# BMI
zbmi_annotate %>% group_by(SubjectID) %>% 
  mutate(early_norm = sum(zbmi < 0.5 & visit <= 24, na.rm = T) >= 5, #max = 9 points
         late_grow = sum(zbmi > 1 & visit > 24, na.rm = T) >= 3, #max = 7 point
         acc = early_norm & late_grow,
         all_high = (sum(zbmi > 0.5 & visit <= 24, na.rm = T) >= 4) & 
           (sum(zbmi > 1.2 & visit > 24, na.rm = T) >= 4)) %>% 
  #ungroup %>% filter(visit == 0) %>% count(acc, all_high)
  ungroup() %>% filter(xor(all_high, acc)) %>% 
  mutate(RACE = recode(mother_ethnicity, 
                       "chinese" = "Chinese", 
                       "malay" = "Malay",
                       "indian" = "Indian")) %>%
  ggplot(aes(x = days, y = bmi)) + geom_point(aes(color = acc)) + #geom_line() +
  #scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7)) +
  scale_x_continuous(breaks = c(0, 365, 731, 1096, 1461, 1826, 2192, 2557),
                     labels = c("birth", "1 year", "2 years", "3 years", 
                                "4 years", "5 years", "6 years", "7 years")) +
  geom_smooth(aes(x = days, y = bmi, color = acc), alpha = 0.5) +
  theme(legend.position = "NONE") + facet_wrap(~ RACE) +
  labs(x = "Age", y = "BMI", title = "Consistently high versus catch-up BMI, by ethnicity.")

###################################
# Export data with growth classifications
###################################
zbmi_classed <- zbmi_annotate %>% group_by(SubjectID) %>% 
  mutate(early_norm = sum(zbmi < 0.5 & visit <= 24, na.rm = T) >= 5, #max = 9 points
         late_grow = sum(zbmi > 1 & visit > 24, na.rm = T) >= 3, #max = 7 point
         acc = early_norm & late_grow,
         all_high = sum(zbmi > 1, na.rm = T) >= 8) %>% ungroup() 

write_dta(zbmi_classed, paste0(base_path,"20181001-zBMI_classes.dta"), version = 12)

###################################
# ZBMI by EDU
###################################
zbmi_annotate %>% filter(mother_ethnicity == "chinese", visit == 78) %>% 
  group_by(mother_highest_education) %>% 
  summarize(mean = mean(zbmi, na.rm = T), sd = sd(zbmi, na.rm = T))
  
zbmi_annotate %>% filter(mother_highest_education != "no_education" & 
                           mother_highest_education != "" &
                           #mother_highest_education == "primary" &
                           mother_ethnicity == "chinese") %>% 
  ggplot(aes(x = days, y = zbmi)) + 
  geom_point(aes(color = mother_highest_education)) +
  geom_smooth(aes(color = mother_highest_education), span = 0.25) + 
  theme(legend.position = "bottom")

###################################
# ZBMI by TX
###################################
zbmi_annotate %>% mutate(GDM_TX = recode(gdm_treatment, 
                                         .default = "non-GDM",
                                         "None" = "untreated",
                                         "Diet" = "diet",
                                         "Insulin" = "insulin")) %>% 
  filter(visit == 78, !is.na(zbmi)) %>% 
  group_by(GDM_TX) %>% 
  summarize(n = n(),
            mean_zbmi = mean(zbmi, na.rm = T), 
            sd_zbmi = sd(zbmi, na.rm = T),
            mean_bmi = mean(cbmi, na.rm = T), 
            sd_bmi = sd(cbmi, na.rm = T),
            mean_fast = mean(ogtt_fasting_pw26, na.rm = T), 
            sd_fast = sd(ogtt_fasting_pw26, na.rm = T),
            min_fast = min(ogtt_fasting_pw26, na.rm = T),
            max_fast = max(ogtt_fasting_pw26, na.rm = T),
            mean_2hr = mean(ogtt_2hour_pw26, na.rm = T), 
            sd_2hr = sd(ogtt_2hour_pw26, na.rm = T),
            mean_ppBMI = mean(ppBMI, na.rm = T), 
            sd_ppBMI = sd(ppBMI, na.rm = T))

zbmi_annotate %>% mutate(GDM_TX = recode(gdm_treatment, 
                                         .default = "non-GDM",
                                         "None" = "untreated",
                                         "Diet" = "diet",
                                         "Insulin" = "insulin")) %>% 
  ggplot(aes(x = days, y = zbmi)) + 
  geom_point(aes(color = GDM_TX)) +
  scale_x_continuous(breaks = c(0, 365, 731, 1096, 1461, 1826, 2192, 2557),
                     labels = c("Birth", "1 year", "2 years", "3 years", 
                                "4 years", "5 years", "6 years", "7 years")) +
  labs(x = "Age", y = "BMI", color = "GDM treatment status: ",
       title = "Smoothed average BMI Z-score, by GDM treatment status.") +
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(aes(color = GDM_TX), span = 0.25) + 
  theme(legend.position = "bottom")


zbmi_annotate %>% mutate(GDM_TX = recode(gdm_treatment, 
                                         .default = "normoglycemic",
                                         "None" = "none reported",
                                         "Diet" = "diet",
                                         "Insulin" = "insulin")) %>% 
  ggplot(aes(x = days, y = cbmi)) + 
  geom_point(aes(color = GDM_TX)) +
  geom_smooth(aes(color = GDM_TX), span = 0.25) + 
  scale_x_continuous(breaks = c(0, 365, 731, 1096, 1461, 1826, 2192, 2557),
                     labels = c("Birth", "1 year", "2 years", "3 years", 
                                "4 years", "5 years", "6 years", "7 years")) +
  labs(x = "Age", y = "BMI", color = "GDM treatment status: ",
       title = "Smoothed average BMI, by GDM treatment status.") +
  scale_color_brewer(palette = "Spectral") +
  theme(legend.position = "bottom")

zbmi_annotate %>% mutate(GDM_TX = recode(gdm_treatment, 
                                         .default = "normoglycemic",
                                         "None" = "none reported",
                                         "Diet" = "diet",
                                         "Insulin" = "insulin"),
                         FAST_CAT = factor((ogtt_fasting_pw26 > 3.5) + 
                                             (ogtt_fasting_pw26 > 4) + 
                                             (ogtt_fasting_pw26 > 4.5) + 
                                             (ogtt_fasting_pw26 >= 5.1) +
                                             (ogtt_fasting_pw26 > 5.5)) ) %>% 
  # filter(visit == 78) %>%
  # group_by(FAST_CAT) %>% 
  # summarize(n = n(), mu_fast = mean(ogtt_fasting_pw26, na.rm = T), 
  #           mu_bmi = mean(bmi, na.rm = T), 
  #           min_ppBMI = min(ppBMI, na.rm = T),
  #           max_ppBMI = max(ppBMI, na.rm = T))
  filter(!is.na(FAST_CAT), visit < 78) %>% 
  ggplot(aes(x = days, y = bmi)) + 
  geom_point(aes(color = FAST_CAT)) +
  geom_smooth(aes(color = FAST_CAT), span = 0.25, alpha = 0.2) + 
  scale_x_continuous(breaks = c(0, 365, 731, 1096, 1461, 1826, 2192, 2557),
                     labels = c("Birth", "1 year", "2 years", "3 years", 
                                "4 years", "5 years", "6 years", "7 years")) +
  labs(x = "Age", y = "BMI", color = "Fasting glucose category: ",
       title = "Smoothed average BMI, by fasting glucose category.") +
  facet_wrap(~GDM_TX) +
  scale_color_brewer(palette = "Spectral", 
                     labels = c(bquote("\u2264 3.5 mmol/l"), 
                                "(3.5, 4] mmol/l", 
                                "(4, 4.5] mmol/l", 
                                "(4.5, 5.1) mmol/l", 
                                "[5.1, 5.5] mmol/l", 
                                "> 5.5 mmol/l"), 
                     guide = guide_legend(title.position = "top", nrow = 1)) +
  theme(legend.position = "bottom")

  ###################################
  # fasting glucose by ppBMI and ethnicity 
  ###################################  
  zbmi_annotate %>% filter(visit == 78, mother_ethnicity != "others", mother_ethnicity != "") %>% 
    ggplot(aes(x = ogtt_fasting_pw26, y = ppBMI, color = mother_ethnicity)) + 
    geom_smooth(alpha = 0.2) + geom_jitter() +
    labs(title = "ppBMI by fasting glucose and ethnicity") +
    theme(legend.position = "bottom")
  
  zbmi_annotate %>% filter(visit == 78, mother_ethnicity != "others", mother_ethnicity != "") %>% 
      ggplot(aes(y = ogtt_fasting_pw26, x = ppBMI, color = mother_ethnicity)) + 
      geom_smooth(alpha = 0.2) + geom_jitter() +
      labs(title = "fasting glucose by ppBMI and ethnicity") +
    theme(legend.position = "bottom")
  
  
  zbmi_annotate %>% filter(visit == 78) %>% 
    ggplot(aes(y = ogtt_fasting_pw26, x = ppBMI, color = gdm_treatment)) + 
    geom_smooth(method = glm, alpha = 0.2) + geom_jitter(aes(size = ogtt_2hour_pw26^3)) +
    labs(title = "F3. fasting glucose by ppBMI and GDM treatment modality",
         size = "scaled 2 hour glucose",
         color = "GDM treatment modality",
         y = "fasting glucose (mmol/l)") +
    scale_color_discrete(label = c("normoglycemic", "diet", "insulin", "untreated")) +
    theme(legend.position = "bottom")
  