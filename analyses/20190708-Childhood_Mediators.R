##############################################
# Sensitivity analyses unrelated to genetics
# To explain associations with growth:
# EFW
# Mode of delivery
# Breastfeeding
# Child infections
##############################################
library(tufte)
library(knitr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(haven)
library(tidyverse)
library(table1)

base_path <- "C:/Users/jhuangyh/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"
merged3 <- read_dta(paste0(base_path,"20181219-IVF-final.dta"))
add_data <- read_xlsx("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Clean dataset/forma428&429_20181224.xlsx", na = c("N/A", "NA"," ", "")) %>% 
  mutate(IVF = substring(SubjectID, 1,3) %in% c("019","029"))

#########################
# PREGNANCY / PERINATAL
#########################
# Age by Parity
merged3 %>% 
  filter(visit == 0, mother_ethnicity %in% c("chinese","indian","malay"), sex != "") %>% 
  group_by(parity) %>% summarize(mu = mean(mother_age_delivery, na.rm = T), 
                                 se = sd(mother_age_delivery, na.rm = T)/n(), 
                                 LL = mu - 1.96*se, 
                                 UL = mu + 1.96*se,
                                 N = n()) #%>% 
  #ggplot() + geom_point(aes(x = parity, y = mu)) + geom_ribbon(aes(x = parity, ymin = LL, ymax = UL), alpha = 0.4)

#########
# Maternal alcohol use - same or less
#########
add_data %>% select(IVF, contains("wine"), contains("beer")) %>% names()

add_data %>% group_by(IVF) %>% select(IVF, contains("wine")) %>% 
  summarize_all(funs(sum(., na.rm = T)/n())) %>% gather(var, val, -IVF)

add_data %>% group_by(IVF) %>% select(IVF, contains("beer")) %>% 
  summarize_all(funs(sum(., na.rm = T)/n())) %>% gather(var, val, -IVF)

add_data %>% group_by(IVF) %>% count(beer_times_per_month_pp, beer_times_per_month_preg)
add_data %>% group_by(IVF) %>% count(trad_wine_times_per_month_pp, trad_wine_times_per_month_preg) %>% print(n = Inf)
add_data %>% group_by(IVF) %>% count(wine_times_per_month_pp, wine_times_per_month_preg) %>% print(n = Inf)
add_data %>% group_by(IVF) %>% count(wine_times_per_month_preg, trad_wine_times_per_month_preg, beer_times_per_month_preg)

########
# EFW
########
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

########
# MODE OF DELIVERY - DOESN'T DIFFER BY IVF
########
merged3 %>% filter(visit == 0, sex != "") %>% filter(SubjectID %in% sub_id) %>%
  #mutate(VD = if_else(substring(mode_of_delivery,1,1) %in% c("1", "2"), 1, 0)) %>% 
  group_by(IVF) %>% mutate(tot = n()) %>% mutate(mod = substr(mode_of_delivery,1,1)) %>% 
  mutate(mod = factor(mod, levels = c("1","2","3","4","5"), 
                      labels = c("Spont. VD", "Spont. Assist. VD", "Induced VD", "Elective CS", "Emerg. CS"))) %>% 
  #ggplot() + geom_area(aes(x = mode_of_delivery, y = stat(scaled), fill = factor(IVF)), stat = "density") # pretty, but wrong
  ggplot() + geom_histogram(aes(mod, fill = factor(IVF)), stat = "count") +
  labs(x = "Mode of Delivery", y = "Frequency") + theme(legend.position = "none")

count(mode_of_delivery, tot) %>% mutate(pct = round(100 * n/tot, 1))

# GESTATIONAL AGE
# 22 < 35 wks: 19 non-IVF, 3 IVF
merged3 %>% filter(mother_ethnicity %in% c("chinese", "indian", "malay"), visit == 0, parity == 0) %>%
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Non-IVF", "IVF"))) %>% 
  ggplot() + geom_histogram(aes(GA, fill = IVF))  + facet_wrap(~mother_ethnicity) +
  theme(legend.position = "bottom") +
  labs(title = "Gestational age at birth, by IVF status",
       x = "Gestational age (weeks)",
       y = "Frequency")

# GESTATIONAL WEIGHT GAIN
# Number of weight obs by IVF
mweight_long %>% group_by(SubjectID) %>% 
  summarize(obs = sum(!is.na(mwt)), IVF = max(IVF)) %>% 
  #ggplot() + geom_histogram(aes(x = obs, fill = as.factor(IVF)))
  group_by(IVF) %>% count(obs)

########
# Maternal GWG by IVF
########
mweight <- merged3 %>% 
  filter(sex != "", visit == 0, !(mother_ethnicity %in% c("", "others"))) %>% #filter(parity == 0) %>% 
  select(SubjectID, mwt_0 = ppWeight, mother_ethnicity, parity, mht = m_height_pw26, mwt_1 = weight1, mwt_2 = weight2, 
         mwt_3 = weight3, mwt_4 = weight4, mwt_5 = weight5, mwt_6 = weight6, mwt_7 = weight7, mwt_8 = weight8, mwt_9 = weight9,
         mwtwk_1 = gestationalWeek1, mwtwk_2 = gestationalWeek3, mwtwk_4 = gestationalWeek4, mwtwk_5 = gestationalWeek5,
         mwtwk_6 = gestationalWeek6, mwtwk_7 = gestationalWeek7, mwtwk_8 = gestationalWeek8, mwtwk_9 = gestationalWeek9) %>% 
  mutate(mwtwk_0 = 0) %>% gather(measure, val, -SubjectID, -mother_ethnicity, -parity, -mht) 

mweight_long <- mweight %>% 
  separate(measure, c("measure", "GA_visit")) %>% 
  spread(measure, val) %>% 
  mutate(IVF = substring(SubjectID,1,3) %in% c("019","029"),
         mbmi = mwt/(mht/100)^2)

# GWG: look for gross discrepancies btw GA and visit #
mweight_long %>% filter(mwtwk != 0) %>% ggplot() + geom_point(aes(x = mwtwk, y = mwt, color = GA_visit))
mweight_long %>% filter(mwtwk != 0) %>% ggplot() + geom_point(aes(x = mwtwk, y = mbmi, color = GA_visit))

# GWG: differences by IVF
mweight_long %>% filter(mwtwk != 0) %>% filter(SubjectID %in% sub_id) %>% #filter(mwt < 120) %>% 
  ggplot(aes(x = mwtwk, y = mbmi, color = IVF)) + geom_point() + geom_smooth() +
  geom_line(aes(x = mwtwk, y = mbmi, group = SubjectID), alpha  = 0.3) +
  labs(title = "Maternal BMI change during pregnancy, by IVF status.", 
       x = "Estimated gestational age at visit (Weeks)", 
       y = expression(paste("Maternal BMI ( ", kg/m^2, ")")), 
       color = "IVF status") +
  scale_color_discrete(labels = c("Spontaneous (subcohort)", "IVF")) +
  theme(legend.position = "bottom")

# All GWG velocities
mweight_long %>% filter(mwtwk != 0) %>% mutate(GA_visit = as.double(GA_visit)) %>% group_by(SubjectID) %>% 
  mutate(wt_dif = mwt - lag(mwt, order_by = GA_visit),
         ga_dif = mwtwk - lag(mwtwk, order_by = GA_visit),
         vel = wt_dif/ga_dif) %>% group_by(IVF) %>% 
  #filter(ga_dif >= 7 & ga_dif <= 14 & mwtwk < 37) %>% ggplot(aes(x = mwtwk, y = vel, color = as.factor(IVF))) + geom_point() + geom_smooth() + labs(title = "Velocity over time, by IVF Status")
  #ggplot(aes(x = ga_dif, y = wt_dif, color = mwtwk)) + geom_point() + scale_color_distiller(palette = "Spectral") + labs(title = "Abs wt change by elapsed weeks")
  ggplot(aes(x = ga_dif, y = vel, color = mwtwk)) + geom_point() + geom_smooth() + scale_color_distiller(palette = "Spectral") + 
  labs(x = "Elapsed weeks",
       y = "Velocity (kg/week)",
       title = "GWG velocity (kg / week), by elapsed weeks",
       color = "GA at calculation") +
  theme(legend.position = "bottom") 
#ggplot(aes(x = mwtwk, y = vel, color = ga_dif)) + geom_point() + geom_smooth() + scale_color_distiller(palette = "Spectral") + labs(title = "Vel by GA") 

# calculate GWG by max_distance (last - first weight) / (last - first GA)
# USING ANY MEASURES:
# Overall:    0.392 (SC) vs. 0.395 (IVF), 28.6 vs. 27.5 wks
# Subcohort:  0.381 (SC) vs. 0.395 (IVF), 29.1 vs. 27.5 wks
mweight_long %>% filter(mwtwk != 0) %>% group_by(SubjectID) %>% 
  summarize(max = max(mwtwk, na.rm = T), min = min(mwtwk, na.rm = T)) %>% 
  left_join(mweight_long, .) %>% 
  group_by(SubjectID) %>% 
  filter(mwtwk == min | mwtwk == max) %>% 
  summarize(wt_change = max(mwt) - min(mwt), 
            ga_change = max(mwtwk) - min(mwtwk),
            vel = wt_change/ga_change) %>% 
  mutate(IVF = substring(SubjectID,1,3) %in% c("019","029"),
         IVF = factor(IVF, levels = c(F,T), labels = c("SC","IVF"))) %>% 
  #filter(SubjectID %in% sub_id) %>% 
  #group_by(IVF) %>% count(is.na(vel))
  #filter(!is.na(vel)) %>% pull(SubjectID)
  #group_by(IVF) %>% filter(!is.na(vel)) %>% summarize(n = n(), mu = sum(vel, na.rm = T)/n(), avg_ga = sum(ga_change, na.rm = T)/n())
  ggplot() + geom_histogram(aes(vel, fill = IVF))
# not much difference in overall weight gain velocity defined by earliest versus latest visit

# calculate 2nd/3rd trimester GWG (13 wks +)
# subjects with >2 measures (1018 SC, 77 IVF)
# Overall:    0.715 (SC) vs. 0.720 (IVF), 15.3 vs. 15.5 wks
# Subcohort:  0.685 (SC) vs. 0.720 (IVF), 16.2 vs. 15.5 wks
mweight_long %>% group_by(SubjectID) %>% summarize(IVF = max(IVF), gt13 = sum(mwtwk >= 13,na.rm = T)) %>% group_by(IVF) %>% count(gt13 > 2) 
# SC = 1018 with 2+ measures; IVF = 77 with 2+ measures
mweight_long %>% filter(mwtwk > 13) %>% group_by(SubjectID) %>% 
  summarize(max = max(mwtwk, na.rm = T), min = min(mwtwk, na.rm = T)) %>% 
  left_join(mweight_long, .) %>% 
  group_by(SubjectID) %>% 
  filter(mwtwk == min | mwtwk == max) %>%
  mutate(last_wt = if_else(mwtwk == max, mwt, 0), 
         first_wt = if_else(mwtwk == min, mwt, 0)) %>% 
  summarize(last_ga = max(mwtwk),
            first_ga = min(mwtwk),
            last_wt = max(last_wt),
            first_wt = max(first_wt),
            wt_change = last_wt - first_wt, 
            ga_change = last_ga - first_ga,
            vel = wt_change/ga_change) %>% 
  mutate(IVF = substring(SubjectID,1,3) %in% c("019","029"),
         IVF = factor(IVF, levels = c(F,T), labels = c("SC","IVF"))) %>% 
  #filter(SubjectID %in% sub_id) %>% 
  #ggplot() + geom_histogram(aes(vel, fill = IVF))
  #filter(vel > 1) %>% arrange(vel)
  #group_by(IVF) %>% filter(!is.na(vel)) %>% summarize(n = n(), mu = sum(vel, na.rm = T)/n(), min(first_ga), max(last_ga), avg_ga = sum(ga_change, na.rm = T)/n())
  filter(!is.na(vel)) %>% ggplot() + geom_histogram(aes(ga_change, fill = IVF), binwidth = 0.5) # distributions appear similar
ggplot() + geom_point(aes(x = last_ga, y = vel, size = ga_change)) #should be okay, the high velocities a product of short/late follow-up
# IVF 2-3 trimester weight gain similarly to overall cohort, higher than subfertile subcohort

mweight_long %>% filter(mwtwk > 0) %>% 
  filter(SubjectID %in% sub_id) %>% 
  #ggplot() + geom_line(aes(x = mwtwk, y = mwt, color = SubjectID)) + theme(legend.position = "none") + facet_wrap(~IVF)
  ggplot() + geom_line(aes(x = mwtwk, y = mbmi, color = SubjectID)) + theme(legend.position = "none") + facet_wrap(~IVF)

#mweight_long %>% ggplot() + geom_point(aes(x = mwtwk, y = mwt, color = IVF)) + geom_smooth(aes(x = mwtwk, y = mwt, color = IVF))
#mweight_long %>% ggplot() + geom_point(aes(x = mwtwk, y = mbmi, color = IVF)) + geom_smooth(aes(x = mwtwk, y = mbmi, color = IVF))
mweight_long %>% filter(parity == 0) %>% 
  ggplot() + geom_point(aes(x = mwtwk, y = mwt, color = IVF)) + 
  geom_smooth(aes(x = mwtwk, y = mwt, color = IVF)) + facet_wrap(~mother_ethnicity)

# Number of weight measures
mweight_long %>% filter(!is.na(mwt)) %>% group_by(SubjectID, IVF) %>% count(!is.na(mwt)) %>% 
  ggplot() + geom_histogram(aes(n, fill = IVF)) + 
  labs(x = "Number of measurements", y = "Count", title = "Number of antepartum weight measurements per subject, by IVF status.") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9, 10)) + theme(legend.position = "bottom")

mweight_long %>% filter(!is.na(mwt)) %>% ggplot() + geom_histogram(aes(mwtwk, fill = IVF))
merged3 %>% filter(visit == 0, sex != "") %>% ggplot() + geom_histogram(aes(last_antenatal_weight_GA, fill = IVF))

# Number of subjects measured per visit
mweight_long %>% filter(!is.na(mwt)) %>% group_by(IVF, GA_visit) %>% summarize(n = sum(!is.na(SubjectID))) %>% 
  ggplot() + geom_bar(aes(x = GA_visit, y = n, fill = IVF), stat = "identity", bins = 10) +
  labs(x = "Visit number", y = "Count", title = "Number of subjects with antepartum weight measurements per visit, by IVF status.") +
  scale_x_discrete(labels = c("pre-preg.","1", "2", "3", "4", "5", "6" , "7", "8", "9")) + theme(legend.position = "bottom")

# Maternal Blood Pressure
merged3 %>% filter(visit == 0) %>% 
  select(IVF, sbp_1 = sbp1, sbp_2 = sbp2, sbp_3 = sbp3, dbp_1 = dbp1, dbp_2 = dbp2, dbp_3 = dbp3) %>% 
  gather(measure, val, - IVF) %>% 
  separate(measure, c("measure", "visit")) %>% 
  ggplot() + 
  geom_jitter(aes(x = as.factor(visit), y = val, color = as.factor(visit))) + facet_wrap(~measure)


merged3 %>% mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Spontaneous", "IVF")),
                   htn = if_else(HTN_clean %in% c("Chronic HTN superimposed PE", "Eclampsia", "PE", "PIH"), 1, 
                                 if_else(is.na(hi_bp), NA_real_, 0)),
                   htn = factor(htn, levels = c(0,1), labels = c("No","Yes"))) %>% 
  #filter(!is.na(hi_bp)) %>% 
  filter(!is.na(htn)) %>% 
  ggplot() + 
  #geom_point(aes(x = days, y = zlen, color = as.factor(hi_bp))) + 
  #geom_smooth(aes(x = days, y = zlen, color = as.factor(hi_bp)), alpha = 0.4, method = "loess") +
  geom_point(aes(x = days, y = zlen, color = as.factor(htn))) + 
  geom_smooth(aes(x = days, y = zlen, color = as.factor(htn)), alpha = 0.4) +
  labs(x = "Age (Days)", y = "Length-for-age (SDS)", color = "Any Elevated Blood Pressure",
       title = "Length-for-age, by maternal pregnancy blood pressure and IVF status.") +
  facet_wrap(~IVF) +
  #facet_wrap(~mother_ethnicity) +
  theme(legend.position = "bottom")

# Antepartum weight and height characteristics by IVF and ethnicity
anthro <- merged3 %>% filter(visit == 0, !(mother_ethnicity %in% c("", "others"))) %>% #filter(parity == 0) %>% 
  ggplot()

anthro + geom_jitter(aes(x = as.factor(IVF), y = ppWeight), width = .1) +
  geom_violin(aes(x = as.factor(IVF), y = ppWeight), draw_quantiles = c(.25,.5,.75), alpha = 0.5) + facet_grid(~mother_ethnicity)

anthro + geom_jitter(aes(x = as.factor(IVF), y = m_height_pw26), width = .1) +
  geom_violin(aes(x = as.factor(IVF), y = m_height_pw26), draw_quantiles = c(.25,.5,.75), alpha = 0.5) + facet_grid(~mother_ethnicity)

anthro + geom_jitter(aes(x = as.factor(IVF), y = weight1), width = .1) +
  geom_boxplot(aes(x = as.factor(IVF), y = weight1), draw_quantiles = c(.25,.5,.75), alpha = 0.5) + facet_grid(~mother_ethnicity)

anthro + geom_jitter(aes(x = as.factor(IVF), y = weight1), width = .1) +
  geom_violin(aes(x = as.factor(IVF), y = weight1), draw_quantiles = c(.25,.5,.75), alpha = 0.5) + facet_grid(~mother_ethnicity)

anthro + geom_jitter(aes(x = as.factor(IVF), y = last_antenatal_weight), width = .1) +
  geom_violin(aes(x = as.factor(IVF), y = last_antenatal_weight), draw_quantiles = c(.25,.5,.75), alpha = 0.5) + facet_grid(~mother_ethnicity)


#########################
# DURATION OF BREASTFEEDING - DOESN'T DIFFER BY IVF
#########################
merged3 %>% filter(visit == 0, sex != "") %>% filter(SubjectID %in% sub_id, Dur_any_BF != "") %>%
  mutate(ANY_BF = factor(Dur_any_BF, levels = c("lt_1M", "1M_to_lt_3M", "3M_to_lt_6M", "6M_to12M", "12M_and_above"),
                         labels = c("<1 mo.", "1 to <3 mo.", "3 to <6 mo.", "6 to 12 mo.", "12+ mo."))) %>% 
  group_by(IVF) %>% mutate(tot = n()) %>% 
  ggplot() + geom_histogram(aes(ANY_BF, fill = factor(IVF)), stat = "count") +
  labs(x = "Duration of Any Breastfeeding", y = "Frequency") + theme(legend.position = "none")


##################################################
# INFECTION AND HOSPITALIZATION EPISODES
##################################################
CHILD_PATH <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Raw LORIS data/GUSTO child health/"

#### 12 MONTHS
child_12 <- read_csv(paste0(CHILD_PATH, "Child_Questionnaire_Month12.csv"))

# DIAGNOSES (ANY = 32.5% SC; 36.8% IVF)
child_12 %>% select(contains("qn_1_27")) %>% names()
child_12 %>% group_by(IVF) %>% count(qn_1_27) %>% print(n = Inf)
child_12 <- child_12 %>%   mutate(IVF = substring(PSCID,3,3) == "9") %>% 
  mutate(n_dx = case_when(!is.na(qn_1_27_yes_SN4_diagnosis) ~ 4,
                          !is.na(qn_1_27_yes_SN3_diagnosis) ~ 3,
                          !is.na(qn_1_27_yes_SN2_diagnosis) ~ 2,
                          !is.na(qn_1_27_yes_SN1_diagnosis) ~ 1,
                          TRUE ~ NA_real_))
child_12 %>% filter(PSCID %in% sub_id) %>%
  #filter(substring(PSCID,3,3) != "7") %>% 
  #group_by(IVF) %>% summarize(mean(n_dx > 0, na.rm = T))
  ggplot() + geom_bar(aes(x = n_dx, fill = IVF), alpha = 0.4, position = "dodge")

# HOSPITALIZATIONS - only 1 in IVF (ANY = 3.3% SC; 1.5% IVF)
child_12 %>% select(contains("qn_1_25")) %>% names()
child_12 %>% group_by(IVF) %>% count(qn_1_25_admission) %>% print(n = Inf)
child_12 %>% group_by(IVF) %>% summarize(mean(qn_1_25_admission == "1_yes", na.rm = T))

# FEVERS (ANY = 37.7% SC; 35.3% IVF)
child_12 %>% select(contains("qn_1_23")) %>% names()
child_12 %>% group_by(IVF) %>% summarize(mean(qn_1_23_episodes, na.rm = T))
child_12 %>% group_by(IVF) %>% summarize(mean(!is.na(qn_1_23_episodes)))
child_12 %>% group_by(IVF) %>% mutate(fev = if_else(is.na(qn_1_23_episodes), 0, qn_1_23_episodes)) %>% 
                                      summarize(mean(fev))

# DIARRHEA (ANY = 16.5% SC; 14.7% IVF)
child_12 %>% group_by(IVF) %>% summarize(mean(qn_1_17_diarrhoea == "1_yes", na.rm = T))

# ANTIBIOTICS
child_12 %>% group_by(IVF) %>% count(qn_1_26)
child_12 %>% group_by(IVF) %>% summarize(mean(qn_1_26 == "1_yes", na.rm = T))
child_12 %>% count(qn_1_26_yes_SN1_duration, qn_1_26_yes_SN2_duration,
                   qn_1_26_yes_SN3_duration, qn_1_26_yes_SN4_duration) %>% print(n = Inf)
child_12 %>% group_by(IVF) %>% summarize(sum(!is.na(qn_1_26_yes_SN2_duration)), sum(is.na(qn_1_26_yes_SN2_duration)))
child_12 %>% mutate(qn_1_26_yes_SN1_duration = 
                      if_else(qn_1_26_yes_SN1_duration %in% c("dnk", "not_answered", NA_character_), 
                              NA_real_, as.double(qn_1_26_yes_SN1_duration))) %>% 
  mutate_at(c("qn_1_26_yes_SN2_duration", "qn_1_26_yes_SN3_duration", "qn_1_26_yes_SN4_duration"),
            list(~ if_else(is.na(.), 0, as.double(.)))) %>% 
  mutate(abx_days_12 = qn_1_26_yes_SN1_duration + 
           qn_1_26_yes_SN2_duration + 
           qn_1_26_yes_SN3_duration + 
           qn_1_26_yes_SN4_duration) %>% 
  group_by(IVF) %>% #summarize(mean(abx_days_12, na.rm = T))
  ggplot() + geom_density(aes(abx_days_12, fill = IVF), alpha = 0.4)

child_12 <- child_12 %>% filter(!is.na(qn_1_13_cough)) %>% 
  mutate(IVF = factor(IVF, levels = c(TRUE, FALSE), labels = c("IVF", "Spontaneous Conception"))) %>% 
  mutate(any_fevers = case_when(qn_1_23_fever_vaccination == "0_no" ~ "No",
                                qn_1_23_fever_vaccination == "1_yes" ~ "Yes",
                                TRUE ~ NA_character_)) %>% 
  mutate_at(c("qn_1_17_diarrhoea", "qn_1_25_admission", "qn_1_26", "qn_1_27"),
            list(~ case_when(
              . == "0_no" ~ "No",
              . == "1_yes" ~ "Yes",
              TRUE ~ NA_character_))) %>% 
  mutate(qn_1_26_yes_SN1_duration = 
           if_else(qn_1_26_yes_SN1_duration %in% c("dnk", "not_answered", NA_character_), 
                   NA_real_, as.double(qn_1_26_yes_SN1_duration))) %>% 
  mutate_at(c("qn_1_26_yes_SN2_duration", "qn_1_26_yes_SN3_duration", "qn_1_26_yes_SN4_duration"),
            list(~ if_else(is.na(.), 0, as.double(.)))) %>% 
  mutate(abx_days_12 = qn_1_26_yes_SN1_duration + 
           qn_1_26_yes_SN2_duration + 
           qn_1_26_yes_SN3_duration + 
           qn_1_26_yes_SN4_duration)

label(child_12$qn_1_17_diarrhoea) <- "Diarrhea lasting 2+ days"
label(child_12$any_fevers) <- "Any fevers > 38 degrees C"
label(child_12$qn_1_23_episodes) <- "Average # fevers (among any fevers)"
label(child_12$qn_1_25_admission) <- "Any hospital admissions"
label(child_12$qn_1_26) <- "Any antibiotics use"
label(child_12$abx_days_12) <- "Average # days of antibiotics (among users)"
label(child_12$qn_1_27) <- "Any other diagnoses"
label(child_12$n_dx) <- "Average # other diagnoses (among any diagnoses)"

table1(~ qn_1_17_diarrhoea + any_fevers + qn_1_23_episodes + qn_1_25_admission + 
         qn_1_26 + abx_days_12 + qn_1_27 + n_dx | IVF, 
       data = child_12, 
       topclass = "Rtable1-zebra",
       render.continuous = my.render.cont) 

#### 15 MONTHS
child_15 <- read_csv(paste0(CHILD_PATH, "Child_Questionnaire_Month15.csv"))
child_15 <- child_15 %>% mutate(IVF = substring(PSCID,3,3) == "9")
child_15 %>% group_by(IVF) %>% summarize(mean(qn_1_17_diarrhoea == "1_yes", na.rm = T))
child_15 %>% group_by(IVF) %>% summarize(mean(!(qn_1_21_episodes %in% c(NA_character_, "dnk"))))
child_15 %>% group_by(IVF) %>% count(qn_1_21_episodes)

merged_tmp <- child_15 %>% 
  mutate(dia_15 = 
           if_else(qn_1_17_diarrhoea == "0_no", 0,
                   if_else(qn_1_17_diarrhoea == "1_yes", 1, 
                           NA_real_))) %>% 
  select(SubjectID = PSCID, dia_15) %>% left_join(merged3, ., by = "SubjectID") %>% 
  filter(visit == 15)

glm(data = merged_tmp, formula = dia_15 ~ IVF + mother_ethnicity + 
      ppBMI + parity + mother_age_delivery, family = binomial) %>% summary()
