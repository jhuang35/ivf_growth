library(haven)
library(tidyverse)
library(ggrepel)
# library(vtree)

plot_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/PLOTS/Sample Size/"
base_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

######################################
# DATA INPUT
######################################
# dating - to estimate GA at scans
us11 <- read_csv("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Raw LORIS data/GUSTO ultrasound data cleaning/Dating_Ultrasound_Scan_CRF_Pregnancy_Week11.csv")
us11 %>% ggplot() + geom_freqpoly(aes(expected_date_delivery), bins=365)
us11 %>% count(!is.na(expected_date_delivery))
tmp <- us11 %>% select(SubjectID = PSCID, EDD = expected_date_delivery)
# %>% left_join(.,ultrasound19)
# tmp %>% mutate(GA_19 = 280 - (EDD - date_of_ultrasound_scan)) %>% summarize(mu = mean(GA_19, na.rm = T))
#   ggplot() + geom_histogram(aes(GA_19/7), binwidth = 1) + coord_cartesian(xlim = c(0,30))
  
#Create EFW using Shepard formula
#-1.7492+ 0.166*BPD +0.046*AC - 2.646*(AC*BPD)/1,000 =log10(bw)

ultrasound19 <- read_csv("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Raw LORIS data/GUSTO ultrasound data cleaning/Fetal_Anomaly_Ultrasound_Scan_CRF_Pregnancy_Week19.csv")
us19 <- ultrasound19 %>% 
  select(SubjectID = PSCID,
         BPD = biparietal_diameter,
         AC = abdominal_circumference,
         DATE = date_of_ultrasound_scan) %>% 
  mutate(BPD = BPD/10, AC = AC/10,
         EFW = 10^(-1.7492 + 0.166*BPD + 0.046*AC - 2.646*(AC*BPD)/1000)) 
us19 <- us19 %>% left_join(tmp, .) %>% 
  mutate(DAYS = 280 - (EDD - DATE), 
         DAYS = if_else(DAYS < 0, as.difftime(NA_real_, unit = "days"), DAYS)) %>% 
  add_column(visit = 19)

ultrasound26 <- read_csv("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Raw LORIS data/GUSTO ultrasound data cleaning/Growth_Ultrasound_Scan_CRF_Pregnancy_Week26.csv")
us26 <- ultrasound26 %>% 
  select(SubjectID = PSCID, 
         BPD = biparietal_diameter, 
         AC = abdominal_circumference, 
         DATE = date_of_ultrasound_scan) %>% 
  mutate(BPD = BPD/10, AC = AC/10,
         EFW = 10^(-1.7492 + 0.166*BPD + 0.046*AC - 2.646*(AC*BPD)/1000))
us26 <- us26 %>% left_join(tmp, .) %>% 
  mutate(DAYS = 280 - (EDD - DATE), 
         DAYS = if_else(DAYS < 0, as.difftime(NA_real_, unit = "days"), DAYS)) %>% 
  add_column(visit = 26)
  
ultrasound32 <- read_csv("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Raw LORIS data/GUSTO ultrasound data cleaning/Growth_and_Doppler_Ultrasound_Scan_CRF_Pregnancy_Week32.csv")
us32 <- ultrasound32 %>% 
  select(SubjectID = PSCID, 
         BPD = biparietal_diameter, 
         AC = abdominal_circumference,
         DATE = date_of_ultrasound_scan) %>% 
  mutate(BPD = BPD/10, AC = AC/10,
         EFW = 10^(-1.7492 + 0.166*BPD + 0.046*AC - 2.646*(AC*BPD)/1000)) 
us32 <- us32 %>% left_join(tmp, .) %>% 
  mutate(DAYS = 280 - (EDD - DATE), 
         DAYS = if_else(DAYS < 0, as.difftime(NA_real_, unit = "days"), DAYS)) %>% 
  add_column(visit = 32)

us_all <- full_join(us19, full_join(us26, us32))
us_all_original <- full_join(us19, full_join(us26, us32))

######################################
# DATA VERIFICATION
######################################
# us visit timing
us_all_original %>% ggplot() + 
  geom_histogram(aes(DAYS, fill = as.factor(visit)), binwidth = 1) + 
  labs(title = "Calculated gestational age, by visit.",
       x = "calculated gestational age (days)",
       fill = "Visit Week") +
  annotate("text", x = 300, y = 5, label = "010-10017") +
  annotate("text", x = 475, y = 5, label = "019-30036") +
  theme(legend.position = "bottom")
us_all_original %>% filter(DAYS < 77 | DAYS > 280) # < 11 wks or > 40 wks 
us_all <- us_all %>% 
  mutate(DATE = if_else(SubjectID == "019-30036" & visit == 26, as.Date("19/2/2010", "%d/%m/%Y"), DATE), # one corrected date - per Der Horng 23/NOV/2018
         DAYS = if_else(SubjectID == "019-30036", 280 - (EDD - DATE), DAYS),
         DAYS = if_else(SubjectID == "010-10017", as.difftime(NA_real_, unit = "days"), DAYS)) # EDD incorrect, individual dropS anyway

# EFW by visit/days
params <- us_all %>% group_by(visit) %>% summarize(mu_BPD = mean(BPD, na.rm = T), mu_AC = mean(AC, na.rm = T), mu_EFW = mean(EFW, na.rm = T),
                                                   sd_BPD = sd(BPD, na.rm = T), sd_AC = sd(AC, na.rm = T), sd_EFW = sd(EFW, na.rm = T),
                                                   min_BPD = mu_BPD - 4*sd_BPD, min_AC = mu_AC - 4*sd_AC, min_EFW = mu_EFW - 4*sd_EFW,
                                                   max_BPD = mu_BPD + 4*sd_BPD, max_AC = mu_AC + 4*sd_AC, max_EFW = mu_EFW + 4*sd_EFW)
us_all_original %>% filter(EFW > 4)
us_all_original %>% filter(DAYS < 300) %>% ggplot() +
  geom_point(aes(x = DAYS, y = EFW, color = as.factor(visit))) +
  labs(title = "Estimated fetal weight, by visit.",
       x = "calculated gestational age (days)",
       y = "EFW (kg)",
       color = "Visit Week") +
  annotate("text", x = 142, y = 17, label = "010-20099") +
  annotate("text", x = 148, y = 22, label = "010-21774") +
  theme(legend.position = "bottom")

us_all <- us_all %>% mutate(BPD = if_else(SubjectID %in% c("010-20099","010-21774") & visit == 19, NA_real_, BPD),
                            EFW = if_else(SubjectID %in% c("010-20099","010-21774") & visit == 19, NA_real_, EFW)) # BPD clearly incorrect (N = 2; >>> 4 SD)
us_all %>% ggplot() +
  geom_point(aes(x = DAYS, y = EFW, color = as.factor(visit))) +
  theme(legend.position = "bottom")

# BPD by AC (N = 8 < 4 SD either AC or BPD)
OOR_TBL <- us_all_original %>% mutate(OOR_ID = if_else(AC < 10.2 | BPD < 3.43, paste0(SubjectID,"(",visit,")"), "")) %>% filter(OOR_ID != "") 

us_all_original %>% filter(DAYS < 300, EFW < 4) %>% ggplot() + 
  geom_point(aes(x = BPD, y = AC, color = as.factor(visit))) + 
  geom_text_repel(data = OOR_TBL, aes(x = BPD, y = AC, label = SubjectID)) +
  labs(title = "BPD versus AC, by visit.",
       x = "BPD (cm)",
       y = "AC (cm)",
       color = "Visit Week") +
  theme(legend.position = "bottom")

OOR_TBL <- us_all %>% mutate(OOR_ID = if_else(AC < 10.2 | BPD < 3.43, paste0(SubjectID,"(",visit,")"), "")) %>% filter(OOR_ID != "") 
us_all <- us_all %>% mutate(BPD = if_else(BPD < 3.43, NA_real_, BPD),
                            AC = if_else(AC < 10.2, NA_real_, AC),
                            EFW = if_else(is.na(BPD) | is.na(AC), NA_real_, EFW)) 
OOR_ID <- c(as.character(OOR_TBL$SubjectID), "010-10017", "010-20099", "010-21774", "019-30036")

###########################################
# FIRST PASS FOR OBVIOUS VISIT AND PARAMETER ERRORS (N = 13 modified):
# 2 CORRECTED DATE OF VISIT: 019-30036, 010-10017
# 2 SET TO MISSING DUE TO BPD > 4 SD: 010-20099, 010-21774
# 9 SET TO MISSING DUE TO AC OR BPD < 4SD: "010-04010" "010-20332" "010-21127" "010-21572" "010-21729" "010-04064" "010-21170" "010-21266" "010-20799"
###########################################

# show points which have discrepancy between visit number and estimated GA -- they seem consistent with a different visit date
us_all %>% filter((visit == 19 & DAYS > 182) | 
                    (visit == 26 & DAYS < 133) |
                    (visit == 26 & DAYS > 224) |
                    (visit == 32 & DAYS < 182)) # if any visit dates less that supposed previous visit, or more than next visit

# generate sub-tables to label flagged points
# 8 flags based on these overlaps
FLAG_TBL <- us_all_original %>% filter(DAYS < 300) %>% 
  mutate(flag = (visit == 19 & DAYS > 182) | 
           (visit == 26 & DAYS < 133) |
           (visit == 26 & DAYS > 224) |
           (visit == 32 & DAYS < 182), 
         WEEKS = as.double(DAYS/7)) %>% filter(flag) %>% arrange(BPD)
FLAG_ID <- as.character(FLAG_TBL$SubjectID)

canvas <- us_all_original %>% filter(!(SubjectID %in% OOR_ID)) %>% 
  ggplot() + 
  geom_point(aes(x = BPD, y = AC, color = as.factor(visit), alpha = 0.8)) + 
  scale_alpha(guide = "none") + 
  theme(legend.position = "bottom")

canvas + geom_point(data = FLAG_TBL, 
                    aes(x = BPD, y = AC, color = as.factor(visit), size = 2)) +
  labs(title = "BPD, AC, and EFW, by visit and discordant visits.",
       x = "BPD (cm)",
       y = "AC (cm)",
       color = "Visit Week",
       size = "Flagged discordant") +
  scale_size(guide = "none") +
  geom_text_repel(data = FLAG_TBL, aes(x = BPD, y = AC, label = SubjectID, color = as.factor(visit)), point.padding = 5)

us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "019-30288") # 1. 26 wk EFW consistent with 26 wk, drop DAYS for 26 wk
us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "010-22055") # 2. 26 wk EFW consistent with 26 wk, drop DAYS for 26 wk
us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "020-08019") # 3. 19 wk EFW consistent with 26 wk, change 19 wk measures to 26 wk (set existing 19 wk to missing)
us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "010-22164") # 4. 26 wk EFW consistent with 26 wk, drop DAYS for 26 wk
us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "010-21198") # 5. 19 wk duplicate entry, DROP entry for 19 wks <<-- remove entry for 19 wks
us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "010-21873") # 6. "26 wk" and "32 wk" entries swapped, swap visit values
us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "010-21927") # 7. 26 wk duplicate entry, delete entry for 26 wks <<-- remove entry for 26 wks
us_all %>% mutate(WKS = as.double(DAYS/7)) %>% filter(SubjectID == "010-21569") # 8. 26 wk duplicate entry, delete entry for 26 wks

us_all <- us_all %>% 
  mutate(DAYS = if_else(SubjectID == "019-30288" & visit == 26, as.difftime(NA_real_, unit = "days"), DAYS),
         DAYS = if_else(SubjectID == "010-22055" & visit == 26, as.difftime(NA_real_, unit = "days"), DAYS),
         visit = if_else(SubjectID == "020-08019" & visit == 19, 26,
                         if_else(SubjectID == "020-08019" & is.na(DATE), 19, visit)),
         DAYS = if_else(SubjectID == "010-22164" & visit == 26, as.difftime(NA_real_, unit = "days"), DAYS),
         
         EDD = if_else(SubjectID == "010-21198" & visit == 19, as.Date("0000-00-00", "%Y-%m-%d"), EDD),
         BPD = if_else(SubjectID == "010-21198" & visit == 19, NA_real_, BPD),
         AC = if_else(SubjectID == "010-21198" & visit == 19, NA_real_, AC),
         DATE = if_else(SubjectID == "010-21198" & visit == 19, as.Date("0000-00-00", "%Y-%m-%d"), DATE),
         EFW = if_else(SubjectID == "010-21198" & visit == 19, NA_real_, EFW),
         DAYS = if_else(SubjectID == "010-21198" & visit == 19, as.difftime(NA_real_, unit = "days"), DAYS),
         
         visit = if_else(SubjectID == "010-21873" & visit == 32, 26, 
                         if_else(SubjectID == "010-21873" & visit == 26, 32, visit)),
           
         EDD = if_else(SubjectID == "010-21927" & visit == 26, as.Date("0000-00-00", "%Y-%m-%d"), EDD),
         BPD = if_else(SubjectID == "010-21927" & visit == 26, NA_real_, BPD),
         AC = if_else(SubjectID == "010-21927" & visit == 26, NA_real_, AC),
         DATE = if_else(SubjectID == "010-21927" & visit == 26, as.Date("0000-00-00", "%Y-%m-%d"), DATE),
         EFW = if_else(SubjectID == "010-21927" & visit == 26, NA_real_, EFW),
         DAYS = if_else(SubjectID == "010-21927" & visit == 26, as.difftime(NA_real_, unit = "days"), DAYS),
         
         EDD = if_else(SubjectID == "010-21569" & visit == 26, as.Date("0000-00-00", "%Y-%m-%d"), EDD),
         BPD = if_else(SubjectID == "010-21569" & visit == 26, NA_real_, BPD),
         AC = if_else(SubjectID == "010-21569" & visit == 26, NA_real_, AC),
         DATE = if_else(SubjectID == "010-21569" & visit == 26, as.Date("0000-00-00", "%Y-%m-%d"), DATE),
         EFW = if_else(SubjectID == "010-21569" & visit == 26, NA_real_, EFW),
         DAYS = if_else(SubjectID == "010-21569" & visit == 26, as.difftime(NA_real_, unit = "days"), DAYS))

###########################################
# SECOND PASS FOR DISCREPANCIES BY VISIT (n = 8 Subject IDs modified)
# "010-21198" "020-08019" "010-21569" "010-21873" "010-21927" "010-22055" "010-22164" "019-30288"
# RUNNING TOTAL = 21 Subject IDs modified
###########################################

ALL_ID <- c(FLAG_ID, OOR_ID) # N = 21 modified Subject IDs

FLAG_TBL_CHK <- us_all %>% 
  mutate(flag = (visit == 19 & DAYS > 182) | 
           (visit == 26 & DAYS < 133) |
           (visit == 26 & DAYS > 224) |
           (visit == 32 & DAYS < 182)) %>% filter(flag)
canvas + geom_point(data = FLAG_TBL_CHK, 
                    aes(x = BPD, y = AC, color = as.factor(visit), size = 2)) +
  geom_text_repel(data = FLAG_TBL_CHK, aes(x = BPD, y = AC, label = SubjectID, color = as.factor(visit)), point.padding = 5)

LBL_TBL <- us_all_original %>% filter(SubjectID == "020-66132", visit == 19) %>% 
  select(BPD, AC, EFW, DAYS, visit, SubjectID) %>% 
  gather(param, value, -DAYS, -visit, -SubjectID) %>% 
  mutate(param = if_else(param == "BPD", "BPD (cm)",
                         if_else(param == "AC", "AC (cm)", 
                                 if_else(param == "EFW", "Estimated Fetal Weight (kg)","")))) %>% 
  mutate(DAYS = DAYS/7)

us_plots <- us_all_original %>% filter(!(SubjectID %in% c(FLAG_ID, OOR_ID))) %>% select(BPD, AC, EFW, DAYS, visit) %>% 
  gather(param, value, -DAYS, -visit) %>% 
  mutate(param = if_else(param == "BPD", "BPD (cm)",
                         if_else(param == "AC", "AC (cm)", 
                                 if_else(param == "EFW", "Estimated Fetal Weight (kg)","")))) %>% 
  ggplot() + 
  geom_point(aes(x = DAYS/7, y = value, color = as.factor(visit))) + 
  geom_vline(xintercept = 19, linetype = "dashed") +
  geom_vline(xintercept = 26, linetype = "dashed") +
  geom_vline(xintercept = 32, linetype = "dashed") +
  labs(title = "Ultrasound parameters and estimated fetal weight, by visit.", 
       y = "", x = "Estimated Gestational Age (weeks)", color = "Visit Week") +
  geom_text_repel(data = LBL_TBL, aes(x = DAYS, y = value, label = SubjectID)) +
  scale_x_continuous(breaks = c(10, 15, 19, 26, 32, 35, 40, 45), limits = c(15,40)) +
  facet_wrap(~param, scales = "free", ncol = 1) +
  theme(legend.position = "bottom")

us_plots

###########################################
# Final check for outliers
###########################################
params_max <- params %>% mutate(smax_BPD = mu_BPD + 8*sd_BPD,
                                smax_AC = mu_AC + 8*sd_AC,
                                smax_EFW = mu_EFW + 8*sd_EFW)

us_all_out <- left_join(us_all, params_max)

us_all_out %>% 
  ggplot() + geom_histogram(aes(BPD),binwidth = .05) + 
  geom_vline(aes(xintercept = min_BPD), linetype = "dashed") +
  geom_vline(aes(xintercept = max_BPD), linetype = "dashed") +
  geom_vline(aes(xintercept = smax_BPD), linetype = "solid") +
  facet_wrap(~visit, ncol = 1) + 
  labs(title = "BPD by visit",
       x = "BPD (cm)",
       caption = "dashed = Mean +/- 4SD; solid = Mean +/- 8SD")

us_all_out %>% 
  ggplot() + geom_histogram(aes(AC),binwidth = .05) + 
  geom_vline(aes(xintercept = min_AC), linetype = "dashed") +
  geom_vline(aes(xintercept = max_AC), linetype = "dashed") +
  geom_vline(aes(xintercept = smax_AC), linetype = "solid") +
  facet_wrap(~visit, ncol = 1)

us_all_out %>% 
  ggplot() + geom_histogram(aes(EFW),binwidth = .01) + 
  geom_vline(aes(xintercept = min_EFW), linetype = "dashed") +
  geom_vline(aes(xintercept = max_EFW), linetype = "dashed") +
  geom_vline(aes(xintercept = smax_EFW), linetype = "solid") +
  facet_wrap(~visit, ncol = 1)

us_all_out %>% filter(BPD < min_BPD) %>% select(SubjectID, BPD, min_BPD) # 2 okay, within 0.1
us_all_out %>% filter(AC < min_AC) %>% select(SubjectID, AC, min_AC) # 0
us_all_out %>% filter(EFW < min_EFW) %>% select(SubjectID, EFW, min_EFW) # 0

us_all_out %>% filter(BPD > smax_BPD) %>% select(SubjectID, BPD, smax_BPD) # 1 -- NO AC MEASURE, REMOVE POINT
us_all_out %>% filter(AC > smax_AC) %>% select(SubjectID, AC, smax_AC) # 0
us_all_out %>% filter(EFW > smax_EFW) #%>% select(SubjectID, EFW, smax_EFW) # 1 -- OKAY B/C "19 wk visit" actually at 25 wks


us_all <- us_all %>% mutate(BPD = if_else(SubjectID == "020-66132" & visit == 19, NA_real_, BPD))

ALL_ID <- c(OOR_ID, FLAG_ID, "020-66132")

write_dta(us_all, paste0(base_path,"20181211-us_efw.dta"))

final <- us_all %>% select(BPD, AC, EFW, DAYS, visit) %>% 
  gather(param, value, -DAYS, -visit) %>% 
  mutate(param = if_else(param == "BPD", "BPD (cm)",
                         if_else(param == "AC", "AC (cm)", 
                                 if_else(param == "EFW", "Estimated Fetal Weight (kg)","")))) %>% 
  ggplot() + 
  geom_point(aes(x = DAYS/7, y = value, color = as.factor(visit))) + 
  geom_vline(xintercept = 19, linetype = "dashed") +
  geom_vline(xintercept = 26, linetype = "dashed") +
  geom_vline(xintercept = 32, linetype = "dashed") +
  labs(title = "Ultrasound parameters and estimated fetal weight, by visit.", 
       y = "", x = "Estimated Gestational Age (weeks)", color = "Visit Week") +
  scale_x_continuous(breaks = c(10, 15, 19, 26, 32, 35, 40, 45), limits = c(15,40)) +
  facet_wrap(~param, scales = "free", ncol = 1) +
  theme(legend.position = "bottom")

final


######################################
# FINAL TALLY: N = 22 Subject IDs modified
######################################






######################################
######################################
######################################
# DEPRECATED
######################################
# raw_data <- read_dta("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/03 - COLLABORATION/Determinants of School Failure/Archive Data/File_for_Jon_2.dta")
# plot individual measures
# us_all %>% ggplot() + geom_histogram(aes(BPD),binwidth = .1) + facet_wrap(~visit, ncol = 1)
# us_all %>% ggplot() + geom_histogram(aes(AC),binwidth = .1) + facet_wrap(~visit, ncol = 1)
# us_all %>% ggplot() + geom_histogram(aes(EFW),binwidth = .1) + facet_wrap(~visit, ncol = 1)
# 
# us_all %>% group_by(visit) %>% summarize(mu_BPD = mean(BPD, na.rm = T), mu_AC = mean(AC, na.rm = T), mu_EFW = mean(EFW, na.rm = T),
#                                          sd_BPD = sd(BPD, na.rm = T), sd_AC = sd(AC, na.rm = T), sd_EFW = sd(EFW, na.rm = T),
#                                          min_BPD = mu_BPD - 4*sd_BPD, min_AC = mu_AC - 4*sd_AC, min_EFW = mu_EFW - 4*sd_EFW,
#                                          max_BPD = mu_BPD + 4*sd_BPD, max_AC = mu_AC + 4*sd_AC, max_EFW = mu_EFW + 4*sd_EFW)
# 
# us_all %>% filter(BPD > 13.2 | AC > 43.9 | EFW > 4.75) 
# us_all <- us_all %>% mutate(BPD = if_else(BPD > 13.2, NA_real_, BPD),
#                             EFW = if_else(EFW > 4.75, NA_real_, EFW))
# 
# params <- us_all %>% group_by(visit) %>% 
#   summarize(mu_AC = mean(AC, na.rm = TRUE), sd_AC = sd(AC, na.rm = TRUE),
#             mu_BPD = mean(BPD, na.rm = TRUE), sd_BPD = sd(BPD, na.rm = TRUE),
#             mu_EFW = mean(EFW, na.rm = TRUE), sd_EFW = sd(EFW, na.rm = TRUE)) %>% 
#   mutate(AC_min = mu_AC-4*sd_AC,
#          BPD_min = mu_BPD-4*sd_BPD,
#          EFW_min = mu_EFW-4*sd_EFW,
#          AC_max = mu_AC+4*sd_AC,
#          BPD_max = mu_BPD+4*sd_BPD,
#          EFW_max = mu_EFW+4*sd_EFW) %>% 
#   select(visit, contains("min"), contains("max"))
# 
# us_all %>% left_join(.,params) %>% group_by(visit) %>% count(AC < AC_min) 
# us_all %>% left_join(.,params) %>% group_by(visit) %>% count(AC > AC_max) 
# 
# us_all %>% left_join(.,params) %>% group_by(visit) %>% count(BPD < BPD_min)
# us_all %>% left_join(.,params) %>% group_by(visit) %>% count(BPD > BPD_max)
# 
# us_all %>% left_join(.,params) %>% group_by(visit) %>% count(EFW < EFW_min)
# us_all %>% left_join(.,params) %>% group_by(visit) %>% count(EFW > EFW_max)
# 
# us_all <- full_join(ultrasound19, full_join(ultrasound26, ultrasound32)) %>% left_join(.,params) %>% 
#   mutate(AC = if_else(AC < AC_min, NA_real_, 
#                       if_else(AC > AC_max, NA_real_, AC)),
#          BPD = if_else(BPD < BPD_min, NA_real_, 
#                        if_else(BPD > BPD_max, NA_real_, BPD)),
#          EFW = if_else(EFW < EFW_min, NA_real_, 
#                        if_else(EFW > EFW_max, NA_real_, EFW))) 
# 
# us_all %>% mutate(IVF = substr(SubjectID,1,3) %in% c("019","029")) %>% 
#   select(-SubjectID, -contains("min"), -contains("max")) %>% 
#   group_by(visit, IVF) %>% summarize_all(funs(mean(., na.rm=T)))
# 
# us_all %>% mutate(IVF = substr(SubjectID,1,3) %in% c("019","029")) %>% 
#   ggplot() + geom_boxplot(aes(x = as.factor(visit), y = EFW, color = IVF)) + 
#   #geom_jitter(aes(x = visit, y = EFW, color = IVF)) + 
#   theme(legend.position = "bottom")
# 
# write_dta(us_all, paste0(base_path,"20181203-us_efw.dta"))
# 
# # Saved as static covariate, i.e. can't plot time series
# params_all <- params %>% gather(measure, val, -visit) %>% 
#   mutate(var = paste0(measure,"_",as.character(visit))) %>% 
#   select(var, val) %>% spread(var, val)
# us_discrete_all <- left_join(us32, left_join(us26, us19)) %>% 
#   mutate(AC_19 = if_else(AC_19 > params_all$AC_max_19, NA_real_, if_else(AC_19 < params_all$AC_min_19, NA_real_, AC_19)),
#          BPD_19 = if_else(BPD_19 > params_all$BPD_max_19, NA_real_, if_else(BPD_19 < params_all$BPD_min_19, NA_real_, BPD_19)),
#          EFW_19 = if_else(EFW_19 > params_all$EFW_max_19, NA_real_, if_else(EFW_19 < params_all$EFW_min_19, NA_real_, EFW_19)),
#          AC_26 = if_else(AC_26 > params_all$AC_max_26, NA_real_, if_else(AC_26 < params_all$AC_min_26, NA_real_, AC_26)),
#          BPD_26 = if_else(BPD_26 > params_all$BPD_max_26, NA_real_, if_else(BPD_26 < params_all$BPD_min_26, NA_real_, BPD_26)),
#          EFW_26 = if_else(EFW_26 > params_all$EFW_max_26, NA_real_, if_else(EFW_26 < params_all$EFW_min_26, NA_real_, EFW_26)),
#          AC_32 = if_else(AC_32 > params_all$AC_max_32, NA_real_, if_else(AC_32 < params_all$AC_min_32, NA_real_, AC_32)),
#          BPD_32 = if_else(BPD_32 > params_all$BPD_max_32, NA_real_, if_else(BPD_32 < params_all$BPD_min_32, NA_real_, BPD_32)),
#          EFW_32 = if_else(EFW_32 > params_all$EFW_max_32, NA_real_, if_else(EFW_32 < params_all$EFW_min_32, NA_real_, EFW_32)))
# 
# us_discrete_all %>% count(AC_19 < params_all$AC_min_19 | AC_19 > params_all$AC_max_19) 
# us_discrete_all %>% count(AC_26 < params_all$AC_min_26 | AC_26 > params_all$AC_max_26) 
# us_discrete_all %>% count(AC_32 < params_all$AC_min_32 | AC_32 > params_all$AC_max_32) 
# 
# us_discrete_all %>% count(BPD_19 < params_all$BPD_min_19 | BPD_19 > params_all$BPD_max_19) 
# us_discrete_all %>% count(BPD_26 < params_all$BPD_min_26 | BPD_26 > params_all$BPD_max_26) 
# us_discrete_all %>% count(BPD_32 < params_all$BPD_min_32 | BPD_32 > params_all$BPD_max_32) 
# 
# us_discrete_all %>% count(EFW_19 < params_all$EFW_min_19 | EFW_19 > params_all$EFW_max_19) 
# us_discrete_all %>% count(EFW_26 < params_all$EFW_min_26 | EFW_26 > params_all$EFW_max_26) 
# us_discrete_all %>% count(EFW_32 < params_all$EFW_min_32 | EFW_32 > params_all$EFW_max_32) 
# 
# write_dta(us_discrete_all, paste0(base_path,"20181205-us_parms_all.dta"))                             
# 
# 
# 
# 
# 
# # plot growth by ethnicity and dropout (3 IVF dropped out prior to birth vs 82 non-IVF)
# # renumber to "fit" visit schedule and associate covar data (19 wk ga -> 0 month; 26 wks ga -> 1 mo; 32 wk ga -> 3 mo) -- quick and dirty
# merge_us <- us_all %>% mutate(visit = if_else(visit == 19, 0,
#                                               if_else(visit == 26, 1,
#                                                       if_else(visit == 32, 3, visit)))) %>% 
#   left_join(merged2, .)
# # count
# merge_us %>% filter(visit %in% c(0,1,3)) %>% group_by(visit) %>% 
#   count(mother_ethnicity, IVF, measured = !is.na(EFW)) %>% filter(measured == T)
# #plot
# merge_us %>% filter(visit %in% c(0,1,3), mother_ethnicity != "others") %>% 
#   mutate(visit = factor(visit, levels = c(0,1,3), labels = c("19 Weeks", "26 Weeks", "32 Weeks"))) %>% 
#   mutate(mother_ethnicity = if_else(mother_ethnicity == "", "Not followed at delivery", 
#                                     if_else(mother_ethnicity == "chinese", "Chinese",
#                                             if_else(mother_ethnicity == "indian", "Indian",
#                                                     if_else(mother_ethnicity == "malay", "Malay", ""))))) %>% 
#   ggplot() + geom_boxplot(aes(x = as.character(visit), y = EFW, color = as.factor(IVF))) + 
#   labs(x = "Gestational Age", 
#        y = "Estimated Fetal Weight (Shepard formula)",
#        color = "IVF status",
#        title = "Fetal size (estimated from ultrasound), by IVF, ethnicity, and inclusion status.") +
#   scale_color_discrete(labels = c("Non-IVF", "IVF")) +
#   facet_wrap(~mother_ethnicity) +
#   theme(legend.position = "bottom")
# 




## DEPRECATED
# fetal <- left_join(us26, us32, by = "SubjectID")
# fetal <- fetal %>% mutate(velocity = EFW_26.y - EFW_26.x)
# 
# fetal_long <- fetal %>% select(SubjectID, EFW_26.x, EFW_26.y, velocity) %>% 
#   mutate(IVF = substring(SubjectID,1,3) %in% c("019", "029"), 
#          decel = EFW_26.y < EFW_26.x) %>% 
#   gather(visit, EFW, -SubjectID, -IVF, -decel, -velocity) %>% 
#   mutate(visit = if_else(visit == "EFW_26.x", 26, 32))
# fetal_long %>% ggplot(aes(x = visit, y = EFW)) + geom_jitter() + facet_wrap(~decel)
# 
# fetal_long %>% ggplot() + geom_histogram(aes(EFW, fill = as.factor(visit)))
# fetal_long %>% filter(visit == 26) %>% ggplot() + geom_histogram(aes(velocity))
# 
# merged <- raw_data %>% rename(SubjectID = subjectid) %>% left_join(fetal)
# merged <- merged %>% rename(BPD_26 = BPD.x, AC_26 = AC.x, EFW_26 = EFW_26.x, BPD_32 = BPD.y, AC_32 = AC.y, EFW_32 = EFW_26.y)
# write_dta(merged, paste0(base_path,"20181128-growth_visit_add.dta"))
# 
# merged %>% group_by(gdm) %>% 
#   filter(!is.na(birthweight_grams)) %>% 
#   summarize(n = n(), mu_bw = mean(birthweight_grams, na.rm = T), mu_vel = mean(velocity, na.rm = T))
# 
# merged %>% #filter(!is.na(gdm)) %>% 
#   ggplot(aes(x = velocity, y = birthweight_grams)) + 
#   geom_point() + geom_smooth() + facet_wrap(~gdm, ncol = 1)


## LOOK AT U/S DATES
# too_small <- fetal %>% filter(velocity < 0.3) %>% 
#   select(SubjectID, velocity) %>% 
#   arrange(velocity)
# 
# date1 <- ultrasound26 %>% filter(PSCID %in% too_small[["SubjectID"]]) %>% select(PSCID, date26 = date_of_ultrasound_scan)
# date2 <- ultrasound32 %>% filter(PSCID %in% too_small[["SubjectID"]]) %>% select(PSCID, date32 = date_of_ultrasound_scan)
# 
# date1 <- ultrasound26 %>% select(PSCID, date26 = date_of_ultrasound_scan)
# date2 <- ultrasound32 %>% select(PSCID, date32 = date_of_ultrasound_scan)
# 
# diffs <- left_join(date1, date2) %>% mutate(time = as.double(date32 - date26)) 
# diffs %>% filter(time <= 0) # Two clearly incorrect dates
# 
# diffs %>% #filter(date32 > date26) %>% 
#   mutate(date32 = if_else(PSCID == "010-20053", as.Date("2010-02-09"), 
#                           if_else(PSCID == "019-30036", as.Date("2011-03-31"), date32)),
#          time = as.double(date32 - date26)) %>% #filter(time <= 0)
#   ggplot() + geom_histogram(aes(time), binwidth = 1)
# 
# merged <- diffs %>% 
#   mutate(date32 = if_else(PSCID == "010-20053", as.Date("2010-02-09"), 
#                           if_else(PSCID == "019-30036", as.Date("2011-03-31"), date32)),
#          time = as.double(date32 - date26)) %>% 
#   select(SubjectID = PSCID, date26, date32, time) %>% right_join(merged)
# 
# # No clear association between ultrasound timing (elapsed date, 26 week visit, or 32 week visit) and birthweight (good thing)
# merged %>% ggplot() + geom_point(aes(x = time, y = birthweight_grams))
# merged %>% ggplot() + geom_point(aes(x = date26, y = birthweight_grams))
# merged %>% ggplot() + geom_point(aes(x = date32, y = birthweight_grams))
# 
# merged %>% ggplot() + geom_point(aes(x = GA, y = birthweight_grams)) + geom_smooth(aes(x = GA, y = birthweight_grams))
# merged %>% ggplot() + geom_histogram(aes(GA))
# merged %>% ggplot(aes(x = velocity, y = birthweight_grams)) + 
#   geom_point(aes(x = velocity, y = birthweight_grams)) + 
#   geom_smooth(aes(color = "red"), formula = y ~ x, method = glm) +
#   geom_smooth() + theme(legend.position = "none")
# 
# quantile(merged$velocity, c(.25, .5, .75), na.rm = T)
# 
# vtree(fetal, "BPD.x BPD.y", check.is.na = T, horiz = F,
#       summary = c("BPD.x \nBPD_26 \nmean = %mean% \nmiss = %mv% %leafonly%",
#       "BPD.y \nBPD_32 \nmean = %mean% \nmiss = %mv% %leafonly%"))
# 
# vtree(merged, c("mother_ethnicity", "mother_highest_education"), 
#       check.is.na = T,
#       showempty = T, horiz = F, vp = F)
# 
# ########################################
# # Wide anthro data
# ########################################
# base_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"
# base_path <- "C:/Users/Jon Huang/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"
# 
# zbmi_annotate <- read_dta(paste0(base_path,"20181001-zBMI.dta"))
# zbmi_grow <- fetal %>% mutate(BPD_vel = BPD.y - BPD.x, AC_vel = AC.y - AC.x, EFW_vel = velocity) %>%
#   select(SubjectID, AC_vel, BPD_vel, EFW_vel) %>% 
#   right_join(zbmi_annotate)
# 
# # identify fetal growth outliers graphically (by plotting; N = 8)
# zbmi_grow %>% filter(visit == 0) %>% ggplot(aes(x = AC_vel, y = weight)) + geom_point() #(N = 7)
# zbmi_grow %>% filter(visit == 0) %>% ggplot(aes(x = BPD_vel, y = weight)) + geom_point() #(N = 4)
# zbmi_grow %>% filter(visit == 0, (AC_vel <= 0 | AC_vel > 15 | BPD_vel > 5)) %>% 
#   select(SubjectID, AC_vel, BPD_vel, EFW_vel) %>% arrange(EFW_vel)
# 
# # label them in a plot of birthweight ~ EFW_vel
# zbmi_grow %>% filter(visit == 0) %>% 
#   ggplot(aes(x = EFW_vel, y = weight, label = SubjectID)) + geom_point() + geom_smooth() + 
#   geom_text(data = subset(zbmi_grow, visit == 0 & (AC_vel <= 0 | AC_vel > 15 | BPD_vel > 5)), nudge_y = 0.11)
# 
# # plot non-outliers @ birth (N = 1090)
# zbmi_grow %>% filter(visit == 0, AC_vel > 0, AC_vel < 15, BPD_vel < 5) %>% 
#   ggplot(aes(x = EFW_vel, y = weight)) + geom_point() + geom_smooth()
# 
# EFW_quart <- quantile(zbmi_grow$EFW_vel, c(0.25, 0.5, 0.75), na.rm = T)
# 
# # BY 3rd trimester fetal growth velocity quartiles 
# zbmi_grow <- zbmi_grow %>% 
#   mutate(EFW_group = if_else(EFW_vel > EFW_quart[3], 4,
#                              if_else(EFW_vel > EFW_quart[2], 3,
#                                      if_else(EFW_vel > EFW_quart[1], 2, 1))))
# zbmi_grow %>% write_dta(paste0(base_path,"20181111-EFW_ZBMI.dta"))
# zbmi_grow <- read_dta(paste0(base_path,"20181111-EFW_ZBMI.dta"))
# 
#   
# # N = 105 are missing 3rd trimester velocity measures
# zbmi_grow %>% filter(visit == 0, !is.na(weight)) %>% 
#   group_by(EFW_group) %>% 
#   summarize(n = n(), mu_vel = mean(EFW_vel, na.rm = T))
# 
# # By ethnicity (total N = 1177; 1072 with EFW)
# zbmi_grow %>% filter(visit == 0, mother_ethnicity %in% c("chinese", "indian", "malay")) %>% group_by(mother_ethnicity) %>% summarize(n = n(), pct = n/1177)
# # Chinese = 662 (56.2%); Indian = 216 (18.4%); Malay = 299 (25.4%)
# zbmi_grow %>% filter(visit == 0, mother_ethnicity %in% c("chinese", "indian", "malay")) %>% group_by(mother_ethnicity, EFW_group) %>% summarize(n = n(), pct = n/1177)
# 
# 
# zbmi_grow %>% filter(visit != 0, mother_ethnicity %in% c("chinese", "indian", "malay")) %>% ggplot(aes(x = days, y = zbmi, group = EFW_group)) + 
#   geom_point(aes(color = as.factor(EFW_group)), alpha = 0.2) +
#   geom_smooth(aes(color = as.factor(EFW_group)), se = T, method = "loess", span = 0.25, alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dotted") + 
#   geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
#   geom_hline(yintercept = c(-2, 2), linetype = "longdash") +
#   labs(title = "Smoothed BMI z-scores (3 wks - 6.5 years), by quartiles of 3rd trimester growth velocity.",
#        color = "Quartiles of EFW velocity 26-32 weeks gestation",
#        x = "") +
#   scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7)) +
#   scale_x_continuous(breaks = c(0, 365, 731, 1096, 1461, 1826, 2192, 2557),
#                      labels = c("Birth", "1 year", "2 years", "3 years", 
#                                 "4 years", "5 years", "6 years", "7 years")) +
#   facet_wrap(~mother_ethnicity, nrow = 1) +
#   theme(legend.position = "bottom")       
# 
# # GDM treatment trajectories
# zbmi_grow %>% filter(visit != 0, mother_ethnicity %in% c("chinese", "indian", "malay")) %>% ggplot(aes(x = days, y = zbmi, group = gdm_treatment)) + 
#   geom_point(aes(color = as.factor(gdm_treatment)), alpha = 0.3) +
#   geom_smooth(aes(color = as.factor(gdm_treatment)), se = T, method = "loess", span = 0.5, alpha = 0.1) +
#   geom_hline(yintercept = 0, linetype = "dotted") + 
#   geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
#   geom_hline(yintercept = c(-2, 2), linetype = "longdash") +
#   labs(title = "Smoothed BMI z-scores (3 wks - 6.5 years), by GDM treatment modality.",
#        color = "GDM treatment modality",
#        x = "") +
#   scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7)) +
#   scale_x_continuous(breaks = c(0, 365, 731, 1096, 1461, 1826, 2192, 2557),
#                      labels = c("Birth", "1 year", "2 years", "3 years", 
#                                 "4 years", "5 years", "6 years", "7 years")) +
#   facet_wrap(~mother_ethnicity, nrow = 1) +
#   theme(legend.position = "bottom")     
