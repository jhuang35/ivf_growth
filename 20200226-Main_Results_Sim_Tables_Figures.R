#######################################################################################
#
# Smaller stature in childhood following assisted reproductive technologies (ART) 
# is not explained by parental or gestational factors, early childhood environment, 
# or fetal DNA methylation
#
# CORRESPONDING AUTHOR: Jon Huang
# CONTACT INFORMATION: jonathan_huang@sics.a-star.edu.sg
#
# VERSION DATE: 2020 February 26
#
# PURPOSE: 
#
# 1. Demonstrate main standard (GLM) and machine learning-based (C-TMLE) results 
#    using anonymized, structure preserving data*
#
# 2. Re-create manuscript tables and figures from model output data-files**
#
# ADDITIONAL DETAILS:
# 
# * All code for analyses will be included separated. However, since release of original data under GUSTO Study 
# data sharing agreements currently requires specific approval, we are including here an example dataset generated 
# by plasmode-like simulation is provided here to demonstrate the functionality of the statistical code.
# The advantage of plasmode simulation is that it preserves covariate and exposure relationships from the observed
# data while inserting a choosen effect size. 
# This allows evaluation of real-world performance of the desired statistical approaches in the actual data setting. 
# For the purposes of this release, missing covariate data were singly-imputed to provide a complete data set. 
# For the manuscript, results are summarized across multiple imputations using full-conditional specification for
# missing covariate data. 
# More details about plasmode simulation can be read in: 
# Franklin JM, Schneeweiss S, Polinski JM, Rassen JA. Plasmode simulation for the evaluation of 
# pharmacoepidemiologic methods in complex healthcare databases. Comput Stat Data Anal. 2014; 72: 219-226.
#
# ** Similarly, all code used to produce tables and figures will be provided, however, 
#    only summary results can be provided so descriptive tables will not be able to be generated
#
# Obtaining original data can be discussed by contacting the corresponding author 
#
#######################################################################################

# REQUIRED LIBRARIES

library(haven)
library(readxl)
library(tidyverse)
library(ctmle)
library(arm)
library(xgboost)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)
library(magick)
library(ggplotify)

#library(tufte)
#library(knitr)
#library(table1)

# SET PATH FOR DATA FILES (ALL INPUT FILES WILL BE RELATIVE TO THIS PATH)
data_path <- "~ROOT/data/" # replace "~ROOT/" with the correct path to the /data/ subdirectory

###################################################################################################
###################################################################################################
# SECTION 1. Demonstrate main standard (GLM) and machine learning-based (C-TMLE) results 
#    using anonymized, structure preserving data*
#
#    Exact code used for analyses will also be included separately
###################################################################################################
###################################################################################################

# Read in simulated data structure:
simdata <- readxl::read_xlsx(paste0(data_path, "simdata.xlsx"), col_types = c("text", "numeric", "numeric", "text", rep("numeric", 15)))

# 1178 subjects x 3 outcomes = 3534 rows
# 19 variables = 19 columns 
# 1 ID variable, 1 outcome indicator (OUTCOME), 1 outcome value (Y), 1 exposure (A1), 15 covariates (L0.a to L0.o)
# Note: because data are simulated, IDs are repeated because they are drawn/simulated more than once; this can be safetly ignored in analyses
#
# A1 = ART status (1 = IVF, 0 = spontaneous conception)
#
# Y = outcome value, corresponding to the standardized anthropometric indicated by the "OUTCOME" column
outcomes <- c("zlen", "zwei", "zbmi")
#
# Covariates (L0.): a = maternal age at delivery, b = maternal education (factor), c = maternal ethnicity (factor), d = household income (factor),
#                   e = maternal height, f = maternal pre-pregnancy BMI, g = parity, h = smoke exposure in pregnancy (factor), i = paternal height,
#                   j = paternal weight, k = child sex, l = polygenic risk score for obesity, m = paternal age at delivery, 
#                   n = paternal hypertension history, o = paternal diabetes history
#

# Correct treatment model equation 
gform <- "A1 ~ L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k + L0.l + L0.m + L0.n + L0.o"

# Correct outcome model equation
Qform <- "Y ~ A1 + L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j + L0.k + L0.l + L0.m + L0.n + L0.o"

#### Run multivariable regressions for effect of IVF (A1) on anthropometrics (Y)
reg_res <- NULL
for(i in 1:length(outcomes)){
  stmp <- filter(sim_final, OUTCOME == outcomes[i])
  pn <- summary(glm(data = stmp, formula = A1 ~ 1, family = "gaussian"))$coefficient[1,1]
  pd <- glm(data = stmp, formula = as.formula(gform), family = "binomial")
  stmp <- stmp %>% add_column(pd = pd$fitted.values) %>% 
    mutate(wt0 = 1, # unweighted
           wt1 = if_else(A1 == 1, pn/pd, (1-pn)/(1-pd))) %>% # Stabilized IPTW
    mutate(wt2 = case_when(wt1 > quantile(wt1, 0.95) ~ quantile(wt1, 0.95),
                           wt1 < quantile(wt1, 0.05) ~ quantile(wt1, 0.05),
                           T ~ wt1)) # Truncated IPTW
  res <- stmp %>% glm(formula = outForm, weight = wt2, family = "gaussian") %>% summary() # run truncated-weighted regression
  reg_res <- rbind(reg_res, cbind(EST = res$coefficients[2,1], 
                                  SE = res$coefficients[2,2],
                                  OUTCOME = outcomes[i]))
}
as_tibble(reg_res) %>% mutate_at(c("EST", "SE"), list(~ as.numeric(.))) %>% 
  mutate(LB = EST - 1.96*SE, UB = EST + 1.96*SE) %>% select(OUTCOME, EST, SE, LB, UB) # DISPLAY RESULTS

# The target effect sizes were set to -0.5 SD to align with the observed magnitude in the study
# Note that because the datasets are simulated stochastically the individual estimate from a given draw/sample will not match exactly: 
# this is only meant to be illustrative of the method and compare efficiency against the machine-learning model

# USING A MULTIVARIABLE REGRESSION MODEL, THE OBSERVED EFFECT OF IVF ON OUTCOMES IN THE SIMULATED DATA SET IS AS FOLLOWS:
# Mean Height-for-age Z-score (zlen): -0.506 (SE = 0.152)
# Mean Weight-for-age Z-score (zwei): -0.689 (SE = 0.151)
# Mean BMI-for-age Z-score (zbmi): -0.589 (SE = 0.140)


#### Run Collaborative-TMLE for effect of IVF (A1) on anthropometrics (Y)

# Specify the learning algorithms to be used by the ensemble learner (SuperLearner)
SL.lib <- c("SL.glm", "SL.glm.interaction", "SL.bayesglm", "SL.mean", "SL.nnet", "SL.xgboost", "SL.glmnet", "SL.polymars")

# Specify the covariates
Lnodes <- c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k", "L0.l", "L0.m", "L0.n", "L0.o")

set.seed(12345)
ctmle_results <- NULL
  for(y in 1:length(outcomes)){
    print(paste0("Estimating effects on ", outcomes[y], "..."))
    subset <- simdata %>% filter(OUTCOME == outcomes[y]) %>% dplyr::select(contains("L0"),A1, Y) 
    W <- subset %>% dplyr::select(Lnodes)
    A <- subset %>% dplyr::select(A1) %>% pull(A1)
    Y <- subset %>% dplyr::select(Y) %>% pull(Y)
    model <- ctmleDiscrete(Y = Y, A = A, W = data.frame(W), V = 5, 
                           SL.library = SL.lib,
                           preOrder = FALSE, detailed = TRUE)
    ATE <- model$est
    SE <- model$var.psi^0.5
    OUT = outcomes[y]
    ctmle_results <- rbind(ctmle_results, cbind(ATE, SE, OUT))
  }
as_tibble(ctmle_results) %>% mutate_at(c("ATE", "SE"), list(~ as.numeric(.))) %>% 
  mutate(LB = ATE - 1.96*SE, UB = ATE + 1.96*SE) %>% select(OUT, ATE, SE, LB, UB) # DISPLAY C-TMLE RESULTS

# USING COLLABORATIVE TARGETED MAXIMUM LIKELIHOOD ESTIMATION, THE OBSERVED EFFECT OF IVF ON OUTCOMES IN THE SIMULATED DATA SET IS AS FOLLOWS:
# Mean Height-for-age Z-score (zlen): -0.405 (SE = 0.015)
# Mean Weight-for-age Z-score (zwei): -0.441 (SE = 0.015)
# Mean BMI-for-age Z-score (zbmi): -0.478 (SE = 0.017)
#
# *Note the increased consistency and overall smaller bias and variance in the C-TMLE models
  
#### ALTERNATIVE: cross-validated-TMLE (CV-TMLE) fitting procedure
# library(tmle3)
# library(sl3)
# 
# SL.lib <- c("SL.glm", "SL.glm.interaction", "SL.bayesglm", "SL.mean", "SL.nnet", "SL.xgboost", "SL.glmnet", "SL.polymars")
# # initialize the structure equation model
# npsem <- list(define_node("Z", 
#                           c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k",
#                             "L0.l", "L0.m", "L0.n", "L0.o")),
#               define_node("A", c("A1"), c("Z")),
#               define_node("Y", c("YT"), c("A", "Z")))
# # initialize the superlearner library
# lrnr_SL <- make_learner(Lrnr_pkg_SuperLearner, SL.lib)
# set.seed(12345)
# cvtmle_results <- NULL
# for(y in 1:length(outcomes)){
#   subset <- simdata %>% filter(OUTCOME == outcomes[y]) %>% dplyr::select(contains("L0"),A1, Y) %>% 
#     mutate(YT = (Y-min(Y))/(max(Y)-min(Y))) %>% data.frame(.) # generate a bounded Y for TMLE
#   tmle_task <- tmle3_Task$new(subset, npsem = npsem)
#   factor_list <- list(define_lf(LF_emp, "Z"), define_lf(LF_fit, "A", lrnr_SL), define_lf(LF_fit, "Y", lrnr_SL))
#   likelihood_def <- Likelihood$new(factor_list)
#   likelihood <- likelihood_def$train(tmle_task)
#   ate_params <- list(Param_ATE$new(likelihood, define_lf(LF_static, "A", value = 1), 
#                                    define_lf(LF_static, "A", value = 0)))
#   updater <- tmle3_Update$new()
#   targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
#   tmle3_fit <- fit_tmle3(tmle_task, targeted_likelihood, ate_params, updater)
#   ATE <- tmle3_fit$summary$tmle_est * (max(tset$Y5)-min(tset$Y5)) # backtransform ATE
#   SE <- tmle3_fit$summary$se * (max(tset$Y5)-min(tset$Y5)) # backtransform SE
#   OUT = outcomes[y]
#   cvtmle_results <- rbind(cvtmle_results, cbind(ATE, SE, OUT))
# }
  
  
###################################################################################################
###################################################################################################
#
# SECTION 2. Re-create manuscript tables and figures from output data-files
#
###################################################################################################
###################################################################################################

# Input core data and format variables
# (Raw data files not currently available)
# merged <- read_dta(paste0(data_path,"20181219-IVF-final.dta"))
# sub_id <- read_dta(paste0(data_path,"20190107-IVF-final-SUB.dta")) %>% filter(visit == 0) %>% pull(SubjectID)
# ewas_sub <- read_xls(paste0(data_path,"20190110-cpg-results-sub.xls"))
# dad_hx <- read_dta(paste0(data_path,"dad_hx.dta"))
  merged <- left_join(merged, dad_hx)
  DAD_COND <- merged %>% mutate(IVF = factor(substring(SubjectID,1,3) %in% c("019","029"))) %>% 
  filter(IVF == F, visit == 0) %>% 
  filter((father_age_delivery > 40 & parity == 0) | 
           (f_weight_m24/(f_height_m24/100)^2 > 35 & parity == 0) |
           father_diabetes == 1 | 
           father_HBP == 1) %>% pull(SubjectID)

# display helper function 
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3, digits.pct = 1), 
       c("","Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

# format data for covariates / parental characteristics
formatted <- merged %>% filter(visit == 0, sex != "") %>% 
  mutate(mom_age = mother_age_delivery,
         mom_edu = factor(mother_highest_education, 
                          levels = c("no_education", "primary", "secondary", "ite_ntc", "gce", "university"),
                          labels = c("No Education", "Primary", "Secondary", "ITE / NTC", "GCE A-Level", "University")),
         ethnic = factor(mother_ethnicity,
                         levels = c("chinese", "malay", "indian", "others"),
                         labels = c("Chinese", "Malay", "Indian", "Other")),
         hh_inc = factor(household_income,
                         levels = c("0_999", "1000_1999", "2000_3999", "4000_5999", "more_than_6000"),
                         labels = c("0 to 999", "1000 to 1999", "2000 to 3999", "4000 to 5999", "6000+")),
         nullip = (parity == 0),
         smk_home = factor(smk_home,
                           levels = c(0,1), labels = c("No", "Yes")),
         dad_dm = factor(father_diabetes,
                         levels = c(0,1), labels = c("No", "Yes")),
         dad_htn = factor(father_HBP,
                         levels = c(0,1), labels = c("No", "Yes")),
         htn = if_else(HTN_clean %in% c("Chronic HTN superimposed PE", "Eclampsia", "PE", "PIH"), 1, if_else(is.na(hi_bp), NA_real_, 0)),
         htn = factor(htn,
                      levels = c(0,1), labels = c("No", "Yes")))

label(formatted$mom_age) <- "Mother's Age at Delivery (Years)"
label(formatted$mom_edu) <- "Mother's Highest Educational Qualification"
label(formatted$ethnic) <- "Mother's Ethnicity"
label(formatted$hh_inc) <- "Monthly Household Income (SGD)"
label(formatted$nullip) <- "Nulliparous"
label(formatted$m_height_pw26) <- "Mother's Height (cm)"
label(formatted$ppBMI) <- "Mother's pre-pregnancy BMI (kg/m^2)"
label(formatted$father_age_delivery) <- "Father's Age at Delivery (Years)"
label(formatted$f_height_m24) <- "Father's Height (cm)"
label(formatted$f_weight_m24) <- "Father's Weight (kg)"
label(formatted$dad_dm) <- "Paternal Diabetes"
label(formatted$dad_htn) <- "Paternal High Blood Pressure History"
label(formatted$smk_home) <- "Any Smoking in Home (During Pregnancy)"
label(formatted$ogtt_fasting_pw26) <- "Maternal Fasting Glucose (mmol/L)"
label(formatted$ogtt_2hour_pw26) <- "Maternal 2-hour post-OGTT Glucose (mmol/L)"
label(formatted$htn) <- "Maternal High Blood Pressure History"

formatted_sub <- formatted %>% filter(SubjectID %in% sub_id)
label(formatted_sub$mom_age) <- "Mother's Age at Delivery (Years)"
label(formatted_sub$mom_edu) <- "Mother's Highest Educational Qualification"
label(formatted_sub$ethnic) <- "Mother's Ethnicity"
label(formatted_sub$hh_inc) <- "Monthly Household Income (SGD)"
label(formatted_sub$nullip) <- "Nulliparous"
label(formatted_sub$m_height_pw26) <- "Mother's Height (cm)"
label(formatted_sub$ppBMI) <- "Mother's pre-pregnancy BMI (kg/m^2)"
label(formatted_sub$father_age_delivery) <- "Father's Age at Delivery (Years)"
label(formatted_sub$f_height_m24) <- "Father's Height (cm)"
label(formatted_sub$f_weight_m24) <- "Father's Weight (kg)"
label(formatted_sub$dad_dm) <- "Paternal Diabetes"
label(formatted_sub$dad_htn) <- "Paternal High Blood Pressure History"
label(formatted_sub$smk_home) <- "Any Smoking in Home (During Pregnancy)"
label(formatted_sub$ogtt_fasting_pw26) <- "Maternal Fasting Glucose (mmol/L)"
label(formatted_sub$ogtt_2hour_pw26) <- "Maternal 2-hour post-OGTT Glucose (mmol/L)"
label(formatted_sub$htn) <- "Maternal High Blood Pressure History"

formatted <- formatted %>% 
  mutate(IVF = factor(IVF, levels = c(0, 1),
                      labels = c("Spontaneous Conception, All", "IVF")))

formatted_sub <- formatted_sub %>% 
  mutate(IVF = factor(IVF, levels = c(0, 1),
                      labels = c("Spontaneous, Possibly subfertile", "IVF")))

# format data for birth / child outcomes
child <- merged %>% filter(sex != "") %>% 
  dplyr::select(weight, height, triceps, subscapular, suprailiac, biceps, SBP, DBP,
         visit, SubjectID, sex, IVF, mode_of_delivery, GA, Dur_full_BF, Dur_any_BF, glucose_1) %>% 
  gather(measure, val, -(visit:glucose_1)) %>% unite(tmp, measure, visit) %>% 
  spread(tmp, val) %>% 
  mutate(mod = substr(mode_of_delivery,1,1),
         mod = if_else(mod %in% c("1","2","3"), "Vaginal", 
                       if_else(mod %in% c("4", "5"), "Caesarean", NA_character_)),
         Full_BF = factor(Dur_full_BF, 
                          levels = c("", "lt_1M", "1M_to_lt_3M", "3M_to_lt_6M", "6M_to12M", "12M_and_above"),
                          labels = c(NA_character_, "< 1", "1 to <3", "3 to <6", "6 to 12", "12+")))

label(child$sex) <- "Child Sex"
label(child$GA) <- "Gestational Age @ Birth (weeks)"
label(child$mod) <- "Mode of Delivery"
label(child$glucose_1) <- "Fasting Glucose @ 6 yr (mmol/L)"
label(child$weight_0) <- "Birth Weight (kg)"
label(child$weight_12) <- "Weight @ 1 yr (kg)"
label(child$weight_24) <- "Weight @ 2 yr (kg)"
label(child$weight_36) <- "Weight @ 3 yr (kg)"
label(child$weight_48) <- "Weight @ 4 yr (kg)"
label(child$weight_60) <- "Weight @ 5 yr (kg)"
label(child$weight_72) <- "Weight @ 6 yr (kg)"
label(child$height_0) <- "Birth Length (cm)"
label(child$height_12) <- "Length @ 1 yr (cm)"
label(child$height_24) <- "Height @ 2 yr (cm)"
label(child$height_36) <- "Height @ 3 yr (cm)"
label(child$height_48) <- "Height @ 4 yr (cm)"
label(child$height_60) <- "Height @ 5 yr (cm)"
label(child$height_72) <- "Height @ 6 yr (cm)"
label(child$SBP_72) <- "Systolic BP @ 6 yr (mmHg)"
label(child$DBP_72) <- "Diastolic BP @ 6 yr (mmHg)"
label(child$Full_BF) <- "Months Exclusively Breastfed"

child <- child %>% 
  mutate(IVF = factor(IVF, levels = c(0, 1),
                      labels = c("Spontaneous Conception, All", "IVF")))

one_cond <- add_data %>% filter(n_cond == 1 | IVF == TRUE) %>% pull(SubjectID)
alt_sub <- add_data %>% filter((num_miss == 1 | IVF == TRUE | SubjectID %in% sub_id) &
                                 num_miss < 2) %>% pull(SubjectID)
inv_sub <- merged %>% filter(visit == 0, !(SubjectID %in% sub_id)) %>% pull(SubjectID)
set.seed(42782)
rand_sub <- tibble(sub_id = inv_sub, choose = rbinom(length(inv_sub), 1, 0.2)) %>% filter(choose == 1) %>% pull(sub_id)

sub_ivf <- child %>% filter(IVF == "IVF") %>% mutate(IVF = 0) # Group = 0
sub_org <- child %>% filter(SubjectID %in% sub_id, IVF != "IVF") %>% mutate(IVF = 1) # Group = 1
sub_both <- child %>% filter(SubjectID %in% c(DAD_COND, sub_id),  IVF != "IVF") %>% mutate(IVF = 2) # Group = 2
sub_dad <- child %>% filter(SubjectID %in% DAD_COND,  IVF != "IVF") %>% mutate(IVF = 3) # Group = 3
sub_inv <- child %>% filter(SubjectID %in% rand_sub, IVF != "IVF") %>% mutate(IVF = 4) # Group = 4

child_sub <- bind_rows(sub_ivf, sub_org, sub_both, sub_dad, sub_inv) %>% 
  mutate(IVF = factor(IVF, levels = c(0, 1, 2, 3, 4),
                      labels = c("IVF",
                                 "Spontaneous (infertility indications only)",
                                 "Spontaneous (Indications + paternal risk factors)",
                                 "Spontaneous (paternal risk factors only)",
                                 "Spontaneous (no indications, 20% random sample)")))

label(child_sub$sex) <- "Child Sex"
label(child_sub$GA) <- "Gestational Age @ Birth (weeks)"
label(child_sub$mod) <- "Mode of Delivery"
label(child_sub$glucose_1) <- "Fasting Glucose @ 6 yr (mmol/L)"
label(child_sub$weight_0) <- "Birth Weight (kg)"
label(child_sub$weight_12) <- "Weight @ 1 yr (kg)"
label(child_sub$weight_24) <- "Weight @ 2 yr (kg)"
label(child_sub$weight_36) <- "Weight @ 3 yr (kg)"
label(child_sub$weight_48) <- "Weight @ 4 yr (kg)"
label(child_sub$weight_60) <- "Weight @ 5 yr (kg)"
label(child_sub$weight_72) <- "Weight @ 6 yr (kg)"
label(child_sub$height_0) <- "Birth Length (cm)"
label(child_sub$height_12) <- "Length @ 1 yr (cm)"
label(child_sub$height_24) <- "Height @ 2 yr (cm)"
label(child_sub$height_36) <- "Height @ 3 yr (cm)"
label(child_sub$height_48) <- "Height @ 4 yr (cm)"
label(child_sub$height_60) <- "Height @ 5 yr (cm)"
label(child_sub$height_72) <- "Height @ 6 yr (cm)"
label(child_sub$SBP_72) <- "Systolic BP @ 6 yr (mmHg)"
label(child_sub$DBP_72) <- "Diastolic BP @ 6 yr (mmHg)"
label(child_sub$Full_BF) <- "Months Exclusively Breastfed"

# child 6 year BP and chemistries
# lab6_clean <- read_csv(paste0(data_path,"lab6_clean.csv"))
child_out <- merged %>% filter(visit == 72) %>% 
  left_join(., lab6_clean, by = "SubjectID") %>% 
  mutate(HOMA_IR = (INS * glucose_1)/22.5,
         HOMA_B = (20 * INS)/(glucose_1 - 3.5)) %>% 
  mutate(HOMA_B = if_else(glucose_1 <= 3.5, NA_real_, HOMA_B))

###############################################
# TABLES
###############################################
# TABLE 1 - Overall parental descriptives
# (Raw data not available for Table 1)
# table1(~ mom_age + nullip + ethnic + mom_edu + hh_inc + 
#          m_height_pw26 + ppBMI + ogtt_fasting_pw26 + ogtt_2hour_pw26 + htn +
#          smk_home + father_age_delivery + f_height_m24 + f_weight_m24 + dad_dm + dad_htn | IVF, 
#        data = formatted, overall  = "Overall",
#        topclass = "Rtable1-zebra",
#        render.continuous = my.render.cont) 

# TABLE 2 - Overall child outcomes @ 6 years, adjusted models
# (Regression output file included - 22 Mar 2019-IVF-size_MI.xlsx)

###############################################
# FIGURES
###############################################
# Figure 1 - Describe
# (Raw data not available for Figure 1)
# merged %>% mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART", "Spontaneous Conception"))) %>% 
#   #filter(SubjectID %in% sub_id) %>% 
#   dplyr::select(IVF, visit, days, height, weight, bmi, zlen, zwei, zbmi) %>% 
#   gather(MEASURE, val, -IVF, -visit, -days) %>% 
#   mutate(TRANS = case_when(IVF == "IVF" ~ 0.5, TRUE ~ 0.08)) %>% 
#   mutate(MEASURE = case_when(MEASURE == "height" ~ "A.~Height~(cm)",
#                              MEASURE == "weight" ~ "B.~Weight~(kg)",
#                              MEASURE == "bmi" ~ "C.~BMI~(kg/m^2)",
#                              MEASURE == "zlen" ~ "D.~'Height-for-age Z-Score'~(SDS)",
#                              MEASURE == "zwei" ~ "E.~'Weight-for-age Z-Score'~(SDS)",
#                              MEASURE == "zbmi" ~ "F.~'BMI Z-Score'~(SDS)")) %>% 
#   # ggplot(aes(x = days/30.5, y = val, color = IVF)) + 
#   # geom_point(aes(alpha = TRANS)) + guides(alpha = FALSE) + scale_alpha_identity() +
#   # geom_smooth() +
#   # scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 36, 48, 60, 72, 78)) +
#   ggplot(aes(x = as.factor(visit), y = val, color = IVF)) +
#   geom_quasirandom(dodge.width = 0.75, alpha = 0.3) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.75) +
#   labs(x = "Visit Month", y = "", color = NULL) + 
#   facet_wrap(~MEASURE, scales = "free", labeller = label_parsed) +
#   theme(legend.position = c(0.11, 0.96), legend.background = element_blank())

# Figure 2 - Anthro
reg_out_new <- read_excel(paste0(data_path,"22 Mar 2019-IVF-size_MI.xls"), sheet = "min")
out_all_visit <- reg_out_new %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei"))
all_visit = rep(c(0, 1, 3, 6, 9, 12, 15, 18, 24, 36, 48, 54, 60, 66, 72, 78), 6)
out_all_visit <- out_all_visit %>% add_column(all_visit)
out_all_visit %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  group_by(outcome) %>% summarize(min = min(N), max = max(N))
out_all_visit %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  mutate(outcome = case_when(outcome == "height" ~ "A.~Height~(cm)",
                             outcome == "log_weight" ~ "B.~Weight~(log(kg))",
                             outcome == "log_bmi" ~ "C.~BMI~(log(kg/m^2))",
                             outcome == "zlen" ~ "D.~'Height-for-age'~'Z-score'~(SDS)",
                             outcome == "zwei" ~ "E.~'Weight-for-age'~'Z-score'~(SDS)",
                             outcome == "zbmi" ~ "F.~BMI~'Z-score'~(SDS)")) %>% 
  ggplot() + 
  geom_point(aes(x = all_visit, y = estimate)) +
  geom_errorbar(aes(x = all_visit, ymin = min95, ymax = max95)) + 
  geom_ribbon(aes(x = all_visit, ymin = min95, ymax = max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Age (Months)", y = "") + theme(legend.position = "none") +
  facet_wrap(~outcome, scale = "free", labeller = label_parsed)

# Figure 3 - bp and skinfolds
BP <- reg_out_new %>% filter(outcome %in% c("SBP", "ZSBP", "DBP", "ZDBP")) %>% 
  add_column(BP_visit = rep(c(36, 48, 60, 72), 4)) %>% 
  mutate(outcome = case_when(outcome == "SBP" ~ "A. Systolic, absolute (mmHg)", 
                             outcome == "DBP" ~ "B. Diastolic, absolute (mmHg)",
                             outcome == "ZSBP" ~ "C. Systolic, standardized (%ile rank)",
                             outcome == "ZDBP" ~ "D. Diastolic, standardized (%ile rank)")) %>% 
  ggplot() + geom_point(aes(x = BP_visit, y = estimate)) +
  geom_errorbar(aes(x = BP_visit, ymin = min95, ymax = max95), width = 1) + 
  geom_ribbon(aes(x = BP_visit, ymin = min95, ymax = max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(limits = c(-20, 15)) +
  scale_x_continuous(breaks = seq(36,72, 12)) +
  facet_wrap(~ outcome) + theme(legend.position = "none") +
  labs(x = "Age (Months)", y = "Absolute (mmHg) or standardized (%ile rank) 
       difference, (IVF - SC)")

SF <- reg_out_new %>% filter(outcome %in% c("log_sub", "log_tri", "log_bi", "log_supra")) %>% 
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
  geom_errorbar(aes(x = skin_visit, ymin = min95, ymax = max95)) + 
  geom_ribbon(aes(x = skin_visit, ymin = min95, ymax = max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-30,10)) +
  scale_x_continuous(breaks = c(0, 18, 24, 36, 48, 54, 60, 66, 72, 78))+
  labs(x = "", y = "% difference in skinfold thickness
       (IVF - SC)") + theme(legend.position = "none") +
  facet_wrap(~outcome)

ggarrange(SF, BP, labels = c("A", "B"), ncol = 1, nrow = 2)

# Figure 4 - DNA methylation plots
# MANHATTAN PLOT
annotate <- read_csv(paste0(data_path, "annotate_results.csv"))
FLAG <- annotate %>% filter(lpval > -log10(0.05))
flag_id <- FLAG %>% arrange(-lpval) %>% pull(IlmnID)
SIG <- annotate %>% filter(lpval > -log10(0.05/281))
sig_id <- SIG %>% arrange(-lpval) %>% pull(IlmnID)
LABELS <- annotate %>% filter((str_detect(UCSC_RefGene_Name,"HIF") |  
                                 str_detect(UCSC_RefGene_Name,"NECAB")) & pval2 < 0.05) %>% 
  dplyr::select(lpval, IlmnID, UCSC_RefGene_Name, CHR) %>% arrange(-lpval) %>% 
  mutate(name = if_else(substring(UCSC_RefGene_Name,1,1) == "N",
                        substring(UCSC_RefGene_Name,1,6),
                        substring(UCSC_RefGene_Name,1,5)))

mp <- annotate %>% ggplot() + 
  geom_point(aes(x = CHR, y = lpval, color = IlmnID %in% flag_id)) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05/281), linetype = "dashed") +
  geom_hline(yintercept = -log10(3.6e-8), linetype = "solid") +
  geom_text_repel(data = LABELS, aes(x = CHR, y = lpval, 
                                     label = paste0(IlmnID, " (",name,")")), point.padding = 1) + 
  labs(x = "", y = "-log(p)") +
  scale_x_continuous(breaks = c(1:22)) +
  theme(legend.position = "none")

p1 <- read_xlsx(paste0(data_path,"cpg_ethnicity.xlsx")) %>% 
  ggplot(aes(x = IVF, y = val)) + 
  geom_quasirandom(alpha = 0.3, aes(fill = IVF, color = IVF)) +
  geom_boxplot(alpha = 0.5, aes(fill = IVF, color = IVF)) +
  #geom_signif(comparisons = list(c("Spontaneous", "IVF")), test = "t.test", map_signif_level = T) +
  facet_grid( ~ IlmnID) +
  labs(group = "IVF Status", fill = "IVF Status", color = "IVF Status", x = "", y = "") +
  theme(legend.position = "none", axis.text = element_blank(), strip.text = element_text(face = "italic"))

p2 <- read_xlsx(paste0(data_path,"cpg_ethnicity.xlsx")) %>% 
  ggplot(aes(x = IVF, y = val)) + 
  geom_quasirandom(alpha = 0.3, aes(fill = IVF, color = IVF)) +
  geom_boxplot(alpha = 0.5, aes(fill = IVF, color = IVF)) +
  #geom_signif(comparisons = list(c("Spontaneous", "IVF")), test = "t.test", map_signif_level = T, vjust = 1.5) +
  facet_grid(mother_ethnicity ~ IlmnID) +
  labs(group = "IVF Status", fill = "IVF Status", color = "IVF Status", x = "", y = "") +
  theme(legend.position = "none", axis.text = element_blank(), strip.text = element_text(face = "italic"))

ggarrange(mp,
          ggarrange(p1, p2, labels = c("B", "C"), ncol = 2, nrow = 1), 
          labels = c("A"), ncol = 1, nrow = 2)

# Figure 5 - decomposition of NECAB3, in subcohort
decomp_anthro_necab <- read_excel(paste0(data_path, "12 Jul 2019-gformula_SUB_NEC.xls")) %>% 
  rename(TCE_b = gcomp1, TCE_se = gcomp2, NDE_b = gcomp3, NDE_se = gcomp4,
         NIE_b = gcomp5, NIE_se = gcomp6, CDE_b = gcomp7, CDE_se = gcomp8, visit = gcomp9) %>% 
  add_column(measure = rep(c("zlen", "height", "zwei", "log_weight", "log_bmi", "zbmi"),16)) %>% 
  gather(var, val, -visit, -measure) %>% separate(var, c("effect", "param")) %>% 
  spread(key = param, val = val) %>% 
  mutate(ub = b + 1.96*se, lb = b - 1.96*se) %>% 
  mutate(te_ub = if_else(effect == "TCE", ub, NA_real_), te_lb = if_else(effect == "TCE", lb, NA_real_),
         nde_ub = if_else(effect == "NDE", ub, NA_real_), nde_lb = if_else(effect == "NDE", lb, NA_real_),
         nie_ub = if_else(effect == "NIE", ub, NA_real_), nie_lb = if_else(effect == "NIE", lb, NA_real_),
         cde_ub = if_else(effect == "CDE", ub, NA_real_), cde_lb = if_else(effect == "CDE", lb, NA_real_)) %>%
  mutate(measure = case_when(measure == "height" ~ "A. Height (cm)",
                             measure == "log_weight" ~ "B. Weight (% change)",
                             measure == "log_bmi" ~ "C. BMI (% change)",
                             measure == "zlen" ~ "D. Height-for-age (SDS)",
                             measure == "zbmi" ~ "F. BMI Z-score (SDS)", 
                             measure == "zwei" ~ "E. Weight-for-age (SDS)")) 
nie_bounds <- decomp_anthro_necab %>% filter(effect == "NIE")
nec <- decomp_anthro_necab %>% ggplot(aes(x = visit)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_point(aes(y = if_else(effect == "NDE", b, NA_real_))) +
  geom_errorbar(aes(ymin = nde_lb, ymax = nde_ub)) +
  geom_point(aes(y = if_else(effect == "NIE", b, NA_real_)), shape = 18) +
  geom_errorbar(aes(ymin = nie_lb, ymax = nie_ub), linetype = "dashed") + 
  geom_ribbon(data = nie_bounds, aes(x = visit, ymin = lb, ymax = ub), alpha = 0.15) +
  labs(x = "Age (Months)", y = "") +
  facet_wrap(~measure, scales = "free") + 
  scale_x_continuous(breaks = seq(0,72,12))
nec

###############################################
# SUP TABLES
###############################################
# SUP TABLE 1 - Subcohort descriptives
# (Raw data not available for )
# table1(~ mom_age + nullip + ethnic + mom_edu + hh_inc + 
#          m_height_pw26 + ppBMI + ogtt_fasting_pw26 + ogtt_2hour_pw26 + htn +
#          smk_home + father_age_delivery + f_height_m24 + f_weight_m24 + dad_dm + dad_htn | IVF, 
#        data = formatted_sub, overall  = "Overall",
#        topclass = "Rtable1-zebra",
#        render.continuous = my.render.cont) 
# 
# # SUP TABLE 2 - Child subcohort descriptives
# table1(~ sex + GA + mod + weight_0 + height_0 + Full_BF +
#          weight_12 + weight_24 + weight_60 + weight_72 + 
#          height_12 + height_24 + height_36 + height_48 + height_60 + height_72 + 
#          glucose_1 + SBP_72 + DBP_72 | IVF, 
#        data = child_sub, overall  = NULL,
#        droplevels = T,
#        topclass = "Rtable1-zebra",
#        render.continuous = my.render.cont)

# SUP TABLE 3 - child QMR and biomarkers at 6 years, adjusted
# (Regression output file included -  11 Jul 2019-IVF-chem_6_MI.xls)

# SUP TABLE 4 - child health @ 12 months
# (Raw data currently unavailable)
# CHILD_PATH <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/Raw LORIS data/GUSTO child health/"
# child_12 <- read_csv(paste0(CHILD_PATH, "Child_Questionnaire_Month12.csv"))
# child_12 <- child_12 %>%   mutate(IVF = substring(PSCID,3,3) == "9") %>% 
#   mutate(n_dx = case_when(!is.na(qn_1_27_yes_SN4_diagnosis) ~ 4,
#                           !is.na(qn_1_27_yes_SN3_diagnosis) ~ 3,
#                           !is.na(qn_1_27_yes_SN2_diagnosis) ~ 2,
#                           !is.na(qn_1_27_yes_SN1_diagnosis) ~ 1,
#                           TRUE ~ NA_real_))
# child_12 <- child_12 %>% filter(!is.na(qn_1_13_cough)) %>% 
#   mutate(IVF = factor(IVF, levels = c(TRUE, FALSE), labels = c("IVF", "Spontaneous Conception"))) %>% 
#   mutate(any_fevers = case_when(qn_1_23_fever_vaccination == "0_no" ~ "No",
#                                 qn_1_23_fever_vaccination == "1_yes" ~ "Yes",
#                                 TRUE ~ NA_character_)) %>% 
#   mutate_at(c("qn_1_17_diarrhoea", "qn_1_25_admission", "qn_1_26", "qn_1_27"),
#             list(~ case_when(
#               . == "0_no" ~ "No",
#               . == "1_yes" ~ "Yes",
#               TRUE ~ NA_character_))) %>% 
#   mutate(qn_1_26_yes_SN1_duration = 
#            if_else(qn_1_26_yes_SN1_duration %in% c("dnk", "not_answered", NA_character_), 
#                    NA_real_, as.double(qn_1_26_yes_SN1_duration))) %>% 
#   mutate_at(c("qn_1_26_yes_SN2_duration", "qn_1_26_yes_SN3_duration", "qn_1_26_yes_SN4_duration"),
#             list(~ if_else(is.na(.), 0, as.double(.)))) %>% 
#   mutate(abx_days_12 = qn_1_26_yes_SN1_duration + 
#            qn_1_26_yes_SN2_duration + 
#            qn_1_26_yes_SN3_duration + 
#            qn_1_26_yes_SN4_duration)
# label(child_12$qn_1_17_diarrhoea) <- "Diarrhea lasting 2+ days"
# label(child_12$any_fevers) <- "Any fevers > 38 degrees C"
# label(child_12$qn_1_23_episodes) <- "Average # fevers (among any fevers)"
# label(child_12$qn_1_25_admission) <- "Any hospital admissions"
# label(child_12$qn_1_26) <- "Any antibiotics use"
# label(child_12$abx_days_12) <- "Average # days of antibiotics (among users)"
# label(child_12$qn_1_27) <- "Any other diagnoses"
# label(child_12$n_dx) <- "Average # other diagnoses (among any diagnoses)"
# 
# table1(~ qn_1_17_diarrhoea + any_fevers + qn_1_23_episodes + qn_1_25_admission + 
#          qn_1_26 + abx_days_12 + qn_1_27 + n_dx | IVF, 
#        data = child_12, 
#        topclass = "Rtable1-zebra",
#        render.continuous = my.render.cont) 

# SUP TABLE 5 - cpg list 
cpg_merge <- read_dta(paste0(base_path,"20190125-genome_data.dta")) %>% select(-SCORE)
annotated_ewas <- read_csv(paste0(out_path, "annotate_results.csv"))
requested <- read_xlsx("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DNA Methylation/20190506-FinalCpGs.xlsx") %>% 
  filter(!is.na(cpg)) %>% rename(IlmnID = cpg) 
# requested %>% count(paper) %>% summarize(sum(n)) # searched for 172 from EWAS + 52 ANRIL/RXRA + 113 from NECAB/HIF
annotated_ewas <- left_join(annotated_ewas, requested) 
annotated_ewas %>% filter(!is.na(paper)) %>% count(IlmnID) # found 143 from EWAS
annotated_ewas %>% filter(is.na(paper) & str_detect(UCSC_RefGene_Name,"CDKN")) %>% count(IlmnID) # found 7 ANRIL 
annotated_ewas %>% filter(is.na(paper) & str_detect(UCSC_RefGene_Name,"RXRA")) %>% count(IlmnID) # found 37 ANRIL 
annotated_ewas %>% filter(is.na(paper) 
                          & !str_detect(UCSC_RefGene_Name,"CDKN") 
                          & !str_detect(UCSC_RefGene_Name,"RXRA")) %>% count(IlmnID) # found 94 NECAB/HIF 
annotated_ewas <- annotated_ewas %>% 
  mutate(paper = case_when(!is.na(paper) ~ paper,
                           str_detect(UCSC_RefGene_Name,"CDKN") ~ "based on Lillycrop K, et al. EBioMedicine. 2017.",
                           str_detect(UCSC_RefGene_Name,"RXRA") ~ "based on Godfrey KM, et al. Diabetes 2017.",
                           TRUE ~ "based on intial cg13403462 finding")) 
write_csv(annotated_ewas, paste0(out_path, "annotate_results_paper.csv"))

###############################################
# SUP FIGURES
###############################################
# SUP FIGURE 1 - EFW
scatter <- read_xlsx(paste0(data_path,"efw.xlsx")) %>% 
  ggplot(aes(x = DAYS/7, y = EFW*1000, color = IVF)) + geom_point() + geom_smooth(alpha = 0.25) +
  labs(color = "", x = "gestational age (weeks)", y = "estimated fetal weight (grams)") +
  scale_x_continuous(limits = c(19,35), breaks = seq(17,37,2)) +
  theme(legend.position = c(0.38, 0.95), legend.background = element_blank())

regress <- reg_out_new %>% filter(outcome %in% c("EFW_19", "EFW_26", "EFW_32")) %>% 
  add_column(EFW_visit = c(19, 26, 32)) %>% 
  mutate_at(c("estimate", "min95", "max95"), list(~ 1000*.)) %>% 
  ggplot() + geom_point(aes(x = EFW_visit, y = estimate)) +
  geom_errorbar(aes(x = EFW_visit, ymin = min95, ymax= max95), width = 1) + 
  geom_ribbon(aes(x = EFW_visit, ymin = min95, ymax= max95), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-200,200), breaks = seq(-200,200, 50)) +
  scale_x_continuous(limits = c(17,33), breaks = seq(17,39,2)) +
  labs(x = "gestational age (weeks)", y = "estimated fetal weight difference (ART - SC), in grams")

ggarrange(scatter, regress, labels = c("A", "B"), ncol = 2, nrow = 1)

# SUP FIGURE 2 - QMR at 5 and 6 
# (Raw data not available)
# merged %>% filter(visit == 0) %>% select(IVF, contains("_5_qmr"), contains("_6_qmr")) %>% 
#   gather(measure, val, -IVF) %>% separate(measure, c("part", "year", "tech")) %>% filter(!is.na(val)) %>% 
#   group_by(year, part, IVF) %>% summarize(sum(!is.na(val)))

# QMR_5 <- merged %>% filter(visit == 0) %>% dplyr::select(IVF, contains("_5_")) %>% 
#   mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART", "Spontaneous \nConception"))) %>%
#   gather(measure, val, -IVF, -age_5_qr) %>% separate(measure, c("part", "year", "tech")) %>% 
#   dplyr::select(-year, -tech) %>% 
#   mutate(part = factor(part, levels = c("wt","lean","fat"), labels = c("Weight (kg)", "Lean mass (kg)", "Fat mass (kg)"))) %>% 
#   ggplot(aes(x = IVF, y = val, color = IVF)) + geom_quasirandom(alpha = 0.75) + geom_boxplot(alpha = 0.25) + facet_wrap(~part) + 
#   labs(x = "Year 5 visit (ART N = 14; SC N = 233)", y = "", color = NULL) + scale_x_discrete(label = NULL) +
#   theme(legend.position = c(0.15, 0.1), legend.background = element_blank())
# 
# QMR_6 <- merged %>% filter(visit == 0) %>% dplyr::select(IVF, contains("_6_")) %>% 
#   mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART", "Spontaneous \nConception"))) %>%
#   gather(measure, val, -IVF, -age_6_qr) %>% separate(measure, c("part", "year", "tech")) %>% 
#   dplyr::select(-year, -tech) %>% 
#   mutate(part = factor(part, levels = c("wt","lean","fat"), labels = c("Weight (kg)", "Lean mass (kg)", "Fat mass (kg)"))) %>% 
#   ggplot(aes(x = IVF, y = val, color = IVF)) + geom_quasirandom(alpha = 0.75) + geom_boxplot(alpha = 0.25) + facet_wrap(~part) + 
#   labs(x = "Year 6 visit (ART N = 23; SC N = 356)", y = "", color = NULL) + scale_x_discrete(label = NULL) +
#   theme(legend.position = "none")
# 
# ggarrange(QMR_5, QMR_6, labels = c("A", "B"), nrow = 1, ncol = 2)

# SUP FIGURE 3 - 6 year chems and BP, by IVF
# (Raw data not currently available)
# child_out %>% 
#   dplyr::select(SubjectID, IVF, glucose_1, HOMA_IR, HOMA_B, SBP, DBP) %>% 
#   left_join(., lab6_clean, by = "SubjectID") %>%
#   dplyr::select(-contains("FLAG"), -twin_id) %>% gather(var, val, -SubjectID, -IVF) %>% 
#   #filter(SubjectID %in% sub_id) %>% 
#   mutate(var = case_when(var == "ALT" ~ "B. ALT (U/L)",
#                          var == "AST" ~ "C. AST (U/L)",
#                          var == "CHDL" ~ "I. total cholesterol/HDL ratio",
#                          var == "CHOL" ~ "G. total cholesterol (mmol/L)",
#                          var == "CRE" ~ "A. creatinine (umol/L)",
#                          var == "GGT" ~ "D. GGT (U/L)",
#                          var == "glucose_1" ~ "J. fasting glucose (mmol/L)",
#                          var == "HDL" ~ "E. HDL (mmol/L)",
#                          var == "LDL" ~ "F. LDL (mmol/L)",
#                          var == "HOMA_B" ~ "L. HOMA-B (%)",
#                          var == "HOMA_IR" ~ "M. HOMA-IR",
#                          var == "HSCRP" ~ "N. hs-CRP (mg/L)",
#                          var == "INS" ~ "K. insulin (mIU/L)",
#                          var == "TG" ~ "H. triglycerides (mmol/L)",
#                          var == "SBP" ~ "O. systolic BP (mmHg)",
#                          var == "DBP" ~ "P. diastolic BP (mmHg)")) %>% 
#   #group_by(var, IVF) %>% summarize(obs = sum(!is.na(val))) %>% arrange(IVF, obs) %>% print(n = Inf)
#   #group_by(var) %>% summarize(obs = sum(!is.na(val))) %>% arrange(obs) %>% print(n = Inf)
# mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART (N = 28)", "SC (N = 384)"))) %>% 
#   ggplot() + 
#   geom_density(aes(x = val, fill = IVF), alpha = 0.4) + 
#   facet_wrap(~var, scales = "free") + 
#   labs(fill = NULL, x = "", y = "") +
#   theme(legend.position = c(0.93,0.95), legend.background = element_blank())

# SUP FIGURE 4 - anthro, skinfolds, bp by C-TMLE
ctmle_res <- read_csv(paste0(data_path, "ivf_ctmle_anthro.csv"))
ANT <- ctmle_res %>% 
  mutate(OUT = case_when(OUT == "height" ~ "A.~Height~(cm)",
                         OUT == "log_weight" ~ "B.~Weight~(log(kg))",
                         OUT == "log_bmi" ~ "C.~BMI~(log(kg/m^2))",
                         OUT == "zlen" ~ "D.~'Height-for-age'~'Z-score'~(SDS)",
                         OUT == "zwei" ~ "E.~'Weight-for-age'~'Z-score'~(SDS)",
                         OUT == "zbmi" ~ "F.~BMI~'Z-score'~(SDS)")) %>% 
  ggplot(aes(x = VIS, y = ATE)) + geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(x = VIS, ymin = ATE-(1.96*SE), ymax = ATE+(1.96*SE))) +
  geom_ribbon(aes(x = VIS, ymin = ATE-(1.96*SE), ymax = ATE+(1.96*SE), fill = OUT), alpha = 0.25) +
  labs(x = "Age (Months)", y = "") + 
  facet_wrap(~OUT, scale = "free", labeller = label_parsed) + theme(legend.position = "NONE")

ctmle_skin_res <- read_csv(paste0(data_path,"ivf_ctmle_skin.csv"))
SF <- ctmle_skin_res %>% 
  mutate(OUT = case_when(OUT == "log_tri" ~ "A. Triceps (% difference)",
                         OUT == "log_sub" ~ "B. Subscapular (% difference)",
                         OUT == "log_bi" ~ "C. Biceps (% difference)",
                         OUT == "log_supra" ~ "D. Suprailiac (% difference)")) %>% 
  mutate(LB = ATE-(1.96*SE), UB = ATE+(1.96*SE)) %>% 
  mutate_at(c("ATE","LB","UB"),funs((exp(.)-1)*100)) %>% 
  ggplot(aes(x = VIS, y = ATE)) + geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(x = VIS, ymin = LB, ymax = UB)) +
  geom_ribbon(aes(x = VIS, ymin = LB, ymax = UB, fill = OUT), alpha = 0.25) +
  scale_y_continuous(limits = c(-25,15)) +
  scale_x_continuous(breaks = c(0, 18, 36, 48, 60, 72)) +
  labs(x = "Age (Months)", y = "") + theme(legend.position = "none") +
  facet_wrap(~OUT) + theme(legend.position = "NONE")

ctmle_bp_res <- read_csv(paste0(data_path, "ivf_ctmle_bp.csv"))
BP <- ctmle_bp_res %>% 
  mutate(OUT = case_when(OUT == "SBP" ~ "A. Systolic, absolute (mmHg)", 
                         OUT == "DBP" ~ "B. Diastolic, absolute (mmHg)",
                         OUT == "ZSBP" ~ "C. Systolic, standardized (%ile rank)",
                         OUT == "ZDBP" ~ "D. Diastolic, standardized (%ile rank)")) %>% 
  mutate(LB = ATE-(1.96*SE), UB = ATE+(1.96*SE)) %>% 
  ggplot(aes(x = VIS, y = ATE)) + geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(x = VIS, ymin = LB, ymax = UB), width = 1) +
  geom_ribbon(aes(x = VIS, ymin = LB, ymax = UB, fill = OUT), alpha = 0.25) +
  scale_y_continuous(limits = c(-25, 10)) +
  scale_x_continuous(breaks = seq(36,72, 12)) +
  labs(x = "Age (Months)", y = "") +
  facet_wrap(~OUT) + theme(legend.position = "NONE")

ggarrange(ANT, ggarrange(SF, BP, labels = c("B", "C"), ncol = 2, nrow = 1), 
          labels = c("A"), ncol = 1, nrow = 2)

# SUP FIGURE 5 - anthro, skinfolds, bp in SUB
reg_out_sub <- read_excel(paste0(data_path,"8 Jan 2019-IVF-size_MI_SUB.xls"), sheet = "min")
out_sub_visit <- reg_out_sub %>% 
  filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  add_column(all_visit = rep(c(0, 1, 3, 6, 9, 12, 15, 18, 24, 36, 48, 54, 60, 66, 72, 78), 6))
out_sub_visit %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  group_by(outcome) %>% summarize(min = min(N), max = max(N))

ANTHRO_SUB <- out_sub_visit %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
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
  labs(x = "Age (Months)", y = "") + theme(legend.position = "NONE") +
  facet_wrap(~outcome, scale = "free", labeller = label_parsed)

BP_SUB <- reg_out_sub %>% filter(outcome %in% c("SBP", "ZSBP", "DBP", "ZDBP")) %>% 
  add_column(BP_visit = rep(c(36, 48, 60, 72), 4)) %>% 
  mutate(outcome = case_when(outcome == "SBP" ~ "A. Systolic, absolute (mmHg)", 
                             outcome == "DBP" ~ "B. Diastolic, absolute (mmHg)",
                             outcome == "ZSBP" ~ "C. Systolic, standardized (%ile rank)",
                             outcome == "ZDBP" ~ "D. Diastolic, standardized (%ile rank)")) %>% 
  ggplot() + geom_point(aes(x = BP_visit, y = estimate)) +
  geom_errorbar(aes(x = BP_visit, ymin = min95, ymax= max95), width = 1) + 
  geom_ribbon(aes(x = BP_visit, ymin = min95, ymax= max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  scale_y_continuous(limits = c(-40, 15)) +
  scale_x_continuous(breaks = seq(36,72, 12)) +
  facet_wrap(~ outcome) + theme(legend.position = "NONE") +
  labs(x = "Age (Months)", y = "")

SF_SUB <- reg_out_sub %>% filter(outcome %in% c("log_sub", "log_tri", "log_bi", "log_supra")) %>% 
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
  scale_y_continuous(limits = c(-40,20)) +
  scale_x_continuous(breaks = c(0, 18, 36, 48, 60, 72))+
  labs(x = "Age (Months)", y = "") + theme(legend.position = "NONE") +
  facet_wrap(~outcome)

ggarrange(ANTHRO_SUB, ggarrange(SF_SUB, BP_SUB, labels = c("B", "C"), ncol = 2, nrow = 1), 
          labels = c("A"), ncol = 1, nrow = 2)

# SUP FIGURE 6 - zlen by alternate subgroupings
# (Raw data not currently available)
# model_sub <- merged %>% filter(SubjectID %in% sub_id | IVF == 1) %>% mutate(model = "A. Infertility indications (N = 93)")
# model_both <- merged %>% filter(SubjectID %in% c(sub_id, DAD_COND) | IVF == 1) %>% mutate(model = "B. Indications + paternal risk factors (N = 200)")
# model_dad <- merged %>% filter(SubjectID %in% DAD_COND | IVF == 1) %>% mutate(model = "C. Paternal risk factors only (N = 121)")
# model_rand_inv <- merged %>% filter(SubjectID %in% rand_sub | IVF == 1) %>% mutate(model = "D. No indications (20% random sample; N = 204)")
# all_models <- bind_rows(model_sub, model_both, model_dad, model_rand_inv)
# all_models %>% mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART", "Spontaneous Conception"))) %>% 
#   ggplot(aes(x = days, y = zlen, color = IVF), alpha = 0.4) +
#   geom_point() + geom_smooth() + facet_wrap(~model) +
#   labs(y = "length-for-age Z-score (SDS)", x = "Age (Days)", color = NULL) +
#   theme(legend.position = c(0.88, 0.08), legend.background = element_blank())

# SUP FIGURE 7 - anthro by IPCW
reg_out_IPCW <- read_excel(paste0(data_path,"12 Feb 2019-IVF-size_wt.xls"), sheet = "min")
out_IPCW_visit <- reg_out_IPCW %>% 
  filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  add_column(all_visit = rep(c(0, 1, 3, 6, 9, 12, 15, 18, 24, 36, 48, 54, 60, 66, 72, 78), 6))
out_IPCW_visit %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  group_by(outcome) %>% summarize(min = min(N), max = max(N))

ANTHRO_IPCW <- out_IPCW_visit %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  mutate(outcome = case_when(outcome == "height" ~ "A.~Height~(cm)",
                             outcome == "log_weight" ~ "B.~Weight~(log(kg))",
                             outcome == "log_bmi" ~ "C.~BMI~(log(kg/m^2))",
                             outcome == "zlen" ~ "D.~'Height-for-age'~'Z-score'~(SDS)",
                             outcome == "zwei" ~ "E.~'Weight-for-age'~'Z-score'~(SDS)",
                             outcome == "zbmi" ~ "F.~BMI~'Z-score'~(SDS)")) %>% 
  ggplot() + 
  geom_point(aes(x = all_visit, y = estimate)) +
  geom_errorbar(aes(x = all_visit, ymin = min95, ymax= max95)) + 
  geom_ribbon(aes(x = all_visit, ymin = min95, ymax= max95), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Age (Months)", y = "") +
  facet_wrap(~outcome, scale = "free", labeller = label_parsed)
ANTHRO_IPCW

# SUP FIG 8 - gestational weight gain, mode of delivery, and duration breastfeeding
# (Raw data not currently available)
# mweight_long <- merged %>% 
#   filter(sex != "", visit == 0, !(mother_ethnicity %in% c("", "others"))) %>% 
#   dplyr::select(SubjectID, IVF, mwt_0 = ppWeight, mother_ethnicity, parity, mht = m_height_pw26, 
#          mwt_1 = weight1, mwt_2 = weight2, mwt_3 = weight3, mwt_4 = weight4, mwt_5 = weight5, 
#          mwt_6 = weight6, mwt_7 = weight7, mwt_8 = weight8, mwt_9 = weight9,
#          mwtwk_1 = gestationalWeek1, mwtwk_2 = gestationalWeek3, mwtwk_4 = gestationalWeek4, 
#          mwtwk_5 = gestationalWeek5, mwtwk_6 = gestationalWeek6, mwtwk_7 = gestationalWeek7, 
#          mwtwk_8 = gestationalWeek8, mwtwk_9 = gestationalWeek9) %>% 
#   mutate(mwtwk_0 = 0) %>% gather(measure, val, -SubjectID, -mother_ethnicity, -parity, -mht, -IVF) %>% 
#   separate(measure, c("measure", "GA_visit")) %>% spread(measure, val) %>% 
#   mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART","Spontaneous Conception"))) %>% 
#   mutate(mbmi = mwt/(mht/100)^2)
# 
# neonate_data <- read_excel("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/GUSTO_neonatal complications.xlsx")
# comp_list <- neonate_data %>% dplyr::select(contains("Neonatal complication")) %>% names()
# severe_cond <- neonate_data %>% dplyr::select(PSCID, comp_list) %>% summarize_all(list(~sum(!is.na(.)))) %>% 
#   gather(var, N) %>% arrange(-N) %>% add_column(n = seq(1,57,1)) %>% 
#   filter(n %in% c(3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 15, 16, 21, 24, 25, 54)) %>% pull(var)
# 
# # plots
# gwg <- mweight_long %>% filter(mwtwk != 0) %>% 
#   filter(SubjectID %in% sub_id) %>% 
#   ggplot(aes(x = mwtwk, y = mbmi, color = IVF)) + #geom_point(alpha = 0.25) + 
#   geom_smooth(alpha = 0.8) +
#   geom_line(aes(x = mwtwk, y = mbmi, group = SubjectID), alpha  = 0.25) +
#   scale_x_continuous(breaks = seq(0,40,7)) +
#   labs(x = "gestational age at visit (weeks)", 
#        y = expression(paste("Maternal BMI ( ", kg/m^2, ")")), 
#        color = NULL) +
#   theme(legend.position = c(0.25, 0.92), legend.background = element_blank())
# 
# complications <- neonate_data %>% dplyr::select(PSCID, severe_cond) %>% 
#   mutate_at(severe_cond, list(~!is.na(.))) %>% mutate(conds = rowSums(dplyr::select(.,severe_cond)), 
#                                                       IVF = substring(PSCID,3,3) == "9") %>% 
#   dplyr::select(PSCID, IVF, conds) %>% filter(PSCID %in% sub_id) %>% 
#   mutate(IVF = factor(IVF, levels = c(T,F), labels = c("ART (mean count = 0.91)","SC (mean count = 0.97)"))) %>% 
#   ggplot() + geom_bar(aes(conds, fill = IVF), position = "dodge", binwidth = 1) +
#   scale_x_continuous(breaks = seq(0,7,1)) +
#   labs(fill = "", y = "", x = "# delivery, neonatal cardiorespiratory, 
#        or neonatal musculoskeletal complications") +
#   theme(legend.position = c(0.73, 0.925), legend.background = element_blank())
# 
# mod <- merged %>% filter(visit == 0, sex != "") %>% filter(SubjectID %in% sub_id) %>%
#   mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART","Spontaneous Conception"))) %>% 
#   group_by(IVF) %>% mutate(tot = n()) %>% mutate(mod = substr(mode_of_delivery,1,1)) %>% 
#   mutate(mod = factor(mod, levels = c("1","2","3","4","5"), 
#                       labels = c("spontaneous, \nvaginal", "assisted, \nvaginal", 
#                                  "induced, \nvaginal", "elective, \nC/S", "emergent, \nC/S"))) %>% 
#   ggplot() + geom_histogram(aes(mod, fill = IVF), stat = "count", position = "dodge") +
#   labs(x = "mode of delivery", y = NULL, fill = NULL) + 
#   theme(legend.position = "none")
#   #theme(legend.position = c(0.875, 0.925), legend.background = element_blank())
# 
# bf <- merged %>% filter(visit == 0, sex != "") %>% 
#   filter(SubjectID %in% sub_id, Dur_any_BF != "") %>%
#   mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART", "Spontaneous Conception"))) %>% 
#   mutate(ANY_BF = factor(Dur_any_BF, levels = c("lt_1M", "1M_to_lt_3M", "3M_to_lt_6M", "6M_to12M", "12M_and_above"),
#                          labels = c("<1 mo.", "1 to <3 mo.", "3 to <6 mo.", "6 to 12 mo.", "12+ mo."))) %>% 
#   group_by(IVF) %>% mutate(tot = n()) %>% 
#   ggplot() + 
#   geom_bar(aes(x = ANY_BF, y = stat(count), fill = IVF), position = "dodge") +
#   labs(x = "duration of any breastfeeding", y = NULL, fill = NULL) + 
#   theme(legend.position = "none")
#   #theme(legend.position = c(0.15, 0.90), legend.background = element_blank())
# 
# ggarrange(
#   ggarrange(gwg, mod, labels = c("A","B"), nrow = 1, ncol = 2),
#   ggarrange(complications, bf, labels = c("C","D"), nrow = 1, ncol = 2),
#   nrow = 2, ncol = 1)

# SUP FIG 9 - decomposition by HIF3A
decomp_anthro_plot <- function(path){
  decomp_data <- read_excel(path) %>%
  rename(TCE_b = gcomp1, TCE_se = gcomp2, NDE_b = gcomp3, NDE_se = gcomp4,
         NIE_b = gcomp5, NIE_se = gcomp6, CDE_b = gcomp7, CDE_se = gcomp8, visit = gcomp9) %>% 
  add_column(measure = rep(c("zlen", "height", "zwei", "log_weight", "log_bmi", "zbmi"),16)) %>% 
  gather(var, val, -visit, -measure) %>% separate(var, c("effect", "param")) %>% 
  spread(key = param, val = val) %>% 
  mutate(ub = b + 1.96*se, lb = b - 1.96*se) %>% 
  mutate(te_ub = if_else(effect == "TCE", ub, NA_real_), te_lb = if_else(effect == "TCE", lb, NA_real_),
         nde_ub = if_else(effect == "NDE", ub, NA_real_), nde_lb = if_else(effect == "NDE", lb, NA_real_),
         nie_ub = if_else(effect == "NIE", ub, NA_real_), nie_lb = if_else(effect == "NIE", lb, NA_real_),
         cde_ub = if_else(effect == "CDE", ub, NA_real_), cde_lb = if_else(effect == "CDE", lb, NA_real_)) %>%
  mutate(measure = case_when(measure == "height" ~ "A. Height (cm)",
                             measure == "log_weight" ~ "B. Weight (% change)",
                             measure == "log_bmi" ~ "C. BMI (% change)",
                             measure == "zlen" ~ "D. Height-for-age (SDS)",
                             measure == "zbmi" ~ "F. BMI Z-score (SDS)", 
                             measure == "zwei" ~ "E. Weight-for-age (SDS)")) 
nie_bounds <- decomp_data %>% filter(effect == "NIE")
plot_obj <- decomp_data %>% ggplot(aes(x = visit)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_point(aes(y = if_else(effect == "NDE", b, NA_real_))) +
  geom_errorbar(aes(ymin = nde_lb, ymax = nde_ub)) +
  geom_point(aes(y = if_else(effect == "NIE", b, NA_real_)), shape = 18) +
  geom_errorbar(aes(ymin = nie_lb, ymax = nie_ub), linetype = "dashed") + 
  geom_ribbon(data = nie_bounds, aes(x = visit, ymin = lb, ymax = ub), alpha = 0.15) +
  labs(x = "Age (Months)", y = "") +
  facet_wrap(~measure, scales = "free") + 
  scale_x_continuous(breaks = seq(0,72,12))
plot_obj
}

hif3a <- decomp_anthro_plot(paste0(data_path, "15 Jul 2019-gformula_HIF.xls"))
hif3a_sub <- decomp_anthro_plot(paste0(data_path, "12 Jul 2019-gformula_SUB_HIF.xls"))

ggarrange(hif3a, hif3a_sub, labels = c("A", "B"), nrow = 2, ncol = 1)

# SUP FIG 10 - lost to follow-up
# (Raw data not currently available)
# merged_status %>% group_by(IVF) %>% filter(visit == 0, STATUS != "Ineligible") %>% summarize(sum(!is.na(mother_ethnicity))) # a
# merged_status %>% group_by(IVF) %>% filter(visit == 0, STATUS != "Ineligible") %>% summarize(sum(!is.na(EFW_19))) 
# merged_status %>% group_by(IVF) %>% filter(visit == 0, STATUS != "Ineligible") %>% summarize(sum(!is.na(EFW_26))) 
# merged_status %>% group_by(IVF) %>% filter(visit == 0, STATUS != "Ineligible") %>% summarize(sum(!is.na(ogtt_fasting_pw26))) 
# merged_status %>% group_by(IVF) %>% filter(visit == 0, STATUS != "Ineligible") %>% summarize(sum(!is.na(EFW_32))) 
# merged_status %>% group_by(IVF) %>% filter(visit == 0, STATUS != "Ineligible") %>% summarize(sum(!is.na(weight))) 
# merged_status %>% group_by(IVF, visit) %>% 
#   summarize(n = sum(!is.na(weight) & STATUS != "Ineligible")) %>% print(n = Inf)
# 
# merged_status <- read_xlsx(paste0(base_path, "FormA428&429_20181001.xlsx"), na = c("N/A", "NA"," ", "")) %>% 
#   dplyr::select(SubjectID, STATUS = participant_status) %>% left_join(merged, .) %>% 
#   mutate(STATUS = case_when(STATUS == "active" ~ "Active",
#                             STATUS == "death" ~ "Death",
#                             STATUS == "dropout" ~ "Dropout",
#                             STATUS == "ineligible" ~ "Ineligible")) %>% 
#   mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART (N = 85)", "Spontaneous Conception (N = 1152)")))

# merged_status %>% group_by(IVF, visit) %>%
#   summarize(active = sum(STATUS == "Active" & !is.na(height)),
#             dropout = sum(STATUS == "Dropout" & !is.na(height)),
#             death = sum(STATUS == "Death" & !is.na(height)),
#             height = sum(!is.na(height)), n = sum(STATUS != "Ineligible"),
#             pct_no_obs = round(100*(n - height)/n, 0)) %>% print(n = Inf)

# ltfu <- merged_status %>% filter(STATUS != "Ineligible") %>% 
#   group_by(IVF, visit) %>% mutate(n = sum(STATUS != "Active" & !is.na(height)),
#                                   obs_ht = sum(!is.na(height)), 
#                                   tot = sum(STATUS != "Ineligible"),
#                                   pct_no_obs = round(100*(tot - obs_ht)/tot, 0)) %>% 
#   ungroup() %>% 
#   ggplot() + geom_boxplot(aes(x = as.factor(visit), y = height, color = STATUS)) + 
#   geom_text(aes(as.factor(visit), 25, label = n)) +
#   geom_text(aes(0.5, 30, label = "N eventual losses observed at visit = "), hjust = "left") +
#   # geom_text(aes(as.factor(visit), 25, label = pct_no_obs)) +
#   # geom_text(aes(0.5, 30, label = "Not unobserved at visit (%): "), hjust = "left") + 
#   theme(legend.position = c(0.1,0.85), legend.background = element_blank()) +
#   scale_y_continuous(breaks = seq(0,150,25)) +
#   facet_wrap(~IVF) +
#   labs(x = "visit month", y = "height (cm)", color = "Last known status:")
# 
# visits <- merged_status %>% count(visit) %>% pull(visit)
# obs <- merged_status %>% 
#   group_by(IVF, visit) %>% 
#   summarize(perc = 100*(sum(!is.na(height))/sum(STATUS != "Ineligible"))) %>% 
#   ggplot(aes(x = visit, y = perc, group = IVF)) + 
#   geom_line(aes(color = IVF, linetype = IVF)) + 
#   scale_x_continuous(breaks = visits) +
#   scale_y_continuous(breaks = seq(40,100,5)) +
#   labs(x = "visit month", y = "% completing visit", color = NULL, linetype = NULL) +
#   theme(legend.position = c(0.825,0.925), legend.background = element_blank())
# 
# ggarrange(obs, ltfu, labels = c("A", "B"), nrow = 2, ncol = 1)

# SUP FIG 11 - margins plots for mixed models
zlen <- image_read(paste0(data_path, "mixed_zlen_margins.png"))     %>% as.ggplot()
height <- image_read(paste0(data_path, "mixed_height_margins.png")) %>% as.ggplot()
zwei <- image_read(paste0(data_path, "mixed_zwei_margins.png"))     %>% as.ggplot()
weight <- image_read(paste0(data_path, "mixed_weight_margins.png")) %>% as.ggplot()
zbmi <- image_read(paste0(data_path, "mixed_zbmi_margins.png"))     %>% as.ggplot()
bmi <- image_read(paste0(data_path, "mixed_bmi_margins.png"))       %>% as.ggplot()
ggarrange(height, zlen,
          weight, zwei, 
          bmi, zbmi, labels = c("A","B","C","D","E","F"), 
          nrow = 3, ncol = 2)
