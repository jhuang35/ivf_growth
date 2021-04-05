#######################################################################################
#
# Analyses of child cardiometabolic phenotype following assisted reproductive technologies using a pragmatic trial emulation approach
#
# CORRESPONDING AUTHOR: Jon Huang
# CONTACT INFORMATION: jonathan_huang@sics.a-star.edu.sg
#
# VERSION DATE: 2021 April 5
# 
# Available at: 
# github.com/jhuang35/ivf_growth/ 
# www.doi.org/10.5281/zenodo.4659508
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
# TABLE 1 - Overall study population parental characteristics, by conception status.
# (Raw data for Table 1 available upon request)
# table1(~ mom_age + nullip + ethnic + mom_edu + hh_inc + 
#          m_height_pw26 + ppBMI + ogtt_fasting_pw26 + ogtt_2hour_pw26 + htn +
#          smk_home + father_age_delivery + f_height_m24 + f_weight_m24 + dad_dm + dad_htn | IVF, 
#        data = formatted, overall  = "Overall",
#        topclass = "Rtable1-zebra",
#        render.continuous = my.render.cont) 

# TABLE 2 - Child anthropometrics, skinfolds, and blood pressure at 6 years, adjusted for pre-pregnancy and pregnancy factors.
# (Regression output file included - 22 Mar 2019-IVF-size_MI.xlsx)

###############################################
# FIGURES
###############################################
# Figure 1 - Associations between ART status and anthropometry, adjusted for pre-pregnancy characteristics.
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
  labs(x = "Age (Months)", y = "") + 
  facet_wrap(~outcome, scale = "free", labeller = label_parsed) +
  theme_bw() + theme(legend.position = "none")

# Figure 2 - Associations between ART status and skinfold thickness (a) and blood pressure (b), adjusted for pre-pregnancy characteristics.
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
  facet_wrap(~ outcome) + 
  labs(x = "Age (Months)", y = "Absolute (mmHg) or standardized (%ile rank) 
       difference, (ART - SC)") + 
  theme_bw() + theme(legend.position = "none")

SF <- reg_out_new %>% filter(outcome %in% c("log_sub", "log_tri", "log_bi", "log_supra")) %>% 
  add_column(skin_visit = c(0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            48, 54, 60, 66, 72, 78)) %>% 
  mutate_at(c("estimate","min95","max95"),funs((exp(.)-1)*100)) %>% 
  mutate(outcome = case_when(outcome == "log_tri" ~ "A. Triceps",
                             outcome == "log_sub" ~ "B. Subscapular",
                             outcome == "log_bi" ~ "C. Biceps",
                             outcome == "log_supra" ~ "D. Suprailiac")) %>% 
  ggplot() + 
  geom_point(aes(x = skin_visit, y = estimate)) +
  geom_errorbar(aes(x = skin_visit, ymin = min95, ymax = max95)) + 
  geom_ribbon(aes(x = skin_visit, ymin = min95, ymax = max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-30,10)) +
  scale_x_continuous(breaks = c(0, 18, 24, 36, 48, 54, 60, 66, 72, 78))+
  labs(x = "", y = "% difference in skinfold thickness
       (ART - SC)") + theme_bw() + theme(legend.position = "none") +
  facet_wrap(~outcome)

ggarrange(SF, BP, labels = c("a", "b"), ncol = 1, nrow = 2)

# FIGURE 3 - Difference in anthropometrics (a), skinfold thickness (b), and blood pressures (c) comparing ART-conceived foetuses versus a putatively subfertile cohort, adjusted for pre-pregnancy characteristics.
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
  labs(x = "Age (Months)", y = "") + theme_bw() + theme(legend.position = "NONE") +
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
  facet_wrap(~ outcome) + theme_bw() + theme(legend.position = "NONE") +
  labs(x = "Age (Months)", y = "")

SF_SUB <- reg_out_sub %>% filter(outcome %in% c("log_sub", "log_tri", "log_bi", "log_supra")) %>% 
  add_column(skin_visit = c(0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            0, 18, 24, 36, 48, 54, 60, 66, 72, 78, 
                            48, 54, 60, 66, 72, 78)) %>% 
  mutate_at(c("estimate","min95","max95"),funs((exp(.)-1)*100)) %>% 
  mutate(outcome = case_when(outcome == "log_tri" ~ "A. Triceps",
                             outcome == "log_sub" ~ "B. Subscapular",
                             outcome == "log_bi" ~ "C. Biceps",
                             outcome == "log_supra" ~ "D. Suprailiac")) %>% 
  ggplot() + 
  geom_point(aes(x = skin_visit, y = estimate)) +
  geom_errorbar(aes(x = skin_visit, ymin = min95, ymax= max95)) + 
  geom_ribbon(aes(x = skin_visit, ymin = min95, ymax= max95, fill = outcome), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-40,20)) +
  scale_x_continuous(breaks = c(0, 18, 36, 48, 60, 72))+
  labs(x = "Age (Months)", y = "") + theme_bw() + theme(legend.position = "NONE") +
  facet_wrap(~outcome)

ggarrange(ANTHRO_SUB, ggarrange(SF_SUB, BP_SUB, labels = c("b", "c"), ncol = 2, nrow = 1), 
          labels = c("a"), ncol = 1, nrow = 2)

# Figure 4 - Associations between ART status and fetal cord tissue DNA methylation and 281 candidate CpGs (a), three top NECAB3 CpGs (b), three top NECAB3 CpGs stratified by ethnicity (c).
# (a) MANHATTAN PLOT
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
  scale_color_discrete() +
  scale_x_continuous(breaks = c(1:22)) + theme_bw() + 
  theme(legend.position = "none")

# (b) NECAB3 CpGs by ART status
p1 <- read_xlsx(paste0(data_path,"cpg_ethnicity.xlsx")) %>% 
  mutate(IVF = case_when(IVF == "IVF" ~ "ART", T ~ IVF)) %>% 
  ggplot(aes(x = IVF, y = val)) + 
  geom_quasirandom(alpha = 0.75, aes(fill = IVF, color = IVF)) +
  #geom_boxplot(alpha = 0.5, aes(fill = IVF, color = IVF)) +
  #geom_signif(comparisons = list(c("Spontaneous", "ART")), test = "t.test", map_signif_level = T) +
  geom_signif(comparisons = list(c("Spontaneous", "ART")), test = "t.test", map_signif_level = F, vjust = 2) +
  scale_color_viridis_d(direction = -1, option = "cividis") +
  scale_fill_viridis_d(direction = -1, option = "cividis") +
  facet_grid( ~ IlmnID) +
  labs(group = "", fill = "", color = "", x = "", y = "") + theme_bw() + 
  theme(legend.position = c(0.85, 0.1), legend.background = element_blank(), axis.text = element_blank(), strip.text = element_text(face = "italic"))

# (c) NECAB3 CpGs by ART status and self-reported ethnicity
p2 <- read_xlsx(paste0(data_path,"cpg_ethnicity.xlsx")) %>% 
  mutate(IVF = case_when(IVF == "IVF" ~ "ART", T ~ IVF)) %>% 
  ggplot(aes(x = IVF, y = val)) + 
  geom_quasirandom(alpha = 0.75, aes(fill = IVF, color = IVF)) +
  #geom_boxplot(alpha = 0.5, aes(fill = IVF, color = IVF)) +
  #geom_signif(comparisons = list(c("Spontaneous", "ART")), test = "t.test", map_signif_level = T) +
  geom_signif(comparisons = list(c("Spontaneous", "ART")), test = "t.test", map_signif_level = F, vjust = 1.5) +
  scale_color_viridis_d(direction = -1, option = "cividis") +
  scale_fill_viridis_d(direction = -1, option = "cividis") +
  facet_grid(mother_ethnicity ~ IlmnID) +
  labs(group = "ART Status", fill = "ART Status", color = "ART Status", x = "", y = "") + theme_bw() + 
  theme(legend.position = "none", axis.text = element_blank(), strip.text = element_text(face = "italic"))

ggarrange(mp,
          ggarrange(p1, p2, labels = c("b", "c"), ncol = 2, nrow = 1), 
          labels = c("a"), ncol = 1, nrow = 2)


# FIGURES 5 & 6 - Combine mediation by NECAB3 and HIF3A into Figures 5 (all) and 6 (target trial)
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
    geom_ribbon(data = nie_bounds, aes(x = visit, ymin = lb, ymax = ub, fill = measure), alpha = 0.15) +
    labs(x = "Age (Months)", y = "") +
    facet_wrap(~measure, scales = "free") + 
    scale_x_continuous(breaks = seq(0,72,12)) +
    theme_bw() + 
    theme(legend.position = "none")
  plot_obj
}

necab <- decomp_anthro_plot(paste0(data_path,"12 Jul 2019-gformula_NEC.xls"))
necab_sub <- decomp_anthro_plot(paste0(data_path,"15 Jul 2019-gformula_SUB_NEC_ALL.xls"))
hif3a <- decomp_anthro_plot(paste0(data_path,"15 Jul 2019-gformula_HIF.xls"))
hif3a_sub <- decomp_anthro_plot(paste0(data_path,"12 Jul 2019-gformula_SUB_HIF.xls"))

# FIG 5
ggarrange(necab, hif3a, labels = c("a", "b"), nrow = 2, ncol = 1)

# FIG 6
ggarrange(necab_sub, hif3a_sub, labels = c("a", "b"), nrow = 2, ncol = 1)



###############################################
# SUPPLEMENTAL MATERIALS
###############################################
#
# SUPPLEMENTAL FILE 1 - Candidate CPG list 
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

# SUPPLEMENTAL TABLE 1 - Target pragmatic trial amongst possible subfertile couples, by conception status.
# (Raw data for Supplemental Table 1 available upon request)
# table1(~ mom_age + nullip + ethnic + mom_edu + hh_inc + 
#          m_height_pw26 + ppBMI + ogtt_fasting_pw26 + ogtt_2hour_pw26 + htn +
#          smk_home + father_age_delivery + f_height_m24 + f_weight_m24 + dad_dm + dad_htn | IVF, 
#        data = formatted_sub, overall  = "Overall",
#        topclass = "Rtable1-zebra",
#        render.continuous = my.render.cont) 
# 
# SUPPLEMENTAL TABLE 2 - Child characteristics and anthropometrics, by conception status subgroups.
# (Raw data for Supplemental Table 2 available upon request)
# table1(~ sex + GA + mod + weight_0 + height_0 + Full_BF +
#          weight_12 + weight_24 + weight_60 + weight_72 + 
#          height_12 + height_24 + height_36 + height_48 + height_60 + height_72 + 
#          glucose_1 + SBP_72 + DBP_72 | IVF, 
#        data = child_sub, overall  = NULL,
#        droplevels = T,
#        topclass = "Rtable1-zebra",
#        render.continuous = my.render.cont)
#
# SUPPLEMENTAL TABLE 3 - Child fat and lean mass at 5 and 6 years and serum cardiometabolic biomarkers at 6 years, adjusted for pre-pregnancy and pregnancy factors.
# (Regression output file included -  11 Jul 2019-IVF-chem_6_MI.xls)
#
# SUPPLEMENTAL TABLE 4 - Child illnesses reported at 12-month visit, in the past three months.
# (Raw data for Supplemental Table 4 available upon request)
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


###############################################
# SUPPLEMENTAL FIGURES
###############################################
# SUPPLEMENTAL FIGURE 1 - Difference in ultrasound-estimated weight of ART-conceived foetuses versus spontaneously-conceived foetuses, 
# compared to the subfertile cohort (A) and adjusted for pre-pregnancy characteristics (B).
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

# SUPPLEMENTAL FIGURE 2 - QMR-measured body composition at years 5 (A) and 6 (B) and distributions of serum cardiometabolic biomarkers and blood pressure at 6 years (C), by ART status.
# (Raw data for Supplemental Figure 2 available upon request)
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


# SUPPLEMENTAL FIGURE 3 - Difference in anthropometrics (A), skinfold thickness (B), and blood pressures (C) 
# comparing ART- and spontaneously-conceived foetuses, estimated by collaborative-targeted maximum likelihood estimation (C-TMLE).
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

# SUPPLEMENTAL FIGURE 4 - Length-for-age Z-score and ART status over time, by comparison cohort.
# (Raw data for Supplemental Figure 4 available upon request)
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

# SUPPLEMENTAL FIGURE 5 - Associations between ART status and anthropometry, accounting for differential follow-up at each study visit.
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

# SUPPLEMENTAL FIGURE 6 - EWAS of ART status.

# Illumina 450k annotation file not provided, but needed to recreate EWAS figure
# sites <- read_csv("HumanMethylation450_15017482_v1-2.csv", skip = 7)

# # Import core exposure and covariate data 
# dta <- readxl::read_xlsx(paste0(data_path, "FormA428&429_20181001.xlsx"))
# new_dta <- haven::read_dta(paste0(data_path,"20190723-merge_CpG_dad.dta")) %>%  
#   mutate(name = paste0("B",as.numeric(substr(SubjectID,5,10)))) %>% filter(visit == 0) %>% 
#   select(name, IVF, mother_age_delivery, mother_highest_education, mother_ethnicity, sex, 
#          ppBMI, parity, m_height_pw26, smk_home, f_height_m24, f_weight_m24, SCORE,
#          father_age_delivery, father_HBP, father_diabetes)
# 
# # Import methylome data (NCBI GEO accession ID: GSE158064)
# cord_tissue_dir <- "Cord-Tissue-Residuals/"
# cord_tissue_file <- "-v9Residuals_GUSTO_IC_Celltype_Adjusted_Residuals-1019samples.txt"
# chr <- paste0("chr",1:22,cord_tissue_file)
# 
# cpg_res_mat3 <- NULL
# for(c in 1:length(chr)){
#   print("_____________________________________________________")
#   print(paste0("IMPORTING CPGs FROM CHROMOSOME: ", c))
# cpg_dta <- data.table::fread(paste0(data_path, cord_tissue_dir,chr[c])) %>% as_tibble() %>% select(-N) %>% 
#   pivot_longer(cols = starts_with("B")) %>% pivot_wider(names_from = CpG, values_from = value)
# cpg_nums <- cpg_dta %>% select(-name) %>% names()
# merge <- left_join(new_dta, cpg_dta, by = "name")
# 
# # run regression
# for(i in 1:length(cpg_nums)){
#   if(i == 1 | i%%500 == 0){
#     print(paste0("ESTIMATING ASSOCIATIONS WITH CPGs FROM CHROMOSOME ", c, " (", round(100*(i/length(cpg_nums)),3), "%)"))
#   }
#   model <- lm(data = merge, 
#   formula = eval(as.name(cpg_nums[i])) ~ IVF + as.factor(mother_ethnicity) + mother_age_delivery + ppBMI + sex) 
#   BETA <- summary(model)$coefficients[2,1]
#   SE <- summary(model)$coefficients[2,2]
#   P <- summary(model)$coefficients[2,4]
#   S <- -log10(P)
#   CPG <- cpg_nums[i]
#   cpg_res_mat3 <- rbind(cpg_res_mat3, cbind(CPG, BETA, SE, P, S))
# }
# }
# writexl::write_xlsx(as_tibble(cpg_res_mat3), paste0(data_path,"ewas_cord_tissue_final_2.xlsx"))

cpg_res_mat3 <- readxl::read_xlsx(paste0(data_path,"ewas_cord_tissue_final_2.xlsx"))

annotate <- as_tibble(rbind(cpg_res_mat3)) %>% 
  mutate(IlmnID = CPG) %>% 
  mutate_at(2:5, list(~as.numeric(.))) %>% 
  left_join(., sites, by = "IlmnID") %>% 
  mutate(CHR = as.double(CHR)) %>% filter(!is.na(CHR)) %>% 
  mutate(Coordinate_36 = case_when(is.na(Coordinate_36) ~ 0, T ~ Coordinate_36)) %>% 
  mutate(loc = as.double(paste0(Chromosome_36,".",Coordinate_36))) %>% 
  arrange(loc) %>% mutate(rank = row_number()) 

BREAK <- annotate %>% filter(Chromosome_36 != "MULTI") %>% group_by(as.double(Chromosome_36)) %>% summarize(cent = median(rank)) %>% pull(cent)

LABELS <- annotate %>% 
  filter(S > -log10(3.6e-8)) %>% 
  select(S, IlmnID, UCSC_RefGene_Name, loc, rank) %>% arrange(-S) %>% 
  mutate(name = substring(UCSC_RefGene_Name,1,6)) 

NECAB <- annotate %>% 
  filter(substring(UCSC_RefGene_Name,1,6) == "NECAB3") %>% 
  select(S, IlmnID, UCSC_RefGene_Name, loc, rank)

hilite <- annotate %>% filter(substring(UCSC_RefGene_Name,1,6) == "NECAB3") %>% pull(IlmnID)

annotate %>% 
  filter(Chromosome_36 != "MULTI") %>% 
  ggplot() + geom_point(aes(x = rank, y = S, color = Chromosome_36)) +
  geom_point(data = NECAB, aes(x = rank, y = S), color = "BLACK") +
  geom_hline(aes(yintercept = -log10(3.6e-8)), linetype = "dashed") +
  geom_hline(aes(yintercept = -log10(0.05/336684)), linetype = "solid") +
  geom_text_repel(data = LABELS, aes(x = rank, y = S, label = paste0(IlmnID, " (",name,")")), point.padding = 1) +
  scale_x_continuous(breaks = BREAK, labels = seq(1,22,1)) +
  scale_y_continuous(breaks = seq(0,10,1)) + 
  labs(x = "Chromosome", y = "-log(P)") +
  theme_minimal() + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# qqman::manhattan(annotate, chr = "CHR", bp = "Coordinate_36", snp = "IlmnID", p = "P", 
#           annotatePval = 3.6e-8, annotateTop = F, highlight = hilite,
#           suggestiveline = -log10(0.05/336684))
# 
# qqman::manhattan(subset(annotate, CHR == 20), chr = "CHR", bp = "Coordinate_36", snp = "IlmnID", p = "P", 
#                  annotatePval = 3.6e-8, annotateTop = F, highlight = hilite, xlim = c(3.168e7, 3.175e7),
#                  suggestiveline = -log10(0.05/336684))
# 
# qqman::qq(annotate$P)


# SUPPLMENTAL FIGURE 7 - Negative control analyses of mediation by cg03904042 (A) and cg27146050 (B) methylation in maternal mid-pregnancy peripheral blood, target trial subcohort.
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
    geom_ribbon(data = nie_bounds, aes(x = visit, ymin = lb, ymax = ub, fill = measure), alpha = 0.15) +
    labs(x = "Age (Months)", y = "") +
    facet_wrap(~measure, scales = "free") + 
    scale_x_continuous(breaks = seq(0,72,12)) +
    theme(legend.position = "none")
  plot_obj
}

data_path <- "C:/Users/JHUANGYH.ARES/My Tresors/SICS/Projects/IVF/Manuscript/SUBMISSIONS/Nature Communications/Revision/code/data/"

mom_necab_sub <- decomp_anthro_plot(paste0(data_path,"26 Aug 2020-gformula_MOM_FIRST_SUB_MOM_NEC.xls"))
mom_hif3a_sub <- decomp_anthro_plot(paste0(data_path,"26 Aug 2020-gformula_MOM_FIRST_SUB_MOM_HIF.xls"))

ggpubr::ggarrange(mom_necab_sub, mom_hif3a_sub, labels = c("A", "B"), ncol = 1)


# SUPPLEMENTAL FIGURE 8. Subject retention (A) and distribution of observed heights (B), by ART and eventual loss-to-follow-up status.
# (Raw data for Supplemental Figure 8 available upon request)
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

# SUPPLMENTAL FIGURE 9 - Predicted average anthropometrics by ART status, under a linear mixed effects model with individual-level random intercepts and slopes.
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

