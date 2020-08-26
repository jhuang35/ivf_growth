##########################
# 
#
# New code / analyses for response to reviews
#
# VERSION DATE: 26 August 2020
#
# 1. cord tissue EWAS for ART (SUPPLEMENTAL FIGURE 6)
# 2. plot mediation by maternal CpGs (SUPPLEMENTAL FIGURE 7)
# 3. Additional figure corrections/revisions
#
##########################

# necessary dependencies
library(tidyverse)
library(haven)
library(ggrepel)
library(readxl)
library(ggpubr)
library(ggbeeswarm)

# SET PATH FOR DATA FILES (ALL INPUT FILES WILL BE RELATIVE TO THIS PATH)
data_path <- "~ROOT/data/" # replace "~ROOT/" with the correct path to the /data/ subdirectory

##########################
# 1. cord tissue EWAS for ART (SUPPLEMENTAL FIGURE 6)
##########################

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
# # Import methylome data (NCBI GEO accession # TBD)
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

##########################
# 2. plot mediation by maternal CpGs (SUPPLEMENTAL FIGURE 7)
##########################
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

##########################
# 3. Additional figure corrections/revisions
##########################

# Combine Supplmentary Figures
# Input baseline data and format variables
library(ggbeeswarm)
library(ggpubr)

base_path <- "C:/Users/JHUANGYH.ARES/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"
merged <- read_dta(paste0(base_path,"20181219-IVF-final.dta"))
dad_hx <- read_dta("C:/Users/JHUANGYH.ARES/Dropbox/00 - Singapore/GUSTO/99 - SUPPORT/Data QC/Paternal SES/dad_hx.dta")
merged <- left_join(merged, dad_hx)
lab6_clean <- read_csv("C:/Users/JHUANGYH.ARES/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/child 6 yr/lab6_clean.csv")
child_out <- merged %>% filter(visit == 72) %>% 
  left_join(., lab6_clean, by = "SubjectID") %>% 
  mutate(HOMA_IR = (INS * glucose_1)/22.5,
         HOMA_B = (20 * INS)/(glucose_1 - 3.5)) %>% 
  mutate(HOMA_B = if_else(glucose_1 <= 3.5, NA_real_, HOMA_B))

# SUP FIGURE 2 - QMR at 5 and 6 
# merged %>% filter(visit == 0) %>% select(IVF, contains("_5_qmr"), contains("_6_qmr")) %>% 
#   gather(measure, val, -IVF) %>% separate(measure, c("part", "year", "tech")) %>% filter(!is.na(val)) %>% 
#   group_by(year, part, IVF) %>% summarize(sum(!is.na(val)))

QMR_5 <- merged %>% filter(visit == 0) %>% dplyr::select(IVF, contains("_5_")) %>% 
  mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART", "Spontaneous \nConception"))) %>%
  gather(measure, val, -IVF, -age_5_qr) %>% separate(measure, c("part", "year", "tech")) %>% 
  dplyr::select(-year, -tech) %>% 
  mutate(part = factor(part, levels = c("wt","lean","fat"), labels = c("Weight (kg)", "Lean mass (kg)", "Fat mass (kg)"))) %>% 
  ggplot(aes(x = IVF, y = val, color = IVF)) + geom_quasirandom(alpha = 0.75) + geom_boxplot(alpha = 0.25) + facet_wrap(~part) + 
  labs(x = "Year 5 visit (ART N = 14; SC N = 233)", y = "", color = NULL) + scale_x_discrete(label = NULL) +
  theme(legend.position = c(0.15, 0.1), legend.background = element_blank())

QMR_6 <- merged %>% filter(visit == 0) %>% dplyr::select(IVF, contains("_6_")) %>% 
  mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART", "Spontaneous \nConception"))) %>%
  gather(measure, val, -IVF, -age_6_qr) %>% separate(measure, c("part", "year", "tech")) %>% 
  dplyr::select(-year, -tech) %>% 
  mutate(part = factor(part, levels = c("wt","lean","fat"), labels = c("Weight (kg)", "Lean mass (kg)", "Fat mass (kg)"))) %>% 
  ggplot(aes(x = IVF, y = val, color = IVF)) + geom_quasirandom(alpha = 0.75) + geom_boxplot(alpha = 0.25) + facet_wrap(~part) + 
  labs(x = "Year 6 visit (ART N = 23; SC N = 356)", y = "", color = NULL) + scale_x_discrete(label = NULL) +
  theme(legend.position = "none")

ggarrange(QMR_5, QMR_6, labels = c("A", "B"), nrow = 1, ncol = 2)

# SUP FIGURE 3 - 6 year chems and BP, by IVF
OUT_6 <- child_out %>% 
  dplyr::select(SubjectID, IVF, glucose_1, HOMA_IR, HOMA_B, SBP, DBP) %>% 
  left_join(., lab6_clean, by = "SubjectID") %>%
  dplyr::select(-contains("FLAG"), -twin_id) %>% gather(var, val, -SubjectID, -IVF) %>% 
  #filter(SubjectID %in% sub_id) %>% 
  mutate(var = case_when(var == "ALT" ~ "B. ALT (U/L)",
                         var == "AST" ~ "C. AST (U/L)",
                         var == "CHDL" ~ "I. total cholesterol/HDL ratio",
                         var == "CHOL" ~ "G. total cholesterol (mmol/L)",
                         var == "CRE" ~ "A. creatinine (umol/L)",
                         var == "GGT" ~ "D. GGT (U/L)",
                         var == "glucose_1" ~ "J. fasting glucose (mmol/L)",
                         var == "HDL" ~ "E. HDL (mmol/L)",
                         var == "LDL" ~ "F. LDL (mmol/L)",
                         var == "HOMA_B" ~ "L. HOMA-B (%)",
                         var == "HOMA_IR" ~ "M. HOMA-IR",
                         var == "HSCRP" ~ "N. hs-CRP (mg/L)",
                         var == "INS" ~ "K. insulin (mIU/L)",
                         var == "TG" ~ "H. triglycerides (mmol/L)",
                         var == "SBP" ~ "O. systolic BP (mmHg)",
                         var == "DBP" ~ "P. diastolic BP (mmHg)")) %>% 
  #group_by(var, IVF) %>% summarize(obs = sum(!is.na(val))) %>% arrange(IVF, obs) %>% print(n = Inf)
  #group_by(var) %>% summarize(obs = sum(!is.na(val))) %>% arrange(obs) %>% print(n = Inf)
  mutate(IVF = factor(IVF, levels = c(1,0), labels = c("ART (N = 28)", "SC (N = 384)"))) %>% 
  ggplot() + 
  geom_density(aes(x = val, fill = IVF), alpha = 0.4) + 
  facet_wrap(~var, scales = "free") + 
  labs(fill = NULL, x = "", y = "") +
  theme(legend.position = c(0.93,0.95), legend.background = element_blank())

ggarrange(ggarrange(QMR_5, QMR_6, labels = c("A", "B"), nrow = 1, ncol = 2),
          OUT_6, labels = c("", "C"), nrow = 2, ncol = 1)


# FIGURES 5 & 6 - Combine mediation by NECAB3 and HIF3A into Figures 5 (all) and 6 (target trial)
necab <- decomp_anthro_plot("C:/Users/JHUANGYH.ARES/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/12 Jul 2019-gformula_NEC.xls")
necab_sub <- decomp_anthro_plot("C:/Users/JHUANGYH.ARES/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/15 Jul 2019-gformula_SUB_NEC_ALL.xls")
hif3a <- decomp_anthro_plot("C:/Users/JHUANGYH.ARES/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/15 Jul 2019-gformula_HIF.xls")
hif3a_sub <- decomp_anthro_plot("C:/Users/JHUANGYH.ARES/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/12 Jul 2019-gformula_SUB_HIF.xls")

# FIG 5
ggarrange(necab, hif3a, labels = c("A", "B"), nrow = 2, ncol = 1)
# FIG 6
ggarrange(necab_sub, hif3a_sub, labels = c("A", "B"), nrow = 2, ncol = 1)


