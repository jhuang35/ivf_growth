library(tidyverse)
library(ggplot2)
library(haven)
library(summarytools)
library(readxl)
library(broom)
library(ggrepel)

options(tibble.print_max = 40, tibble.print_min = 30)

# 1. Set paths
out_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/"
data_path <- "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/ANALYTIC/"

# 2. Gather CpG data
PRS <- read.table("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/CpGs/EWAS-BW-987infants-scaled-allethnic-PRS-for-Jonathan.txt",
                  sep="\t", header = TRUE)
PRS <- PRS %>% rename(SubjectID = FID) %>% select(-IID)
CpG <- read.table("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/CpGs/v9Residuals_GUSTO_IC_Celltype_Adjusted_Residuals-1019samples-Jonathan.txt",
                  sep="\t", header = TRUE)
CpG_tbl <- as.tibble(CpG)
CpG_Key <- read_csv("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/CpGs/GUSTO-ICT-covariates-and-celltypePCs-1019samples-Jonathan.csv")
CpG_Key <- CpG_Key %>% select(IID, PSCID)
CpG_clean <- CpG_tbl %>% gather(SubjectID, resid, 3:1021) %>% select(-N) %>% spread(CpG, resid)
CpG_merge <- CpG_clean %>% rename(IID = SubjectID) %>% left_join(CpG_Key) %>% select(-IID) %>% rename(SubjectID = PSCID)
CpG_merge <- left_join(CpG_merge, PRS)

CpG2 <- read.table("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/CpGs/v9Residuals_GUSTO_IC_Celltype_Adjusted_Residuals-1019samples-Jonathan-16012019.txt",
                  sep="\t", header = TRUE)
CpG2_clean <- as_tibble(CpG2)
CpG2_clean <- CpG2_clean %>% gather(SubjectID, resid, 3:1021) %>% select(-N) %>% spread(CpG, resid)
CpG2_clean <- CpG2_clean %>% rename(IID = SubjectID) %>% left_join(CpG_Key) %>% select(-IID) %>% rename(SubjectID = PSCID)

pt2 <- read.table("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/DATA/CpGs/v9Residuals_GUSTO_IC_Celltype_Adjusted_Residuals-part2-1019samples.txt",
                   sep="\t", header = TRUE)
pt2_clean <- as_tibble(pt2)
pt2_clean <- pt2_clean %>% gather(SubjectID, resid, 3:1021) %>% select(-N) %>% spread(CpG, resid)
pt2_clean <- pt2_clean %>% rename(IID = SubjectID) %>% left_join(CpG_Key) %>% select(-IID) %>% rename(SubjectID = PSCID)

CpG_merge <- left_join(CpG_merge, left_join(pt2_clean, CpG2_clean))
#CpG_nums <- CPG %>% select(contains("cg")) %>% gather(name, val) %>% group_by(name) %>% summarize(mean(val, na.rm = T)) %>% pull(name)
#writeClipboard(CpG_nums)

write_dta(CpG_merge, paste0(data_path,"20190125-genome_data.dta"))
#write_dta(CpG_merge, paste0(data_path,"20181219-genome_data.dta"))

# 3. Merge with main data set
merged <- read_dta(paste0(data_path,"20181219-IVF-final.dta"))
CPG <- read_dta(paste0(data_path,"20190125-genome_data.dta"))
#CPG <- read_dta(paste0(data_path,"20181219-genome_data.dta"))
merged3 <- merged %>% left_join(., CPG)
write_dta(merged3, paste0(data_path,"20190125-merge_CpG.dta"))
#write_dta(merged3, paste0(data_path,"20181219-merge_CpG.dta"))

# 4. Re-load merged dataset as necessary
merged3 <- read_dta(paste0(data_path,"20190125-merge_CpG.dta"))
#merged3 <- read_dta(paste0(data_path,"20181219-merge_CpG.dta"))

# 5. Plot CPG distributions by IVF status
# all cpgs
by_cpg <- merged3 %>% 
  filter(visit == 0) %>% select(SubjectID, IVF, contains("cg")) %>% 
  #filter(SubjectID %in% sub_id) %>% 
  gather(site, val, -SubjectID, -IVF)

cpg_plot <- by_cpg %>% #filter(substring(site,1,4) == "cg00" | site == "cg13403462") %>% 
  ggplot(aes(x = val, y = stat(scaled), group = as.factor(IVF), fill = as.factor(IVF), color = as.factor(IVF))) +
  geom_area(alpha = 0.6, stat = "density") + facet_wrap(~site) + 
  labs(x = "CpG Methylation (Standardized Residuals)", y = "Scaled density", caption = "Red = SC; Blue = IVF") +
  theme(legend.position = "off", axis.text.x = element_blank(), axis.text.y = element_blank()) 
cpg_plot

sig_cpg_plot <- by_cpg %>% filter(site %in% c("cg13403462", "cg03904042", "cg14921437")) %>% 
  group_by(IVF, site) %>% summarize(mu = mean(val, na.rm = T)) %>% left_join(by_cpg, ., by = c("IVF","site")) %>% filter(!is.na(mu)) %>% 
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Non-IVF", "IVF"))) %>% 
  ggplot(aes(x = val, group = IVF, fill = IVF, color = IVF)) +
  geom_histogram() +
  #ggplot(aes(x = val, y = stat(scaled), group = IVF, fill = IVF, color = IVF)) +
  #geom_area(alpha = 0.6, stat = "density") + 
  geom_vline(aes(xintercept = mu, linetype = IVF)) + facet_wrap(~site) + 
  labs(x = "CpG Methylation (Standardized Residuals)", y = "Scaled density", group = "IVF status") +
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.text.y = element_blank()) 
sig_cpg_plot

# attempt to highlight
# cpg_plot + geom_rect(data = subset(by_cpg, site == "cg13403462"),
#                      aes(fill = site), 
#                      xmin = -Inf, xmax = Inf,
#                      ymin = -Inf, ymax = Inf, alpha = 0.3)
# Annotation Data
#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", version = "3.8")

# top hit
merged3 %>% filter(visit == 0) %>% 
  ggplot() + 
  geom_violin(aes(y = cg13403462, x = as.factor(IVF)), draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(aes(y = cg13403462, x = as.factor(IVF), color = as.factor(IVF)), position = position_jitterdodge(.15)) +
  labs(title = "Fetal cord tissue DNA methylation distributions at cg13403462, by IVF status.",
       y = "Normalized DNA Methylation (residuals)",
       x = "") +
  scale_x_discrete(labels = c("non-IVF", "IVF")) +
  theme(legend.position = "none")

# 6. Load and display EWAS association results
ewas <- read_xlsx("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/20181205-ewas.xlsx")
ewas %>% arrange(pval) %>% ggplot() + geom_histogram(aes(-log(pval)), bins = 200)
ewas %>% ggplot() + geom_point(aes(x = cpg, y = -log(pval)))

merged3 %>% filter(visit == 0) %>% ggplot() + 
  geom_histogram(aes(cg13403462, fill = as.factor(IVF)), binwidth = 0.005) +
  #geom_density(aes(cg13403462, fill = as.factor(IVF)), alpha = .75) +
  labs(title = "Fetal cord tissue DNA methylation distributions at cg13403462, by IVF status.",
       x = "Normalized DNA Methylation residuals",
       y = "Frequency",
       fill = "IVF status") +
  scale_fill_discrete(labels = c("non-IVF", "IVF")) +
  theme(legend.position = "bottom")

EWAS <- merged3 %>% filter(visit == 0) %>% select(IVF, contains("cg")) %>% gather(cpg, val, -IVF)
EWAS %>%
  ggplot() +
  geom_density(aes(x = val, fill = as.factor(IVF)), alpha = 0.75) +
  #geom_violin(aes(y = val, x = as.factor(IVF)), draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_wrap(~cpg, scales = "free") +
  labs(title = "Distributions of cord tissue candidate CpG methylation, by IVF status.", 
       x = "Normalized DNA methylation (residuals)", y = "") +
  theme(legend.position = "none", 
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        axis.text = element_blank())

#ewas_sub <- read_xls("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/20190110-cpg-results-sub.xls")
ewas_sub <- read_xls(paste0(out_path,"20190125-cpg-results.xls"))

sites <- read_csv("C:/Users/JHUANGYH/Google Drive/GUSTO Archive/HumanMethylation450_15017482_v1-2.csv", skip = 7)

annotate <- ewas_sub %>% mutate(IlmnID = if_else(str_length(pval3) == 5, paste0("cg000",pval3), 
                                                 if_else(str_length(pval3) == 6, paste0("cg00",pval3), 
                                                         if_else(str_length(pval3) == 7, paste0("cg0",pval3),
                                                                 if_else(str_length(pval3) == 8, paste0("cg",pval3), ""))))) %>% 
  mutate(lpval = -log10(pval2)) %>% 
  left_join(., sites) %>% mutate(CHR = as.double(CHR)) %>% filter(!is.na(CHR)) %>% 
  mutate(loc = as.double(paste0(Chromosome_36,".",Coordinate_36)))
# annotate %>% filter(is.na(CHR)) %>% arrange(-lpval) %>% select(lpval, IlmnID) #check which ones have no annotation
annotate %>% filter(str_detect(UCSC_RefGene_Name, "NECAB")) %>% select(pval2, IlmnID, CHR, UCSC_RefGene_Name) %>% arrange(pval2)

FLAG <- annotate %>% filter(lpval > -log10(0.05))
flag_id <- FLAG %>% arrange(-lpval) %>% pull(IlmnID)
SIG <- annotate %>% filter(lpval > -log10(0.05/281))
sig_id <- SIG %>% arrange(-lpval) %>% pull(IlmnID)
#annotate %>% arrange(lpval) %>% ggplot() + geom_histogram(aes(lpval), bins = 200)

CAND_LIST <- annotate %>% filter(str_detect(UCSC_RefGene_Name,"HIF") |  str_detect(UCSC_RefGene_Name,"NECAB")) %>% pull(IlmnID)
writeClipboard(CAND_LIST)

NOTE <- annotate %>% filter((str_detect(UCSC_RefGene_Name,"HIF") |  str_detect(UCSC_RefGene_Name,"NECAB")) & pval2 < 0.05) %>% 
  select(lpval, IlmnID, UCSC_RefGene_Name, CHR) %>% arrange(-lpval) %>% mutate(name = if_else(substring(UCSC_RefGene_Name,1,1) == "N",
                                                                                              substring(UCSC_RefGene_Name,1,6),
                                                                                              substring(UCSC_RefGene_Name,1,5)))

annotate %>% ggplot() + 
  geom_point(aes(x = CHR, y = lpval, color = IlmnID %in% flag_id)) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05/281), linetype = "dashed") +
  geom_hline(yintercept = -log10(1e-8), linetype = "solid") +
  #geom_text_repel(data = SIG, aes(x = CHR, y = lpval, label = paste0("(NECAB3) ",IlmnID), point.padding = 1) + 
  geom_text_repel(data = NOTE, aes(x = CHR, y = lpval, label = paste0(IlmnID, " (",name,")")), point.padding = 1) + 
  labs(x = "Chromosome", y = "-log(p)", title = "Cord tissue candidate CpG-IVF status associations (N = 934)",
       caption = "Adjusted for ethnicity, maternal age, parity, ppBMI, child sex") +
  scale_x_continuous(breaks = c(1:22)) +
  theme(legend.position = "none")

annotate %>% mutate(loc = as.double(paste0(Chromosome_36,".",Coordinate_36))) %>% 
  ggplot() + geom_point(aes(x = loc, y = lpval, color = factor(CHR))) +
  geom_hline(aes(yintercept = -log(0.05/281)), linetype = "dashed") +
  geom_hline(aes(yintercept = -log(0.05)), linetype = "dotted") +
  geom_text_repel(data = SIG, aes(x = loc, y = lpval, label = IlmnID), point.padding = 1) + 
  theme(legend.position = "none")

annotate %>% select(IlmnID, UCSC_RefGene_Name, Chromosome_36, Coordinate_36, lpval) %>% 
  filter(str_detect(UCSC_RefGene_Name,"NEC")) %>% 
  #filter(str_detect(UCSC_RefGene_Name,"HIF")) %>% 
  arrange(-lpval)

# Manhattan plot of candidate CpGs
library(qqman)
manhattan(annotate, chr = "CHR", bp = "Coordinate_36", snp = "IlmnID", p = "pval2", 
          annotatePval = 0.0001779359, annotateTop = F,
          suggestiveline = -log10(0.05/281))

# Investigate direction of associations
merged3 %>% filter(visit == 0) %>% select(flag_id, IVF) %>% 
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Spontaneous", "IVF"))) %>% 
  #group_by(IVF) %>% summarize_all(funs(mean(.,na.rm=T))) %>% 
  gather(IlmnID, val, -IVF) %>% 
  #left_join(., annotate, by = "IlmnID")  %>% 
  ggplot(aes(x = IVF, y = val, fill = IVF, color = IVF)) + 
  #geom_jitter(width = 0.2) +
  geom_violin(alpha = 0.2) +
  geom_boxplot(alpha = 0.4) +
  facet_wrap(~IlmnID) + 
  labs(title = "Distributions of differentially methylated CpGs",
    group = "IVF Status", fill = "IVF Status", color = "IVF Status", x = "", y = "") +
  theme(legend.position = "bottom", axis.text = element_blank())

merged3 %>% select(CAND_LIST, IVF) %>% ggplot() + geom_violin(aes(x = IVF, y = ))
merged3 %>% group_by(IVF) %>% summarize(NECAB1_1 = mean(cg14750836, na.rm = T), 
                                        HIF1A_1 = mean(cg04948941, na.rm = T), 
                                        HIF1A_2 = mean(cg13259118, na.rm = T),
                                        HIF3A_1 = mean(cg27146050, na.rm = T),
                                        HIF3A_2 = mean(cg16672562, na.rm = T),
                                        HIF3A_3 = mean(cg07022477, na.rm = T),
                                        HIF3A_4 = mean(cg11253785, na.rm = T))

# Height vs. CPG by visit and IVF status
merged3 %>% #filter(mother_ethnicity == "chinese") %>% 
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Spontaneous", "IVF"))) %>% 
  #ggplot(aes(x = cg03904042, y = height, group = IVF, color = IVF)) + 
  ggplot(aes(x = cg03904042, y = height)) + 
  geom_point() + geom_smooth(method = loess) + facet_wrap(~visit, scales = "free") +
  labs(title = expression(paste("Associations between cg03904042 (", italic(NECAB3), ") methylation, by visit month"), parse = T),
       x = "", y = "Height (cm)", color = "IVF Status") +
  theme(legend.position = "bottom")


#Q-Q Plots
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}
set.seed(42782)
pvalue <- runif(1000, min=0, max=1)
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

gg_qqplot(annotate$pval2) +
  theme_bw(base_size = 24) +
  labs(title = "QQ plot of IVF-methylation p-values") +
  annotate(geom = "text",
           x = -Inf,
           y = Inf,
           hjust = -0.15,
           vjust = 1 + 0.15 * 3,
           label = sprintf("?? = %.2f", inflation(annotate$pval2)),
           size = 8) +
  theme(axis.ticks = element_line(size = 0.5),
        #panel.grid = element_blank()
        panel.grid = element_line(size = 0.5, color = "grey80"))


# 7. POST-HOC: BASED ON INITIAL NECAB3 FINDING, FIND RELATED GENE CPGs
# look for more candidate CpGs
lin_cpgs <- c("cg00510507", "cg08390209", "cg23671997", "cg14300531", "cg25685359", "cg22383874", "cg02729344", "cg25487405")
sites %>% rename(cg_id = IlmnID) %>% select(cg_id, UCSC_RefGene_Name) %>% filter(cg_id %in% lin_cpgs)

sites %>% rename(cg_id = IlmnID) %>% select(cg_id, UCSC_RefGene_Name) %>% filter(cg_id %in% c("cg03904042","cg13403462", "cg14921437"))

sites %>% rename(cg_id = IlmnID) %>% select(cg_id, UCSC_RefGene_Name) %>% filter(str_detect(UCSC_RefGene_Name,"NECAB"))
sites %>% rename(cg_id = IlmnID) %>% select(cg_id, UCSC_RefGene_Name) %>% filter(str_detect(UCSC_RefGene_Name,"HIF"))

more_sites <- sites %>% rename(cg_id = IlmnID) %>% 
  select(cg_id, UCSC_RefGene_Name) %>% 
  filter(str_detect(UCSC_RefGene_Name,"HIF") | str_detect(UCSC_RefGene_Name,"NECAB"))
write_csv(more_sites, "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/20190115-NECAB_HIF_cg.csv")






########################################
# ADDITIONAL ANALYSES:
# 1. Contribution of PRS
# 2. Associations with IVF in subset
# 3. Assocations with pregnancy characterstics (mediation by)
# 4. Associations with outcomes (mediation of)
# 5. Associations with QMR measures
########################################
sub <- read_dta(paste0(data_path,"20190107-IVF-final-SUB.dta"))
sub_id <- sub %>% filter(visit == 0) %>% pull(SubjectID)
merged3_sub <- merged3 %>% filter(SubjectID %in% sub_id)
merged_MI <- read_dta(paste0(data_path,"20190115-IVF_Gen_MI_Drop.dta"))
miss_CPG_flag <- merged_MI %>% filter(`_mi_m` == 0, visit == 0, is.na(CPG)) %>% pull(SubjectID)
miss_PRS_flag <- merged_MI %>% filter(`_mi_m` == 0, visit == 0, is.na(SCORE)) %>% pull(SubjectID)
merged_CC <- merged_MI %>% filter(`_mi_m` == 1) %>% 
  mutate(miss_CPG = SubjectID %in% miss_CPG_flag,
         miss_PRS = SubjectID %in% miss_PRS_flag)

# Dx was CPG/PRS imputed with bias?
merged_CC %>% ggplot() + geom_histogram(aes(CPG, fill = miss_CPG)) + facet_wrap(~IVF)   # appears reasonable
merged_CC %>% ggplot() + geom_histogram(aes(SCORE, fill = miss_PRS)) + facet_wrap(~IVF) # may not be well imputed for IVF

# 1. Contribution of PRS
# Diagnostics of PRS
# histogram PRs by IVF, overall and subcohort
merged3_sub %>% filter(visit == 0) %>% ggplot() + geom_histogram(aes(SCORE, fill = as.factor(IVF)), binwidth = 0.05)
merged3 %>% filter(visit == 0) %>% ggplot() + geom_histogram(aes(SCORE, fill = as.factor(IVF)), binwidth = 0.05)

merged3 %>% filter(visit == 0) %>% tidy(lm(formula = SCORE ~ IVF + mother_ethnicity + parity + ppBMI))

# mean PRS by IVF, overall and subcohort (subcohort is closer: mean diff = 0.0129 vs. 0.02058 in overall cohort)
merged3_sub %>% filter(visit == 0, !is.na(SCORE)) %>% group_by(IVF) %>% summarize(n = n(), mean = sum(SCORE, na.rm = T)/n())
merged3 %>% filter(visit == 0, !is.na(SCORE)) %>% group_by(IVF) %>% summarize(n = n(), mean = sum(SCORE, na.rm = T)/n())

merged3 %>% filter(!is.na(SCORE)) %>% ggplot(aes(x = days, y = zwei, color = SCORE > 0)) + geom_point() + geom_smooth()

# Dx imputed PRS and CPG relationships
merged_CC %>% filter(visit == 0) %>% count(miss_CPG, miss_PRS)
merged_CC %>% ggplot(aes(x = SCORE, y = CPG, color = as.factor(miss_CPG | miss_PRS))) + geom_point() + geom_smooth() +
  labs(color = "Any Imputed Values") + theme(legend.position = "bottom") # note: among imputed, CPG and PRS perfectly uncorrelated

# 2. Associations with IVF in subset
# DEMONSTRATED IN THE REGRESSION FILE -- ASSOCIATION IS IDENTICAL


# 3. Assocations with pregnancy characterstics (mediation by)
# slightly lower in chinese
merged_CC %>% filter(visit == 0) %>% ggplot() + geom_histogram(aes(CPG, fill = mother_ethnicity)) + theme(legend.position = "bottom") #+ facet_wrap(~miss_CPG)

#relationship with IVF is same or stronger within ethnicity
merged_CC %>% filter(mother_ethnicity %in% c("chinese", "indian", "malay"), visit == 0) %>% 
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Non-IVF", "IVF")),
         ethnic = factor(mother_ethnicity, levels = c("chinese", "indian", "malay"), 
                         labels = c("Chinese (N = 662)", "Indian (N = 216)", "Malay (N = 299)"))) %>% 
  ggplot() + geom_area(aes(cg13403462, fill = IVF), stat = "density") + facet_wrap(~ethnic) + 
  labs(x = "cg13403462 methylation (standardized residuals)",
       y = "Density", title = "cg13403462 methylation distributions, by ethnicity and IVF status") + 
  theme(legend.position = "bottom")

# marginally higher with greater parity (driven predominantly by parity = 5)
merged_CC %>% filter(visit == 0, parity < 7) %>% ggplot(aes(x = parity, y = CPG)) + 
  geom_jitter() + geom_smooth(method = "loess") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))

merged_CC %>% filter(visit == 0, IVF == F) %>% 
  ggplot(aes(x = parity, y = CPG, group = parity)) + 
  geom_boxplot() + geom_smooth(method = "loess") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7))

# Higher in males (mean = 0.008 vs. -0.013)
merged3 %>% filter(visit == 0) %>% ggplot() + geom_histogram(aes(cg13403462, fill = sex)) + theme(legend.position = "bottom")
merged_CC %>% filter(visit == 0) %>% group_by(sex) %>% summarize(mean(CPG))

# unrelated to breastfeeding
merged3 %>% filter(visit == 0) %>% ggplot() + geom_histogram(aes(cg13403462, fill = Dur_any_BF)) + theme(legend.position = "bottom")

# CPG by BP; not obviously related
merged3 %>% filter(visit == 0) %>% 
  mutate(HTN_Hx = if_else(HTN_clean == "chronic HTN", 1, 
                          if_else(HTN_clean == "PIH", 2, 
                                  if_else(HTN_clean == "", 0, 3)))) %>% 
  ggplot() + geom_violin(aes(y = cg13403462, x = factor(HTN_Hx)), draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(alpha = 0.5, aes(y = cg13403462, x = factor(HTN_Hx)))

merged3 %>% filter(visit == 0) %>% ggplot() + geom_point(aes(x = sbp3, y = cg13403462))+ geom_smooth(aes(x = sbp3, y = cg13403462))

merged_CC %>% filter(visit == 0, parity == 0) %>% 
  select(IVF, CPG, sbp_1 = sbp1, sbp_2 = sbp2, sbp_3 = sbp3, dbp_1 = dbp1, dbp_2 = dbp2, dbp_3 = dbp3) %>% 
  gather(measure, val, -IVF, -CPG) %>% #filter(!is.na(val)) %>% group_by(measure) %>% count(IVF) # good coverage across measures by IVF status
  #separate(measure, c("measure", "visit")) %>% 
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Spontaneous", "IVF")),
         measure = factor(measure, levels = c("sbp_1","sbp_2","sbp_3","dbp_1","dbp_2","dbp_3"),
                          labels = c("Systolic (Visit 1)","Systolic (Visit 2)","Systolic (Visit 3)",
                                     "Diastolic (Visit 1)","Diastolic (Visit 2)","Diastolic (Visit 3)"))) %>% 
  ggplot(aes(x = val, y = CPG, color = IVF)) +
  labs(x = "Blood pressure (mmHg)",
       y = "Methylation (standardized residuals)",
       color = "IVF status",
       title = "cg13403462 methylation, by blood pressure measure and IVF status") +
  geom_point(alpha = 0.3) + geom_smooth() +
  facet_wrap(~measure) +
  theme(legend.position = "bottom")

# CPG by gestational weight; no clear relation, possible interaction by IVF, but post-hoc
merged_CC %>% filter(visit == 0) %>% 
  select(IVF, CPG, ppWeight, weight1:weight7) %>% 
  gather(measure, val, -IVF, -CPG) %>% filter(!is.na(val)) %>% #group_by(measure) %>% count(IVF) # good coverage across measures by IVF status
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Spontaneous", "IVF"))) %>% 
  ggplot(aes(x = val, y = CPG, color = IVF)) +
  labs(x = "Maternal weight (kg)",
       y = "Methylation (standardized residuals)",
       color = "IVF status",
       title = "cg13403462 methylation, by maternal weight measure and IVF status") +
  geom_point(alpha = 0.3) + geom_smooth(span = 0.2) +
  facet_wrap(~measure) +
  theme(legend.position = "bottom")
  
# 4. Association with outcomes
# Regression results CPG vs. Anthro
reg_out_new <- read_excel("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/16 Jan 2019-CPG-size_MI.xls", sheet = "adjusted")
out_all_visit <- reg_out_new %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei"))
all_visit = rep(c(0, 1, 3, 6, 9, 12, 15, 18, 24, 36, 48, 54, 60, 66, 72, 78), 6)
out_all_visit <- out_all_visit %>% add_column(all_visit)
reg_out_new %>% group_by(outcome) %>% summarize(max_N = max(N), min_N = min(N))
max_adj <- "Adjusted for: maternal age, education, ethnicity, household income, height, ppBMI, parity, 26 wk OGTT (fasting and 2-hour); 
paternal height and weight; smoke exposure in pregnancy; child sex, genetic risk score for adiposity,
GA at delivery (weeks), and age at measurement (days). Multiple imputation by chained equations for missing covariate values."

out_all_visit %>% filter(outcome %in% c("height", "log_bmi", "log_weight", "zbmi", "zlen", "zwei")) %>% 
  mutate(outcome = if_else(outcome == "zlen", "D. Height-for-age (SDS; N = 821 to 1174)",
                           if_else(outcome == "log_weight", "B. Weight (log(kg); N = 821 to 1178)",
                                   if_else(outcome == "height", "A. Height (cm; N = 821 to 1174)",
                                           if_else(outcome == "log_bmi", "C. BMI (log(kg/m^2); N = 821 to 1174)",
                                                   if_else(outcome == "zbmi", "F. BMI Z-score (SDS; N = 821 to 1174)", 
                                                           if_else(outcome == "zwei", "E. Weight-for-age (SDS; N = 821 to 1174)", ""))))))) %>% 
  ggplot() + 
  geom_point(aes(x = all_visit, y = estimate)) +
  geom_errorbar(aes(x = all_visit, ymin = min95, ymax= max95)) + 
  geom_ribbon(aes(x = all_visit, ymin = min95, ymax= max95), alpha = 0.15) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Associations between cg13403462 methylation and anthropometry.",
       x = "Age in Months",
       y = "",
       caption = max_adj) +
  facet_wrap(~outcome, scale = "free", labeller = label_bquote(.(outcome)))


# 5. Assocations with QMR measures
# merged3 <- read_dta(paste0(data_path,"20190125-merge_CpG.dta"))

merged3 %>% filter(visit == 0) %>% 
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("SC", "IVF"))) %>% 
  ggplot(aes(x = cg04948941, y = fat_6_qmr, group = IVF, color = IVF)) + geom_point() + geom_smooth(method = glm)

merged3 %>% filter(visit == 0) %>% select(IVF, contains("qmr")) %>% count(!is.na(wt_5_qmr), !is.na(wt_6_qmr))

merged3 %>% filter(visit == 0) %>% #filter(mother_ethnicity == "chinese") %>% 
  select(IVF, contains("_6_")) %>% 
  mutate(IVF = factor(IVF, levels = c(0,1), labels = c("Spontaneous", "IVF"))) %>%
  gather(measure, val, -IVF, -age_6_qr) %>% separate(measure, c("part", "year", "tech")) %>% select(-year, -tech) %>% 
  mutate(part = factor(part, levels = c("wt","lean","fat"), labels = c("Weight (kg)", "Lean mass (kg)", "Fat mass (kg)"))) %>% 
  ggplot(aes(x = IVF, y = val, color = IVF)) + geom_quasirandom(alpha = 0.75) + geom_boxplot(alpha = 0.25) + facet_wrap(~part) + 
  labs(x = "", y = "kg") +
  theme(legend.position = "none")

# HIF3A  = cg27146050
# NECAB3 = cg03904042
# HIF1A  = cg04948941


##############################################
# Sensitivity analyses unrelated to genetics
##############################################
# Age by Parity
merged3 %>% 
  filter(visit == 0, mother_ethnicity %in% c("chinese","indian","malay"), sex != "") %>% 
  group_by(parity) %>% summarize(mu = mean(mother_age_delivery, na.rm = T), 
                                 se = sd(mother_age_delivery, na.rm = T)/n(), 
                                 LL = mu - 1.96*se, 
                                 UL = mu + 1.96*se)

# MODE OF DELIVERY - DOESN'T DIFFER BY IVF
merged3 %>% filter(visit == 0, sex != "") %>% filter(SubjectID %in% sub_id) %>%
  #mutate(VD = if_else(substring(mode_of_delivery,1,1) %in% c("1", "2"), 1, 0)) %>% 
  group_by(IVF) %>% mutate(tot = n()) %>% 
  #ggplot() + geom_area(aes(x = mode_of_delivery, y = stat(scaled), fill = factor(IVF)), stat = "density") # pretty, but wrong
  ggplot() + geom_histogram(aes(mode_of_delivery, fill = factor(IVF)), stat = "count")
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
  

# Maternal GWG by IVF
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


