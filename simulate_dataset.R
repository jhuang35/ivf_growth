### CREATE SIMPLE PLASMODE FILE FOR PAPER

# Use PlasmodeContFixed function from 20200226-Plasmode_Correctedv4.R

{
  tictoc::tic()
  mi_data <- read_dta("D:/SICS/Data and Instruments/Archive/20190723-IVF_Dad_MI.dta")
  tictoc::toc()
}
mi2 <- mi_data %>% filter(`_mi_m` == 2, visit == 0)

outcomes <- c("height", "zlen", "log_weight", "zwei", "zbmi", "log_bmi")
anthro <- mi2 %>% filter(visit == 0) %>% 
  mutate(L0.b = as.factor(MOM_EDU), L0.c = as.factor(ETHNIC), 
         L0.d = as.factor(HH_INC), L0.h = as.factor(smk_home)) %>% 
  dplyr::select(A1 = IVF,
                L0.a = mother_age_delivery, L0.b, L0.c, L0.d,
                L0.e = m_height_pw26, L0.f = ppBMI, L0.g = parity, 
                L0.h, L0.i = f_height_m24, L0.j = f_weight_m24,
                L0.k = MALE, L0.l = SCORE, L0.m = father_age_delivery,
                L0.n = father_HBP, L0.o = father_diabetes,
                outcomes, visit) %>% 
  dplyr::select(contains("L0"), A1, outcomes, visit) %>% 
  mutate(ID = paste0("ID",row_number())) %>% data.frame(.)

Lnodes <- c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k", "L0.l", "L0.m", "L0.n", "L0.o")
gform <- "A1 ~ L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j+ L0.k+ L0.l+ L0.m + L0.n + L0.o"
Qpreds <- " ~ A1 + L0.a + L0.b + L0.c + L0.d + L0.e + L0.f + L0.g + L0.h + L0.i + L0.j+ L0.k+ L0.l+ L0.m + L0.n + L0.o"

N_sims <- 50
Effect_Size <- -0.5 # simulated risk difference = -0.5 (SD) for standardized scores

set.seed(42782)
sim_res <- NULL
for(i in 1:length(outcomes)){
  Qform <- paste0(outcomes[i], Qpreds)
  sim_df <- PlasmodeContFixed(formulaOut = as.formula(Qform), objectOut = NULL,
                              formulaExp = as.formula(gform), objectExp = NULL,
                              data = anthro, idVar = "ID", 
                              effectOR = Effect_Size, MMOut = 1, MMExp = 1, 
                              nsim = N_sims, size = nrow(anthro),
                              exposedPrev = NULL)
  sim_res <- rbind(sim_res, cbind(ID = sim_df$Sim_Data$ID1, A1 = sim_df$Sim_Data$EXPOSURE1,
                                  Y = sim_df$Sim_Data$OUTCOME1, OUTCOME = outcomes[i]))
}

#outForm <- as.formula(paste0("Y", Qpreds))

# for outcome only
# sim_final <- as_tibble(sim_res) %>% 
#   mutate_at(c("Y"), list(~ as.numeric(.))) %>%
#   left_join(., select(anthro, -outcomes), by = c("ID"))
  
# for outcome and exposure  
sim_final <- as_tibble(sim_res) %>% 
 mutate_at(c("A1", "Y"), list(~ as.numeric(.))) %>% 
 left_join(., select(anthro, -A1, -outcomes, -visit), by = c("ID"))

writexl::write_xlsx(data.frame(filter(sim_final, OUTCOME %in% c("zbmi", "zlen", "zwei"))),
                    "D:/SICS/Projects/IVF/Manuscript/SUBMISSIONS/Nature Communications/code/simdata.xlsx")
