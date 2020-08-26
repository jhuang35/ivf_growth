merged3 <- read_dta("C:/Users/JHUANGYH.ARES/My Tresors/SICS/Projects/IVF/merged_MI.dta") 

p1 <- merged3 %>% filter(visit == 0, mother_ethnicity %in% c("chinese", "indian", "malay")) %>% 
  mutate(mother_ethnicity = factor(mother_ethnicity, levels = c("chinese", "indian", "malay"),
                                   labels = c("Chinese", "Indian", "Malay"))) %>% 
  dplyr::select(sig_id, IVF, mother_ethnicity) %>%
  mutate(IVF = factor(IVF, levels = c(1,0), labels = c("IVF", "Spontaneous"))) %>% 
  gather(IlmnID, val, -IVF, -mother_ethnicity) %>% 
  mutate(IlmnID = factor(IlmnID, levels = c("cg03904042", "cg14921437", "cg13403462"), 
                         labels = c("cg03904042", "cg14921437", "cg13403462"), 
                         ordered = T)) 

writexl::write_xlsx(p1, paste0(data_path,"cpg_ethnicity.xlsx"))

efw <- merged3 %>% filter(visit == 0, sex != "") %>% filter(SubjectID %in% sub_id) %>% 
  dplyr::select(SubjectID, IVF, contains("EFW"), contains("DAYS_"), -contains("miss")) %>% 
  gather(var, val, -SubjectID, -IVF) %>% separate(var, c("MEASURE", "VISIT")) %>% 
  spread(MEASURE, val) %>% filter(!is.na(VISIT)) %>% 
  dplyr::select(-SubjectID) %>% 
  mutate(IVF = factor(IVF, levels = c(1,0), 
                      labels = c("Assisted Reproductive Technologies (ART)", "Spontaneous Conception (SC)")))

writexl::write_xlsx(efw, paste0(data_path,"efw.xlsx"))
