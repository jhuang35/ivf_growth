library(ctmle)
library(ggpubr)
library(SuperLearner)

Q.lib <- c("SL.glm", "SL.glm.interaction", "SL.bayesglm", "SL.mean", "SL.nnet", "SL.xgboost")

### SIZE
outcomes <- c("height", "zlen", "log_weight", "zwei", "zbmi", "log_bmi")
visits <- ivf %>% count(visit) %>% pull(visit)
N_visits = length(visits)
N_outcomes = length(outcomes)

anthro <- ivf %>% 
  mutate(L0.b = as.factor(MOM_EDU), L0.c = as.factor(ETHNIC), 
         L0.d = as.factor(HH_INC), L0.h = as.factor(smk_home)) %>% 
  dplyr::select(A1 = IVF,
                L0.a = mother_age_delivery, L0.b, L0.c, L0.d,
                L0.e = m_height_pw26, L0.f = ppBMI, L0.g = parity, 
                L0.h, L0.i = f_height_m24, L0.j = f_weight_m24,
                L0.k = MALE, L0.l = SCORE, L0.m = father_age_delivery,
                L0.n = father_HBP, L0.o = father_diabetes,
                outcomes, visit) %>% 
  dplyr::select(contains("L0"), A1, outcomes, visit)

Lnodes <- c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k", "L0.l", "L0.m", "L0.n", "L0.o")

# W <- anthro %>% dplyr::select(Lnodes)
# A <- anthro %>% dplyr::select(A1) %>% pull(A1)
# Y <- anthro %>% dplyr::select(Y2) %>% pull(Y2)
# 
# ctmle_discrete_fit1 <- ctmleDiscrete(Y = Y, A = A, W = data.frame(W), V = 5, preOrder = FALSE, detailed = TRUE)

set.seed(12345)
ctmle_results <- NULL
for(i in 1:N_visits){
  for(y in 1:N_outcomes){
    paste0("Estimating ATE on: ", outcomes[y], " at time: ", visits[i]) %>% print()
    subset <- anthro %>% filter(visit == visits[i])
    subset <- subset %>% dplyr::select(contains("L0"),A1, Y2 = outcomes[y]) %>% filter(!is.na(Y2))
    W <- subset %>% dplyr::select(Lnodes)
    A <- subset %>% dplyr::select(A1) %>% pull(A1)
    Y <- subset %>% dplyr::select(Y2) %>% pull(Y2)
    model <- ctmleDiscrete(Y = Y, A = A, W = data.frame(W), V = 5, 
                           SL.library = Q.lib,
                           preOrder = FALSE, detailed = TRUE)
    ATE <- model$est
    SE <- model$var.psi^0.5
    VIS = visits[i]
    OUT = outcomes[y]
    ctmle_results <- rbind(ctmle_results, cbind(ATE, SE, VIS, OUT))
  }
}

ctmle_res <- as_tibble(ctmle_results) %>% mutate_at(c("ATE","SE","VIS"),list(~as.double(.)))

write_excel_csv(ctmle_res, "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/ivf_ctmle_anthro.csv")
# ctmle_res <- read_csv("C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/ivf_ctmle_anthro.csv")

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


### SKINFOLDS
Lnodes <- c("L0.a", "L0.b", "L0.c", "L0.d", "L0.e", "L0.f", "L0.g", "L0.h", "L0.i", "L0.j", "L0.k", "L0.l", "L0.m", "L0.n", "L0.o")

skin_outcomes <- c("log_tri", "log_bi", "log_sub", "log_supra")
N_skin_outcomes <- length(skin_outcomes)

skin <- ivf %>% 
  mutate(L0.b = as.factor(MOM_EDU), L0.c = as.factor(ETHNIC), 
         L0.d = as.factor(HH_INC), L0.h = as.factor(smk_home)) %>% 
  dplyr::select(A1 = IVF,
                L0.a = mother_age_delivery, L0.b, L0.c, L0.d,
                L0.e = m_height_pw26, L0.f = ppBMI, L0.g = parity, 
                L0.h, L0.i = f_height_m24, L0.j = f_weight_m24,
                L0.k = MALE, L0.l = SCORE, L0.m = father_age_delivery,
                L0.n = father_HBP, L0.o = father_diabetes,
                skin_outcomes, visit) %>% 
  dplyr::select(contains("L0"), A1, skin_outcomes, visit)

count_skin_visits <- function(var){
  skin %>% group_by(visit) %>% summarize(N = sum(!is.na(!!sym(var)))) %>% 
    filter(N != 0) %>% pull(visit)
 }
skin_mat <- lapply(skin_outcomes, count_skin_visits)

set.seed(12345)
ctmle_skin_results <- NULL
for(y in 1:N_skin_outcomes){
  visits <- skin_mat[[y]]
  N_skin_visits <- length(visits)
  for(i in 1:N_skin_visits){
    paste0("Estimating ATE on: ", skin_outcomes[y], " at time: ", visits[i]) %>% print()
    subset <- skin %>% filter(visit == visits[i])
    subset <- subset %>% dplyr::select(contains("L0"),A1, Y2 = skin_outcomes[y]) %>% filter(!is.na(Y2))
    W <- subset %>% dplyr::select(Lnodes)
    A <- subset %>% dplyr::select(A1) %>% pull(A1)
    Y <- subset %>% dplyr::select(Y2) %>% pull(Y2)
    model <- ctmleDiscrete(Y = Y, A = A, W = data.frame(W), V = 5, 
                           SL.library = Q.lib,
                           preOrder = FALSE, detailed = TRUE)
    ATE <- model$est
    SE <- model$var.psi^0.5
    VIS = visits[i]
    OUT = skin_outcomes[y]
    ctmle_skin_results <- rbind(ctmle_skin_results, cbind(ATE, SE, VIS, OUT))
  }
}

ctmle_skin_res <- as_tibble(ctmle_skin_results) %>% mutate_at(c("ATE","SE","VIS"),list(~as.double(.)))

write_excel_csv(ctmle_skin_res, "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/ivf_ctmle_skin.csv")

SF <- ctmle_skin_res %>% 
  mutate(OUT = case_when(OUT == "log_tri" ~ "A. Triceps (N = 800 to 1119)",
                         OUT == "log_sub" ~ "B. Subscapular (N = 765 to 1118)",
                         OUT == "log_bi" ~ "C. Biceps (N = 765 to 892)",
                         OUT == "log_supra" ~ "D. Suprailiac (N = 785 to 862)")) %>% 
  mutate(LB = ATE-(1.96*SE), UB = ATE+(1.96*SE)) %>% 
  mutate_at(c("ATE","LB","UB"),funs((exp(.)-1)*100)) %>% 
  ggplot(aes(x = VIS, y = ATE)) + geom_point() + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(x = VIS, ymin = LB, ymax = UB)) +
  geom_ribbon(aes(x = VIS, ymin = LB, ymax = UB, fill = OUT), alpha = 0.25) +
  scale_y_continuous(limits = c(-25,15)) +
  scale_x_continuous(breaks = c(0, 18, 36, 48, 60, 72)) +
  labs(x = "", y = "% difference in skinfold thickness
       (IVF - SC)") + theme(legend.position = "none") +
  facet_wrap(~OUT) + theme(legend.position = "NONE")


### BLOOD PRESSURE

bp_outcomes <- c("SBP", "DBP", "ZSBP", "ZDBP")
N_bp_outcomes <- length(bp_outcomes)
bp_visits <- c(36, 48, 60, 72)
N_bp_visits <- length(bp_visits)

bp <- ivf %>% 
  mutate(L0.b = as.factor(MOM_EDU), L0.c = as.factor(ETHNIC), 
         L0.d = as.factor(HH_INC), L0.h = as.factor(smk_home)) %>% 
  dplyr::select(A1 = IVF,
                L0.a = mother_age_delivery, L0.b, L0.c, L0.d,
                L0.e = m_height_pw26, L0.f = ppBMI, L0.g = parity, 
                L0.h, L0.i = f_height_m24, L0.j = f_weight_m24,
                L0.k = MALE, L0.l = SCORE, L0.m = father_age_delivery,
                L0.n = father_HBP, L0.o = father_diabetes,
                bp_outcomes, visit) %>% 
  dplyr::select(contains("L0"), A1, bp_outcomes, visit)

set.seed(12345)
ctmle_bp_results <- NULL
for(y in 1:N_bp_outcomes){
  for(i in 1:N_bp_visits){
    paste0("Estimating ATE on: ", bp_outcomes[y], " at time: ", bp_visits[i]) %>% print()
    subset <- bp %>% filter(visit == bp_visits[i])
    subset <- subset %>% dplyr::select(contains("L0"),A1, Y2 = bp_outcomes[y]) %>% filter(!is.na(Y2))
    W <- subset %>% dplyr::select(Lnodes)
    A <- subset %>% dplyr::select(A1) %>% pull(A1)
    Y <- subset %>% dplyr::select(Y2) %>% pull(Y2)
    model <- ctmleDiscrete(Y = Y, A = A, W = data.frame(W), V = 5, 
                           SL.library = Q.lib,
                           preOrder = FALSE, detailed = TRUE)
    ATE <- model$est
    SE <- model$var.psi^0.5
    VIS = bp_visits[i]
    OUT = bp_outcomes[y]
    ctmle_bp_results <- rbind(ctmle_bp_results, cbind(ATE, SE, VIS, OUT))
  }
}

ctmle_bp_res <- as_tibble(ctmle_bp_results) %>% mutate_at(c("ATE","SE","VIS"),list(~as.double(.)))

write_excel_csv(ctmle_bp_res, "C:/Users/JHUANGYH/Dropbox/00 - Singapore/GUSTO/01 - PROJECTS/IVF and Growth/Output/ivf_ctmle_bp.csv")

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
  labs(x = "Age (Months)", y = "Absolute (mmHg) or standardized (%ile rank) 
       difference, (IVF - SC)") +
  facet_wrap(~OUT) + theme(legend.position = "NONE")

ggarrange(ANT, ggarrange(SF, BP, labels = c("B", "C"), ncol = 2, nrow = 1), 
          labels = c("A"), ncol = 1, nrow = 2)











##### C-TMLE GENERAL 5-FOLD CV

N <- dim(tmle_72)[1]
V <- 5
W <- tmle_72 %>% dplyr::select(Lnodes) %>% as.matrix(.)
A <- tmle_72 %>% dplyr::select(A1) %>% pull(A1)
Y <- tmle_72 %>% dplyr::select(Y2) %>% pull(Y2)

folds <- by(sample(1:N,N), rep(1:V, length=N), list)

lasso_fit <- cv.glmnet(x = as(W, "dgCMatrix"), y = A, alpha = 1, nlambda = 100, nfolds = 10)
lasso_lambdas <- lasso_fit$lambda[lasso_fit$lambda <= lasso_fit$lambda.min][1:5]

gn_seq <- build_gn_seq(A = A, W = data.frame(W), SL.library = SL.library, folds = folds)

Q <- cbind(rep(mean(Y[A == 1]), N), rep(mean(Y[A == 0]), N))
ctmle_fit1 <- ctmleGeneral(Y = Y, A = A, W = data.frame(W), folds = folds, 
                           ctmletype = 1, Q = Q,
                           gn_candidates = gn_seq$gn_candidates,
                           gn_candidates_cv = gn_seq$gn_candidates_cv)

# # Build template for glmnet
# SL.glmnet_new <- function (Y, X, newX, family, obsWeights, id, alpha = 1,
#                            nlambda = 100, lambda = 0,...){
#   if (!is.matrix(X)) {
#     X <- model.matrix(~-1 + ., X)
#     newX <- model.matrix(~-1 + ., newX)
#   }
#   fit <- glmnet::glmnet(x = X, y = Y,
#                         lambda = lambda,
#                         family = family$family, alpha = alpha)
#   pred <- predict(fit, newx = newX, type = "response")
#   fit <- list(object = fit)
#   class(fit) <- "SL.glmnet"
#   out <- list(pred = pred, fit = fit)
#   return(out)
# }
# 
# # Use a sequence of estimator to build gn sequence:
# SL.cv1lasso <- function (... , alpha = 1, lambda = lasso_lambdas[1]){
#   SL.glmnet_new(... , alpha = alpha, lambda = lambda)
# }
# SL.cv2lasso <- function (... , alpha = 1, lambda = lasso_lambdas[2]){
#   SL.glmnet_new(... , alpha = alpha, lambda = lambda)
# }
# SL.cv3lasso <- function (... , alpha = 1, lambda = lasso_lambdas[3]){
#   SL.glmnet_new(... , alpha = alpha, lambda = lambda)
# }
# SL.cv4lasso <- function (... , alpha = 1, lambda = lasso_lambdas[4]){
#   SL.glmnet_new(... , alpha = alpha, lambda = lambda)
# }
# SL.library = c('SL.cv1lasso', 'SL.cv2lasso', 'SL.cv3lasso', 'SL.cv4lasso', 'SL.glm')