# libraries for data analysis
########################################################################
library(readstata13); library(cem) # reading NSW datasets
library(rlist); library(locfit); library(plyr); library(dplyr); library(data.table)
library(nnet); library(xtable); library(rlang);library(gridExtra); library(tableone)
library(ggplot2); library(rockchalk); library(nnet); library(stats); library(mgsub); library(reshape2)
library(Matching); library(sandwich); library(clubSandwich); library(lmtest); library(splitstackshape)

library(caret); library(PerformanceAnalytics); library(Hmisc)
library(DOS); library(rmutil); library(rlist); library(glue); library(tidyr)
########################################################################

########################################################################
# source for data NSW_data_analysis
setwd("~/A matching framework for truncation by death problems")
source("Data_analysis/data_processing_and_eda_funcs.R")
source("Simulations/naive_estimation.R")
source("Data_analysis/data_matching.R")
source("Data_analysis/data_regression_estimators.R")
source("Data_analysis/data_aligned_ranktest.R")
source("Data_analysis/data_sensitivity_analyses/data_SA_regression_funcs.R")
source("EM/EM_seq.R")
#source("Simulations/PS_M_weighting.R")
#source("Simulations/PS_M_weighting_SA_CPSR")
source("Data_analysis/DL_SE_boot.R")
########################################################################

# data files
########################################################################
nsw <- read.dta13("NSW_data_analysis/data_files/LL_DW_datasets/nsw_dw.dta") # dehejia and wahba dataset
data(LL, package = "cem") # LaLonde dataset 
########################################################################

set.seed(101)
########################################################################
data_bool = "LL" # "DW" for dehejia and wahba dataset # "LL" for LaLonde dataset 
# EM parameters
# EM-seq or EM-multi
EM_est_seq = TRUE
# two_log_est_EM=FALSE: S(0)=1, is estimated within A=0, with label according to S, before the EM process.
two_log_est_EM = FALSE
iterations_EM = 500; epsilon_EM = 1e-06

covariates_PS = c("age", "black", "hispanic", "married", "re75", "emp75") # "re75_square",
cont_cov_mahal = c("age", "education", "re75")
reg_after_match = c("age", "education", "black", "hispanic", "married", "re75")
reg_BC =          c("age", "education", "black", "hispanic", "married", "re75") 

# parameters and variables for matching and regressions ####
match_on = "e_1_as"  # "e_1_as" # EMest_p_as 
caliper = 0.4 # 0.3 # 0.4
######################################################################## 

######################################################################## 
# adjust data  ####
data = LL
data = adjust_data(data, 1000, data_bool=data_bool) #DW #LL
#data$re75_square = data$re75^2
variables = setdiff(colnames(data), c("id", "A", "S", "Y", "OBS", "emp_74_75"))
######################################################################## 

######################################################################## 
# pi estimation ####
pis_est = pis_est_func(data, xi_est=0)
# effect of A on S
S_on_A_test = prop.test(x = c( sum(filter(data, A==1)$S), sum(filter(data, A==0)$S)),
                        n = c( nrow(filter(data, A==1)), nrow(filter(data, A==0)) ), p = NULL, alternative = "greater", correct = FALSE)
######################################################################## 

######################################################################## 
# naive estimators ####
# composite naive
# naive estimators
naive_sace_estimation = naive_sace_estimation_func(data)
naive_estimators = c(composite_naive = naive_sace_estimation$composite_naive_est, surv_naive = naive_sace_estimation$sur_naive_est)
naive_estimators_SE = c(composite_naive = naive_sace_estimation$composite_naive_se, surv_naive = naive_sace_estimation$sur_naive_se)
naive_estimators_CI = naive_sace_estimation$CI_naive_before_matching
CI_naives_before_matching = data.frame(naive_estimators_CI$composite_naive, naive_estimators_CI$composite_naive)
colnames(CI_naives_before_matching) = c("naive_without_matching", "survivors_naive_without_matching")
######################################################################## 

######################################################################## 
# EM algorithm ####

if(EM_est_seq == TRUE){ # EM-seq
  if(two_log_est_EM == FALSE){
    #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
    fit_S0_in_A0 =
      glm(as.formula(paste0("S ~ ",paste(covariates_PS, collapse="+"))), data=filter(data, A==0), family="binomial")
    beta_S0 = fit_S0_in_A0$coefficients
  }else{beta_S0=NULL}
  
  est_ding_lst = xi_2log_PSPS_M_weighting(Z=data$A, D=data$S,
                           X=as.matrix(subset(data, select = covariates_PS)), Y=data$Y,
                           xi_est=0, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL,
                           iter.max=iterations_EM, error0=epsilon_EM)
}else{ # EM-multi
  # EM
  # est_ding_lst = PSPS_M_weighting(Z=data$A, D=data$S,
  #   X=as.matrix(subset(data, select = covariates_PS)), Y=data$Y, trc = TRUE, ep1 = 1, ep0 = 1, 
  #   beta.a = NULL, beta.n = NULL, iter.max = iterations_EM, error0 = epsilon_EM)
  
  # est_ding_lst = xi_PSPS_M_weighting_SA(Z=data_for_EM$A, D=data_for_EM$S,
  #   X=as.matrix(subset(data_for_EM, select = grep(paste(paste0("^",X_sub_cols[-1], "$"), collapse="|"), colnames(data_for_EM)))),
  #   Y=data_for_EM$Y, eta=xi_assm, iter.max=iterations_EM, error0=epsilon_EM)  
}

#EM error
error = est_ding_lst$error
#EM coeffs
coeff_as = est_ding_lst$beta.a
coeff_pro = est_ding_lst$beta.c
EM_coeffs = rbind(coeff_as, coeff_pro)
colnames(EM_coeffs)[-1] = sub(".*]", "", colnames(EM_coeffs)[-1])


# adjust the cols the same order as in myEM: my order is: as, ns, pro. Ding order: c(prob.c, prob.a, prob.n)
PS_est = data.frame(est_ding_lst$ps.score)
# add the principal scores to the data
data_with_PS = data.table(data, PS_est)
data_with_PS$e_1_as = data_with_PS$EMest_p_as / (data_with_PS$EMest_p_as + data_with_PS$EMest_p_pro)
######################################################################## 

######################################################################## 
# Ding and Lu estimators ####
DL_est = c(DL_est = est_ding_lst$AACE, DL_MA_est = est_ding_lst$AACE.reg)
# bootstrap for DL estimator
#boosting_results = run_boosting(data, BS=500, seed=19, iter.max=iterations, error0=epsilon_EM)
######################################################################## 

######################################################################## 
# matching procedure ####
# match only among the employed
data_list = list(data_with_PS[S==1])
# MATCHING
lst_matching_estimators = list()
replace_vec = c(FALSE, TRUE)
for(j in c(1:length(replace_vec))){
  lst_matching_estimators[[j]] =
    lapply(1:length(data_list), function(l){
      #set.seed(101)
      matching_func_multiple_data(match_on=match_on,
            cont_cov_mahal=cont_cov_mahal, reg_cov=reg_after_match, X_sub_cols=variables, 
            reg_BC=reg_BC, m_data=data_list[[l]], 
            w_mat_bool="NON-INFO", M=1, replace=replace_vec[j], estimand="ATC", mahal_match=2, caliper=caliper, 
            boost_HL=FALSE, vertical_table=TRUE, rnd=1)
    })
}
######################################################################## 

######################################################################## 
# aligned_ranktets ####
data_pairs_lst = lst_matching_estimators[[2]][[1]]$data_pairs_lst
aligned_ranktets_lst = list()
for (measure in names(data_pairs_lst)) {
  data_new_grp = adjust_pairs_to_new_grp(data_pairs_lst[[measure]])
  aligned_ranktets_lst[[measure]] = alignedranktest(outcome=data_new_grp$Y, matchedset=data_new_grp$trt_grp, treatment=data_new_grp$A)
}
######################################################################## 

######################################################################## 
# print to LateX
# EM_coeffs ###
print(EM_coeffs %>% xtable(), size="\\fontsize{9pt}{9pt}\\selectfont", include.rownames=F)

# balance  ####
# balance in the full dataset
balance_full_data = covarites_descriptive_table_cont_disc(dat = data_with_PS, cov_descr = variables)
# balance in the survivors (employed)
balance_employed = covarites_descriptive_table_cont_disc(dat = data_with_PS[S==1], cov_descr = variables)


# balance in the employed and in the matched dataset, using 3 distance measures
BALANCE_TABLE = rbind(lst_matching_estimators[[1]][[1]]$balance_table, lst_matching_estimators[[2]][[1]]$balance_table) 

BALANCE_TABLE_with = filter(BALANCE_TABLE, Replacements==TRUE, Variable != "N") %>% 
  subset(select = -grep(".1|.2|Replacements", colnames(BALANCE_TABLE)))
BALANCE_TABLE_wout = filter(BALANCE_TABLE, Replacements==FALSE, Variable != "N") %>% 
  subset(select = -grep(".1|.2|Replacements", colnames(BALANCE_TABLE)))

# balance in the full dataset, employed and matched dataset using mahalanobis with caliper
BALANCE_TABLE_with = cbind(filter(balance_full_data, !Variable %in% c("N", "S")), 
                           subset(BALANCE_TABLE_with, select = -Variable))
BALANCE_TABLE_wout = cbind(filter(balance_full_data, !Variable %in% c("N", "S")), 
                           subset(BALANCE_TABLE_wout, select = -Variable))

colnames(BALANCE_TABLE_with) <- colnames(BALANCE_TABLE_wout) <- gsub("\\..*","", colnames(BALANCE_TABLE_with))
print(BALANCE_TABLE_with[-c(4:6, 9, 15),] %>% xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")), size="\\fontsize{6pt}{6pt}\\selectfont", include.rownames=F)
print(BALANCE_TABLE_wout[-c(4:6),] %>% xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")), size="\\fontsize{6pt}{6pt}\\selectfont", include.rownames=F)


# matching estimators ####
ESTIMATORS_TABLE = rbind(lst_matching_estimators[[1]][[1]]$summary_table, 
                         lst_matching_estimators[[2]][[1]]$summary_table) %>% data.frame() 
######################################################################## 




