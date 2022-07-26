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
source("Simulations/sim_matching.R")
source("Simulations/sim_post_matching_analysis.R")
source("Simulations/sim_regression_estimators.R")
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
# matching newWF ####
m_data=data_with_PS[S==1]
m_data$id = c(1:nrow(m_data)); m_data$g="unknown"
matching_datasets_lst = list()
replace_vec = c(FALSE, TRUE)
for(j in c(1:length(replace_vec))){
  matching_datasets_lst[[j]] = 
    matching_all_measures_func(m_data=m_data, match_on=match_on, X_sub_cols=variables, 
                               M=1, replace=replace_vec[j], estimand="ATC", mahal_match=2, caliper=caliper)
}

matching_measures = c("PS", "mahal", "mahal_cal")
post_matching_analysis_lst = list()
matching_estimators = matching_estimators_SE = matching_estimators_CI = c()
for(j in c(1:length(replace_vec))){
  replace=replace_vec[j]; all_measures_matched_lst=matching_datasets_lst[[j]] # j=1: wout replacement, j=2: with replacement
  post_matching_analysis_lst[[j]] =
    post_matching_analysis_func(m_data=m_data, replace=replace, all_measures_matched_lst=all_measures_matched_lst,  X_sub_cols=X_sub_cols)
  # extract estimators and SE + CI of crude/BC/HL matching estimators of all distance measures 
  for (l in 1:length(matching_measures)){
    # extract estimators, SE and CI of crude/BC/HL/regression matching estimators
    est_tmp = unlist(lapply(post_matching_analysis_lst[[j]][[l]][1:3], "[[", "SACE_matching_est"))
    SE_tmp = unlist(lapply(post_matching_analysis_lst[[j]][[l]][1:3], "[[", "SACE_matching_SE"))
    CI_tmp = unlist(lapply(post_matching_analysis_lst[[j]][[l]][1:3], "[[", "CI"))       
    names(est_tmp) = names(SE_tmp) = names(CI_tmp) = 
      paste0(matching_measures[l], c("_crude", "_BC", "_wilcox"), c("_No", "_Yes")[j], "_rep") # replace_vec[j]
    
    reg_estimator_tmp_lst = reg_estimator_per_measure(lst_one_measure=post_matching_analysis_lst[[j]][[l]], 
                                                      measure_name=matching_measures[l], replace=replace)
    
    matching_estimators = c(matching_estimators, est_tmp, reg_estimator_tmp_lst$reg_matching_estimators)
    matching_estimators_SE = c(matching_estimators_SE, SE_tmp, reg_estimator_tmp_lst$reg_matching_estimators_SE)
    matching_estimators_CI = c(matching_estimators_CI, CI_tmp, reg_estimator_tmp_lst$reg_matching_estimators_CI)
  }
}

# check ties in BC caliper (only after matching with replacement [[2]])
#unlist(lapply(lapply(post_matching_analysis_lst[[2]], "[[", "BC_inference_lst"), "[[", "BC_ties_multiple_treated"))
BC_ties = c(unlist(post_matching_analysis_lst[[2]]$ps_estimators$BC_inference_lst$BC_ties_multiple_treated), 
            unlist(post_matching_analysis_lst[[2]]$mahal_estimators$BC_inference_lst$BC_ties_multiple_treated),
            unlist(post_matching_analysis_lst[[2]]$mahal_cal_estimators$BC_inference_lst$BC_ties_multiple_treated))
######################################################################## 


######################################################################## 
# TODO adjust to newWF
#aligned_ranktets ####
data_pairs_lst = matching_datasets_lst[[2]][[1]]$data_pairs_lst
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

# TODO adjust to newWF
# balance in the employed and in the matched dataset, using 3 distance measures
#matching_datasets_lst[[1]]$balance_all_measures$mean_by_subset_mahal_cal
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




