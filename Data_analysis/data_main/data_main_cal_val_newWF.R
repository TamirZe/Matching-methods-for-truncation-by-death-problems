# libraries for data analysis
########################################################################
library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); library(rlang);library(glue); library(gridExtra)
library(ggplot2); library(rockchalk); library(nnet); library(stats); library(rlist); library(mgsub); library(reshape2); library(gridExtra)
library(optmatch); library(DOS); library(Matching); library(sandwich); library(rmutil); library(clubSandwich); library(Hmisc)
library(sandwich); library(rmutil);  library(caret); library(splitstackshape); library(MatchIt); library(PerformanceAnalytics)
library(tidyr); library(dplyr); library(data.table); library(tidyr); library(tableone); library(lmtest)
library(readstata13); library(cem) # reading NSW datasets
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

covariates_PS =    c("age", "black", "hispanic", "married", "re75", "emp75") # "re75_square",
# adding intercept is for keeping the format of vars_names[-1] as in the simulations, since X1 in the simulations is the intercept
covariates_mahal = c("intercept", "age", "education", "re75")
reg_after_match =  c("intercept", "age", "education", "black", "hispanic", "married", "re75")
reg_BC =           c("intercept", "age", "education", "black", "hispanic", "married", "re75") 

# parameters and variables for matching and regressions ####
caliper_variable = "pi_tilde_as1"  
######################################################################## 

######################################################################## 
# adjust data  ####
data = LL
data = adjust_data(data, 1000, data_bool=data_bool) #DW #LL
#data$re75_square = data$re75^2
variables = setdiff(colnames(data), c("id", "A", "S", "Y", "OBS", "emp_74_75", "g"))
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
data_with_PS$pi_tilde_as1 = data_with_PS$EMest_p_as / (data_with_PS$EMest_p_as + data_with_PS$EMest_p_pro)
######################################################################## 

######################################################################## 
# Ding and Lu estimators ####
DL_est = c(DL_est = est_ding_lst$AACE, DL_MA_est = est_ding_lst$AACE.reg)
# bootstrap for DL estimator
#boosting_results = run_boosting(data, BS=500, seed=19, iter.max=iterations, error0=epsilon_EM)
######################################################################## 

######################################################################## 
# balance in the full dataset
balance_full_data = covarites_descriptive_table_cont_disc(dat=data_with_PS, cov_descr=variables)
# balance in the survivors (employed)
balance_employed = covarites_descriptive_table_cont_disc(dat=data_with_PS[S==1], cov_descr=variables, metric="Employed")
# combine full dataset and survivors (employed)
variables_remove = c("N", "N_match","N_unq", "EMest_p_as")
variables_names_balance = setdiff(balance_full_data$Variable, variables_remove)
balance_before_matching = merge(balance_full_data, balance_employed ,by="Variable", all.x = F)
balance_before_matching = balance_before_matching[match(variables_names_balance, balance_before_matching$Variable), ]
######################################################################## 

######################################################################## 
list_results_by_caliper = list()
caliper_values = seq(0.05,0.5,0.05)
for (i in 1:length(caliper_values)){
  
  m_data=data_with_PS[S==1]
  matching_datasets_lst = list()
  replace_vec = c(FALSE, TRUE)
  for(j in c(1:length(replace_vec))){
    matching_datasets_lst[[j]] = 
      matching_all_measures_func(m_data=m_data, match_on=caliper_variable, 
                                 covariates_mahal=covariates_mahal, reg_BC=reg_BC, X_sub_cols=variables, 
                                 M=1, replace=replace_vec[j], estimand="ATC", caliper=caliper)
  }
  
  matching_measures = c("PS", "mahal", "mahal_cal")
  post_matching_analysis_lst = list()
  matching_estimators = matching_estimators_SE = matching_estimators_CI = c()
  for(j in c(1:length(replace_vec))){
    replace=replace_vec[j]; all_measures_matched_lst=matching_datasets_lst[[j]] # j=1: wout replacement, j=2: with replacement
    post_matching_analysis_lst[[j]] =
      post_matching_analysis_func(m_data=m_data, replace=replace, all_measures_matched_lst=all_measures_matched_lst, reg_covariates=reg_after_match)
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
  ######################################################################## 
  
  ######################################################################## 
  estimators = c(naive_estimators, DL_est, matching_estimators)
  # estimated SE of matching estimators
  DL_na = -101 # DL estimators do not include SE and CI
  SE = c(naive_estimators_SE, c(DL_est=DL_na, DL_MA_est=DL_na), matching_estimators_SE)
  # CI of matching estimators
  CI = c(unlist(naive_estimators_CI), c(DL_est=DL_na, DL_MA_est=DL_na), matching_estimators_CI)
  ######################################################################## 
  
  ######################################################################## 
  # balance  ####
  # balance in the matched dataset, using 3 distance measures, with and wout replacement
  variables_balance_match = 
    c("Metric","N","N_match","N_unq","EMest_p_as","age","education","re74","re75","black","hispanic","married","nodegree","emp74","emp75")
  balance_match_wout = balance_match_with = data.frame(Variable = variables_balance_match)
  print(names(matching_datasets_lst[[2]])[1:length(matching_measures)]) # check this is the same order as matching_measures
  for (l in 1:length(matching_measures)){
    matching_lst_measure_wout = matching_datasets_lst[[1]][[l]]
    matching_lst_measure_with = matching_datasets_lst[[2]][[l]]
    
    balance_match_wout = merge(balance_match_wout, balance_after_matching_newWF(m_data=m_data, 
          match_obj=matching_lst_measure_wout$match_obj, dt_match=matching_lst_measure_wout$dt_match_S1,
          X_sub_cols=variables, metric=matching_measures[l]), by="Variable")
    
    balance_match_with = merge(balance_match_with, balance_after_matching_newWF(m_data=m_data, 
          match_obj=matching_lst_measure_with$match_obj, dt_match=matching_lst_measure_with$dt_match_S1,
          X_sub_cols=variables, metric=matching_measures[l]), by="Variable")
  }
  
  BALANCE_TABLE_wout = arrange_balance_table(balance_before_matching=balance_before_matching, balance_match=balance_match_wout,
                 variables_balance_match=variables_balance_match, matching_measures=matching_measures)
  BALANCE_TABLE_with = arrange_balance_table(balance_before_matching=balance_before_matching, balance_match=balance_match_with,
                 variables_balance_match=variables_balance_match, matching_measures=matching_measures)
  ######################################################################## 
  
  list_results_by_caliper[[i]] = list(EM_coeffs=EM_coeffs, estimators=estimators,
      BALANCE_TABLE_with=BALANCE_TABLE_with, BALANCE_TABLE_wout=BALANCE_TABLE_wout)
  names(list_results_by_caliper)[i] = paste0("cliper_", caliper_values[i])
}

######################################################################## 
# arrange in a table ####
BALANCE_TABLE_with_caliper_values = list_results_by_caliper$cliper_0.05$BALANCE_TABLE_with
BALANCE_TABLE_with_caliper_values = BALANCE_TABLE_with_caliper_values[,BALANCE_TABLE_with_caliper_values[1,]!="mahal_cal"]

tmp_res3 = NULL
for (i in 1:length(caliper_values)){
  tmp_res = list_results_by_caliper[[i]]$BALANCE_TABLE_with
  variables_col = c("", tmp_res$Variable)
  tmp_res = tmp_res[,tmp_res[1,]=="mahal_cal"]
  tmp_res = rbind(paste0("cal=", caliper_values[i]), tmp_res)
  if(i %% 2 == 1){tmp_res3 = tmp_res}else{tmp_res3 = cbind(tmp_res3, tmp_res)}
  if(i %% 2 == 0){
    tmp_res3 = data.frame(Variable = variables_col, tmp_res3)
    colnames(tmp_res3) <- gsub("\\..*","", colnames(tmp_res3))
    BALANCE_TABLE_with_caliper_values = rbind(BALANCE_TABLE_with_caliper_values, tmp_res3[!tmp_res3$Variable=="Metric",])
    tmp_res3 = NULL
  }
}
######################################################################## 

print(BALANCE_TABLE_with_caliper_values %>% xtable(), size="\\fontsize{9pt}{9pt}\\selectfont", include.rownames=F)
compare_PS_and_mahal = balance_match_with[match(variables_balance_match, balance_match_with$Variable), 
                   balance_match_with[balance_match_with$Variable=="Metric",]!="mahal_cal"]

