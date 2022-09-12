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
source("Data_analysis/data_matching.R")
source("Data_analysis/data_regression_estimators.R")
source("Data_analysis/data_aligned_ranktest.R")
source("Data_analysis/data_sensitivity_analyses/data_SA_regression_funcs.R")
source("Simulations/EM_seq.R")
#source("Simulations/PS_M_weighting.R")
#source("Simulations/PS_M_weighting_SA_CPSR")
source("Data_analysis/DL_SE_boot.R")
########################################################################

# data files
########################################################################
library(readstata13); library(cem)
nsw <- read.dta13("NSW_data_analysis/data_files/LL_DW_datasets/nsw_dw.dta") # dehejia and wahba dataset
data(LL, package = "cem") # LaLonde dataset 
########################################################################

set.seed(101)
########################################################################
data_bool = "LL" # "DW" for dehejia and wahba dataset # "LL" for LaLonde dataset 
# EM parameters
two_log_est_EM = FALSE
iterations_EM = 400; epsilon_EM = 1e-06

# parameters and variables for matching and regressions ####
match_on = "e_1_as"  # "e_1_as" # EMest_p_as 
caliper = 0.3

covariates_PS = c("age", "black", "hispanic", "married", "re75", "emp75")
cont_cov_mahal = c("age", "education", "re75")
reg_after_match = c("age", "education", "black", "hispanic", "married", "re75")
reg_BC =          c("age", "education", "black", "hispanic", "married", "re75") 
######################################################################## 

######################################################################## 
# adjust data  ####
data = LL
data = adjust_data(data, 1000, data_bool=data_bool) #DW #LL
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
# xi_2log_PSPS_M_weighting
est_ding_lst = PSPS_M_weighting(Z=data$A, D=data$S,
        X=as.matrix(subset(data, select = covariates_PS)), Y=data$Y, trc = TRUE, ep1 = 1, ep0 = 1, 
        beta.a = NULL, beta.n = NULL, iter.max = iterations_EM, error0 = epsilon_EM) 

if(two_log_est_EM == FALSE){
#S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
  fit_S0_in_A0 =
  glm(as.formula(paste0("S ~ ",paste(covariates_PS, collapse="+"))), data=filter(data, A==0), family="binomial")
  beta_S0 = fit_S0_in_A0$coefficients
}else{beta_S0=NULL}

est_ding_lst2 = xi_2log_PSPS_M_weighting(Z=data$A, D=data$S,
     X=as.matrix(subset(data, select = covariates_PS)), Y=data$Y,
     xi_est=0, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL,
     iter.max=iterations_EM, error0=epsilon_EM)

est_ding_lst3 = xi_2log_PSPS_M_weighting(Z=data$A, D=data$S,
        X=as.matrix(subset(data, select = covariates_PS)), Y=data$Y, 
        xi_est=0, beta.S0=NULL, beta.ah=NULL, beta.c=NULL, 
        iter.max=iterations_EM, error0=epsilon_EM)

#EM error
error = est_ding_lst$error
#EM coeffs
coeff_as = est_ding_lst$beta.a; coeff_ns = est_ding_lst$beta.n
EM_coeffs = rbind(coeff_as, coeff_ns)
colnames(EM_coeffs)[-1] = sub(".*]", "", colnames(EM_coeffs)[-1])

######################################################################## 
EM_coeffs2 = rbind(fit_S0_in_A0$coefficients, est_ding_lst2$beta.ah, est_ding_lst2$beta.c)
EM_coeffs3 = rbind(fit_S0_in_A0$coefficients, est_ding_lst3$beta.ah, est_ding_lst3$beta.c)
######################################################################## 

# adjust the cols the same order as in myEM: my order is: as, ns, pro. Ding order: c(prob.c, prob.a, prob.n)
PS_est = data.frame(est_ding_lst$PROB[,2], est_ding_lst$PROB[,3], est_ding_lst$PROB[,1])
colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
# add the principal scores to the data
data_with_PS = data.table(data, PS_est)
data_with_PS$e_1_as = data_with_PS$EMest_p_as / (data_with_PS$EMest_p_as + data_with_PS$EMest_p_pro)
######################################################################## 

######################################################################## 
PS2 = data.frame(est_ding_lst2$ps.score)
PS3 = data.frame(est_ding_lst3$ps.score)
diff_PS23 = PS2[,-1] - PS3[,-1]
apply(diff_PS23, 2, mean); apply(diff_PS23, 2, sd); apply(abs(diff_PS23), 2, mean)

diff_PS2 = PS2[,-1] - PS_est
diff_PS3 = PS3[,-1] - PS_est
apply(diff_PS2, 2, mean); apply(diff_PS2, 2, sd); apply(abs(diff_PS2), 2, mean)
apply(diff_PS3, 2, mean); apply(diff_PS3, 2, sd); apply(abs(diff_PS3), 2, mean)
apply(diff_PS2, 2, hist)
apply(diff_PS3, 2, hist)
######################################################################## 

######################################################################## 
# Ding and Lu estimators ####
DL_est = c(DL_est = est_ding_lst$AACE, DL_MA_est = est_ding_lst$AACE.reg)

# bootstrap for DL estimator
'''boosting_results = run_boosting(data, BS=500, seed=19, iter.max=iterations, error0=epsilon_EM)'''
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
      matching_func_multiple_data(match_on = match_on,
          cont_cov_mahal=cont_cov_mahal,  reg_cov=reg_after_match, X_sub_cols=variables, 
          reg_BC=reg_BC, m_data=data_list[[l]], 
          w_mat_bool="NON-INFO", M=1, replace=replace_vec[j], estimand="ATC", mahal_match=2, caliper=caliper 
          #,OBS_table = descrip_all_data_OBS$OBS_table
          ,change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units=FALSE, vertical_table=TRUE, rnd=1)
    })
  }
######################################################################## 

######################################################################## 
'''# matching newWF ####
# perform matching with all units with S=1
# matching_datasets_lst[[1]] - wout replacement, matching_datasets_lst[[2]] - with replacement
m_data=data_with_PS[S==1]
m_data$id = c(1:nrow(m_data))
#Y1_as_mean = mean(filter(m_data, g=="as")$Y1); Y0_as_mean = mean(filter(m_data, g=="as")$Y0)
#SACE
matching_datasets_lst = list()
replace_vec = c(FALSE, TRUE)
for(j in c(1:length(replace_vec))){
  matching_datasets_lst[[j]] = 
    matching_all_measures_func(m_data=m_data, match_on=match_on, X_sub_cols=X_sub_cols, 
                               M=1, replace=replace_vec[j], estimand="ATC", mahal_match=2, caliper=caliper)
}
#dim(matching_datasets_lst[[2]]$only_ps_lst$matched_data)
#length(unique(matching_datasets_lst[[2]]$only_ps_lst$matched_data$id))

# post-matching analysis for 2 (wout/with replacement) datasets
matching_measures = c("PS", "maha", "maha_cal", "BC", "BC_cal", "wilcox")
post_matching_analysis_lst = list()
matching_estimators = matching_estimators_SE = matching_estimators_CI = c()
for(j in c(1:length(replace_vec))){
  post_matching_analysis_lst[[j]] =
    post_matching_analysis_func(m_data=m_data, replace=replace_vec[j], all_measures_matched_lst=matching_datasets_lst[[j]])
  # extract estimators and SE + CI of crude/BC/HL matching estimators of all distance measures (for crude)
  for (l in 1:length(matching_measures)){
    # extract estimators, SE and CI of crude/BC/HL matching estimators
    
    est_tmp = post_matching_analysis_lst[[j]][[l]]$SACE_matching_est
    SE_tmp = post_matching_analysis_lst[[j]][[l]]$SACE_matching_SE
    CI_tmp = post_matching_analysis_lst[[j]][[l]]$CI         
    names(est_tmp) = names(SE_tmp) = names(CI_tmp) = paste0(matching_measures[l], "_rep_", replace_vec[j])
    matching_estimators = c(matching_estimators, est_tmp)
    matching_estimators_SE = c(matching_estimators_SE, SE_tmp)
    matching_estimators_CI = c(matching_estimators_CI, CI_tmp)
  }
}'''
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

# balance in the employed and in the matched dataset, using 3 distance measures
BALANCE_TABLE = rbind(lst_matching_estimators[[1]][[1]]$balance_table, lst_matching_estimators[[2]][[1]]$balance_table) 

BALANCE_TABLE_with = filter(BALANCE_TABLE, Replacements==TRUE, Variable != "N") %>% subset(select = -grep(".1|.2|Replacements", colnames(BALANCE_TABLE)))
BALANCE_TABLE_wout = filter(BALANCE_TABLE, Replacements==FALSE, Variable != "N") %>% subset(select = -grep(".1|.2|Replacements", colnames(BALANCE_TABLE)))
# balance in the full dataset, employed and matched dataset using mahalanobis with caliper
BALANCE_TABLE_with = cbind(filter(balance_full_data, !Variable %in% c("N", "S")), subset(BALANCE_TABLE_with, select = -Variable))
BALANCE_TABLE_wout = cbind(filter(balance_full_data, !Variable %in% c("N", "S")), subset(BALANCE_TABLE_wout, select = -Variable))

colnames(BALANCE_TABLE_with) <- colnames(BALANCE_TABLE_wout) <- gsub("\\..*","", colnames(BALANCE_TABLE_with))
print(BALANCE_TABLE_with[-c(4:6),] %>% xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")), size="\\fontsize{6pt}{6pt}\\selectfont", include.rownames=F)
print(BALANCE_TABLE_wout[-c(4:6),] %>% xtable(caption = paste0("Matched data-set means, ", data_bool ," Sample.")), size="\\fontsize{6pt}{6pt}\\selectfont", include.rownames=F)

# matching estimators ####
ESTIMATORS_TABLE = rbind(lst_matching_estimators[[1]][[1]]$summary_table, 
                         lst_matching_estimators[[2]][[1]]$summary_table) %>% data.frame() 
######################################################################## 




