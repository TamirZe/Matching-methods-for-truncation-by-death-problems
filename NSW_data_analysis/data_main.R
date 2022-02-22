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
source("NSW_data_analysis/data_processing/data_processing_and_eda_funcs.R")
source("NSW_data_analysis/data_matching/data_matching_multiple.R")
source("NSW_data_analysis/data_post_matching_analysis/data_regression_estimators.R")
#source("NSW_data_analysis/data_post_matching_analysis/old/data_OLS_WLS_estimator_old.R")
source("NSW_data_analysis/data_post_matching_analysis/data_aligned_ranktest.R")
source("NSW_data_analysis/data_SA/data_SA_regression_funcs.R")
source("NSW_data_analysis/data_sensitivity_analyses/data_SA_regression_funcs.R")
#source("Simulations_studies/sim_DGM_and_simulations/simulation_run.R")
#source("Simulations_studies/sim_TABLES/table_design_multiple_func.R")
source("Ding_Lu/PS_M_weighting.R")
source("Ding_Lu/PS_M_weighting_SA.R")
source("Ding_Lu/xi_PS_M_weighting_SA.R")
source("Ding_Lu/DL_SE_boot.R")
########################################################################

# data files
########################################################################
library(readstata13); library(cem)
nsw <- read.dta13("NSW_data_analysis/data_files/nsw_dw.dta") # dehejia and wahba datase
data(LL, package = "cem") # LaLonde dataset 
########################################################################

set.seed(101)
########################################################################
data_bool = "LL" # "DW" for dehejia and wahba dataset # "LL" for LaLonde dataset 
# EM convergence parameters
iterations = 400; epsilon_EM = 1e-06

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
pi_as_est = mean(filter(data, A==0)$S)
pi_pro_est = mean(filter(data, A==1)$S) - mean(filter(data, A==0)$S)
pi_ns_est = 1 - mean(filter(data, A==1)$S)
pis_est = c(pi_as_est = pi_as_est, pi_pro_est = pi_pro_est, pi_ns_est = pi_ns_est)
# effect of A on S
S_on_A_test = prop.test(x = c( sum(filter(data, A==1)$S), sum(filter(data, A==0)$S)),
  n = c( nrow(filter(data, A==1)), nrow(filter(data, A==0)) ), p = NULL, alternative = "greater", correct = FALSE)
######################################################################## 

######################################################################## 
# naive estimators ####
# composite naive
most_naive_est = mean(data[A==1, Y]) - mean(data[A==0, Y]) #1794
most_naive_est_se = sqrt(  ( var(data[A==1, Y])  / nrow(data[A==1, ]) ) + 
                             ( var(data[A==0, Y])  / nrow(data[A==0, ]) )  )  
CI_by_SE_and_Z_val_most_naive = round(most_naive_est + c(-1,1) * 1.96 * most_naive_est_se, 3)
CI_by_SE_and_Z_val_most_naive = paste(CI_by_SE_and_Z_val_most_naive, sep = ' ', collapse = " , ")

# survivors naive
sur_naive_est = mean(data[A==1 & S == 1, Y]) - mean(data[A==0 & S == 1, Y]) #1340
sur_naive_est_se = sqrt(  ( var(data[A==1 & S==1, Y])  / nrow(data[A==1 & S==1, ]) ) + 
                            ( var(data[A==0 & S==1, Y])  / nrow(data[A==0 & S==1, ]) )  )
CI_by_SE_and_Z_val_sur_naive = round(sur_naive_est + c(-1,1) * 1.96 * sur_naive_est_se, 3)
CI_by_SE_and_Z_val_sur_naive = paste(CI_by_SE_and_Z_val_sur_naive, sep = ' ', collapse = " , ")

CI_naives_before_matching = data.frame(CI_by_SE_and_Z_val_most_naive, CI_by_SE_and_Z_val_sur_naive)
colnames(CI_naives_before_matching) = c("naive_without_matching", "survivors_naive_without_matching")
######################################################################## 

######################################################################## 
# EM algorithm ####
est_ding_lst = PSPS_M_weighting(Z=data$A, D=data$S,
        X=as.matrix(subset(data, select = covariates_PS)), Y=data$Y, trc = TRUE, ep1 = 1, ep0 = 1, 
        beta.a = NULL, beta.n = NULL, iter.max = iterations, error0 = epsilon_EM) 
#EM error
error = est_ding_lst$error
#EM coeffs
coeff_as = est_ding_lst$beta.a; coeff_ns = est_ding_lst$beta.n
EM_coeffs = rbind(coeff_as, coeff_ns)
colnames(EM_coeffs)[-1] = sub(".*]", "", colnames(EM_coeffs)[-1])

# adjust the cols the same order as in myEM: my order is: as, ns, pro. Ding order: c(prob.c, prob.a, prob.n)
PS_est = data.frame(est_ding_lst$PROB[,2], est_ding_lst$PROB[,3], est_ding_lst$PROB[,1])
colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
# add the principal scores to the data
data_with_PS = data.table(data, PS_est)
data_with_PS$e_1_as = data_with_PS$EMest_p_as / (data_with_PS$EMest_p_as + data_with_PS$EMest_p_pro)
######################################################################## 

######################################################################## 
# Ding and Lu estimators ####
DING_est = est_ding_lst$AACE
DING_model_assisted_est_ps = est_ding_lst$AACE.reg

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




