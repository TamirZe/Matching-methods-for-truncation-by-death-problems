library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist)
library(nnet); library(locfit); library(splitstackshape); library(MASS)
library(Matching); library(sandwich); library(clubSandwich); library(lmtest); library(mgsub)

########################################################################
# source for Simulations_studies
setwd("~/A matching framework for truncation by death problems")
source("Simulations_studies/sim_parameters_and_pis/sim_set_parameters.R")
source("Simulations_studies/sim_parameters_and_pis/sim_check_pis_and_covariates.R")
source("Simulations_studies/sim_DGM_and_simulations/DGM_CPSR.R")
source("Simulations_studies/sim_naive_estimation/naive_estimation.R")
#source("Ding_Lu/PS_M_weighting.R")  # EM with one multinomial regression 
source("Ding_Lu_EM/Sequencial_logistic_regressions/EM_2log_CPSR.R") 
source("Simulations_studies/sim_matching_procedure/sim_matching.R")
source("Simulations_studies/sim_post_matching_analysis/sim_post_matching_analysis.R")
source("Simulations_studies/sim_post_matching_analysis/sim_regression_estimators.R")
#############################################################################################

set.seed(101)

#############################################################################################
# treatment probability
prob_A = 0.5

# parameters for simulating X
# @@@@@@@@@@@@ dim_x includes an intercept @@@@@@@@@@@@@@@
dim_x = 4; cont_x = dim_x - 1
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)
X_sub_cols = paste0("X", c(1:(dim_x)))
#############################################################################################

#############################################################################################
# CPSR parameter ####
xi = 0
xi_est = xi # xi
two_log_models_DGM = TRUE # DGM includes sequential 2 logistic regressions
#############################################################################################

#############################################################################################
# misspec parameters (for PS model and Y models) ####
misspec_PS = 2 # 0: no misspec of PS model # 2: PS functional form misspecification
funcform_factor_sqr=5; funcform_factor_log=-5 # funcform_factor_sqr=-3; funcform_factor_log=3
transform_x = 0
misspec_outcome = 2 # 0: no misspec of Y model # 2: Y functional form misspecification
#############################################################################################

#############################################################################################
# EM  parameters ####
# two_log_est_EM: S(0)=1, is estimated within A=0, with label according to S, before the EM process.
# In both cases, sequencial logistic models are being estimated. For one multinational model, use PS_M_weighting (and seq_log_or_multinomial_EM=FALSE)
two_log_est_EM = FALSE 
iterations_EM = 200; epsilon_EM = 10^-6
# seq_log_or_multinomial_EM = TRUE
#############################################################################################

#################################################################################################################
# beta_and_gamma ####
# scenario is determined according to mat_gamma (and its row, k) and betas_GPI

beta_and_gamma = set_parameters_func(dim_x=dim_x, high_pi_pro=T, AX_interactions=T) # high_pi_pro = T/F
# beta
betas_GPI = beta_and_gamma$betas_GPI
# gamma
mat_gamma = beta_and_gamma$mat_gamma

k=2 # k=1 (correct PS specification: pi as = 0.5) # k=2 (correct PS specification: pi as = 0.75)
gamma_ns = rep(0, dim_x)
gamma_ah = as.numeric(mat_gamma[k, c(1:dim_x)])
gamma_pro =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])  
##################################################################################################################

##################################################################################################################
# extract pis, wout and with PS misspecification ####
extract_pis_lst = extract_pis_from_scenarios(nn=500000, mat_gamma=mat_gamma, xi=xi, misspec_PS=0, two_log_models_DGM=T)
mat_pis_per_gamma = extract_pis_lst$mat_pis
mat_pis_per_gamma
extract_pis_lst_mis_PS = extract_pis_from_scenarios(nn=500000, mat_gamma=mat_gamma, xi=xi, misspec_PS=2, two_log_models_DGM=T)
mat_pis_per_gamma_mis_PS = extract_pis_lst_mis_PS$mat_pis
mat_pis_per_gamma_mis_PS
##################################################################################################################

###############################################################################################
# correlation structure between PO's
var_GPI = as.matrix(rbind(1, 1))
rownames(var_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4 
###############################################################################################

###############################################################################################
# matching arguments
# caliper is in sd
caliper = 0.25; match_on = "O11_posterior_ratio" 
###############################################################################################

###############################################################################################
param_n = 2000; param_n_sim = 50 # param_n = 2000; param_n_sim = 1000
mu_x_fixed = FALSE
###############################################################################################

###############################################################################################
# true SACE parameter from one large simulation
one_large_simulation = simulate_data_func(
  gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns,
  xi=xi, two_log_models_DGM=two_log_models_DGM, param_n=1000000, 
  misspec_PS=0, misspec_outcome=0, transform_x=transform_x,
  funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
  betas_GPI=betas_GPI, var_GPI=var_GPI, rho_GPI_PO=rho_GPI_PO)
true_SACE = one_large_simulation$true_SACE
true_SACE

# SACE parameter from mean of multiple simulations of sample size param_n
SACE_vec = vector(length = 500)
for (i in 1:length(SACE_vec)) {
  one_small_simulation = simulate_data_func(
    gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns,
    xi=xi, two_log_models_DGM=two_log_models_DGM, param_n=param_n, 
    misspec_PS=0, misspec_outcome=0, transform_x=transform_x,
    funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
    betas_GPI=betas_GPI, var_GPI=var_GPI, rho_GPI_PO=rho_GPI_PO)
  SACE_vec[i] = one_small_simulation$true_SACE
}
mean(SACE_vec)

###############################################################################################

###############################################################################################
# matrices and lists to retain results from all iterations ####
list_EM_not_conv <- list_mean_by_g <- balance_wout_rep_lst <- balance_with_rep_lst <- list()
pis_pis_est_obs_mat = NULL
coeff_ah_mat <- coeff_pro_mat <- beta_S0_mat <-  matrix(nrow = param_n_sim, ncol = dim_x) 
#balance_wout_rep <- balance_with_rep <- matrix(nrow = param_n_sim, ncol = (dim_x-1))
matching_estimators_mat <- matrix(nrow = param_n_sim, ncol = 21)
matching_estimators_SE_mat <- matrix(nrow = param_n_sim, ncol = 21)
CI_mat <- matrix(nrow = param_n_sim, ncol = 21)
BC_ties_multiple_treated_mat <- matrix(nrow = param_n_sim, ncol = 12)
OLS_NOint_mat<- OLS_YESint_mat <- matrix(nrow = param_n_sim, ncol =  3)
WLS_NOint_mat <- WLS_YESint_mat <- matrix(nrow = param_n_sim, ncol =  5)
scen_parameter_lst = list(true_SACE=true_SACE, param_n=param_n, param_n_sim=param_n_sim,
    gamma_ah=gamma_ah, gamma_pro=gamma_pro,
    xi=xi, xi_est=xi_est, two_log_models_DGM=two_log_models_DGM, 
    misspec_PS=misspec_PS, misspec_outcome=misspec_outcome, transform_x=transform_x,
    funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
    betas_GPI=betas_GPI, var_GPI=var_GPI, rho_GPI_PO=rho_GPI_PO)
###############################################################################################


# run over param_n_sim different iterations, in each iteration, sample contains param_n observations
index_EM_not_conv = 0; i_EM_not_conv = c() # check proportion of iteration when EM has not converged
#while (i <= param_n_sim){
for (i in 1:param_n_sim){
  print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match. ",
               "index_EM_not_conv: ", index_EM_not_conv, "."))
  start_time1 <- Sys.time()
  
  # DGM: simulate data
  list_data_DGM = simulate_data_func(seed_num=NULL, 
      gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, 
      xi=xi, two_log_models_DGM=two_log_models_DGM, param_n=param_n,
      misspec_PS=misspec_PS, misspec_outcome=misspec_outcome, transform_x=transform_x,
      funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
      betas_GPI=betas_GPI, var_GPI=var_GPI, rho_GPI_PO=rho_GPI_PO)
  
  data_for_EM = list_data_DGM$dt
  list_mean_by_g[[i]] = apply(list_data_DGM$mean_by_g, 2, as.numeric) # x_obs
  x_obs = list_data_DGM$x_obs
  x_PS = data.frame(list_data_DGM$x_PS)
  x_outcome = data.frame(list_data_DGM$x_outcome)
  pis = list_data_DGM$pis
  vec_OBS_table = list_data_DGM$vec_OBS_table 
  # "real" SACE parameter from current iteration
  SACE = list_data_DGM$true_SACE
  pis_est = pis_est_func(data_for_EM=data_for_EM, xi_est=xi_est)
  pis_pis_est_obs_mat = rbind(pis_pis_est_obs_mat, c(pis, t(pis_est), vec_OBS_table))
  colnames(pis_pis_est_obs_mat) = c(colnames(pis), names(pis_est), colnames(vec_OBS_table))
  
  # naive estimators
  naive_sace_estimation = naive_sace_estimation_func(data_for_EM)
  naive_estimators = c(composite_naive = naive_sace_estimation$composite_naive_est, surv_naive = naive_sace_estimation$sur_naive_est)
  naive_estimators_SE = c(composite_naive = naive_sace_estimation$composite_naive_se, surv_naive = naive_sace_estimation$sur_naive_se)
  naive_estimators_CI = naive_sace_estimation$CI_naive_before_matching
  
  # EM and PS estimation
  # If beta_S0=NULL, employ two logistic regressions during the EM
  if(two_log_est_EM == FALSE){
    #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
    fit_S0_in_A0 = glm(as.formula(paste0("S ~ ",paste(X_sub_cols[-1], collapse="+"))), data=filter(data_for_EM, A==0), family="binomial")
    beta_S0 = fit_S0_in_A0$coefficients
  }else{beta_S0=NULL}
  
  # EM
  start_timeDing <- Sys.time()
  est_ding_lst = xi_2log_PSPS_M_weighting(Z=data_for_EM$A, D=data_for_EM$S,
        X=as.matrix(subset(data_for_EM, select = 
        grep(paste(X_sub_cols[-1], collapse="|"), colnames(data_for_EM)))), Y=data_for_EM$Y, 
        xi_est=xi_est, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL, 
        iter.max=iterations_EM, error0=epsilon_EM)
  end_timeDing <- Sys.time()
  print(paste0("Ding EM lasts ", difftime(end_timeDing, start_timeDing)))
  
  # EM coefficient estimators and PS estimators
  PS_est = est_ding_lst$ps.score
  # add PS's to the data
  data_with_PS = data.table(data_for_EM, PS_est)
  
  # Ding and Lu's estimator: DL's plain estimator and DL's model assisted estimator
  DL_est = c(DL_est = est_ding_lst$AACE, DL_MA_est = est_ding_lst$AACE.reg)
  
  # if PS_est contains NAS, it probably implies that the EM process has not converged, so skip this iteration and go to the next
  if( sum(is.na(PS_est)) > 0 | (est_ding_lst$iter == iterations_EM + 1 & est_ding_lst$error >= epsilon_EM) ){ 
    print("EM has not convereged")
    i_EM_not_conv = c(i_EM_not_conv, i)
    index_EM_not_conv = index_EM_not_conv + 1
    list_EM_not_conv$probs[[index_EM_not_conv]] = PS_est
    list_EM_not_conv$coeffs[[index_EM_not_conv]] = data.frame(rbind(coeff_ah=est_ding_lst$beta.ah, coeff_pro=est_ding_lst$beta.c))
    colnames(list_EM_not_conv$coeffs[[index_EM_not_conv]]) = X_sub_cols
    list_EM_not_conv$probs_nas[[index_EM_not_conv]] = c(total_na = sum(is.na(PS_est)), 
              prop_na = sum(is.na(PS_est)) / ( nrow(PS_est) * ncol(PS_est) ) ) %>% round(3)
    next()
  }
  
  # EM coefficient estimators and PS estimators
  beta_S0_mat[i,] = beta_S0
  coeff_ah_mat[i,] = est_ding_lst$beta.ah
  coeff_pro_mat[i,] = est_ding_lst$beta.c
  
  # calculate weights: O11_prior_ratio, O11_posterior_ratio W_1_as, and W_1_as_true
  weights_lst = add_PS_weights_func(data_with_PS=data_with_PS, pis=pis, pis_est=pis_est)
  O11_prior_ratio = weights_lst$O11_prior_ratio
  O11_prior_ratio_true = weights_lst$O11_prior_ratio_true
  data_with_PS = weights_lst$data_with_PS 
  
  # matching
  # perform matching with all units with S=1
  # matching_datasets_lst[[1]] - wout replacement, matching_datasets_lst[[2]] - with replacement
  matching_datasets_lst = list()
  replace_vec = c(FALSE, TRUE)
  for(j in c(1:length(replace_vec))){
    matching_datasets_lst[[j]] = 
        matching_all_measures_func(m_data=data_with_PS[S==1], match_on=match_on, X_sub_cols=X_sub_cols, 
         M=1, replace=replace_vec[j], estimand="ATC", mahal_match=2, caliper = caliper)
  }
  
  # post-matching analysis for 2 (wout/with replacement) datasets
  matching_measures = c("PS", "maha", "maha_cal", "BC", "BC_cal", "wilcox")
  post_matching_analysis_lst = list()
  matching_estimators = matching_estimators_SE = matching_estimators_CI = c()
  for(j in c(1:length(replace_vec))){
    post_matching_analysis_lst[[j]] =
      post_matching_analysis_func(m_data=data_with_PS[S==1], replace=replace_vec[j], all_measures_matched_lst=matching_datasets_lst[[j]])
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
  }
  
  # check ties in BC caliper (only after matching with replacement [[2]])
  BC_ties = c(BC = unlist(post_matching_analysis_lst[[2]]$BC_inference_lst$BC_ties_multiple_treated),
    BC_cal = unlist(post_matching_analysis_lst[[2]]$BC_cal_inference_lst$BC_ties_multiple_treated))
  BC_ties_multiple_treated_mat[i,] = BC_ties
  colnames(BC_ties_multiple_treated_mat) = names(BC_ties)
  
  # regression estimators after matching on mahalanobis + PS caliper
  # est and SE
  OLS_NOint_mat  = c(unlist(post_matching_analysis_lst[[1]]$reg_wout_interactions_inference_lst$estimator_and_se))
  OLS_YESint_mat = c(unlist(post_matching_analysis_lst[[1]]$reg_with_interactions_inference_lst$estimator_and_se))
  WLS_NOint_mat  = c(unlist(post_matching_analysis_lst[[2]]$reg_wout_interactions_inference_lst$estimator_and_se))
  WLS_YESint_mat = c(unlist(post_matching_analysis_lst[[2]]$reg_with_interactions_inference_lst$estimator_and_se))
  regression_matching_estimators = c(OLS_NOint_mat[1], OLS_YESint_mat[1], WLS_NOint_mat[1], WLS_YESint_mat[1])
  regression_matching_estimators_SE = c(OLS_NOint_mat["OLS_se"], OLS_YESint_mat["OLS_se"], WLS_NOint_mat["WLS_clstr_se"], WLS_YESint_mat["WLS_clstr_se"])
  names(regression_matching_estimators) = names(regression_matching_estimators_SE) = c("OLS", "OLS_int", "WLS", "WLS_int")
  # CI of regression estimators
  regression_matching_estimators_CI = c(OLS = unlist(post_matching_analysis_lst[[1]]$reg_wout_interactions_inference_lst$CI_LS)
  ,OLS_int = unlist(post_matching_analysis_lst[[1]]$reg_with_interactions_inference_lst$CI_LS)
  ,WLS = unlist(post_matching_analysis_lst[[2]]$reg_wout_interactions_inference_lst$CI_LS)
  ,WLS_int = unlist(post_matching_analysis_lst[[2]]$reg_with_interactions_inference_lst$CI_LS))
  
  
  # add the current iteration (i.e. in param_n_sim) to the big matrices that contain info re all iterations
  estimators = c(SACE=SACE, naive_estimators, DL_est, matching_estimators, regression_matching_estimators)
  matching_estimators_mat[i,] = estimators
  colnames(matching_estimators_mat) = names(estimators)
  
  DL_mis = -101 # DL estimators do not include SE and CI
  SE = c(SACE=SACE, naive_estimators_SE, c(DL_est=DL_mis, DL_MA_est=DL_mis), matching_estimators_SE, regression_matching_estimators_SE)
  matching_estimators_SE_mat[i,] = SE
  colnames(matching_estimators_SE_mat) = names(SE)
  
  # CI of matching crude, BC and regression estimators after matching
  CI_mat[i,] = c(SACE=SACE, unlist(naive_estimators_CI), c(DL_est=DL_mis, DL_MA_est=DL_mis), matching_estimators_CI, regression_matching_estimators_CI)
  colnames(CI_mat) = c("SACE", names(unlist(naive_estimators_CI)), names(DL_est),
                       names(matching_estimators_CI), names(regression_matching_estimators_CI))
  
  # balance (mean_by_subset) of x_obs (original covariates)
  balance_wout_rep_lst[[i]] = matching_datasets_lst[[1]]$balance_all_measures$mean_by_subset_maha_cal
  balance_with_rep_lst[[i]] = matching_datasets_lst[[2]]$balance_all_measures$mean_by_subset_maha_cal
  
  end_time1 <- Sys.time()
  print(paste0("one iteration lasts ", difftime(end_time1, start_time1)))
  print(difftime(end_time1, start_time1))
} # out of for loop for all the iterations (param_n_sim iterations in total)

#TODO summaries of all iterations of one scenario (according to mat_gamma and its row, k) from the simulations in summaries_newWF.R
#source("Simulations_studies/sim_DGM_and_simulations/summaries_newWF.R")



