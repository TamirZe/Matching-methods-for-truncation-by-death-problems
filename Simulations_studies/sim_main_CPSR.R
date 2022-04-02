library(rlist); library(locfit); library(nnet); library(xtable); library(rlang);library(glue)
# library(multilevelMatching); library(PerformanceAnalytics); library(lmtest); library(caret);
library(matrixStats); library(data.table); library(dplyr); library(plyr); library(reshape); library(MASS); library(Hmisc); 
library(ggplot2); library(rockchalk); library(stats); library(rlist); library(mgsub); library(reshape2); library(gridExtra)
library(optmatch); library(DOS); library(Matching); library(sandwich); library(rmutil); library(clubSandwich); library(tableone)
library(sandwich); library(lmtest); library(rmutil); library(splitstackshape); library(PerformanceAnalytics)

########################################################################
# source for Simulations_studies
setwd("~/A matching framework for truncation by death problems")
source("Simulations_studies/sim_DGM_and_simulations/simulation_run_CPSR.R")
source("Simulations_studies/sim_matching_procedure/matching_PS_multiple.R")
source("Simulations_studies/sim_post_matching_analysis/sim_regression_estimators.R")
source("Simulations_studies/sim_tables_and_figures/table_design_multiple_func.R")
source("Simulations_studies/sim_tables_and_figures/coverage_naive_est.R")
source("Ding_Lu/PS_M_weighting.R")
source("Ding_Lu/PS_M_weighting_SA.R")
#############################################################################################

#############################################################################################
# treatment probability
prob_A = 0.5

# parameters for simulating X
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 6; cont_x = 5; categ_x = 0; vec_p_categ = rep(0.5, categ_x); dim_x_misspec = 2
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)

# misspec parameters (for PS model and Y model:
# misspec_PS: 0 <- NO, 1:only PS model, 2: PS model (possibly also Y)
funcform_mis_out = FALSE; match_and_reg_watch_true_X = FALSE
funcform_factor_sqr=-3; funcform_factor_log=3
mean_x_misspec = rep(0.5, dim_x_misspec)
misspec_PS = 0 # 0: no misspec of PS model # 2: PS functional form misspecification

# CPSR parameter
xi = 0

# EM convergence parameters
iterations = 200; epsilon_EM = 10^-6
#############################################################################################

##########################################################
# beta ####
# with interactions between A and X:
betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0))) # cont_x=3
betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3))) # cont_x=5
betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(3,3,0,1,3),2)))) # cont_x=10 

# without interactions between A and X: "simple effect" is 2
betas_GPI = as.matrix(rbind(c(22,3,4,5), c(20,3,4,5))) # cont_x=3
betas_GPI = as.matrix(rbind(c(22,3,4,5,1,3), c(20,3,4,5,1,3))) # cont_x=5
betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(5,2,1,3,5),2)))) # cont_x=10

rownames(betas_GPI) = c("beta_treatment", "beta_control")
###############################################################################################

##########################################################
# correlation structure between PO'
var_GPI = as.matrix(rbind(1, 1))
rownames(var_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4 #0.4 #1
##########################################################

# mat_gamma ####
#############################################################################################
# values for gamma's 3X ####
# Large pi pro 3X
mat_gamma = matrix(c(
  c(-0.1, rep(0.27, dim_x-1)), c(-0.52, rep(-0.6, dim_x-1))
  ,c(0.6, rep(1.325, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))) ,nrow = 2, byrow = T) 
# small pi pro 3X 
mat_gamma = matrix(c(
  c(-0.07, rep(1.24, dim_x-1)), c(1.25, rep(0.24, dim_x-1)) 
  ,c(0.84, rep(1.2, dim_x-1)), c(0.32, rep(-0.2, dim_x-1))) ,nrow = 2, byrow = T)

# values for gamma's 5x ####
# large pi pro 5X 
mat_gamma = matrix(c(
  c(-0.05, rep(0.16, dim_x-1)), c(-0.4, rep(-0.25, dim_x-1)) 
  ,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))) ,nrow = 2, byrow = T)
# small pi pro 5X 
mat_gamma = matrix(c(
  c(0.45, rep(0.75, dim_x-1)), c(0.62, rep(0.6, dim_x-1))
  ,c(-0.12, rep(1.25, dim_x-1)), c(0.45, rep(-0.2, dim_x-1))) ,nrow = 2, byrow = T)

# values for gamma's 10X #### 
# large pi pro 10X:
mat_gamma = matrix(c(
  c(-0.954, rep(0.25, dim_x-1)), c(-0.31, rep(-0.16, dim_x-1))
  ,c(-1.01, rep(0.8, dim_x-1)), c(-1.3, rep(0.45, dim_x-1))) ,nrow = 2, byrow = T)

# small pi pro 10X 
mat_gamma = matrix(c(
  c(0.024, rep(0.41, dim_x-1)), c(0.26, rep(0.32, dim_x-1))
  ,c(-0.45, rep(0.61, dim_x-1)), c(0.4, rep(-0.02, dim_x-1))) ,nrow = 2, byrow = T)

# assign 0's to gamma_pro and add coefficients names ####
gamma_pro = rep(0, dim_x)
colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("as", "ns"), each = dim_x) )
#############################################################################################

#############################################################################################
# extract strata proportion and mean covariates ####
extract_pis_from_scenarios = function(nn=250000){
  big_lst = list(); mat_x_as <- mat_pis <- mat_x_by_g_A <- NULL
  for( k in c(1 : nrow(mat_gamma)) ){
    gamma_as = as.numeric(mat_gamma[k, c(1:dim_x)])
    gamma_ns =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
    lst_mean_x_and_pi = simulate_data_function(gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, xi,
                                               param_n=nn, misspec_PS=2, funcform_mis_out=FALSE, 
                                               funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, only_mean_x_bool=TRUE)
    big_lst[[k]] = lst_mean_x_and_pi
    mat_x_as = rbind(mat_x_as, lst_mean_x_and_pi$x_as)
    mat_pis = rbind(mat_pis, lst_mean_x_and_pi$pi)
    mat_x_by_g_A = rbind(mat_x_by_g_A, data.frame(Scenar = k, lst_mean_x_and_pi$mean_by_A_g))
  }
  mat_pis = data.frame(pi_as=mat_pis[,1], pi_pro=mat_pis[,3], pi_ns=mat_pis[,2])
  round(mat_pis,3)
  return(list(mat_x_by_g_A=mat_x_by_g_A, big_lst=big_lst, mat_x_as=mat_x_as))
}
mat_gamma[,c(1,2,7,8)]

big_lst = list(); big_mat_x_by_g_A=NULL
for(i in 1:100){
  print(i)
  big_lst[[i]] = extract_pis_from_scenarios(nn=2000)
  big_mat_x_by_g_A = rbind(big_mat_x_by_g_A, big_lst[[i]]$mat_x_by_g_A)
}
big_mat_x_by_g_A = subset(big_mat_x_by_g_A, select = c(Scenar,A,g, grep("X", colnames(big_mat_x_by_g_A))))
big_mat_x_by_g_A = data.table(big_mat_x_by_g_A)[, lapply(.SD, mean), by=c("Scenar", "A", "g")] %>% arrange(Scenar, g, A)
#############################################################################################

#############################################################################################
# extract mean covariates of as ####
mat_x_as = NULL
for( k in c(1 : nrow(mat_gamma)) ){
  gamma_as = as.numeric(mat_gamma[k, c(1:dim_x)]); gamma_ns =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
  lst_mean_x_and_pi = simulate_data_function(gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, xi=xi,
                                             param_n=250000, misspec_PS=0, funcform_mis_out=FALSE, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, only_mean_x_bool=TRUE)
  mat_x_as = rbind(mat_x_as, lst_mean_x_and_pi$x_as)
}
#############################################################################################

param_n = 2000; param_n_sim = 3 # param_n = 2000; param_n_sim = 1000
caliper = 0.25; match_on = "O11_posterior_ratio" 
mu_x_fixed = FALSE; mat_x_as; x_as = mat_x_as[1,]

param_measures = c("mean","med","sd","MSE"); num_of_param_measures_per_param_set = length(param_measures)
list_all_mat_SACE_estimators <- list_all_WLS_NOint_regression_estimators <- list_all_WLS_YESint_regression_estimators <-
  list_all_OLS_NOint_regression_estimators <- list_all_OLS_YESint_regression_estimators <- list_all_CI <- 
  list_all_EM_coeffs <- list_all_excluded_included_matching <-
  list_all_repeated_as_and_pro <- list_all_diff_distance_aspr_asas <- list_all_matched_units <-
  list_all_std_mean_diff <- list_all_means_by_subset  <- list_all_EM_not_conv <- list_all_BCclpr <- list()

# run over different values of gamma's: 1:nrow(mat_gamma)
# param_n_sim * time per run * nrow(mat_gamma)
for ( k in c(1 : nrow(mat_gamma)) ){
  print(paste0("in the outer for loop ", k))
  gamma_as=as.numeric(mat_gamma[k, c(1:dim_x)])
  gamma_ns=as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
  gamma_pro=gamma_pro
  start_time <- Sys.time()
  
  EM_and_matching = simulate_data_run_EM_and_match(return_EM_PS=FALSE, index_set_of_params=k,
                                                   gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, xi=xi,
                                                   misspec_PS=misspec_PS, funcform_mis_out=FALSE, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, 
                                                   match_and_reg_watch_true_X=FALSE, param_n=param_n, param_n_sim=param_n_sim,
                                                   iterations=iterations, epsilon_EM=epsilon_EM, caliper=caliper,
                                                   match_on=match_on, mu_x_fixed=mu_x_fixed, x_as=mat_x_as[k,])
  
  mat_SACE_estimators = EM_and_matching[[1]]
  df_parameters = matrix(rep(as.numeric(mat_gamma[k,])
                             , each = nrow(mat_SACE_estimators))
                         , nrow = nrow(mat_SACE_estimators))
  colnames(df_parameters) = colnames(mat_gamma)
  mat_SACE_estimators = data.frame(mat_SACE_estimators, df_parameters)
  
  list_all_mat_SACE_estimators[[k]] = mat_SACE_estimators
  list_all_WLS_NOint_regression_estimators[[k]] = EM_and_matching[["WLS_NOint_mat_reg_estimators"]]
  list_all_WLS_YESint_regression_estimators[[k]] = EM_and_matching[["WLS_YESint_mat_reg_estimators"]]
  list_all_OLS_NOint_regression_estimators[[k]] = EM_and_matching[["OLS_NOint_mat_reg_estimators"]]
  list_all_OLS_YESint_regression_estimators[[k]] = EM_and_matching[["OLS_YESint_mat_reg_estimators"]]
  list_all_CI[[k]] = EM_and_matching[["CI_mat"]]
  list_all_EM_coeffs[[k]] = EM_and_matching[["coeffs_df"]]
  list_all_excluded_included_matching[[k]] = EM_and_matching[["mat_excluded_included_matching"]]
  list_all_repeated_as_and_pro[[k]] = EM_and_matching[["mean_list_repeated_as_and_pro"]]
  rownames(list_all_repeated_as_and_pro[[k]]) = paste0("s", k, rownames(list_all_repeated_as_and_pro[[k]]))
  list_all_matched_units[[k]] = EM_and_matching[["mean_list_matched_units"]]
  rownames(list_all_matched_units[[k]]) = paste0("s", k, rownames(list_all_matched_units[[k]]))
  
  list_all_diff_distance_aspr_asas[[k]] = EM_and_matching[["mat_diff_distance_aspr_asas"]]
  rownames(list_all_diff_distance_aspr_asas[[k]]) = paste0("s", k, rownames(list_all_diff_distance_aspr_asas[[k]]))
  list_all_std_mean_diff[[k]] = EM_and_matching[["mean_list_std_mean_diff"]]
  rownames(list_all_std_mean_diff[[k]]) = paste0("s", k, rownames(list_all_std_mean_diff[[k]]))
  list_all_means_by_subset[[k]] = EM_and_matching[["mean_list_means_by_subset"]]
  rownames(list_all_means_by_subset[[k]]) = paste0("s", k, rownames(list_all_means_by_subset[[k]]))
  list_all_EM_not_conv[[k]] = EM_and_matching[["list_EM_not_conv"]]
  if(! is_empty(list_all_EM_not_conv[[k]])){names(list_all_EM_not_conv[[k]]) = paste0("s", k, names(list_all_EM_not_conv[[k]]))}
  list_all_BCclpr[[k]] = EM_and_matching[["list_BCclpr"]]
  
  end_time <- Sys.time()
  print(paste0("in the end of outer for loop ", k, ", ", difftime(end_time, start_time)))
}
########################################################################

########################################################################
# check ties in BC with caliper ####
sum(abs(list.rbind(list_all_BCclpr[[1]])$trt_added_by_ties))
sum(abs(list.rbind(list_all_BCclpr[[2]])$trt_added_by_ties))
########################################################################

########################################################################
# function that changes rownames to LETTERS per each parameter (gammas) set ####
change_rownames_to_LETTERS_by_param_set = function(df, mat_params=mat_gamma,
                                                   num_of_param_measures=num_of_param_measures_per_param_set){
  rownames(df) = paste0(rep(LETTERS[1:nrow(mat_params)],each = num_of_param_measures), "_",
                        rep(rownames(df)[1:num_of_param_measures], times = nrow(mat_params)))
  return(df)
}
########################################################################

########################################################################
# summary of SACE estimators ####
mat_all_estimators = list.rbind(lapply(list_all_mat_SACE_estimators, tail, length(param_measures)))
#num_of_param_measures_per_param_set = nrow(mat_all_estimators) / nrow(mat_gamma) # = length(param_measures)
mat_all_estimators = data.frame(subset(mat_all_estimators,
                                       select = grep("gamma", colnames(mat_all_estimators))),
                                subset(mat_all_estimators,
                                       select = -grep("gamma", colnames(mat_all_estimators))))
mat_all_estimators = change_rownames_to_LETTERS_by_param_set(mat_all_estimators)
pi_from_mat_all_estimators = subset(mat_all_estimators, select = grep("pi_", colnames(mat_all_estimators))) %>% round(3)

# list_all_regression_estimators
WLS_NOint_mat_regression_estimators = list.rbind(lapply(list_all_WLS_NOint_regression_estimators, tail, length(param_measures)))
WLS_NOint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(WLS_NOint_mat_regression_estimators)
WLS_YESint_mat_regression_estimators = list.rbind(lapply(list_all_WLS_YESint_regression_estimators, tail, length(param_measures)))
WLS_YESint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(WLS_YESint_mat_regression_estimators)
OLS_NOint_mat_regression_estimators = list.rbind(lapply(list_all_OLS_NOint_regression_estimators, tail, length(param_measures)))
OLS_NOint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(OLS_NOint_mat_regression_estimators)
OLS_YESint_mat_regression_estimators = list.rbind(lapply(list_all_OLS_YESint_regression_estimators, tail, length(param_measures)))
OLS_YESint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(OLS_YESint_mat_regression_estimators)


# summary of EM estimators for the gamms's- the PS coefficient- logistic reg of stratum on X
summary_EM_coeffs = list.rbind(lapply(list_all_EM_coeffs, function(x) x[c((param_n_sim + 1):(param_n_sim + 5))]))
rownames(summary_EM_coeffs) = paste0(rep(LETTERS[1:nrow(mat_gamma)],each = ncol(mat_gamma)), "_",
                                     rep(rownames(summary_EM_coeffs)[1:ncol(mat_gamma)],times = nrow(mat_gamma)))

# covariates means by A and S, before and after matching (+ mean of as)
mat_all_means_by_subset = NULL
first_3_rows = c("mean_as", "mean_A0_S1", "mean_A1_S1_as")
for(i in c(1:length(list_all_means_by_subset))){
  first_3_rows_ind = grep(paste0(first_3_rows, collapse = "|"), rownames(list_all_means_by_subset[[i]]))
  mat_all_means_by_subset = rbind(mat_all_means_by_subset,
                                  list_all_means_by_subset[[i]][-first_3_rows_ind[-c(1:3)], ])
}
########################################################################

########################################################################
# initial summaries
means_by_subset_sum = mat_all_means_by_subset[grep("S1|mean_as", 
                                                   rownames(mat_all_means_by_subset), ignore.case = F) , ] %>% round(3)
means_by_subset_sum = means_by_subset_sum[-grep("_approx", rownames(means_by_subset_sum)),]
#pi_from_mat_all_estimators VS 
pis = mat_all_estimators[grep("_mean",rownames(mat_all_estimators)),grep("pi",colnames(mat_all_estimators))] %>% round(3)
pis = data.frame(pi_as=pis$pi_as, pi_pro=pis$pi_pro, pi_ns=pis$pi_ns)
########################################################################

########################################################################
# TABLE DESIGN ####
list_all_CI_temp = list_all_CI # list_all_CI from main
data_set_names_vec = c("_all", "_wout_O_0_0", "_S1")

lst_final_tables_ALL_est_ALL_dataset = TABLES_before_coverage_and_wout_naive_and_DING(data_set_names_vec,
                                                                                      estimators_names=c("PS Crude", "mahal Crude", "Crude", "HL", "BC", "BC caliper", "BC inter", "BC caliper inter"))
lst_final_tables_ALL_est_ALL_dataset_PLUS_naives = 
  TABLES_add_coverage(data_set_names_vec, lst_final_tables_ALL_est_ALL_dataset, list_all_CI)
lst_final_tables_ALL_est_ALL_dataset = lst_final_tables_ALL_est_ALL_dataset_PLUS_naives$lst_final_tables_ALL_est_ALL_dataset
naives_before_matching_coverage = lst_final_tables_ALL_est_ALL_dataset_PLUS_naives$naives_before_matching_coverage
final_tables = TABLES_add_naive_and_ding(data_set_names_vec, lst_final_tables_ALL_est_ALL_dataset, naives_before_matching_coverage)
final_tables_general = adjustments_for_final_tables(final_tables)
final_tables_crude = adjustments_for_final_tables_crude_est(final_tables)
########################################################################
