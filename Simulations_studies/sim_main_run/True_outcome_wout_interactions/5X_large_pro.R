args = commandArgs(trailingOnly = TRUE)
job_id = NULL
if(length(args)>0){
  job_id = as.numeric(args[1])
}else{
  stop(' ERROR no job id given!!!')
}

library(data.table); library(plyr); library(dplyr); library(rlang); library(rlist)
library(nnet); library(locfit)
library(Matching); library(sandwich); library(clubSandwich); library(lmtest)
library(mgsub)

########################################################################
# source for Simulations_studies
main_path = paste0("/a/home/cc/stud_math/tamirzehavi/MatchingSACE/Simulation_studies/")
path_code = paste0(main_path, "Code/") 
source(paste0(path_code, "Functions/simulation_run_CPSR.R"))
source(paste0(path_code, "Functions/EM_2log_CPSR.R"))
source(paste0(path_code, "Functions/matching_multiple.R"))
source(paste0(path_code, "Functions/sim_regression_estimators.R"))
source(paste0(path_code, "Functions/table_design_multiple_func.R"))
source(paste0(path_code, "Functions/coverage_naive_est.R"))
#############################################################################################

#############################################################################################
# treatment probability
prob_A = 0.5

# parameters for simulating X
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 6; cont_x = 5; categ_x = 0
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)

# misspec parameters (for PS model and Y model:
# misspec_PS: 0 <- NO, 1:only PS model, 2: PS model (possibly also Y)
misspec_PS = 0 # 0: no misspec of PS model # 2: PS functional form misspecification
funcform_factor_sqr=-3; funcform_factor_log=3
misspec_outcome = FALSE

# EM convergence parameters
iterations = 200; epsilon_EM = 10^-6
#############################################################################################

###############################################################################################
# beta ####
# wout interactions between A and X:
betas_GPI = as.matrix(rbind(c(22,3,4,5,1,3), c(20,3,4,5,1,3))) # cont_x=5
rownames(betas_GPI) = c("beta_treatment", "beta_control")
###############################################################################################

###############################################################################################
# correlation structure between PO'
var_GPI = as.matrix(rbind(1, 1))
rownames(var_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4 
###############################################################################################

###############################################################################################
# Large pi pro 5X ####
mat_gamma = matrix(c(
  c(-0.07, rep(0.03, dim_x-1)), c(0.41, rep(0.2, dim_x-1)) 
  ,c(0.39, rep(0.33, dim_x-1)), c(0.85, rep(0.85, dim_x-1))) ,nrow = 2, byrow = T) 
# assign 0's to gamma_ns and add coefficients names ####
gamma_ns = rep(0, dim_x)
colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("ah", "pro"), each = dim_x) )

# mat_gamma[,c(1,2,dim_x+1,dim_x+2)]
# extract_pis_lst = extract_pis_from_scenarios(nn=1000000, xi=xi, misspec_PS=0); mat_pis_per_gamma = extract_pis_lst$mat_pis
# mat_pis_per_gamma
###############################################################################################

param_n = 2000; param_n_sim = 1000 # param_n = 2000; param_n_sim = 1000
caliper = 0.25; match_on = "O11_posterior_ratio" 
mu_x_fixed = FALSE # mat_x_as; x_as = mat_x_as[1,]


if(job_id >=0 & job_id <=3){
  set.seed(101)
  xi_values = c(0, 0.05, 0.1, 0.2)
  # CPSR parameter 
  xi = xi_values[job_id+1]
  xi_est = xi
  print(paste0('xi_',xi))
  
  param_measures = c("mean","med","sd","MSE"); num_of_param_measures_per_param_set = length(param_measures)
  list_all_mat_SACE_estimators <- list_all_WLS_NOint_regression_estimators <- list_all_WLS_YESint_regression_estimators <-
    list_all_OLS_NOint_regression_estimators <- list_all_OLS_YESint_regression_estimators <- list_all_CI <- 
    list_all_EM_coeffs <- list_all_means_by_subset  <- list_all_EM_not_conv <- list_all_BCclpr <- list()
  
  ########################################################################
  # run over different values of gamma's: 1:nrow(mat_gamma) ####
  # param_n_sim * time per run * nrow(mat_gamma)
  for ( k in c(1 : nrow(mat_gamma)) ){
    print(paste0("in the outer for loop ", k))
    gamma_ah=as.numeric(mat_gamma[k, c(1:dim_x)])
    gamma_pro=as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
    gamma_ns=gamma_ns
    start_time <- Sys.time()
    EM_and_matching = simulate_data_run_EM_and_match(only_EM_bool=FALSE, return_EM_PS=FALSE, index_set_of_params=k,
                                                     gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, xi=xi, xi_est=xi_est, two_log_models=TRUE, two_log_est_EM=FALSE,
                                                     misspec_PS=misspec_PS, misspec_outcome=misspec_outcome, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, 
                                                     param_n=param_n, param_n_sim=param_n_sim, iterations=iterations, epsilon_EM=epsilon_EM, caliper=caliper,
                                                     match_on=match_on, mu_x_fixed=mu_x_fixed, x_as=NULL)
    
    mat_SACE_estimators = EM_and_matching[["mat_param_estimators"]]
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
  ties_setA = sum(abs(list.rbind(list_all_BCclpr[[1]])$trt_added_by_ties))
  ties_setB = sum(abs(list.rbind(list_all_BCclpr[[2]])$trt_added_by_ties))
  ties = c(ties_setA=ties_setA, ties_setB=ties_setB)
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
  mat_all_estimators = data.frame(subset(mat_all_estimators, select = grep("gamma", colnames(mat_all_estimators))),
                                  subset(mat_all_estimators, select = -grep("gamma", colnames(mat_all_estimators))))
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
  
  
  # summary of EM estimators for the gamms's- the PS coefficient- logistic reg of stratum on X ####
  ########################################################################
  summary_EM_coeffs = list.rbind(lapply(list_all_EM_coeffs, function(x) x[c((param_n_sim + 1):(param_n_sim + 5))]))
  rownames(summary_EM_coeffs) = paste0(rep(LETTERS[1:nrow(mat_gamma)],each = ncol(mat_gamma)), "_",
                                       rep(rownames(summary_EM_coeffs)[1:ncol(mat_gamma)],times = nrow(mat_gamma)))
  ########################################################################
  
  # covariates means by A and S, before and after matching (+ mean of as) ####
  ########################################################################
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
  #final_tables_general$'_S1'
  ########################################################################
  
  ########################################################################
  # save ####
  Large_pi_pro = TRUE
  path_data = paste0(main_path, "Data/") 
  path = paste0(path_data, "True_outcome_with_interactions/",
                ifelse(misspec_outcome==0, "Correct_spec_outcome/", "Mis_spec_outcome/"),
                ifelse(misspec_PS==0, "Correct_spec_PS/", "Mis_spec_PS/"), paste0(cont_x, "X/"),
                ifelse(Large_pi_pro, "Large_pi_pro/", "Low_pi_pro/"), "xi = ", xi, "/")
  
  save(final_tables_general, file = paste0(path, 'final_tables_general_',job_id,'.Rdata'))
  save(final_tables_crude, file = paste0(path, 'final_tables_crude_',job_id,'.Rdata'))
  save(list_all_CI_temp, file = paste0(path, 'list_all_CI_temp_',job_id,'.Rdata'))
  save(mat_all_means_by_subset, file = paste0(path, 'mat_all_means_by_subset_',job_id,'.Rdata'))
  save(pis, file = paste0(path, 'pis_',job_id,'.Rdata'))
  save(ties, file = paste0(path, 'ties_',job_id,'.Rdata'))
  ########################################################################
  
}else{
  stop(' Error job id not in 0-3')
}