# installing/loading the package:
if(!require(installr)) { install.packages("installr"); require(installr)} #load / install+load installr
# step by step functions:
check.for.updates.R() # tells you if there is a new version of R or not.
install.R() # download and run the latest R installer
copy.packages.between.libraries() # copy your packages to the newest R installation from the one version before it (if ask=T, it will ask you between which two versions to perform the copying)


devtools::install_github("shuyang1987/multilevelMatching")
install.packages("rlist"); install.packages("locfit"); install.packages("plyr")
install.packages("nnet"); install.packages("xtable"); install.packages("doParallel")
install.packages("matrixStats"); install.packages("data.table"); install.packages("rlang")
install.packages("dplyr"); install.packages("reshape"); install.packages("MASS"); install.packages("mgsub")
install.packages("ggplot2"); install.packages("rockchalk"); install.packages("nnet")
installed.packages("stats"); install.packages("optmatch"); install.packages("DOS"); install.packages("PerformanceAnalytics")
install.packages("Matching"); install.packages("sandwich"); install.packages("rmutil")
install.packages("sandwich"); install.packages("rmutil"); install.packages("caret")

library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); library(rlang)
# library(multilevelMatching); library(PerformanceAnalytics); library(lmtest)
library(matrixStats); library(data.table); library(dplyr); library(reshape); library(MASS); library(Hmisc)
library(ggplot2); library(rockchalk); library(nnet); library(stats); library(rlist); library(mgsub)
library(optmatch); library(DOS); library(Matching); library(sandwich); library(rmutil); library(clubSandwich)
library(sandwich); library(rmutil);  library(caret); library(splitstackshape); library(MatchIt); library(PerformanceAnalytics)

source("matching_scripts/matching_PS_basic.R")
source("matching_scripts/matching_PS_multiple.R")
source("matching_scripts/more/matching_matchit.R")
source("matching_scripts/more/matching_PS_caliper_pairmatch.R")
#source("matching_scripts/matching_from_real_data.R")
#source("my_EM_V2.R")
source("EM_V3_eps_stop.R"); source("OLS_WLS_estimator.R")
#source("sim1.R")
source("simulations_scripts/sim2.R")
source("simulations_scripts/sim3_true_differs_same_estim.R")
source("DING_model_assisted_estimator.R")
source("Extra code/TABLES/table_design_multiple_func.R")
#source("Extra code/TABLES/table_design_func.R"); source("Extra code/TABLES/table_design_w_BCest_func.R")


library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)
#foreach(j = seq_along(numeric(length(c(1:n_sim)))), .combine=c) %dopar%{
stopCluster(cl)

#all.equal(pm1, pm2, check.attributes = FALSE)
#identical()

#############################################################################################
# # simulate dependency between errors
# rho <- cbind(c(1, 1), c(1, 1))
# sds = 1
# x <- mvrnorm(10, mu=mean_x, Sigma = diag(sd_x, cont_x))
# var(x)
#############################################################################################


##########################################################
# parameters for the linear regression Y(a,g) on X when GPI is viloated
# shape: beta_ # of coefficient_stratum_treatment
# beta0_as1 = 1; beta1_as1 = 4; sigma_square_as1 = 1; beta_as1 = c(beta0_as1, beta1_as1)
# beta0_as0 = 0; beta1_as0 = 2; sigma_square_as0 = 1; beta_as0 = c(beta0_as0, beta1_as0)
# beta0_pro1 = 1; beta1_pro1 = 2.5; sigma_square_pro1 = 1; beta_pro1 = c(beta0_pro1, beta1_pro1)
# beta0_har0 = 0.5; beta1_har0 = 3; sigma_square_har0 = 1; beta_har0 = c(beta0_har0, beta1_har0)
# betas = as.matrix(rbind(beta_as1, beta_as0, beta_pro1, beta0_har0))
# sigma_square = as.matrix(rbind(sigma_square_as1, sigma_square_as0, sigma_square_pro1, sigma_square_har0))
# ##########################################################

#############################################################################################
# parameters for the functions
prob_A = 0.5; iterations = 12; epsilon_EM = 0.001 # epsilon_EM is for the EM convergence
monotonicity_assumption = "mono"; PI_assum = "strong"
# parameters for simulating x
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 6; cont_x = 5; categ_x = 0; vec_p_categ = rep(0.5, categ_x); dim_x_misspec = 2
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)

# misspec PARAMETERS (for PS model and Y model:
# TODO misspec_PS: 0 <- NO, 1:only PS model, 2: PS model, and matching (mahalanobis, eucl...) and regression
#misspec_PS = 2
misspec_outcome_funcform = FALSE; match_and_reg_watch_true_X = FALSE
U_factor=1.5; funcform_factor_sqr=1; funcform_factor_log=2.5
mean_x_misspec = rep(0.5, dim_x_misspec)


# for covariance between X 
# TODO Sigma here is vcov not sd cor. need to change in the PO sim in sim1 rows ~80
'''corr_x = 0.3
sqrt(var_x)
cov_x <- diag(sqrt(5),cont_x) + corr_x - diag(1,cont_x) * corr_x # cov_x is vcov and corr
X = mvrnorm(param_n, mu = mean_x, Sigma =  cov_x)
apply(X, 2, mean); apply(X, 2, var); cov(X); cor(X)'''

# sensitivity paramters
epsilon_1_GPI = 1
#############################################################################################

##########################################################
# betas under GPI
# TODO YES interactions between A and X:
# c(6,5,2,1)
# betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0)))
betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3)))

# TODO NO interactions between A and X: "simple effect" is 2
#betas_GPI = as.matrix(rbind(c(22,3,2,1), c(20,3,2,1)))
# betas_GPI = as.matrix(rbind(c(22,3,2,1,1,3), c(20,3,2,1,1,3)))

rownames(betas_GPI) = c("beta_treatment", "beta_control")

var_GPI = as.matrix(rbind(1, 1))
rownames(var_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4
#sds_GPI_PO = 1
##########################################################

# add the harmed stratum parameters if needed
# if (monotonicity_assumption == "nothing"){
#   betas = rbind(betas, beta_har0)
#   sigma_square = c(sigma_square, sigma_square_har0)
# }
#############################################################################################

#############################################################################################
# different values of gamma's

# mat_full_gamma = expand.grid(gamma0_as = seq(-0.25,1,0.25), gamma1_as = seq(-0.5,1,0.25),
#                         gamma0_ns = seq(-1,1,0.25), gamma1_ns = seq(-1,1,0.25))

# mat_full_gamma = expand.grid.df( data.frame(gamma0_as = seq(0, 1, 0.5), gamma1_as = seq(0, 1, 0.5)),
#                                  data.frame(gamma0_ns =  -seq(0, 1, 0.5),  gamma1_ns = -seq(0, 1, 0.5)) )
# colnames(mat_gamma) = colnames(mat_full_gamma)
# 
# mat_gamma = mat_full_gamma[-c(2,5),]
# mat_gamma = rbind(mat_full_gamma[-c(2,5),], c(-0.5,-0.5,0.5,0.5))
# mat_gamma = mat_full_gamma[c(1, 5, 9) , ]

###########################################################################


################ new values for gamma's ################ 
gamma_pro = rep(0, dim_x)
mat_gamma = matrix(c(rep(0.25, dim_x),rep(-0.25, dim_x)
                     ,rep(0.5, dim_x), rep(0, dim_x)
                     ,rep(0.1, dim_x),rep(1.25, dim_x))
                   ,nrow = 3, byrow = T)

mat_gamma = matrix(c(rep(0.05, dim_x),rep(-0.05, dim_x)
                     ,rep(-0.25, dim_x), rep(0.25, dim_x)
                     ,rep(1, dim_x),rep(0.25, dim_x))
                   ,nrow = 3, byrow = T)

# ~0.14, 0.39, 0.6, <0.81 befroe PS misspec
# take this 22.10!
mat_gamma = matrix(c(
  rep(-0.05, dim_x), rep(0, dim_x)
  ,rep(0.05, dim_x),rep(-0.05, dim_x)
  ,rep(0.25, dim_x),rep(-0.25, dim_x)
  #,rep(1, dim_x),rep(0.25, dim_x)
)
,nrow = 3, byrow = T)

mat_gamma = matrix(c(
  rep(0, dim_x), rep(0.046, dim_x)
  ,rep(0.024, dim_x),rep(-0.06, dim_x)
  ,rep(0.0746, dim_x),rep(-0.25, dim_x)
  #,rep(1, dim_x),rep(0.25, dim_x)
)
,nrow = 3, byrow = T)

colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("as", "ns"), each = dim_x) )

# if monotonicity != TRUE, times = 3
mat_gamma = mat_gamma[c(1,2), ]
#mat_gamma = matrix(mat_gamma[c(1), ], nrow = 1, ncol = 8)

# TODO calculate with function:
mat_x_as = NULL; pis_mat = NULL
for( k in c(1 : nrow(mat_gamma)) ){
  gamma_as = as.numeric(mat_gamma[k, c(1:dim_x)])
  gamma_ns =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
  lst_mean_x_and_pi = simulate_data_function(gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, param_n=10000, 
                                             misspec_PS=2, misspec_outcome_funcform=FALSE,
                                             U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                                             epsilon_1_GPI = 1, only_mean_x_bool=TRUE)
  mat_x_as = rbind(mat_x_as, lst_mean_x_and_pi$x_as)
  pis_mat = rbind(pis_mat, lst_mean_x_and_pi$pi)
}

param_n = 2000; param_n_sim = 200
#param_n = 1000; param_n_sim = 3
caliper = 0.05; pairmatch_bool = FALSE
match_on = "EMest_p_as" # NULL # Feller and Mialli: "EMest_p_as" # Ding Lu appendix: "O11_posterior_ratio"
param_measures = c("mean","med","sd","MSE"); num_of_param_measures_per_param_set = length(param_measures)
mu_x_fixed = FALSE; mat_x_as; x_as = mat_x_as[1,]
sim=3; misspec_PS_if_sim_not3 = 1


if(sim==3){
  lstlist_all_mat_SACE_estimators <- lstlist_all_WLS_NOint_regression_estimators <- lstlist_all_WLS_YESint_regression_estimators <-
    lstlist_all_OLS_NOint_regression_estimators <- lstlist_all_OLS_YESint_regression_estimators <- lstlist_all_CI <- 
    lstlist_all_EM_coeffs <- lstlist_all_excluded_included_matching <-
    lstlist_all_repeated_as_and_pro <- lstlist_all_diff_distance_aspr_asas <- lstlist_all_matched_units <-
    lstlist_all_std_mean_diff <- lstlist_all_means_by_subset <- list(list(),list(),list())
  set.seed(101); seed_vec = sample(c(1:10000), param_n_sim, replace = FALSE)
}else{  
  lstlist_all_mat_SACE_estimators <- lstlist_all_WLS_NOint_regression_estimators <- lstlist_all_WLS_YESint_regression_estimators <-
    lstlist_all_OLS_NOint_regression_estimators <- lstlist_all_OLS_YESint_regression_estimators <- lstlist_all_CI <- 
    lstlist_all_EM_coeffs <- lstlist_all_excluded_included_matching <-
    lstlist_all_repeated_as_and_pro <- lstlist_all_diff_distance_aspr_asas <- lstlist_all_matched_units <-
    lstlist_all_std_mean_diff <- lstlist_all_means_by_subset  <- list()
  # if not sim3, we rutn per each misspec_PS separately
  misspec_PS = misspec_PS_if_sim_not3
  seed_vec = NULL
}


# run over different values of gamma's: 1:nrow(mat_gamma)
# param_n_sim * time per run * nrow(mat_gamma)
for ( k in c(1 : nrow(mat_gamma)) ){
  print(paste0("in the outer for loop ", k))
  gamma_as=as.numeric(mat_gamma[k, c(1:dim_x)])
  gamma_ns=as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
  gamma_pro=gamma_pro
  start_time <- Sys.time()
  # sim1
  if(sim == 1){
    EM_and_matching = simulate_data_run_EM_and_match(index_set_of_params=k, 
                                                     gamma_as, gamma_ns, gamma_pro, 
                                                     misspec_PS, param_n, param_n_sim, iterations = iterations, epsilon_EM = epsilon_EM,
                                                     caliper, epsilon_1_GPI, match_on = match_on,
                                                     mu_x_fixed=mu_x_fixed, x_as=mat_x_as[k,])
    size_EM_and_matching=1
  }
  # sim2
  if(sim == 2){
    EM_and_matching = simulate_data_run_EM_and_match(return_EM_PS = FALSE, index_set_of_params=k,
                                                     gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, misspec_PS=misspec_PS,
                                                     misspec_outcome_funcform=FALSE, U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, 
                                                     match_and_reg_watch_true_X=FALSE, param_n=param_n, param_n_sim=param_n_sim,
                                                     iterations=iterations, epsilon_EM = epsilon_EM, caliper=caliper, epsilon_1_GPI=epsilon_1_GPI,
                                                     match_on = match_on, mu_x_fixed=mu_x_fixed, x_as=mat_x_as[k,])
    size_EM_and_matching=1
  }
  # EM_and_matching = simulate_data_run_EM_and_match()
  # sim3
  if(sim == 3){
    EM_and_matching = simulate_multi_misspec_by_calculate_EM_forall_together(seed_vec=seed_vec, index_set_of_params=k,
                                                                             gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro,
                                                                             param_n_sim=param_n_sim, first_misspec=0, last_misspec=2,
                                                                             U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, match_on=match_on)
    size_EM_and_matching = length(0:2)      
  }
  
  if(sim == 3){
    for (misspec_ind in c(1:size_EM_and_matching)){
      mat_SACE_estimators = EM_and_matching[[misspec_ind]][[1]]
      # nrow = nrow(mat_gamma) 
      df_parameters = matrix(rep(as.numeric(mat_gamma[k,])
                                 #, each = (param_n_sim + length(param_measures)))
                                 , each = nrow(mat_SACE_estimators))
                             #, nrow = (param_n_sim + length(param_measures))
                             , nrow = nrow(mat_SACE_estimators))
      colnames(df_parameters) = colnames(mat_gamma)
      mat_SACE_estimators = data.frame(mat_SACE_estimators, df_parameters)
      
      lstlist_all_mat_SACE_estimators[[misspec_ind]][[k]] = mat_SACE_estimators
      lstlist_all_WLS_NOint_regression_estimators[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["WLS_NOint_mat_reg_estimators"]]
      lstlist_all_WLS_YESint_regression_estimators[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["WLS_YESint_mat_reg_estimators"]]
      lstlist_all_OLS_NOint_regression_estimators[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["OLS_NOint_mat_reg_estimators"]]
      lstlist_all_OLS_YESint_regression_estimators[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["OLS_YESint_mat_reg_estimators"]]
      lstlist_all_CI[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["CI_mat"]]
      lstlist_all_EM_coeffs[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["coeffs_df"]]
      lstlist_all_excluded_included_matching[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["mat_excluded_included_matching"]]
      lstlist_all_repeated_as_and_pro[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["mean_list_repeated_as_and_pro"]]
      rownames(lstlist_all_repeated_as_and_pro[[misspec_ind]][[k]]) = paste0("s", k, rownames(lstlist_all_repeated_as_and_pro[[misspec_ind]][[k]]))
      lstlist_all_matched_units[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["mean_list_matched_units"]]
      rownames(lstlist_all_matched_units[[misspec_ind]][[k]]) = paste0("s", k, rownames(lstlist_all_matched_units[[misspec_ind]][[k]]))
      
      # TODO in list_all_diff_distance_aspr_asas and list_all_std_mean_diff   
      # when diff_distance_aspr_asas is positive, the matches between as to as are closer in X, 
      #then matches between as to protected, since its the abs diff between as to pro matches - abs diff between as to as matches 
      lstlist_all_diff_distance_aspr_asas[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["mat_diff_distance_aspr_asas"]]
      rownames(lstlist_all_diff_distance_aspr_asas[[misspec_ind]][[k]]) = 
        paste0("s", k, rownames(lstlist_all_diff_distance_aspr_asas[[misspec_ind]][[k]]))
      lstlist_all_std_mean_diff[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["mean_list_std_mean_diff"]]
      rownames(lstlist_all_std_mean_diff[[misspec_ind]][[k]]) = paste0("s", k, rownames(lstlist_all_std_mean_diff[[misspec_ind]][[k]]))
      lstlist_all_means_by_subset[[misspec_ind]][[k]] = EM_and_matching[[misspec_ind]][["mean_list_means_by_subset"]]
      rownames(lstlist_all_means_by_subset[[misspec_ind]][[k]]) = 
        paste0("s", k, rownames(lstlist_all_means_by_subset[[misspec_ind]][[k]]))
      print(paste0("sim is ", sim))
    }
  }else{
    mat_SACE_estimators = EM_and_matching[[1]]
    # nrow = nrow(mat_gamma) 
    df_parameters = matrix(rep(as.numeric(mat_gamma[k,])
                               #, each = (param_n_sim + length(param_measures)))
                               , each = nrow(mat_SACE_estimators))
                           #, nrow = (param_n_sim + length(param_measures))
                           , nrow = nrow(mat_SACE_estimators))
    colnames(df_parameters) = colnames(mat_gamma)
    mat_SACE_estimators = data.frame(mat_SACE_estimators, df_parameters)
    
    lstlist_all_mat_SACE_estimators[[k]] = mat_SACE_estimators
    lstlist_all_WLS_NOint_regression_estimators[[k]] = EM_and_matching[["WLS_NOint_mat_reg_estimators"]]
    lstlist_all_WLS_YESint_regression_estimators[[k]] = EM_and_matching[["WLS_YESint_mat_reg_estimators"]]
    lstlist_all_OLS_NOint_regression_estimators[[k]] = EM_and_matching[["OLS_NOint_mat_reg_estimators"]]
    lstlist_all_OLS_YESint_regression_estimators[[k]] = EM_and_matching[["OLS_YESint_mat_reg_estimators"]]
    lstlist_all_CI[[k]] = EM_and_matching[["CI_mat"]]
    lstlist_all_EM_coeffs[[k]] = EM_and_matching[["coeffs_df"]]
    lstlist_all_excluded_included_matching[[k]] = EM_and_matching[["mat_excluded_included_matching"]]
    lstlist_all_repeated_as_and_pro[[k]] = EM_and_matching[["mean_list_repeated_as_and_pro"]]
    rownames(lstlist_all_repeated_as_and_pro[[k]]) = paste0("s", k, rownames(lstlist_all_repeated_as_and_pro[[k]]))
    lstlist_all_matched_units[[k]] = EM_and_matching[["mean_list_matched_units"]]
    rownames(lstlist_all_matched_units[[k]]) = paste0("s", k, rownames(lstlist_all_matched_units[[k]]))
    
    # TODO in list_all_diff_distance_aspr_asas and list_all_std_mean_diff   
    # when diff_distance_aspr_asas is positive, the matches between as to as are closer in X, 
    #then matches between as to protected, since its the abs diff between as to pro matches - abs diff between as to as matches 
    lstlist_all_diff_distance_aspr_asas[[k]] = EM_and_matching[["mat_diff_distance_aspr_asas"]]
    rownames(lstlist_all_diff_distance_aspr_asas[[k]]) = paste0("s", k, rownames(lstlist_all_diff_distance_aspr_asas[[k]]))
    lstlist_all_std_mean_diff[[k]] = EM_and_matching[["mean_list_std_mean_diff"]]
    rownames(lstlist_all_std_mean_diff[[k]]) = paste0("s", k, rownames(lstlist_all_std_mean_diff[[k]]))
    lstlist_all_means_by_subset[[k]] = EM_and_matching[["mean_list_means_by_subset"]]
    rownames(lstlist_all_means_by_subset[[k]]) = paste0("s", k, rownames(lstlist_all_means_by_subset[[k]]))
    print(paste0("sim is ", sim))
  }
  end_time <- Sys.time()
  print(paste0("in the end of outer for loop ", k, ", ", difftime(end_time, start_time)))
}
########################################################################

########################################################################
# save lists from simulation
save(lstlist_all_mat_SACE_estimators, file = "lstlist_all_mat_SACE_estimators.RData")
save(lstlist_all_WLS_NOint_regression_estimators, file = "lstlist_all_WLS_NOint_regression_estimators.RData")
save(lstlist_all_WLS_YESint_regression_estimators, file = "lstlist_all_WLS_YESint_regression_estimators.RData")
save(lstlist_all_OLS_NOint_regression_estimators, file = "lstlist_all_OLS_NOint_regression_estimators.RData")
save(lstlist_all_OLS_YESint_regression_estimators, file = "lstlist_all_OLS_YESint_regression_estimators.RData")
save(lstlist_all_CI, file = "lstlist_all_CI.RData")
save(lstlist_all_EM_coeffs, file = "lstlist_all_EM_coeffs.RData")
save(lstlist_all_excluded_included_matching, file = "lstlist_all_excluded_included_matching.RData")
save(lstlist_all_repeated_as_and_pro, file = "lstlist_all_repeated_as_and_pro.RData")
save(lstlist_all_matched_units, file = "lstlist_all_matched_units.RData")
# TODO list_all_diff_distance_aspr_asas is already in mat_SACE_estimators,
save(lstlist_all_diff_distance_aspr_asas, file = "lstlist_all_diff_distance_aspr_asas.RData")
save(lstlist_all_std_mean_diff, file = "lstlist_all_std_mean_diff.RData")
save(lstlist_all_means_by_subset, file = "lstlist_all_means_by_subset.RData")


# when sim == 3, divide each misspec (in 1:3) list
list_all_mat_SACE_estimators = lstlist_all_mat_SACE_estimators[[1]]
list_all_WLS_NOint_regression_estimators = lstlist_all_WLS_NOint_regression_estimators[[1]]
list_all_WLS_YESint_regression_estimators = lstlist_all_WLS_YESint_regression_estimators[[1]]
list_all_OLS_NOint_regression_estimators = lstlist_all_OLS_NOint_regression_estimators[[1]]
list_all_OLS_YESint_regression_estimators = lstlist_all_OLS_YESint_regression_estimators[[1]]
list_all_CI = lstlist_all_CI[[1]]
list_all_EM_coeffs = lstlist_all_EM_coeffs[[1]]
list_all_excluded_included_matching = lstlist_all_excluded_included_matching[[1]]
list_all_repeated_as_and_pro = lstlist_all_repeated_as_and_pro[[1]]
list_all_matched_units = lstlist_all_matched_units[[1]]
list_all_diff_distance_aspr_asas = lstlist_all_diff_distance_aspr_asas[[1]]
list_all_std_mean_diff = lstlist_all_std_mean_diff[[1]]
list_all_means_by_subset = lstlist_all_means_by_subset[[1]]

########################################################################


########################################################################
param_measures_computing = function(df, abs_bool = FALSE){
  if(abs_bool == TRUE){df = abs(df)}
  df[c("mean", "med", "SD"), ] = rbind(
    apply(df, 2, mean), apply(df, 2, median), apply(df, 2, sd)                                    
  )
  return(df)
}
########################################################################

########################################################################
# TODO calculate mean med and sd for pairmatch with small amount of na's
# TODO if I want to come back to the original list list_all_mat_SACE_estimators, just load it from rdata
param_measures_excluding_na = function(df, param_n_sim){
  df[c("mean", "med", "SD"), ] = rbind(
    apply(df[c(1:param_n_sim), ], 2, mean, na.rm=TRUE),
    apply(df[c(1:param_n_sim), ], 2, median, na.rm=TRUE),
    apply(df[c(1:param_n_sim), ], 2, sd, na.rm=TRUE)                                    
  )
  #tail(rownames(df),3) = c("mean", "med", "SD")
  return(df)
}

df = list_all_mat_SACE_estimators[[1]]
df = is.na(df)
colSums(is.na(df))
mat_nans = cbind(colSums(is.na(list_all_mat_SACE_estimators[[1]])), 
                 colSums(is.na(list_all_mat_SACE_estimators[[2]])))

list_all_mat_SACE_estimators_excluding_na = 
  lapply(1:length(list_all_mat_SACE_estimators), function(l){
    param_measures_excluding_na(list_all_mat_SACE_estimators[[l]], param_n_sim)
  })
list_all_mat_SACE_estimators = list_all_mat_SACE_estimators_excluding_na
########################################################################

########################################################################
# function that changes rownames to A,B,C per each parameter (gammas) set
#df = mat_all_estimators
change_rownames_to_LETTERS_by_param_set = function(df, mat_params=mat_gamma, num_of_param_measures=num_of_param_measures_per_param_set){
  rownames(df) = paste0(rep(LETTERS[1:nrow(mat_params)],each = num_of_param_measures), "_",
                        rep(rownames(df)[1:num_of_param_measures], times = nrow(mat_params)))
  return(df)
}
####################################LETTERS####################################

# TODO calcultae MSE and coverage
# a = list_all_mat_SACE_estimators[[1]]
# a = a[1:6,]
# RMSE <- function (x) sqrt(mean((x-mean(x))^2))
# x <- list.cbind(lapply(a, FUN = RMSE)) %>% data.frame()
# sqrt(mean((a$DING_est - mean(a$DING_est))^2))

# TODO :) start calculating

# summary of SACE estimators 
mat_all_estimators = list.rbind(lapply(list_all_mat_SACE_estimators, tail, length(param_measures)))
#num_of_param_measures_per_param_set = nrow(mat_all_estimators) / nrow(mat_gamma) # = length(param_measures)
mat_all_estimators = data.frame(subset(mat_all_estimators,
                                       select = grep("gamma", colnames(mat_all_estimators))),
                                subset(mat_all_estimators,
                                       select = -grep("gamma", colnames(mat_all_estimators))))
mat_all_estimators = change_rownames_to_LETTERS_by_param_set(mat_all_estimators)
# if I want to delete weights est
# mat_all_estimators = subset(mat_all_estimators,
#                     select = -grep("MATCH_w_", colnames(mat_all_estimators)))

mean_all_estimators = mat_all_estimators[grep("mean", rownames(mat_all_estimators)) , ]
mean_estimators_over_all_parameters = 
  rbind(apply(mean_all_estimators, 2, mean), apply(mean_all_estimators, 2, sd))

# list_all_regression_estimators
#mat_regression_estimators = list.rbind(lapply(list_all_regression_estimators, tail, 3))
WLS_NOint_mat_regression_estimators = list.rbind(lapply(list_all_WLS_NOint_regression_estimators, tail, length(param_measures)))
WLS_NOint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(WLS_NOint_mat_regression_estimators)
WLS_YESint_mat_regression_estimators = list.rbind(lapply(list_all_WLS_YESint_regression_estimators, tail, length(param_measures)))
WLS_YESint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(WLS_YESint_mat_regression_estimators)
OLS_NOint_mat_regression_estimators = list.rbind(lapply(list_all_OLS_NOint_regression_estimators, tail, length(param_measures)))
OLS_NOint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(OLS_NOint_mat_regression_estimators)
OLS_YESint_mat_regression_estimators = list.rbind(lapply(list_all_OLS_YESint_regression_estimators, tail, length(param_measures)))
OLS_YESint_mat_regression_estimators = change_rownames_to_LETTERS_by_param_set(OLS_YESint_mat_regression_estimators)



# summary of EM estimators for the gamms's- the PS coefficient- logistic reg of stratum on X
#sapply(list_all_EM_coeffs, `[[`, 1001)
summary_EM_coeffs = list.rbind(lapply(list_all_EM_coeffs, function(x) x[c((param_n_sim + 1):(param_n_sim + 5))]))
rownames(summary_EM_coeffs) = paste0(rep(LETTERS[1:nrow(mat_gamma)],each = ncol(mat_gamma)), "_",
                                     rep(rownames(summary_EM_coeffs)[1:ncol(mat_gamma)],times = nrow(mat_gamma)))

# summary of list_all_excluded_included_matching
mean_excluded_included_matching = list.rbind(list_all_excluded_included_matching)
mean_excluded_included_matching = mean_excluded_included_matching[grep(c("mean|sd"), rownames(mean_excluded_included_matching)), ]

# mean_excluded_included_matching = mean_excluded_included_matching[
# (as.numeric(rownames(mean_excluded_included_matching)) %% (param_n_sim + 2)) %in% c((param_n_sim + 1), 0) , ]

# rownames(mean_excluded_included_matching) =
#   paste0(rep(LETTERS[1:nrow(mat_gamma)], each = 2), "_", c("mean", "sd"))

#colnames(mean_excluded_included_matching) = colnames(mat_excluded_included_matching)


# list_all_repeated_as_and_pro
mat_all_repeated_as_and_pro = list.rbind(list_all_repeated_as_and_pro)
# hist
repeated_histogram(mat_all_repeated_as_and_pro)
mat_all_repeated_as_and_pro_t = t(mat_all_repeated_as_and_pro)
reT_mat_all_repeated_as_and_pro = subset(mat_all_repeated_as_and_pro_t,
                                         select = -grep("repF_", colnames(mat_all_repeated_as_and_pro_t)))


# list_all_matched_units
mat_all_matched_units = list.rbind(list_all_matched_units)

# list_all_std_mean_diff
mat_all_std_mean_diff = list.rbind(list_all_std_mean_diff)

#list_all_std_mean_diff
# list_abs_std_mean_diff = lapply(1:length(list_all_std_mean_diff), function(l){
#   param_measures_computing(list_all_std_mean_diff[[l]], abs_bool = TRUE)
# })
# mat_abs_std_mean_diff = list.rbind(lapply(list_abs_std_mean_diff, tail, 3))


# list_all_means_by_subset
#mat_all_means_by_subset = list.rbind(list_all_means_by_subset)
mat_all_means_by_subset = NULL
first_3_rows = c("mean_as", "mean_A0_S1", "mean_A1_S1_as")
for(i in c(1:length(list_all_means_by_subset))){
  first_3_rows_ind = grep(paste0(first_3_rows, collapse = "|"), rownames(list_all_means_by_subset[[i]]))
  mat_all_means_by_subset = rbind(mat_all_means_by_subset,
                                  list_all_means_by_subset[[i]][-first_3_rows_ind[-c(1:3)], ])
}
########################################################################

########################################################################
# save summaries
save(mat_all_estimators, file = "mat_all_estimators.RData")
save(mean_estimators_over_all_parameters, file = "mean_estimators_over_all_parameters.RData")
save(WLS_NOint_mat_regression_estimators, file = "WLS_NOint_mat_regression_estimators.RData")
save(WLS_YESint_mat_regression_estimators, file = "WLS_YESint_mat_regression_estimators.RData")
save(OLS_NOint_mat_regression_estimators, file = "OLS_NOint_mat_regression_estimators.RData")
save(OLS_YESint_mat_regression_estimators, file = "OLS_YESint_mat_regression_estimators.RData")

save(summary_EM_coeffs, file = "summary_EM_coeffs.RData")
save(mean_excluded_included_matching, file = "mean_excluded_included_matching.RData")
#save(mat_abs_std_mean_diff, file = "mat_abs_std_mean_diff.RData")
save(mat_all_repeated_as_and_pro, file = "mat_all_repeated_as_and_pro.RData")
save(reT_mat_all_repeated_as_and_pro, file = "reT_mat_all_repeated_as_and_pro.RData")
save(mat_all_matched_units, file = "mat_all_matched_units.RData")
save(mat_all_std_mean_diff, file = "mat_all_std_mean_diff.RData")
save(mat_all_means_by_subset, file = "mat_all_means_by_subset.RData")

########################################################################

########################################################################
# TABLES for LaTeX
# mat_nans %>% xtable(caption = 
#                       "5 continous covariates, matching on mahal, ")
# 


# print(t(mat_all_estimators) %>% xtable(caption = 
#  "nsim = 1000, 5 continous covariates, PM matching on PS, PA = 0.5"), include.rownames = FALSE)
# print(t(mean_estimators_over_all_parameters) %>% xtable(), include.rownames = FALSE)
# print(t(mat_regression_estimators) %>% xtable(), include.rownames = FALSE)
# print(summary_EM_coeffs %>% xtable(), include.rownames = FALSE)
# print(t(mean_excluded_included_matching) %>% xtable(), include.rownames = FALSE)
# print(t(mat_all_repeated_as_and_pro) %>% xtable(), include.rownames = FALSE)
# print(t(mat_all_std_mean_diff) %>% xtable(), include.rownames = FALSE)


t(mat_all_estimators) %>% xtable(digits=c(3), caption = paste0("True Model with interactions. ",
                                                               cont_x, " continous covariates, nsim=", param_n_sim, " n=", param_n, 
                                                               ", 1 over var matching on X with caliper ", caliper, " sd PS, PA = ", prob_A))
t(mean_estimators_over_all_parameters) %>% xtable(digits=c(3))
print(mat_all_means_by_subset %>% xtable(digits=c(4)), size="\\fontsize{9pt}{11pt}\\selectfont")
t(WLS_NOint_mat_regression_estimators) %>% xtable(digits=c(3), caption = "WLS wout interactions")
t(WLS_YESint_mat_regression_estimators) %>% xtable(digits=c(3), caption = "WLS with interactions")
t(OLS_NOint_mat_regression_estimators) %>% xtable(digits=c(3), caption = "OLS wout interactions")
t(OLS_YESint_mat_regression_estimators) %>% xtable(digits=c(3), caption = "OLS with interactions")

summary_EM_coeffs %>% xtable(digits=c(3))
t(mean_excluded_included_matching) %>% xtable(digits=c(3))
# t(mat_abs_std_mean_diff) %>% xtable(digits=3,
#                         caption = "mean abs std diff of covariates")

a = t(mat_all_matched_units)
colnames(a) = substring(colnames(a), 1,8)
print(xtable(a, digits=c(0),caption = 
               paste0("True Model with interactions. ",
                      cont_x, " continous covariates, nsim=", param_n_sim, " n=", param_n, 
                      ", 1 over var matching on X with caliper ", caliper, " sd PS, PA = ", prob_A)), 
      size="\\fontsize{5pt}{16pt}\\selectfont")
# size="\\tiny",
print(reT_mat_all_repeated_as_and_pro %>% xtable(digits=c(3)), size="\\fontsize{7pt}{10pt}\\selectfont")
print(t(mat_all_std_mean_diff) %>% xtable(digits=c(3)), size="\\fontsize{5pt}{16pt}\\selectfont")
#t(mat_all_std_mean_diff) %>% xtable(digits=c(3))
########################################################################

########################################################################
# TABLE DESIGN
list_all_CI_temp = list_all_CI # list_all_CI from main
data_set_names_vec = c("_all", "_wout_O_0_0", "_S1")

lst_final_tables_ALL_est_ALL_dataset = TABLES_before_coverage_and_wout_naive_and_DING(data_set_names_vec,
                                                                                      estimators_names=c("PS Crude", "mahal Crude", "Crude", "HL", "BC", "BC caliper"))
lst_final_tables_ALL_est_ALL_dataset_PLUS_naives = 
  TABLES_add_coverage(data_set_names_vec, lst_final_tables_ALL_est_ALL_dataset, list_all_CI)
lst_final_tables_ALL_est_ALL_dataset = lst_final_tables_ALL_est_ALL_dataset_PLUS_naives$lst_final_tables_ALL_est_ALL_dataset
naives_before_matching_coverage = lst_final_tables_ALL_est_ALL_dataset_PLUS_naives$naives_before_matching_coverage
final_tables = TABLES_add_naive_and_ding(data_set_names_vec, lst_final_tables_ALL_est_ALL_dataset, naives_before_matching_coverage)
final_tables_general = adjustments_for_final_tables(final_tables)
final_tables_crude = adjustments_for_final_tables_crude_est(final_tables)
save(final_tables_general, file = "final_tables_general.RData")
save(final_tables_crude, file = "final_tables_crude.RData")

pis = mat_all_estimators[grep("_mean",rownames(mat_all_estimators)),grep("pi",colnames(mat_all_estimators))] %>% round(3)
pis = data.frame(pi_as=pis$pi_as, pi_pro=pis$pi_pro, pi_ns=pis$pi_ns)

print(final_tables_general$`_S1` %>% xtable(digits=c(3),
                                            caption = paste0("True model with interactions. ", ifelse(misspec_PS==1, "ONLY PS model misspecification", ifelse(misspec_PS==2, "PS, maha and regression model misspecification", "No misspecification")), 
                                                             ". ", (dim_x-1) ," X's, Delta Method, ",
                                                             param_n_sim, " replications, ", "N=", param_n, ", pA = ", prob_A, ". OLS clustered SE")),
      size="\\fontsize{11pt}{11pt}\\selectfont")
print(pis %>% xtable, size="\\fontsize{15pt}{15pt}\\selectfont")
print(final_tables_crude$`_S1` %>% xtable(digits=c(3),
                                          caption = paste0("True model with interactions. ", ifelse(misspec_PS==1, "ONLY PS model misspecification", ifelse(misspec_PS==2, "PS, maha and regression model misspecification", "No misspecification")), 
                                                           ". ", (dim_x-1) ," X's, Delta Method, ",
                                                           param_n_sim, " replications, ", "N=", param_n, ", pA = ", prob_A, ". OLS clustered SE")),
      size="\\fontsize{11pt}{11pt}\\selectfont")
print(mat_gamma %>% xtable, size="\\fontsize{6pt}{6pt}\\selectfont")
print(betas_GPI %>% xtable, size="\\fontsize{12pt}{12pt}\\selectfont")
########################################################################

########################################################################
# for the matrix with many parameters
round_table = round(mat_all_estimators, 2)
View(round_table)
pdf(file = "q.pdf")
grid.table(t(round_table), theme=ttheme_minimal(base_size = 3.2))
# ttheme_default(base_size = 12, base_colour = "black", base_family = "",
# parse = FALSE, padding = unit(c(4, 4), "mm"), ...)
dev.off()
########################################################################

#repeated_histogram(mat_all_repeated_as_and_pro)
repeated_histogram = function(mat_all_repeated_as_and_pro){
  mat_all_repeated_as_and_pro_t = t(mat_all_repeated_as_and_pro)
  reT_mat_all_repeated_as_and_pro = subset(mat_all_repeated_as_and_pro_t,
                                           select = -grep("repF_", colnames(mat_all_repeated_as_and_pro_t)))
  
  reT_mat_all_repeated_hist = reT_mat_all_repeated_as_and_pro[-grep("mean", 
                                                                    rownames(reT_mat_all_repeated_as_and_pro)),] 
  len_each_strat = nrow(reT_mat_all_repeated_hist) / 2
  title_vec = colnames(reT_mat_all_repeated_hist)
  #par(mfrow = c(2,3))
  pdf(file = "repeated_by_str_hist.pdf")
  for(i in 1 : length(title_vec) ){
    #temp = t(data.frame(reT_mat_all_repeated_hist[,i]))
    #as = data.frame(re = c(1 : len_each_strat), count = temp[,grep("as_", colnames(temp))])
    #pro = data.frame(re = c(1 : len_each_strat), count = temp[,grep("pro", colnames(temp))])
    #ggplot(data = as) + geom_bar(aes(x = re, y = count), stat = "identity")
    
    temp = data.frame(reT_mat_all_repeated_hist[,i])
    temp = data.frame( count = temp[,1], re = rep(c(1 : len_each_strat), times=2),
                       stratum = rep(c("as","pro"), each=len_each_strat) )
    
    print( ggplot(data = temp, aes(x = re, fill = stratum)) +
             stat_identity(data = temp, aes(x = re, y = count), geom = "bar", alpha = 1) + 
             ggtitle(paste0("after matching and excluding pairs with non-surv ", 
                            title_vec[i])) )
    
    # print( ggplot(data = temp, aes(x = re, fill = stratum)) + 
    #          geom_bar(aes(x = re, y = count), stat = "identity") + geom_density(alpha = 0.5) + 
    #          ggtitle(paste0("after matching and excluding pairs with non-surv ", 
    #                         title_vec[i])) )
    # 
    
    
  } 
  dev.off()
}


data = data.frame(x=c(1:10), y=10*c(1:10))
ggplot() +
  stat_identity(data = data, aes(x, y), geom = "bar", alpha = 0.5)


####################################################################
df1 = data.frame(id=rep(c(1:10), each=20), day=rep(c(1:20), time=10),
                 x1=rnorm(200,0,1), x2=rnorm(200,0,1), x3=rnorm(200,0,1))
df1[,-c(1,2)][df1[,-c(1,2)]>1.25] = Inf
df2 = df1; df2[sapply(df2, is.infinite)] <- -1000
data.table(df2)[,lapply(.SD, max), by=id]

# x<- df2[,2:7]
# colMax <- function(x)apply(x, 2, max)
# mymax<- t(sapply(unique(id), function(i)colMax(x[which(id==i), ])))

replace_inf_with_max = function(x){
  x[is.infinite(x)] = max(x[!is.infinite(x)])
  return(x)
}
df_final= data.table(df1)[, lapply(.SD, replace_inf_with_max), by=id]
####################################################################
