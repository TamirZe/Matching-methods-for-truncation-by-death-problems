# installing/loading the package:
'''if(!require(installr)) { install.packages("installr"); require(installr)} #load / install+load installr
# step by step functions:
check.for.updates.R() # tells you if there is a new version of R or not.
install.R() # download and run the latest R installer
copy.packages.between.libraries()''' # copy your packages to the newest R installation from the one version before it (if ask=T, it will ask you between which two versions to perform the copying)

'''devtools::install_github("shuyang1987/multilevelMatching")
install.packages("rlist"); install.packages("locfit"); install.packages("plyr")
install.packages("nnet"); install.packages("xtable"); install.packages("doParallel")
install.packages("matrixStats"); install.packages("data.table"); install.packages("rlang")
install.packages("dplyr"); install.packages("reshape"); install.packages("MASS"); install.packages("mgsub")
install.packages("ggplot2"); install.packages("rockchalk"); install.packages("nnet")
installed.packages("stats"); install.packages("optmatch"); install.packages("DOS"); install.packages("PerformanceAnalytics")
install.packages("Matching"); install.packages("sandwich"); install.packages("rmutil")
install.packages("sandwich"); install.packages("rmutil"); install.packages("caret")'''

install.packages("BiocManager")
library("BiocManager")
BiocManager::install("genefilter")

library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); library(rlang);library(glue)
# library(multilevelMatching); library(PerformanceAnalytics); library(lmtest)
library(matrixStats); library(data.table); library(dplyr); library(reshape); library(MASS); library(Hmisc); library(lmtest)
library(ggplot2); library(rockchalk); library(stats); library(rlist); library(mgsub); library(reshape2); library(gridExtra)
library(optmatch); library(DOS); library(Matching); library(sandwich); library(rmutil); library(clubSandwich); library(tableone)
library(sandwich); library(rmutil);  library(caret); library(splitstackshape); library(MatchIt); library(PerformanceAnalytics)

source("matching_scripts/matching_PS_multiple.R")
source("matching_scripts/matching_PS_basic.R")
#source("matching_scripts/matching_from_real_data.R")
#source("my_EM_V2.R")
source("EM_V3_eps_stop.R"); source("OLS_WLS_estimator.R")
#source("sim1.R")
#source("simulations_scripts/sim2.R")
source("simulations_scripts/sim2DingLuEst.R")
source("DING_model_assisted_estimator.R")
source("TABLES/table_design_multiple_func.R")
source("TABLES/coverage_naive_est.R")
#source("Extra code/TABLES/table_design_func.R"); source("Extra code/TABLES/table_design_w_BCest_func.R")
#swog_path = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/new papers data/B1191Ding/"
source(paste0("B1191Ding/", "PS_M_weighting.R"))
source(paste0("B1191Ding/", "PS_M_weighting_SA.R"))

'''library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)
#foreach(j = seq_along(numeric(length(c(1:n_sim)))), .combine=c) %dopar%{
stopCluster(cl)'''

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
# parameters 

# treatment probability
prob_A = 0.5

# parameters for simulating x
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 6; cont_x = 5; categ_x = 0; vec_p_categ = rep(0.5, categ_x); dim_x_misspec = 2
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)

# EM convergence parameters
# iterations = 1000; epsilon_EM = 10^-6 
iterations = 200; epsilon_EM = 10^-6 # epsilon_EM is for the EM convergence

# monotonicity and PI assumptions
monotonicity_assumption = "mono"; PI_assum = "strong"

# misspec PARAMETERS (for PS model and Y model:
# TODO misspec_PS: 0 <- NO, 1:only PS model, 2: PS model, and matching (mahalanobis, eucl...) and regression
misspec_outcome_funcform = FALSE; match_and_reg_watch_true_X = FALSE
U_factor=1; funcform_factor_sqr=-3; funcform_factor_log=3
mean_x_misspec = rep(0.5, dim_x_misspec)
sim=2
misspec_PS = 0 # 0: no misspec # 2: ff misspec
# sensitivity paramters
epsilon_1_GPI = 1

# for covariance between X 
# TODO Sigma here is vcov not sd cor. need to change in the PO sim in sim1 rows ~80
'''corr_x = 0.3
sqrt(var_x)
cov_x <- diag(sqrt(5),cont_x) + corr_x - diag(1,cont_x) * corr_x # cov_x is vcov and corr
X = mvrnorm(param_n, mu = mean_x, Sigma =  cov_x)
apply(X, 2, mean); apply(X, 2, var); cov(X); cor(X)'''


#############################################################################################

##########################################################
# betas under GPI # betas_GP
# TODO YES interactions between A and X:
betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0)))
betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3)))
betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(3,3,0,1,3),2))))

# TODO NO interactions between A and X: "simple effect" is 2
betas_GPI = as.matrix(rbind(c(22,3,4,5), c(20,3,4,5)))
betas_GPI = as.matrix(rbind(c(22,3,4,5,1,3), c(20,3,4,5,1,3)))
betas_GPI = as.matrix(rbind(c(22,rep(c(5,2,1,3,5),2)), c(20,rep(c(5,2,1,3,5),2))))

rownames(betas_GPI) = c("beta_treatment", "beta_control")

# print beta
###############################################################################################
'''beta_print = rbind.fill(data.frame(betas_GPIwo1), data.frame(betas_GPIwo2), data.frame(betas_GPIwo3),
                        data.frame(betas_GPIw1) , data.frame(betas_GPIw2), data.frame(betas_GPIw3)) 
beta_print = apply(beta_print, 2, as.character)
print(beta_print %>% xtable(), size="\\fontsize{15pt}{15pt}\\selectfont", include.rownames=F)'''
###############################################################################################

var_GPI = as.matrix(rbind(1, 1))
rownames(var_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4 #0.4 #1
#sds_GPI_PO = 1
##########################################################

# add the harmed stratum parameters if needed
# if (monotonicity_assumption == "nothing"){
#   betas = rbind(betas, beta_har0)
#   sigma_square = c(sigma_square, sigma_square_har0)
# }
#############################################################################################


#############################################################################################
#TODO mat_gamma ####
#############################################################################################
################ new values for gamma's 3X ####
#TODO 25,50,75: large pi pro 3X
mat_gamma = matrix(c(
  #c(0.11, rep(-0.34, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
  #,c(-0.1, rep(0.315, dim_x-1)), c(-0.27, rep(-0.5, dim_x-1))
  c(0.08, rep(-0.355, dim_x-1)), c(0.2, rep(-0.11, dim_x-1))
  ,c(-0.1, rep(0.27, dim_x-1)), c(-0.52, rep(-0.6, dim_x-1))
  ,c(0.6, rep(1.325, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))
)
,nrow = 3, byrow = T) 

#TODO 25,50,75: small pi pro 3X:
mat_gamma = matrix(c(
  c(1.05, rep(-0.4, dim_x-1)), c(0.8, rep(0.72, dim_x-1))
  ,c(0.92, rep(1, dim_x-1)), c(0.9, rep(0.87, dim_x-1))
  ,c(0.84, rep(1.2, dim_x-1)), c(0.32, rep(-0.2, dim_x-1))
)
,nrow = 3, byrow = T) 

#TODO 25,50,75: small pi pro 3X for smaller prop in misspec: @@@@
mat_gamma = matrix(c(
  # c(1, rep(0.84, dim_x-1)), c(0.7, rep(0.85, dim_x-1))
  # ,c(1.5, rep(0.6, dim_x-1)), c(0.84, rep(0.8, dim_x-1))
  c(1.05, rep(-0.39, dim_x-1)), c(0.8, rep(0.72, dim_x-1))
  ,c(-0.07, rep(1.24, dim_x-1)), c(1.25, rep(0.24, dim_x-1)) # @@@@
  ,c(0.84, rep(1.2, dim_x-1)), c(0.32, rep(-0.2, dim_x-1))
  #,c(1.05, rep(-0.2, dim_x-1)), c(0.5, rep(-0.3, dim_x-1))
) ,nrow = 3, byrow = T)
#############################################################################################

#############################################################################################
################  new values for gamma's 5x ####
#TODO 25,50,75: large pi pro 5X:
mat_gamma = matrix(c(
      c(0.3, rep(-0.3, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
     ,c(-0.6, rep(0.35, dim_x-1)), c(-0.27, rep(-0.5, dim_x-1))
      ,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))
    )
    ,nrow = 3, byrow = T)


#TODO 25,50,75: large pi pro 5X for larger prop in misspec: @@@@:
mat_gamma = matrix(c(
  c(0.3, rep(-0.3, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
  ,c(-0.05, rep(0.16, dim_x-1)), c(-0.4, rep(-0.25, dim_x-1)) # @@@@
  ,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))
)
,nrow = 3, byrow = T)


#TODO 25,50,75: small pi pro 5X:
mat_gamma = matrix(c(
                    # TODO 25,50,75 (take 2,5,7)
                     # c(0.1, rep(-0.2, dim_x-1)), c(0.3, rep(-0.1, dim_x-1))
                     #c(0.3, rep(-0.3, dim_x-1)), c(0.275, rep(-0.1, dim_x-1))
                     # ,c(0.2, rep(0.2, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
                     # ,c(-0.12, rep(0.25, dim_x-1)), c(-0.25, rep(-0.1, dim_x-1))
                     #  c(-0.25, rep(0.5, dim_x-1)), c(0.25, rep(0.1, dim_x-1))

                     # TODO 0.5 as, 0.35 pro
                    #,c(-0.6, rep(0.35, dim_x-1)), c(-0.27, rep(-0.5, dim_x-1))

                     # ,c(-1.3, rep(1.25, dim_x-1)), c(-0.4, rep(0.75, dim_x-1))
                     #,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))

                     # TODO 20,40,60 small pi pro
                     # ,c(0.8, rep(-0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
                     # ,c(0.15, rep(0.5, dim_x-1)), c(0.25, rep(0.5, dim_x-1))
                     # ,c(-0.675, rep(0.675, dim_x-1)), c(-0.2, rep(-0.2, dim_x-1))

                     # TODO 25,50,75: big pi pro:
                     #  c(0.2, rep(-0.425, dim_x-1)), c(-0.3, rep(-0.2, dim_x-1)),
                     #    c(-0.6, rep(0.375, dim_x-1)), c(0, rep(-0.5, dim_x-1)),
                     #    c(-0.42, rep(1.25, dim_x-1)), c(-0.45, rep(-0.1, dim_x-1))

                     # TODO 25,50,75: small pi pro:
                      c(0.75, rep(0.125, dim_x-1)), c(0.6, rep(0.6, dim_x-1))
                     ,c(0.45, rep(0.75, dim_x-1)), c(0.62, rep(0.6, dim_x-1))
                     ,c(-0.12, rep(1.25, dim_x-1)), c(0.45, rep(-0.2, dim_x-1))
                     )
                   ,nrow = 3, byrow = T)

# COVERAGE RATE proporions 50 20 
mat_gamma = matrix(c(
  c(-0.275, rep(0.5, dim_x-1)), c(0.235, rep(0.1, dim_x-1))
  ,c(-1.31, rep(1.24, dim_x-1)), c(-0.42, rep(0.75, dim_x-1))
  ,c(0.45, rep(0.75, dim_x-1)), c(0.625, rep(0.6, dim_x-1))
)
,nrow = 3, byrow = T)
###################################################################


################ new values for gamma's 10X ################ 
#TODO 25,50,75: large pi pro 10X:
mat_gamma = matrix(c(
  c(0.58, rep(-0.2, dim_x-1)), c(0.54, rep(-0.1, dim_x-1))
  ,c(-0.954, rep(0.25, dim_x-1)), c(-0.31, rep(-0.16, dim_x-1))
  #,c(-0.19, rep(0.1, dim_x-1)), c(0.58, rep(-0.41, dim_x-1))
  #,c(-0.3, rep(0.475, dim_x-1)), c(0, rep(-0.1, dim_x-1))
  ,c(-1.01, rep(0.8, dim_x-1)), c(-1.3, rep(0.45, dim_x-1))
  #,c(-0.78, rep(0.6, dim_x-1)), c(0.14, rep(-0.2, dim_x-1))
)
,nrow = 3, byrow = T)

#TODO 25,50,75: small pi pro 10X:
mat_gamma = matrix(c(
  c(0.042, rep(0.29, dim_x-1)), c(-0.1, rep(0.51, dim_x-1))
  ,c(0.024, rep(0.41, dim_x-1)), c(0.26, rep(0.32, dim_x-1))
  ,c(-0.43, rep(0.7, dim_x-1)), c(-0.275, rep(0.3, dim_x-1))
)
,nrow = 3, byrow = T)


#TODO 25,50,75: small pi pro 10X for smaller prop in misspec: @@@@::
mat_gamma = matrix(c(
  c(0.042, rep(0.29, dim_x-1)), c(-0.1, rep(0.51, dim_x-1))
  ,c(0.024, rep(0.41, dim_x-1)), c(0.26, rep(0.32, dim_x-1))
  ,c(-0.45, rep(0.61, dim_x-1)), c(0.4, rep(-0.02, dim_x-1))
)
,nrow = 3, byrow = T)


#############################################################################################

gamma_pro = rep(0, dim_x)
colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("as", "ns"), each = dim_x) )

# TODO calculate with function: ####
extract_pis_from_scenarios = function(nn=250000){
  big_lst = list(); mat_x_as <- mat_pis <- mat_x_by_g_A <- NULL
  for( k in c(1 : nrow(mat_gamma)) ){
    gamma_as = as.numeric(mat_gamma[k, c(1:dim_x)])
    gamma_ns =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
    lst_mean_x_and_pi = simulate_data_function(gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, 
                   param_n=nn, misspec_PS=2, misspec_outcome_funcform=FALSE, 
                   U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                   epsilon_1_GPI = 1, only_mean_x_bool=TRUE)
    big_lst[[k]] = lst_mean_x_and_pi
    mat_x_as = rbind(mat_x_as, lst_mean_x_and_pi$x_as)
    mat_pis = rbind(mat_pis, lst_mean_x_and_pi$pi)
    mat_x_by_g_A = rbind(mat_x_by_g_A, data.frame(Scenar = k, lst_mean_x_and_pi$mean_by_A_g))# %>%
      #subset(select = c(Scenar,A,g, grep("X", colnames(lst_mean_x_and_pi$mean_by_A_g))))
    
  }
  mat_pis = data.frame(pi_as=mat_pis[,1], pi_pro=mat_pis[,3], pi_ns=mat_pis[,2])
  round(mat_pis,3)
  return(list(mat_x_by_g_A=mat_x_by_g_A, big_lst=big_lst))
}
mat_gamma[,c(1,2,7,8)]

big_lst = list(); big_mat_x_by_g_A=NULL
for(i in 1:100){
  print(i)
  big_lst[[i]] = extract_pis_from_scenarios(param_n=2000)
  big_mat_x_by_g_A = rbind(big_mat_x_by_g_A, big_lst[[i]]$mat_x_by_g_A)
}
big_mat_x_by_g_A = subset(big_mat_x_by_g_A, select = c(Scenar,A,g, grep("X", colnames(big_mat_x_by_g_A))))
big_mat_x_by_g_A = data.table(big_mat_x_by_g_A)[, lapply(.SD, mean), by=c("Scenar", "A", "g")] %>% arrange(Scenar, g, A)
#############################################################################################

param_n = 2000; param_n_sim = 500 # param_n = 2000; param_n_sim = 1000
caliper = 0.25; pairmatch_bool = FALSE
match_on = "O11_posterior_ratio" # NULL # Feller and Mialli: "EMest_p_as" # Ding Lu appendix: "O11_posterior_ratio"
mu_x_fixed = FALSE; mat_x_as; x_as = mat_x_as[1,]

param_measures = c("mean","med","sd","MSE"); num_of_param_measures_per_param_set = length(param_measures)
list_all_mat_SACE_estimators <- list_all_WLS_NOint_regression_estimators <- list_all_WLS_YESint_regression_estimators <-
list_all_OLS_NOint_regression_estimators <- list_all_OLS_YESint_regression_estimators <- list_all_CI <- 
list_all_EM_coeffs <- list_all_excluded_included_matching <-
list_all_repeated_as_and_pro <- list_all_diff_distance_aspr_asas <- list_all_matched_units <-
list_all_std_mean_diff <- list_all_means_by_subset  <- list_all_EM_not_conv <- list_all_BCclpr <- list()
# for naive estimators and their SE and CI
only_naive_bool = F

#mat_gamma = mat_gamma[c(2,3),]
# run over different values of gamma's: 1:nrow(mat_gamma)
# param_n_sim * time per run * nrow(mat_gamma)
for ( k in c(1 : nrow(mat_gamma)) ){
  print(paste0("in the outer for loop ", k))
  gamma_as=as.numeric(mat_gamma[k, c(1:dim_x)])
  gamma_ns=as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
  gamma_pro=gamma_pro
  start_time <- Sys.time()
  

  if(only_naive_bool == TRUE){
    EM_and_matching = simulate_data_run_EM_and_match(return_EM_PS = FALSE, index_set_of_params=k,
       gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, misspec_PS=misspec_PS,
       misspec_outcome_funcform=FALSE, U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, 
       match_and_reg_watch_true_X=FALSE, param_n=param_n, param_n_sim=param_n_sim,
       iterations=iterations, epsilon_EM = epsilon_EM, caliper=caliper, epsilon_1_GPI=epsilon_1_GPI,
       match_on = match_on, mu_x_fixed=mu_x_fixed, x_as=mat_x_as[k,], only_naive_bool=only_naive_bool)
    
    mat_SACE_estimators = EM_and_matching$mat_param_estimators
    # nrow = nrow(mat_gamma) 
    df_parameters = matrix(rep(as.numeric(mat_gamma[k,]), each=nrow(mat_SACE_estimators)), nrow=nrow(mat_SACE_estimators))
    colnames(df_parameters) = colnames(mat_gamma)
    mat_SACE_estimators = data.frame(mat_SACE_estimators, df_parameters)

    list_all_mat_SACE_estimators[[k]] = mat_SACE_estimators
    list_all_CI[[k]] = EM_and_matching$CI_mat
    next()
  }
  
  # sim1
  if(sim == 1){
    EM_and_matching = simulate_data_run_EM_and_match(index_set_of_params=k, 
       gamma_as, gamma_ns, gamma_pro, 
       misspec_PS, param_n, param_n_sim, iterations = iterations, epsilon_EM = epsilon_EM,
       caliper, epsilon_1_GPI, match_on = match_on,
       mu_x_fixed=mu_x_fixed, x_as=mat_x_as[k,])
    #size_EM_and_matching=1
  }
  # sim2
  if(sim == 2){
    EM_and_matching = simulate_data_run_EM_and_match(return_EM_PS = FALSE, index_set_of_params=k,
        gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, misspec_PS=misspec_PS,
        misspec_outcome_funcform=FALSE, U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, 
        match_and_reg_watch_true_X=FALSE, param_n=param_n, param_n_sim=param_n_sim,
        iterations=iterations, epsilon_EM = epsilon_EM, caliper=caliper, epsilon_1_GPI=epsilon_1_GPI,
        match_on = match_on, mu_x_fixed=mu_x_fixed, x_as=mat_x_as[k,], only_naive_bool=only_naive_bool)
    #size_EM_and_matching=1
  }
  
  mat_SACE_estimators = EM_and_matching[[1]]
  # nrow = nrow(mat_gamma) 
  df_parameters = matrix(rep(as.numeric(mat_gamma[k,])
                         #, each = (param_n_sim + length(param_measures)))
                          , each = nrow(mat_SACE_estimators))
                         #, nrow = (param_n_sim + length(param_measures))
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
  
  # TODO in list_all_diff_distance_aspr_asas and list_all_std_mean_diff   
  # when diff_distance_aspr_asas is positive, the matches between as to as are closer in X, 
  #then matches between as to protected, since its the abs diff between as to pro matches - abs diff between as to as matches 
  list_all_diff_distance_aspr_asas[[k]] = EM_and_matching[["mat_diff_distance_aspr_asas"]]
  rownames(list_all_diff_distance_aspr_asas[[k]]) = paste0("s", k, rownames(list_all_diff_distance_aspr_asas[[k]]))
  list_all_std_mean_diff[[k]] = EM_and_matching[["mean_list_std_mean_diff"]]
  rownames(list_all_std_mean_diff[[k]]) = paste0("s", k, rownames(list_all_std_mean_diff[[k]]))
  list_all_means_by_subset[[k]] = EM_and_matching[["mean_list_means_by_subset"]]
  rownames(list_all_means_by_subset[[k]]) = paste0("s", k, rownames(list_all_means_by_subset[[k]]))
  list_all_EM_not_conv[[k]] = EM_and_matching[["list_EM_not_conv"]]
  if(! is_empty(list_all_EM_not_conv[[k]])){names(list_all_EM_not_conv[[k]]) = paste0("s", k, names(list_all_EM_not_conv[[k]]))}
  list_all_BCclpr[[k]] = EM_and_matching[["list_BCclpr"]]
  #rownames(list_all_EM_not_conv[[k]]) = paste0("s", k, rownames(list_all_EM_not_conv[[k]]))
  print(paste0("sim is ", sim))

  end_time <- Sys.time()
  print(paste0("in the end of outer for loop ", k, ", ", difftime(end_time, start_time)))
}
########################################################################

sum(abs(list.rbind(list_all_BCclpr[[1]])$trt_added_by_ties))
sum(abs(list.rbind(list_all_BCclpr[[2]])$trt_added_by_ties))
save(list_all_BCclpr, file = "list_all_BCclpr.RData")

dsfsdfsd# save lists from simulation ####
########################################################################
save(list_all_mat_SACE_estimators, file = "list_all_mat_SACE_estimators.RData")
save(list_all_WLS_NOint_regression_estimators, file = "list_all_WLS_NOint_regression_estimators.RData")
save(list_all_WLS_YESint_regression_estimators, file = "list_all_WLS_YESint_regression_estimators.RData")
save(list_all_OLS_NOint_regression_estimators, file = "list_all_OLS_NOint_regression_estimators.RData")
save(list_all_OLS_YESint_regression_estimators, file = "list_all_OLS_YESint_regression_estimators.RData")
save(list_all_CI, file = "list_all_CI.RData")
save(list_all_EM_coeffs, file = "list_all_EM_coeffs.RData")
save(list_all_excluded_included_matching, file = "list_all_excluded_included_matching.RData")
save(list_all_repeated_as_and_pro, file = "list_all_repeated_as_and_pro.RData")
save(list_all_matched_units, file = "list_all_matched_units.RData")
# TODO list_all_diff_distance_aspr_asas is already in mat_SACE_estimators,
save(list_all_diff_distance_aspr_asas, file = "list_all_diff_distance_aspr_asas.RData")
save(list_all_std_mean_diff, file = "list_all_std_mean_diff.RData")
save(list_all_means_by_subset, file = "list_all_means_by_subset.RData")
save(list_all_EM_not_conv, file = "list_all_EM_not_conv.RData")
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
change_rownames_to_LETTERS_by_param_set = function(df, mat_params=mat_gamma,
                         num_of_param_measures=num_of_param_measures_per_param_set){
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
pi_from_mat_all_estimators = subset(mat_all_estimators, select = grep("pi_", colnames(mat_all_estimators))) %>% round(3)
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
# initial summaries
means_by_subset_sum = mat_all_means_by_subset[grep("S1|mean_as", 
     rownames(mat_all_means_by_subset), ignore.case = F) , ] %>% round(3)
means_by_subset_sum = means_by_subset_sum[-grep("_approx", rownames(means_by_subset_sum)),]
#pi_from_mat_all_estimators
ctr_as_matched_after = sum(mat_all_repeated_as_and_pro["s1repT_S1",-c(29,30)] * c(c(1:14), c(1:14)))
mat_all_matched_units_S1 = data.frame(t(mat_all_matched_units[grep("F_S1|T_S1", rownames(mat_all_matched_units)),]))
ctr_as_matched_after - mat_all_matched_units_S1["ctr_asAfS1","s1repT_S1"]

pis = mat_all_estimators[grep("_mean",rownames(mat_all_estimators)),grep("pi",colnames(mat_all_estimators))] %>% round(3)
pis = data.frame(pi_as=pis$pi_as, pi_pro=pis$pi_pro, pi_ns=pis$pi_ns)
#mat_pis

saveRDS(pi_from_mat_all_estimators, file = "pi_from_mat_all_estimators.rds")
saveRDS(means_by_subset_sum, file = "means_by_subset_sum.rds")
saveRDS(mat_all_matched_units_S1, file = "mat_all_matched_units_S1.rds")
########################################################################

#TODO NAIVE estimators ####
########################################################################
if(only_naive_bool == TRUE){
  scenarios = rownames(mat_all_estimators) %>% substr(1,1) %>% unique
  length(list_all_CI) == length(scenarios)
  list_coverage = list(); coverage_mat = NULL
  tmp = mat_all_estimators
  for (i in 1:length(scenarios)) {
    # i is the index for the scenario, A, B and C
    list_coverage[[i]] = calculate_coverage(list_all_CI[[i]])
    coverage = list_coverage[[i]]$coverage
    rownames(coverage) = paste0(scenarios[i], "_coverage")
    coverage_mat = rbind(coverage_mat, coverage)
  }
  if(only_naive_bool == TRUE){
  est_mean = tmp[-grep("_med", rownames(tmp)), setdiff(grep("*\\_est$|SACE",colnames(tmp)), grep("^pi_",colnames(tmp)))]
  est_se = tmp[,setdiff(grep("*\\_est_se$",colnames(tmp)), grep("^pi_",colnames(tmp)))]
  est_se = est_se[grep("mean", rownames(est_se)),]
  colnames(est_se) = setdiff(colnames(est_mean), "SACE")
  rownames(est_se) = paste0(scenarios, "_se_est")
  #scenario_ind = which(substr(rownames(est_mean),1,1)==scenarios[i])
  colnames(coverage_mat) = mgsub(colnames(coverage_mat), c("Coverage_naive_without_matching",
       "Coverage_survivors_naive_without_matching"), setdiff(colnames(est_mean), "SACE")) 
  naive_estimators_sum = rbind(est_mean, data.frame(SACE = 101, est_se), 
        data.frame(SACE = 101, coverage_mat)) %>% round(3) # rbind.fill
  naive_estimators_sum[-grep("MSE", rownames(naive_estimators_sum)),] = round(naive_estimators_sum[-grep("MSE", rownames(naive_estimators_sum)),], 3)
  naive_estimators_sum = naive_estimators_sum[order(rownames(naive_estimators_sum)),]
  rownames(naive_estimators_sum) = mgsub(rownames(naive_estimators_sum), c("_mean", "_sd"), c("_est", "_sd_emp"))
  naive_estimators_sum = data.frame(Scenario = substr(rownames(naive_estimators_sum),1,1), 
      Value = substr(rownames(naive_estimators_sum),3,100), naive_estimators_sum)
  target <- paste0(rep(scenarios, each = length(c("est", "sd_emp", "se_est", "MSE", "coverage"))),
                   "_", rep(c("est", "sd_emp", "se_est", "MSE", "coverage"), length(scenarios)))
  naive_estimators_sum = naive_estimators_sum[match(target, rownames(naive_estimators_sum)),]
  save(naive_estimators_sum, file = "naive_estimators_sum.RData")
  }
  
  # apply this even when only_naive_bool = FALSE if we need to compute coverage and MSE for naives after actuall simulation.
  if(only_naive_bool == FALSE){
    coverage_mat = subset(coverage_mat, select = grep("Naive", colnames(coverage_mat), ignore.case = T))
    est_se = tmp[,setdiff(grep("*\\_est_se$",colnames(tmp)), grep("^pi_",colnames(tmp)))]
    est_se = est_se[grep("mean", rownames(est_se)),]
    rownames(est_se) = paste0(scenarios, "_se_est")
    est_mean = tmp[grep("mean|sd|MSE", rownames(tmp)), 
                   setdiff( grep("SACE|naive", colnames(tmp), ignore.case = T), grep(paste(colnames(est_se), collapse="|") , colnames(tmp)) )] 
    #scenario_ind = which(substr(rownames(est_mean),1,1)==scenarios[i])
    colnames(coverage_mat) = mgsub(colnames(coverage_mat), c("Coverage_naive_without_matching", "Coverage_survivors_naive_without_matching"), 
                         colnames(est_mean)[-grep("SACE", colnames(est_mean), ignore.case = T)]) 
    colnames(est_se) = colnames(coverage_mat)
    naive_estimators_sum = rbind(est_mean[,-grep("SACE", colnames(est_mean), ignore.case = T)], 
                                 est_se, coverage_mat) %>% round(3) # rbind.fill
    #naive_estimators_sum[-grep("MSE", rownames(naive_estimators_sum)),] = round(naive_estimators_sum[-grep("MSE", rownames(naive_estimators_sum)),], 2)
    naive_estimators_sum = naive_estimators_sum[order(rownames(naive_estimators_sum)),]
    naive_estimators_sum = merge(naive_estimators_sum, round(est_mean,3), by=0, all.x = TRUE) %>% subset(select = -SACE_conditional)
    sum(naive_estimators_sum$most_naive_est.x - naive_estimators_sum$most_naive_est.y, na.rm = T)
    naive_estimators_sum$Row.names <- mgsub( naive_estimators_sum$Row.names, c("_mean", "_sd"), c("_est", "_sd_emp"))
    colnames(naive_estimators_sum) = mgsub(colnames(naive_estimators_sum), ".x", "")
    naive_estimators_sum = subset(naive_estimators_sum, select = -grep("*\\.y$", colnames(naive_estimators_sum)))
    save(naive_estimators_sum, file = "naive_estimators_sum.RData")
  }
}
########################################################################


########################################################################
# TABLE DESIGN
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
save(final_tables_general, file = "final_tables_general.RData")
save(final_tables_crude, file = "final_tables_crude.RData")

print(means_by_subset_sum %>% xtable(digits=c(3),
        caption = "no misspec, small pi pro, 3X"), size="\\fontsize{11pt}{11pt}\\selectfont")

# TODO print
print(adjustments_for_tables_before_present(final_tables_crude$`_S1`) %>% xtable(digits=c(2),
    caption = paste0("True model with interactions. ", ifelse(misspec_PS==1, "misspecefication: Unobserved in PS and Y",
        ifelse(misspec_PS==2, "misspecefication: functional form of PS, no misspecification in Y", "No misspecification")), 
                     ". ", (dim_x-1) ," X's, Delta Method, ",
                     param_n_sim, " replications, ", "N=", param_n, ", pA = ", prob_A, ". OLS clustered SE")),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
print(pis %>% xtable(digits=2), size="\\fontsize{15pt}{15pt}\\selectfont", include.rownames=F)
if(misspec_PS != 0){
  if(misspec_PS == 1){
    # beta interactions of U in Y ~ X + U is: # c(5,2), c(2,1)) # c(1,1), c(1,1)
    print(cbind(beta_U1=c(1,1), beta_U2= c(1,1)) %>% 
            xtable, size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)
  }
  print(cbind(U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log) %>% 
          xtable, size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F) 
}

print(adjustments_for_tables_before_present(final_tables_crude$`_S1`) %>% xtable(digits=c(2),
  caption = paste0("True model with interactions. ", ifelse(misspec_PS==1, "misspecefication: Unobserved in PS and Y",
      ifelse(misspec_PS==2, "misspecefication: functional form of PS, no misspecification in Y", "No misspecification")), 
                   ". ", (dim_x-1) ," X's, Delta Method, ",
                   param_n_sim, " replications, ", "N=", param_n, ", pA = ", prob_A, ". OLS clustered SE")),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
print(pis %>% xtable(digits=2), size="\\fontsize{15pt}{15pt}\\selectfont", include.rownames=F)
print(mat_gamma[,c(1,2,dim_x+1,dim_x+2)] %>% xtable(digits=3), size="\\fontsize{10pt}{10pt}\\selectfont", include.rownames=F)
print(betas_GPI %>% xtable, size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)


# for summary of each scenario, including EM coeffs
#print(mat_gamma[,c(1,2,7,8)] %>% xtable(digits=3), size="\\fontsize{10pt}{10pt}\\selectfont", include.rownames=F)
#print(pis %>% xtable(digits=2), size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)
print(round(mat_x_as,3) %>% xtable(digits=3), size="\\fontsize{10pt}{10pt}\\selectfont", include.rownames=F)
print(subset(summary_EM_coeffs, select = c(mean, parameter, diff, perc)) %>% xtable(digits=3,
  caption = "Check EM coefficients in order to understand why Ding is poorly perform when pi pro is small (pi pro is 0.1)."), 
      size="\\fontsize{12pt}{12pt}\\selectfont")

###################################################
print(adjustments_for_tables_before_present(final_tables_general$`_S1`) %>% xtable(digits=c(2)),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
no_mis = adjustments_for_tables_before_present(final_tables_general$`_S1`)
ff_mis = adjustments_for_tables_before_present(final_tables_general$`_S1`) %>% subset(select=-c(1:2))
colnames(ff_mis) = paste0("mis_", colnames(ff_mis))
print(cbind(no_mis, ff_mis) %>% 
  filter(!Estimator %in% c("HL No", "BC Yes", "BC caliper inter Yes", "BC inter Yes", "HL Yes")) %>% 
  xtable(digits=c(2)), size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)

print(adjustments_for_tables_before_present(final_tables_crude$`_S1`) %>% xtable(digits=c(2)),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
print(naive_estimators_sum %>% xtable(digits=c(2)),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
###################################################
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
# for the matrix with many parameters
round_table = round(mat_all_estimators, 2)
View(round_table)
pdf(file = "q.pdf")
grid.table(t(round_table), theme=ttheme_minimal(base_size = 3.2))
# ttheme_default(base_size = 12, base_colour = "black", base_family = "",
# parse = FALSE, padding = unit(c(4, 4), "mm"), ...)
dev.off()
########################################################################


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
df_final = data.table(df1)[, lapply(.SD, replace_inf_with_max), by=id]
####################################################################


####################################################################
n.sample <- 100000
x <- rnorm(n.sample)
beta0 <- -5
beta1 <- log(1.5)
pr <- expit(beta0 + beta1*x)
hist(pr)
mean(pr)
y <- rbinom(n.sample, 1, pr)
mean(y)

#all.equal(pm1, pm2, check.attributes = FALSE)
#identical()
####################################################################