library(data.table); library(plyr); library(dplyr); library(rlist); library(nnet); library(MASS); library(rockchalk); library(locfit)
setwd("~/A matching framework for truncation by death problems")
source("Simulations_studies/sim_gamma_and_pi_CPSR/sim_check_pis_and_covariates.R")
source("Simulations_studies/sim_DGM_and_simulations/simulation_run_CPSR.R")
source("Ding_Lu_EM/Sequencial_logistic_regressions/EM_2log_CPSR.R") 

#############################################################################################
# treatment probability
prob_A = 0.5
# parameters for simulating X
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 4; cont_x = 3; categ_x = 0; vec_p_categ = rep(0.5, categ_x); dim_x_misspec = 2
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)
#############################################################################################

#############################################################################################
# misspec parameters (for PS model and Y model):
# misspec_PS: 0 <- NO, 1:only PS model, 2: PS model (possibly also Y)
misspec_PS = 2 # 0: no misspec of PS model # 2: PS functional form misspecification
funcform_factor_sqr=5; funcform_factor_log=-5 # funcform_factor_sqr=-3; funcform_factor_log=3
mean_x_misspec = rep(0.5, dim_x_misspec)
misspec_outcome = 0
#############################################################################################

#############################################################################################
# CPSR parameter 
xi = 0
xi_est = xi

# EM convergence parameters
iterations = 200; epsilon_EM = 10^-5
#############################################################################################

###############################################################################################
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

#################################################################################################################
# mat_gamma[,c(1,2,dim_x+1,dim_x+2)]
extract_pis_lst = extract_pis_from_scenarios(nn=500000, xi=xi, misspec_PS=0, two_log_models=T); mat_pis_per_gamma = extract_pis_lst$mat_pis
mat_pis_per_gamma
##################################################################################################################

##################################################################################################################
caliper = 0.25; match_on = "O11_posterior_ratio" 
mu_x_fixed = FALSE; mat_x_as; x_as = mat_x_as[1,]
##################################################################################################################

##################################################################################################################
# one_log_true_ah # two_log_EM_initial_ah # one_log_EM # two_log_EM 
log_EM = simulate_data_run_EM_and_match(only_EM_bool=TRUE, return_EM_PS=FALSE, index_set_of_params=1,
                                            gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, xi=xi, xi_est=xi_est,
                                            two_log_models=TRUE, two_log_est=FALSE, 
                                            misspec_PS=misspec_PS, misspec_outcome=0,
                                            funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                                            param_n=2000, param_n_sim=150,
                                            iterations=iterations, epsilon_EM=epsilon_EM, caliper=0.25,
                                            match_on=match_on, mu_x_fixed=mu_x_fixed, x_as=NULL)

apply(list.rbind(log_EM$list_beta_S0), 2, mean)
coeff_ah = list.rbind(log_EM$list_coeff_ah)
apply(coeff_ah, 2, mean)
coeff_pro = list.rbind(log_EM$list_coeff_pro)
apply(coeff_pro, 2, mean)
#################################################################################################################

##################################################################################################################
# summary of EM after one run, and DL estimators

res_lst = list(); res_mat = matrix(nrow = 25,ncol = 2)
for ( k in c(1 : nrow(mat_gamma)) ){  
  for(i in 1:25){
    print(k); print(i)
    # assign 0's to gamma_ns and add coefficients names 
    gamma_ns = rep(0, dim_x); gamma_ah = as.numeric(mat_gamma[k, c(1:dim_x)]); gamma_pro =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
    
    sum_EM = simulate_data_run_EM_and_match(only_EM_bool=FALSE, return_EM_PS=TRUE, index_set_of_params=1,
                                            gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, xi=xi, xi_est=xi_est,
                                            two_log_models=TRUE, two_log_est=FALSE, 
                                            misspec_PS=misspec_PS, misspec_outcome=0,
                                            funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                                            param_n=2000, param_n_sim=1,
                                            iterations=iterations, epsilon_EM=epsilon_EM, caliper=0.25,
                                            match_on=match_on, mu_x_fixed=mu_x_fixed, x_as=NULL)
    
    res_mat[i,] = c(sum_EM$SACE, sum_EM$DL_MA)
    #sum_EM$PS_true_EM_compr %>% View
    #mean(sum_EM$true_x_PS$X_sqr); mean(sum_EM$true_x_PS$X_log); sum_EM$mean_by_g
    #mean(sum_EM$PS_true_EM_compr$prob_as); mean(sum_EM$PS_true_EM_compr$EMest_p_as);
    #mean(sum_EM$PS_true_EM_compr$diff); mean(abs(sum_EM$PS_true_EM_compr$diff)); sd(sum_EM$PS_true_EM_compr$diff)
    #sum_EM$pis; sum_EM$pis_est; sum_EM$EM_coeffs
  }
  res_lst[[k]] = res_mat
}
for ( k in c(1 : nrow(mat_gamma)) ){ print(apply(res_lst[[k]], 2, mean)) }
#################################################################################################################