setwd("~/A matching framework for truncation by death problems")
source("Simulations_studies/sim_DGM_and_simulations/simulation_run_CPSR.R")
source("Ding_Lu_EM/Sequencial_logistic_regressions/EM_2log_CPSR.R") 

#############################################################################################
# treatment probability
prob_A = 0.5
# parameters for simulating X
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 6; cont_x = 5; categ_x = 0; vec_p_categ = rep(0.5, categ_x); dim_x_misspec = 2
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x)
#############################################################################################

#############################################################################################
# misspec parameters (for PS model and Y model:
# misspec_PS: 0 <- NO, 1:only PS model, 2: PS model (possibly also Y)
funcform_mis_out = FALSE
funcform_factor_sqr=-3; funcform_factor_log=3
mean_x_misspec = rep(0.5, dim_x_misspec)
misspec_PS = 0 # 0: no misspec of PS model # 2: PS functional form misspecification
#############################################################################################

#############################################################################################
# CPSR parameter 
xi = 0.5
# for now, xi for DGM and xi for estimation are the same
# xi_est = (or call THIS eta, as in DL EM)

# EM convergence parameters
iterations = 200; epsilon_EM = 10^-6
#############################################################################################

###############################################################################################
betas_GPI = as.matrix(rbind(c(22,5,2,1,3,5), c(20,3,3,0,1,3))) # cont_x=5
rownames(betas_GPI) = c("beta_treatment", "beta_control")
###############################################################################################

##########################################################
# correlation structure between PO'
var_GPI = as.matrix(rbind(1, 1))
rownames(var_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4 #0.4 #1
##########################################################

#################################################################################################################
mat_gamma = matrix(c(
  c(-0.05, rep(0.16, dim_x-1)), c(-0.4, rep(-0.25, dim_x-1)) 
  ,c(-0.5, rep(1.5, dim_x-1)), c(-0.25, rep(0.25, dim_x-1))) ,nrow = 2, byrow = T)
gamma_ns = rep(0, dim_x)
gamma_ah = as.numeric(mat_gamma[1, c(1:dim_x)])
gamma_pro =  as.numeric(mat_gamma[1, (dim_x+1): (2*dim_x)])
##################################################################################################################

##################################################################################################################
caliper = 0.25; match_on = "O11_posterior_ratio" 
mu_x_fixed = FALSE; mat_x_as; x_as = mat_x_as[1,]
##################################################################################################################

##################################################################################################################
# one_log_true_ah # two_log_EM_initial_ah # one_log_EM # two_log_EM 
two_log_EM = simulate_data_run_EM_and_match(only_EM_bool=TRUE, return_EM_PS=FALSE, index_set_of_params=1,
                                            gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, xi=xi,
                                            two_log_models=TRUE, two_log_est=TRUE, 
                                            misspec_PS=misspec_PS, funcform_mis_out=FALSE,
                                            funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                                            param_n=100, param_n_sim=2,
                                            iterations=iterations, epsilon_EM=epsilon_EM, caliper=0.25,
                                            match_on=match_on, mu_x_fixed=mu_x_fixed, x_as=NULL)
save(two_log_EM, file = "two_log_EM.RData")

apply(list.rbind(two_log_EM$list_beta_S0), 2, mean)
coeff_ah = list.rbind(two_log_EM$list_coeff_ah)
apply(coeff_ah, 2, mean)
coeff_pro = list.rbind(two_log_EM$list_coeff_pro)
apply(coeff_pro, 2, mean)
#################################################################################################################