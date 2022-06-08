library(data.table); library(plyr); library(dplyr); library(rlist); library(nnet); library(MASS); library(rockchalk); library(locfit)
setwd("~/A matching framework for truncation by death problems")
source("Simulations_studies/sim_parameters_and_pis/sim_set_parameters.R")
source("Simulations_studies/sim_parameters_and_pis/sim_check_pis_and_covariates.R")
source("Simulations_studies/sim_verifications/sim_EM_and_DL_one_dataset.R")
source("Simulations_studies/sim_DGM_and_simulations/DGM_CPSR.R")
source("Ding_Lu_EM/Sequencial_logistic_regressions/EM_2log_CPSR.R") 

#############################################################################################
# treatment probability
prob_A = 0.5
# parameters for simulating X
#@@@@@@@@@@@@ dim_x includes intercept @@@@@@@@@@@@@@@
dim_x = 2; cont_x = 1
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x) # var_x = rep(1, cont_x)
#############################################################################################

#############################################################################################
# misspec parameters (for PS model and Y model):
# misspec_PS: 0 <- NO, 1:only PS model, 2: PS model (possibly also Y)
misspec_PS = 2 # 0: no misspec of PS model # 2: PS functional form misspecification
funcform_factor_sqr=5; funcform_factor_log=-5 # thesis: funcform_factor_sqr=-3; funcform_factor_log=3
misspec_outcome = 0
#############################################################################################

#############################################################################################
# CPSR parameter 
xi = 0
xi_est = 0  # xi
# EM convergence parameters
iterations = 200; epsilon_EM = 10^-5
#############################################################################################

##########################################################
# correlation structure between POs
var_GPI = as.matrix(rbind(1, 1))
rownames(var_GPI) = c("var_treatment", "var_control")
rho_GPI_PO = 0.4 #0.4 
##########################################################

#################################################################################################################
beta_and_gamma = set_parameters_func(dim_x=dim_x, AX_interactions=T)
# beta
betas_GPI = beta_and_gamma$betas_GPI
# gamma
mat_gamma = beta_and_gamma$mat_gamma
# extract pis
extract_pis_lst = extract_pis_from_scenarios(nn=300000, mat_gamma=mat_gamma, xi=xi, misspec_PS=2, two_log_models=T)
mat_pis_per_gamma = extract_pis_lst$mat_pis; mat_pis_per_gamma
##################################################################################################################

##################################################################################################################
# summary of EM after one run, and DL estimators ####
set.seed(101)
k=2
gamma_ns = rep(0, dim_x)
gamma_ah = as.numeric(mat_gamma[k, c(1:dim_x)])
gamma_pro =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])

sum_EM = simulate_data_EM_and_DL(gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, xi=xi, xi_est=xi_est,
                                  two_log_models=TRUE, two_log_est_EM=FALSE,
                                  misspec_PS=2, misspec_outcome=misspec_outcome, transform_x=0,
                                  funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, 
                                  param_n=10000, iterations=iterations, epsilon_EM=epsilon_EM,
                                  mu_x_fixed=FALSE, x_as=NULL)

SACE_end_DL_est = c(sum_EM$SACE, sum_EM$DL_MA, sum_EM$DL)
    
# O(1,1) 
dfrm_A1_S1 = filter(sum_EM$data_with_PS, A==1&S==1)

# comparison ####
max(abs(sum_EM$w1a- dfrm_A1_S1$W_1_as))
max(abs(sum_EM$w1a_all- sum_EM$data_with_PS$W_1_as))
mean(sum_EM$w1a); mean(dfrm_A1_S1$W_1_as)
unique(dfrm_A1_S1$O11_prior_ratio); mean(dfrm_A1_S1$O11_posterior_ratio)
plot(sum_EM$w1a, dfrm_A1_S1$W_1_as)
mean(sum_EM$w1a_all); mean(sum_EM$data_with_PS$W_1_as)
unique(sum_EM$w0a); unique(sum_EM$w0a_all)

# plots- scatter and histograms ####
#pdf(file= "~/A matching framework for truncation by death problems/Plots2.pdf" )

# weights ####
# O(1,1)
plot(dfrm_A1_S1$W_1_as_true, dfrm_A1_S1$W_1_as, cex.main=0.7, main = paste0("W_1_as, ", " O(1,1)", ". set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
res = dfrm_A1_S1$W_1_as-dfrm_A1_S1$W_1_as_true
plot(dfrm_A1_S1$W_1_as_true, res, cex.main=0.7, main = paste0("W_1_as, res:est-true", " O(1,1)", ". set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
# principal scores
plot(dfrm_A1_S1$prob_as, dfrm_A1_S1$EMest_p_as, cex.main=0.7, main = paste0("set", k, ",  O(1,1), as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
plot(dfrm_A1_S1$prob_pro, dfrm_A1_S1$EMest_p_pro, cex.main=0.7, main = paste0("set", k, ",  O(1,1), as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
# weights
hist(dfrm_A1_S1$W_1_as_true, xlim=c(0,1), col="blue", main = paste0("W_1_as, O(1,1)", ". blue-true, green-est"), breaks = 50)
hist(dfrm_A1_S1$W_1_as, add=T, col=rgb(0, 1, 0, 0.5), breaks = 50)

# full dataset
datfrm = sum_EM$PS_true_EM_compr
# weights
plot(datfrm$W_1_as_true, datfrm$W_1_as_est, main = paste0("W_1_as, set", k, ". as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ",")))
abline(a=0,b=1, lwd=3)
# principal scores
plot(datfrm$prob_as, datfrm$EMest_p_as, main = paste0("set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
plot(datfrm$prob_pro, datfrm$EMest_p_pro, main = paste0("set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)

# by stratum ####
g_options = unique(datfrm$g)
for (j in 1:length(g_options)){
  datfrm_g = subset(datfrm, g==g_options[j])
  
  # weights
  plot(datfrm_g$W_1_as_true, datfrm_g$W_1_as_est, main = paste0("W_1_as, ", "g=", g_options[j], ". set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ",")))
  abline(a=0,b=1, lwd=3)
  hist(datfrm_g$W_1_as_true, xlim=c(0,1), col="blue", main = paste0("W_1_as, ", "g=", g_options[j], ". blue-true, green-est"), breaks = 30)
  hist(datfrm_g$W_1_as_est, add=T, col=rgb(0, 1, 0, 0.5), breaks = 30)
  
  #plot(datfrm_g$O11_posterior_ratio_true, datfrm_g$O11_posterior_ratio_est, main = paste0("O11_posterior_ratio, set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ",")))
  #abline(a=0,b=1, lwd=3)
  #abline(a=datfrm_g$O11_prior_ratio_true[1], b=0); abline(a=datfrm_g$O11_prior_ratio_est[1], b=0)
  #hist(datfrm_g$O11_posterior_ratio_true, xlim=c(0,1), col="blue", main = "O11_posterior_ratio, blue-true, green-est")
  #hist(datfrm_g$O11_posterior_ratio_est, add=T, col=rgb(0, 1, 0, 0.5) )
  
  # principal scores
  # plot(datfrm_g$prob_as, datfrm_g$EMest_p_as, main = paste0("set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
  # hist(datfrm_g$prob_as, xlim=c(0,1), col="blue", main = "prob_as, blue-true, green-est"); hist(datfrm_g$EMest_p_as, add=T, col=rgb(0, 1, 0, 0.5) )
  # plot(datfrm_g$prob_pro, datfrm_g$EMest_p_pro, main = paste0("set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
  # hist(datfrm_g$prob_pro, xlim=c(0,1), col="blue",  main = "prob_pro, blue-true, green-est"); hist(datfrm_g$EMest_p_pro, add=T, col=rgb(0, 1, 0, 0.5) )
  # plot(datfrm_g$prob_ns, datfrm_g$EMest_p_ns, main = paste0("set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
  # hist(datfrm_g$prob_ns, xlim=c(0,1), col="blue",  main = "prob_ns, blue-true, green-est"); hist(datfrm_g$EMest_p_ns, add=T, col=rgb(0, 1, 0, 0.5) )
}
#dev.off()
#################################################################################################################



