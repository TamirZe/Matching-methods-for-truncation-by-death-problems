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
dim_x = 6; cont_x = 5; categ_x = 0; vec_p_categ = rep(0.5, categ_x); dim_x_misspec = 2
mean_x = rep(0.5, cont_x); var_x = rep(1, cont_x) # var_x = rep(1, cont_x)
#############################################################################################

#############################################################################################
# misspec parameters (for PS model and Y model):
# misspec_PS: 0 <- NO, 1:only PS model, 2: PS model (possibly also Y)
misspec_PS = 2 # 0: no misspec of PS model # 2: PS functional form misspecification
funcform_factor_sqr=5; funcform_factor_log=-5 # funcform_factor_sqr=-3; funcform_factor_log=3
mean_x_misspec = rep(2, dim_x_misspec) # mean_x_misspec = rep(0.5, dim_x_misspec)
misspec_outcome = 0
#############################################################################################

#############################################################################################
# CPSR parameter 
xi = 0
xi_est = 0  # xi

# EM convergence parameters
iterations = 200; epsilon_EM = 10^-5
#############################################################################################

###############################################################################################
# beta ####
# with interactions between A and X:
betas_GPI = as.matrix(rbind(c(22,5,2,1), c(20,3,3,0))) # cont_x=3
betas_GPI = as.matrix(rbind(c(22,10,-2,10), c(20,3,3,0))) # cont_x=3
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
extract_pis_lst = extract_pis_from_scenarios(nn=300000, xi=xi, misspec_PS=0, two_log_models=T); mat_pis_per_gamma = extract_pis_lst$mat_pis
mat_pis_per_gamma
##################################################################################################################

##################################################################################################################
caliper = 0.25; match_on = "O11_posterior_ratio" 
mu_x_fixed = FALSE
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
set.seed(101)
iters_n_sim=50; res_mat = matrix(nrow = iters_n_sim, ncol = 2); NOT_conv_vec = vector(length=iters_n_sim)
res_lst <- NOT_conv_lst <- all_em_list <- list()
for ( k in c(1 : nrow(mat_gamma)) ){ 
  set_em_list = list()
  gamma_ns = rep(0, dim_x); gamma_ah = as.numeric(mat_gamma[k, c(1:dim_x)]); gamma_pro =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
  for(i in 1:iters_n_sim){
    print(k); print(i)
    #sum_EM_tmp=sum_EM
    sum_EM = simulate_data_run_EM_and_match(only_EM_bool=FALSE, return_EM_PS=TRUE, index_set_of_params=1,
                                            gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, xi=xi, xi_est=xi_est,
                                            two_log_models=TRUE, two_log_est_EM=FALSE, 
                                            misspec_PS=2, misspec_outcome=0, transform_x=0,
                                            funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                                            param_n=10000, param_n_sim=1,
                                            iterations=iterations, epsilon_EM=epsilon_EM, caliper=0.25,
                                            match_on=match_on, mu_x_fixed=mu_x_fixed, x_as=NULL)
    
    res_mat[i,] = c(sum_EM$SACE, sum_EM$DL_MA, sum_EM$DL)
    NOT_conv_vec[i] = sum_EM$real_iter_ind
    set_em_list[[i]] = sum_EM

    #sum_EM$PS_true_EM_compr %>% View
    #mean(sum_EM$true_x_PS$X_sqr); mean(sum_EM$true_x_PS$X_log); sum_EM$mean_by_g
    #mean(sum_EM$PS_true_EM_compr$prob_as); mean(sum_EM$PS_true_EM_compr$EMest_p_as);
    #mean(sum_EM$PS_true_EM_compr$diff); mean(abs(sum_EM$PS_true_EM_compr$diff)); sd(sum_EM$PS_true_EM_compr$diff)
    #sum_EM$pis; sum_EM$pis_est; sum_EM$EM_coeffs
    
    max(abs(sum_EM$w1a- sum_EM$data_with_PS[A==1&S==1,]$W_1_as))
    max(abs(sum_EM$w1a_all- sum_EM$data_with_PS$W_1_as))
    mean(sum_EM$w1a); mean(sum_EM$data_with_PS[A==1&S==1,]$W_1_as)
    plot(sum_EM$w1a, sum_EM$data_with_PS[A==1&S==1,]$W_1_as)
    mean(sum_EM$w1a_all); mean(sum_EM$data_with_PS$W_1_as)
    unique(sum_EM$w0a); unique(sum_EM$w0a_all)
    
    # plots- scatter and histograms ####
    pdf(file= "~/A matching framework for truncation by death problems/Plots2.pdf" )
    
    # O(1,1)
    dfrm_A1_S1 = filter(sum_EM$data_with_PS, A==1&S==1)
    # weights
    par(mfrow=c(1,1))
    plot(dfrm_A1_S1$W_1_as_true, dfrm_A1_S1$W_1_as, cex.main=0.7, main = paste0("W_1_as, ", " O(1,1)", ". set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
    res = dfrm_A1_S1$W_1_as-dfrm_A1_S1$W_1_as_true
    plot(dfrm_A1_S1$W_1_as_true, res, cex.main=0.7, main = paste0("W_1_as, res:est-true", " O(1,1)", ". set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
    # principal scores
    plot(dfrm_A1_S1$prob_as, dfrm_A1_S1$EMest_p_as, cex.main=0.7, main = paste0("set", k, ",  O(1,1), as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
    plot(dfrm_A1_S1$prob_pro, dfrm_A1_S1$EMest_p_pro, cex.main=0.7, main = paste0("set", k, ",  O(1,1), as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
    # weights
    hist(dfrm_A1_S1$W_1_as_true, xlim=c(0,1), col="blue", main = paste0("W_1_as, O(1,1)", ". blue-true, green-est"), breaks = 30)
    hist(dfrm_A1_S1$W_1_as, add=T, col=rgb(0, 1, 0, 0.5), breaks = 30)
    
    # full dataset
    datfrm = sum_EM$PS_true_EM_compr
    # weights
    plot(datfrm$W_1_as_true, datfrm$W_1_as_est, main = paste0("W_1_as, set", k, ". as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ",")))
    abline(a=0,b=1, lwd=3)
    # principal scores
    plot(datfrm$prob_as, datfrm$EMest_p_as, main = paste0("set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
    plot(datfrm$prob_pro, datfrm$EMest_p_pro, main = paste0("set", k, ", as: ", paste(mat_gamma[k,c(1,2)], collapse = ","), ", pro: " , paste(mat_gamma[k,c(dim_x+1,dim_x+2)], collapse = ","))); abline(a=0,b=1, lwd=3)
    
    # by stratum
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
      dev.off()
  }
  all_em_list[[k]] = set_em_list
  res_lst[[k]] = res_mat; NOT_conv_lst[[k]] = NOT_conv_vec
  res_mat = matrix(nrow = iters_n_sim, ncol = 2); NOT_conv_vec = vector(length=iters_n_sim)
}
for ( k in c(1 : nrow(mat_gamma)) ){ print(apply(res_lst[[k]], 2, mean)); print(mean(NOT_conv_lst[[k]])) }
#################################################################################################################

save(sum_EM_S1V, file = "sum_EM_S1V.Rdata")

#w = abs(rnorm(20000, mean = 1))
w = runif(n = 2000, min = 0.5, max = 1.5)
weights = w / sum(w)
smpl_weights = sample(w)
plot(weights, smpl_weights)
Y = rnorm(2000, mean = 10, sd=2)
wgt_mean1 = sum(weights*Y)
wgt_mean11 = sum(sample(weights)*Y)
wgt_mean2 = mean(w*Y) 
wgt_mean22 = mean(smpl_weights*Y)
