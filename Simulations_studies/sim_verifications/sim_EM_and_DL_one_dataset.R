simulate_data_EM_and_DL = function(gamma_ah, gamma_pro, gamma_ns, xi, xi_est,
                                   two_log_models_DGM=TRUE, two_log_est_EM=FALSE,
                                   misspec_PS, misspec_outcome=0, transform_x=0, 
                                   funcform_factor_sqr, funcform_factor_log, 
                                   param_n, iterations, epsilon_EM = 0.001,
                                   mu_x_fixed=FALSE, x_as,
                                   betas_GPI, var_GPI, rho_GPI_PO, only_mean_x_bool=FALSE){
  
  X_sub_cols = paste0("X", c(1:(dim_x)))
  
  # simulate data
  list_data_for_EM_and_X = simulate_data_func(gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, 
                                xi=xi, two_log_models_DGM=two_log_models_DGM, param_n=param_n,
                                misspec_PS=misspec_PS, misspec_outcome=misspec_outcome, transform_x=transform_x,
                                funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                                betas_GPI=betas_GPI, var_GPI=var_GPI, rho_GPI_PO=rho_GPI_PO, only_mean_x_bool=FALSE)

  data_for_EM = list_data_for_EM_and_X$dt
  mean_by_g = list_data_for_EM_and_X$mean_by_g
  x = list_data_for_EM_and_X$x_obs; x_PS = data.frame(list_data_for_EM_and_X$x_PS)
  x_outcome = data.frame(list_data_for_EM_and_X$x_outcome)
  OBS_table = list_data_for_EM_and_X$OBS_table
  pis = list_data_for_EM_and_X$pis; pis_est = list_data_for_EM_and_X$pis_est
  vec_OBS_table = t(c(OBS_table)); colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
  # "real" SACE parameter
  SACE = list_data_for_EM_and_X$true_SACE
  
  # EM and PS estiamtion
  # If beta_S0=NULL, employ two logistic regressions during the EM
  if(two_log_est_EM == FALSE){
    #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
    fit_S0_in_A0 = glm(as.formula(paste0("S ~ ",paste(X_sub_cols[-1], collapse="+"))), data=filter(data_for_EM, A==0), family="binomial")
    beta_S0 = fit_S0_in_A0$coefficients
  }else{beta_S0=NULL}
  
  start_timeDing <- Sys.time()
  est_ding_lst = xi_2log_PSPS_M_weighting(Z=data_for_EM$A, D=data_for_EM$S,
                                          X=as.matrix(subset(data_for_EM, select = 
                                                               grep(paste(X_sub_cols[-1], collapse="|"), colnames(data_for_EM)))), Y=data_for_EM$Y, 
                                          xi_est=xi_est, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL, 
                                          iter.max=iterations, error0=epsilon_EM)
  end_timeDing <- Sys.time()
  print(paste0("Ding EM lasts ", difftime(end_timeDing, start_timeDing)))
  
  coeff_ah = est_ding_lst$beta.ah ; coeff_pro = est_ding_lst$beta.c
  EM_coeffs = rbind(coeff_ah, coeff_pro)
  PS_est = est_ding_lst$ps.score
  data_with_PS = data.table(data_for_EM, PS_est)
  
  # Ding estimator DL plain estimator and model assisted estimator
  DING_est = est_ding_lst$AACE
  DING_model_assisted_est_ps = est_ding_lst$AACE.reg
  
  # calculate weights: O11_prior_ratio, O11_posterior_ratio W_1_as, and W_1_as_true
  weights_lst = add_PS_weights_func(data_with_PS=data_with_PS, pis=pis, pis_est=pis_est)
  O11_prior_ratio = weights_lst$O11_prior_ratio
  O11_prior_ratio_true = weights_lst$O11_prior_ratio_true
  data_with_PS = weights_lst$data_with_PS 
  
  # EM summary
  pis = data.frame(pis)
  PS_true_EM_compr = subset( data_with_PS, select = grep("^id$|^g$|prob.|EM|O11_posterior_ratio|O11_prior_ratio|W_1_as",colnames(data_with_PS)) )
  PS_true_EM_compr = rapply(object = PS_true_EM_compr, f = round, classes = "numeric", how = "replace", digits = 3)
  
  PS_true_EM_compr = data.frame( id = PS_true_EM_compr$id, g = PS_true_EM_compr$g,
                                 prob_as = PS_true_EM_compr$prob_as, EMest_p_as=PS_true_EM_compr$EMest_p_as, diff = PS_true_EM_compr$prob_as - PS_true_EM_compr$EMest_p_as,
                                 prob_har = PS_true_EM_compr$prob_har, EMest_p_har=PS_true_EM_compr$EMest_p_har, 
                                 prob_ns = PS_true_EM_compr$prob_ns, EMest_p_ns=PS_true_EM_compr$EMest_p_ns,
                                 prob_pro = PS_true_EM_compr$prob_pro, EMest_p_pro=PS_true_EM_compr$EMest_p_pro,
                                 O11_prior_ratio_true = PS_true_EM_compr$O11_prior_ratio_true, O11_prior_ratio_est = PS_true_EM_compr$O11_prior_ratio,
                                 O11_posterior_ratio_true = PS_true_EM_compr$O11_posterior_ratio_true, O11_posterior_ratio_est = PS_true_EM_compr$O11_posterior_ratio,
                                 W_1_as_true = PS_true_EM_compr$W_1_as_true, W_1_as_est = PS_true_EM_compr$W_1_as)
  
  return(list(data_with_PS=data_with_PS, PS_true_EM_compr=PS_true_EM_compr, true_x_PS=x_PS,
          pis=pis, pis_est=pis_est, EM_coeffs=EM_coeffs, gamma = c(gamma_ah=gamma_ah,gamma_ah=gamma_ah),
          O11_prior_ratio_true=O11_prior_ratio_true, O11_prior_ratio=O11_prior_ratio, OBS_table=OBS_table, 
          beta_S0=beta_S0, error=est_ding_lst$error, mean_by_g=mean_by_g,
          SACE=SACE, DL=DING_est, DL_MA=DING_model_assisted_est_ps,
          w1a=est_ding_lst$w1a, w0a=est_ding_lst$w0a, w1a_all=est_ding_lst$w1a_all, w0a_all=est_ding_lst$w0a_all,
          weighted.Y1 = est_ding_lst$weighted.Y.a1, weighted.Y1.adj=est_ding_lst$weighted.Y1a, weighted.ra=est_ding_lst$weighted.ra,
          weighted.Y1.unb=est_ding_lst$weighted.Y.a1.unb, weighted.Y1.adj.unb=est_ding_lst$weighted.Y1a.unb, weighted.ra.unb=est_ding_lst$weighted.ra.unb,
          em_NOT_conv = est_ding_lst$iter == iterations+1 & est_ding_lst$error>=epsilon_EM))

}



