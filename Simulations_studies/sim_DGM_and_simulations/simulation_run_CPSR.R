# misspec_PS: 0 <- NO mis, 2: add transformations to PS model, and remain original X's in outcome model.
simulate_data_function = function(seed_num=NULL, gamma_ah, gamma_pro, gamma_ns, xi, xi_est, two_log_models=TRUE, param_n, 
                                  misspec_PS, misspec_outcome=0,
                                  funcform_factor_sqr=0, funcform_factor_log=0, only_mean_x_bool=FALSE){
  if(!is.null(seed_num)){set.seed(seed_num)}
  
  # draw covariate matrix
  x_obs <- matrix( c( rep(1,param_n), 
                      mvrnorm(param_n, mu=mean_x, Sigma = diag(var_x, cont_x))), nrow = param_n )
  colnames(x_obs) = paste0("X", c(1:dim_x))
  # x are the outcome (Y) covariates 
  x = x_obs
  # x_misspec for PS misspec
  x_misspec = as.matrix(data.frame(X_sqr = x_obs[,(ncol(x_obs) - 1)]^2,
                                   X_log = log(x_obs[,ncol(x_obs)] - (min(x_obs[,ncol(x_obs)]) - 0.1)))) %>% as.data.frame
  
  # no PS misspecification
  if(misspec_PS == 0){
    x_PS = x_obs
    gamma_ah_adj = gamma_ah; gamma_pro_adj = gamma_pro; gamma_ns_adj = gamma_ns; betas_GPI_adj = betas_GPI 
  }
  
  # misspec2: replace 2 X's with x^2 and ~log(X), to PS model and possibly to outcome model (if misspec_outcome != 0) 
  if(misspec_PS == 2){
    # PS true model covariates
    x_PS = as.matrix( data.frame( x_obs[,-c((ncol(x_obs) - 1) ,ncol(x_obs))], x_misspec ) )
    
    gamma_ah_adj = c(gamma_ah[-c((ncol(x_obs) - 1) ,ncol(x_obs))], funcform_factor_sqr*gamma_ah[2], funcform_factor_log*gamma_ah[2])
    gamma_pro_adj = c(gamma_pro[-c((ncol(x_obs) - 1) ,ncol(x_obs))], funcform_factor_sqr*gamma_pro[2], funcform_factor_log*gamma_pro[2])
    gamma_ns_adj = rep(0, length(gamma_ah_adj)) 
    betas_GPI_adj = betas_GPI
    colnames(betas_GPI_adj) = rep("", ncol(betas_GPI_adj))
  }
  
  # Y true model covariates:
  # if misspec_outcome == 0 (default), Y on original (obs) X 
  # if misspec_outcome == 2, Y on the transformation of X, as in the misspec in the PS model
  if(misspec_outcome == 2){ #TODO change it so Y misspec is not necessarily the same as PS misspec 
    x = as.matrix( data.frame( x_obs[,-c((ncol(x_obs) - 1) ,ncol(x_obs))], x_misspec ) )
    } 
  
  if(two_log_models==TRUE){ # two logistic models for s(0) and S(1) given S(0)=1
    #1) log reg of S(0)
    prob_S0 = exp(x_PS%*%gamma_ah_adj) / ( 1 + exp(x_PS%*%gamma_ah_adj) )
    S0_vec = rbinom( length(prob_S0), 1, prob_S0 ) # ah - 1, pro and ns - 0
    
    #2.a) if S(0)==1, assign to as with p = 1/1+xi and to har, with p = xi/x+xi (log reg with a constant only)
    g_vec_num = S0_vec
    g_vec_num[g_vec_num == 0] = -1 # pro and ns
    g_vec_num[g_vec_num == 1] = rbinom( length(g_vec_num[g_vec_num == 1]), 1, (1 / (1+xi)) ) # as - 1, har - 0
    #2.b) if S(0)==0, log reg for S(1)
    prob_S1 = exp(x_PS%*%gamma_pro_adj) / ( 1 + exp(x_PS%*%gamma_pro_adj) ) # prob_S1=1 given S(0)=0
    # +2 for converting ns to 2 and pro to 3
    g_vec_num[g_vec_num == -1] = rbinom(length(prob_S1[g_vec_num == -1]), 1, prob_S1[g_vec_num == -1]) + 2 # ns - 2 and pro - 3 
    
    g_vec = mapvalues(g_vec_num, from = c(0:3), to = c("har", "as", "ns", "pro"))
    prob = data.frame(prob_har = prob_S0*(xi/(1+xi)), prob_as = prob_S0*(1/(1+xi)), 
                      prob_ns = (1-prob_S0)*(1-prob_S1), prob_pro = (1-prob_S0)*prob_S1)
    }else{ # single multinomial model for G
    # vector of probabilities
    vProb = cbind(exp(x_PS%*%gamma_ah_adj), exp(x_PS%*%gamma_ns_adj), exp(x_PS%*%gamma_pro_adj)) 
    prob = vProb / apply(vProb, 1, sum) 
    probs_mean = apply(vProb, 2, mean) / sum(apply(vProb, 2, mean))
    # multinomial draws
    mChoices = t(apply(prob, 1, rmultinom, n = 1, size = 1))
    # 1-ah, 2-pro, 3-ns
    g_vec_num = apply(mChoices, 1, function(z) which(z==1))
    # within ah, randomize to as or har, according to xi
    g_vec_num[g_vec_num==1] = rbinom( length(g_vec_num[g_vec_num==1]), 1, (1 / (1+xi)) )
    # 0-har, 1-ah, 2-pro, 3-ns
    g_vec = mapvalues(g_vec_num, from = c(0:3), to = c("har", "as", "ns", "pro"))
    prob = data.frame(prob_har = prob[,1]*(xi/(1+xi)), prob_as = prob[,1]*(1/(1+xi)), prob_ns = prob[,2], prob_pro = prob[,3])
  }
  
  # descriptive of the principal scores
  pis = table(g_vec) / param_n
  pis = t(c(pis)); colnames(pis) = paste0("pi_", colnames(pis))
  
  # generate data ####
  # data is going to be used in the EM first, in simulate_data_run_EM_and_match. Thus, data contains the "obs" X (x_obs).
  data = data.frame(prob, x_obs, g = g_vec, g_num = g_vec_num,
                    A = rbinom(param_n, 1, prob_A))
  data$S = ifelse((data$g == "as") | (data$g == "pro" & data$A == 1) | (data$g == "har" & data$A == 0), 1, 0)
  mean_by_g = data.table(data)[, lapply(.SD, mean), by="g"]
  mean_by_g$g = mapvalues(mean_by_g$g, from = c("har", "as", "ns", "pro"), to = c(0:3))
  mean_by_A_g = data.table(data)[, lapply(.SD, mean), by=c("A", "g")] %>% arrange(g,A)
  x_har = filter(mean_by_g, g=="har") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_as = filter(mean_by_g, g=="as") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_pro = filter(mean_by_g, g=="pro") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_ns = filter(mean_by_g, g=="ns") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  if(only_mean_x_bool==TRUE){
    return(list(x_har=x_har, x_as=x_as, x_pro=x_pro, x_ns=x_ns, pis=pis, mean_by_A_g=mean_by_A_g))
  }
  
  # models for Y(1) & y(0) 
  PO_by_treatment_and_stratum = function(n, x_outcome, dim_x, coeffs, sigma_square_param){
    return(rnorm(n, mean = x_outcome %*% matrix(coeffs,nrow = dim_x, ncol = 1), sd = sqrt(sigma_square_param)))
  }
  
  # TODO model with PI with pair dependent errors 
  if(rho_GPI_PO != 0){
    print(paste0("rho_GPI_PO", " is ", rho_GPI_PO))
    # simulate dependency between errors
    cov_GPI_PO = rho_GPI_PO * sqrt(var_GPI[1]) * sqrt(var_GPI[2])
    cov_mat <- cbind(c(var_GPI[1], cov_GPI_PO), c(cov_GPI_PO, var_GPI[2]))
    mu_x_beta_Y1 = x %*% matrix(betas_GPI_adj[1,], ncol = 1)
    mu_x_beta_Y0 = x %*% matrix(betas_GPI_adj[2,], ncol = 1)
    two_PO = lapply(1:param_n, function(l){
      mvrnorm(1, mu = c(mu_x_beta_Y1[l], mu_x_beta_Y0[l]), cov_mat)
    })
    two_PO = data.frame(list.rbind(two_PO))
  }
  # model with PI with pair independent errors
  if(rho_GPI_PO == 0){
    print(paste0("rho_GPI_PO", " is ", 0))
    # wout dependency
    two_PO = lapply(1 : nrow(betas_GPI_adj), function(l){
      PO_by_treatment_and_stratum(n=param_n, x_outcome=x, dim_x=dim_x, coeffs=betas_GPI_adj[l,], sigma_square_param=sigma_square_ding[l])
    })
    two_PO = data.frame(list.cbind(two_PO))
  }
  
  colnames(two_PO) = c("Y1", "Y0")
  mean(two_PO$Y1, na.rm = T); mean(two_PO$Y0, na.rm = T)
  dt = data.frame(data, two_PO)
  # generate Y with SUTVA
  dt$Y = (dt$A * dt$Y1 + (1 - dt$A) * dt$Y0) * dt$S
  dt = data.table(id = c(1:param_n), dt)
  dt$OBS = paste0("O(", dt$A, ",", dt$S, ")")
  
  #OBS table
  obs_table = table(dt$OBS)
  OBS_table = matrix(c(obs_table[1], obs_table[3], obs_table[2], obs_table[4]), nrow = 2, ncol = 2)
  OBS_table = OBS_table/ param_n
  rownames(OBS_table) = c("A=0", "A=1"); colnames(OBS_table) = c("S=0", "S=1")
  
  # estimation of strata proportions (under CPSR) 
  p1 = mean(filter(dt, A==1)$S); p0 = mean(filter(dt, A==0)$S)
  pi_har_est = (xi_est/(1+xi_est))*p0
  pi_as_est = (1/(1 + xi_est))*p0
  pi_ns_est = 1 - p1 - (xi_est/(1+xi_est))*p0
  pi_pro_est = p1 - (1/(1+xi_est))*p0
  pis_est = c(pi_har_est = pi_har_est, pi_as_est = pi_as_est, pi_ns_est = pi_ns_est, pi_pro_est = pi_pro_est)
  
  return(list(dt=dt, x_obs=x_obs, x_PS=x_PS, x_outcome=x, mean_by_g=mean_by_g,
              OBS_table=OBS_table, pis=pis, pis_est=pis_est))
}


simulate_data_run_EM_and_match = function(only_EM_bool=FALSE, return_EM_PS=FALSE, index_set_of_params, gamma_ah, gamma_pro, gamma_ns, xi, xi_est, 
                                          two_log_models=TRUE, two_log_est_EM=FALSE,
                                          misspec_PS, misspec_outcome=0, funcform_factor_sqr=0, funcform_factor_log=0, 
                                          param_n, param_n_sim, iterations, epsilon_EM = 0.001,
                                          caliper, match_on = NULL, mu_x_fixed=FALSE, x_as, only_naive_bool=FALSE){
  
  X_sub_cols = paste0("X", c(1:(dim_x)))
  list_dat_EM <- list_coeff_ah <- list_coeff_pro <- list_beta_S0 <- list()
  # run over param_n_sim different samples, each with param_n observations
  WLS_NOint_mat_reg_estimators <- WLS_YESint_mat_reg_estimators <-
    OLS_NOint_mat_reg_estimators <- OLS_YESint_mat_reg_estimators <-
    mat_param_estimators <- mat_excluded_included_matching <- CI_mat <- mat_std_mean_diff <-  NULL
 list_std_mean_diff <- list_means_by_subset <- list_EM_not_conv <- list_BCclpr <- list_mean_by_g <- list()
  
  
  # run over param_n_sim different samples, each with param_n observations
  # i is the index. we dont consider iteration when EM does not converge
  i = 1; real_iter_ind = 1;  index_EM_not_conv = 0
  while (i <= param_n_sim) {
    #for (i in 1:param_n_sim)
    print(paste0("this is index_set_of_params ", index_set_of_params))
    print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match. ",
                 "index_EM_not_conv: ", index_EM_not_conv, ". real number of iterations: "  , real_iter_ind, "."))
    start_time1 <- Sys.time()
    list_data_for_EM_and_X = simulate_data_function(seed_num=NULL, gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, 
            xi=xi, xi_est=xi_est, two_log_models=two_log_models, param_n=param_n,
            misspec_PS=misspec_PS, misspec_outcome=misspec_outcome, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log)
    
    data_for_EM = list_data_for_EM_and_X$dt
    mean_by_g = list_data_for_EM_and_X$mean_by_g
    x = list_data_for_EM_and_X$x_obs; x_PS = data.frame(list_data_for_EM_and_X$x_PS)
    x_outcome = data.frame(list_data_for_EM_and_X$x_outcome)
    OBS_table = list_data_for_EM_and_X$OBS_table
    pis = list_data_for_EM_and_X$pis; pis_est = list_data_for_EM_and_X$pis_est
    vec_OBS_table = t(c(OBS_table)); colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
    
    # real parameter
    SACE = mean(data_for_EM[g=="as" , Y1]) - mean(data_for_EM[g=="as", Y0])
    SACE_conditional = mean(data_for_EM[A==1 & g=="as" , Y]) - mean(data_for_EM[A==0 & g=="as", Y])
    
    # naive estimators
    # naive
    most_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y]) 
    most_naive_est_se = sqrt(  ( var(data_for_EM[A==1, Y])  / nrow(data_for_EM[A==1, ]) ) + 
                                 ( var(data_for_EM[A==0, Y])  / nrow(data_for_EM[A==0, ]) )  )  
    CI_by_SE_and_Z_val_most_naive = round(most_naive_est + c(-1,1) * 1.96 * most_naive_est_se, 3)
    CI_by_SE_and_Z_val_most_naive = paste(CI_by_SE_and_Z_val_most_naive, sep = ' ', collapse = " , ")
    
    # survivors naive
    sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
    sur_naive_est_se = sqrt(  ( var(data_for_EM[A==1 & S==1, Y])  / nrow(data_for_EM[A==1 & S==1, ]) ) + 
                                ( var(data_for_EM[A==0 & S==1, Y])  / nrow(data_for_EM[A==0 & S==1, ]) )  )
    CI_by_SE_and_Z_val_sur_naive = round(sur_naive_est + c(-1,1) * 1.96 * sur_naive_est_se, 3)
    CI_by_SE_and_Z_val_sur_naive = paste(CI_by_SE_and_Z_val_sur_naive, sep = ' ', collapse = " , ")
    
    CI_naives_before_matching = data.frame(CI_by_SE_and_Z_val_most_naive, CI_by_SE_and_Z_val_sur_naive)
    colnames(CI_naives_before_matching) = c("naive_without_matching", "survivors_naive_without_matching")
    
    if(only_naive_bool==TRUE){
      # TODO all naive estimators together in the current row of mat_param_estimators
      mat_param_estimators = rbind( mat_param_estimators,
            data.frame(SACE, most_naive_est, most_naive_est_se, sur_naive_est, sur_naive_est_se, pis, t(pis_est) ))
      CI_mat = rbind( CI_mat, data.frame(SACE, CI_naives_before_matching) )
      i = i + 1
      next()
    }
    
    # Ding estimator ####
    
    if(two_log_est_EM == FALSE){
      #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
      fit_S0_in_A0 = glm(as.formula(paste0("S ~ ",paste(X_sub_cols[-1], collapse="+"))), data=filter(data_for_EM, A==0), family="binomial")
      beta_S0 = fit_S0_in_A0$coefficients
      #P_S0 = predict(fit_S0_in_A0, newdata=data_for_EM, type = "response") 
      # P_S0[i] # expit(predict(fit_S0_in_A0, newdata=data_for_EM)[i]) # expit(t(beta_S0)%*%X[i, ]) # fit_S0_in_A0$fitted.values
    }else{beta_S0=NULL}
    #S(1)=1: logistic regression S(1)=1 on X, using S|A=1 (1 for pro) with weights 1-P_S0, so ah get less weight
    #fit_pseudo_S1_given_A1 = glm(as.formula(paste0("S ~ ",paste(X_sub_cols[-1], collapse="+"))), data=filter(data_for_EM, A==1),
    #                            weights = 1-P_S0[data_for_EM$A==1], family="quasibinomial")
    #beta_S1 = fit_pseudo_S1_given_A1$coefficients
    
    # EM
    start_timeDing <- Sys.time()
    est_ding_lst = xi_2log_PSPS_M_weighting(Z=data_for_EM$A, D=data_for_EM$S,
                    X=as.matrix(subset(data_for_EM, select = 
                    grep(paste(X_sub_cols[-1], collapse="|"), colnames(data_for_EM)))), Y=data_for_EM$Y, 
                    xi_est=xi_est, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL, # beta.S0=beta_S0 # beta.S0=NULL 
                    iter.max=iterations, error0=epsilon_EM)
    coeff_ah = est_ding_lst$beta.ah ; coeff_pro = est_ding_lst$beta.c
    list_beta_S0[[i]] = beta_S0; list_coeff_ah[[i]] = coeff_ah; list_coeff_pro[[i]] = coeff_pro
    EM_coeffs = rbind(est_ding_lst$beta.ah, est_ding_lst$beta.c)
    end_timeDing <- Sys.time()
    print(paste0("Ding EM lasts ", difftime(end_timeDing, start_timeDing)))
    PS_est = est_ding_lst$ps.score
    data_with_PS = data.table(data_for_EM, PS_est)
    
    # if PS_est contains NAS, it probably implies that the EM process has not converged, so skip this iteration and go to the next
    if( sum(is.na(PS_est)) > 0 ){ # or if(est_ding_lst$iter == iterations)
      index_EM_not_conv = index_EM_not_conv + 1
      list_EM_not_conv$probs[[index_EM_not_conv]] = PS_est
      #coeff_ah = est_ding_lst$beta.ah ; coeff_pro = est_ding_lst$coeff_pro
      list_EM_not_conv$coeffs[[index_EM_not_conv]] = data.frame(rbind(coeff_ah=coeff_ah, coeff_pro=coeff_pro))
      colnames(list_EM_not_conv$coeffs[[index_EM_not_conv]]) = X_sub_cols
      list_EM_not_conv$probs_nas[[index_EM_not_conv]] = c(total_na = sum(is.na(PS_est)), 
                                                          prop_na = sum(is.na(PS_est)) / ( nrow(PS_est) * ncol(PS_est) ) ) %>% round(3)
      real_iter_ind = real_iter_ind + 1
      next()
    }
    
    # calculate O11_prior_ratio, O11_posterior_ratio and W_1_as
    O11_prior_ratio = pis_est["pi_as_est"] / (pis_est["pi_as_est"] + pis_est["pi_pro_est"])
    data_with_PS[, `:=` ( O11_posterior_ratio = EMest_p_as / (EMest_p_as + EMest_p_pro), O11_prior_ratio = O11_prior_ratio )]
    data_with_PS$W_1_as = data_with_PS$O11_posterior_ratio / O11_prior_ratio #data_with_PS$O11_posterior_ratio / O11_prior_ratio
    data_with_PS$O11_posterior_ratio_true = data_with_PS$prob_as / (data_with_PS$prob_as + data_with_PS$prob_pro)
    # DL plain estimator and model assisted estimator
    DING_est = est_ding_lst$AACE
    DING_model_assisted_est_ps = est_ding_lst$AACE.reg
    
    # return only EM coefficients
    if(only_EM_bool){
      i = i + 1
      next()
    }
    
    # EM summary
    if(return_EM_PS == TRUE){
      pis = data.frame(pis)
      O11_prior_ratio_true = pis["pi_as"] / (pis["pi_as"] + pis["pi_pro"])
      PS_true_EM_compr = subset( data_with_PS, select = grep("^id$|^g$|prob.|EM|O11_posterior_ratio|W_1_as",colnames(data_with_PS)) )
      PS_true_EM_compr = rapply(object = PS_true_EM_compr, f = round, classes = "numeric", how = "replace", digits = 3)
      PS_true_EM_compr = data.frame( id = PS_true_EM_compr$id, g = PS_true_EM_compr$g,
        prob_as = PS_true_EM_compr$prob_as, EMest_p_as=PS_true_EM_compr$EMest_p_as, diff = PS_true_EM_compr$prob_as - PS_true_EM_compr$EMest_p_as,
        prob_har = PS_true_EM_compr$prob_har, EMest_p_har=PS_true_EM_compr$EMest_p_har, 
        prob_ns = PS_true_EM_compr$prob_ns, EMest_p_ns=PS_true_EM_compr$EMest_p_ns,
        prob_pro = PS_true_EM_compr$prob_pro, EMest_p_pro=PS_true_EM_compr$EMest_p_pro,
        O11_posterior_ratio_true = PS_true_EM_compr$O11_posterior_ratio_true, O11_posterior_ratio = PS_true_EM_compr$O11_posterior_ratio,
        W_1_as_true = PS_true_EM_compr$O11_posterior_ratio_true / O11_prior_ratio, W_1_as = PS_true_EM_compr$W_1_as)
      return(list(data_with_PS=data_with_PS, PS_true_EM_compr=PS_true_EM_compr, true_x_PS=x_PS,
                  pis=pis, pis_est=pis_est, EM_coeffs=EM_coeffs, 
                  O11_prior_ratio_true=O11_prior_ratio_true, O11_prior_ratio=O11_prior_ratio, OBS_table=OBS_table, 
                  beta_S0=beta_S0, error=est_ding_lst$error, mean_by_g=mean_by_g,
                  SACE=SACE, DL=DING_est, DL_MA=DING_model_assisted_est_ps, 
                  list_EM_not_conv=list_EM_not_conv, real_iter_ind=real_iter_ind))
    }
    
    # run for all options (3 options - full dataset, wout A=0,S=0, only S=1)
    data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 
    lst_matching_estimators = list()
    replace_vec = c(FALSE, TRUE)
    for(j in c(1:length(replace_vec))){
      lst_matching_estimators[[j]] =
        lapply(1:length(data_list), function(l){
          my_matching_func_multiple(match_on = match_on, X_sub_cols, data_list[[l]],
                                    weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2,
                                    min_PS = min_PS, min_diff_PS = min_diff_PS,
                                    caliper = caliper, OBS_table, change_id = TRUE, mu_x_fixed=mu_x_fixed, x_as=x_as, pass_tables_matched_units=FALSE)
        })
    }
    

    # matching_estimators 
    matching_estimators = lapply(1:length(replace_vec), function(j){
      data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators[[j]], head, 8))))) 
    })
    
    # MULTIPLE MATCHING WITH BC
    matching_estimators = list.cbind(matching_estimators)
    colnames(matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(matching_estimators)/2)),
                                           "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), each=2, times=3*2*2),
                                           rep(c("_PS", "_maha", "","_HL","_BC","_BCclpr", "_BC_inter","_BCclpr_inter"), each=6, times=2),
                                           rep(c("_est","_SE"), times=length(matching_estimators)/2)) 
    
    CI_matching_estimators = lapply(1:length(replace_vec), function(j){
      as.vector(t(t(unlist(list.rbind(lapply(lst_matching_estimators[[j]],
                                             "[[", "CI_crude_HL_BC"))))))
    })
    CI_matching_estimators = data.frame((unlist(CI_matching_estimators))) %>% t
    rownames(CI_matching_estimators) = ""
    
    # MULTIPLE MATCHING WITH BC
    # with and wout BC inter
    colnames(CI_matching_estimators) = paste0(paste0("rep", 
                                                     rep(substr(replace_vec,1,1), each=length(CI_matching_estimators)/2)),
                                              "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), times=3*2),
                                              rep(c("_PS", "_maha","","_HL","_BC","_BCclpr"
                                                    , "_BC_inter","_BCclpr_inter" # remove if we dont have BC with interactions
                                              ), each=3, times=2))
    
    
    # rep_bool_false_true = 2 for WLS and 1 for OLS
    arrange_lin_reg_estimators = function(rep_bool_false_true, estimator_str,
                                          CI_or_TABLE_EST_SE="estimator_and_se", name=""){
      LS_lin_reg_estimators = lapply(rep_bool_false_true:rep_bool_false_true, function(j){
        lapply(lst_matching_estimators[[j]], "[[",
               estimator_str)})
      LS_lin_reg_estimators = lapply(LS_lin_reg_estimators[[1]], "[[", CI_or_TABLE_EST_SE)
      if(CI_or_TABLE_EST_SE=="estimator_and_se"){
        colnas = paste0(rep(colnames(LS_lin_reg_estimators[[1]]), times=3), "_",
                        rep(c("all", "wout_O_0_0", "S1"), each=length(LS_lin_reg_estimators[[1]])))
        LS_lin_reg_estimators = list.cbind(LS_lin_reg_estimators)
        colnames(LS_lin_reg_estimators) = colnas
        
      }else{
        colnas = name
        LS_lin_reg_estimators = list.cbind(LS_lin_reg_estimators)
        colnames(LS_lin_reg_estimators) =  paste0(colnas, "_", c("all", "wout_O_0_0", "S1"))
      }
      
      return(LS_lin_reg_estimators)
    }
    
    # WLS ONLY FOR replace_vec == TRUE
    # WLS wout interactions
    WLS_NOint_matching_reg_estimators = arrange_lin_reg_estimators(2, "WLS_NOinteractions_reg_adj_estimators_and_se")
    WLS_NOint_matching_reg_estimators_CI = arrange_lin_reg_estimators(2, "WLS_NOinteractions_reg_adj_estimators_and_se", "CI_LS", name="WLS_NOint")
    # WLS with interactions
    WLS_YESint_matching_reg_estimators = arrange_lin_reg_estimators(2, "WLS_YESinteractions_reg_adj_estimators_and_se")
    WLS_YESint_matching_reg_estimators_CI = arrange_lin_reg_estimators(2, "WLS_YESinteractions_reg_adj_estimators_and_se", "CI_LS", name="WLS_YESint")
    # OLS ONLY FOR replace_vec == FALSE
    # OLS wout interactions
    OLS_NOint_matching_reg_estimators = arrange_lin_reg_estimators(1, "OLS_NOinteractions_reg_adj_estimators_and_se")
    OLS_NOint_matching_reg_estimators_CI = arrange_lin_reg_estimators(1, "OLS_NOinteractions_reg_adj_estimators_and_se", "CI_LS", name="OLS_NOint")
    # OLS with interactions
    OLS_YESint_matching_reg_estimators = arrange_lin_reg_estimators(1, "OLS_YESinteractions_reg_adj_estimators_and_se")
    OLS_YESint_matching_reg_estimators_CI = arrange_lin_reg_estimators(1, "OLS_YESinteractions_reg_adj_estimators_and_se", "CI_LS", name="OLS_YESint")
    
    # TODO Coverage CI
    
    
    # means_by_subset
    means_by_subset_lst = lapply(1:length(replace_vec), function(j){
      lapply(lst_matching_estimators[[j]], "[[", "means_by_subset")})
    means_by_subset_mat = list.rbind(unlist(means_by_subset_lst, recursive = FALSE))
    rownames(means_by_subset_mat) = paste0(rep(c("reF_", "reT_"), each=18), 
                                           rep(c("all", "wout_0_0", "S1"), each=6), "_", rownames(means_by_subset_mat))
    
    
    # check ties in BC caliper
    BCclpr_untrt_surv_matched_untrt_matched_trt = lapply(1:length(replace_vec), function(j){
      list.rbind(lapply(lst_matching_estimators[[j]],
                        "[[", "BCclpr_untrt_surv_matched_untrt_matched_trt"))
    })
    
    # SUMMARIES
    # put all results together in the current row of mat_param_estimators
    mat_param_estimators = rbind( mat_param_estimators,
                                  data.frame(SACE, SACE_conditional,
                                             DING_est, DING_model_assisted_est_ps
                                             ,matching_estimators
                                             ,most_naive_est, most_naive_est_se, sur_naive_est, sur_naive_est_se
                                             ,pis, t(pis_est), vec_OBS_table
                                  ))
    
    # regression estimators
    WLS_NOint_mat_reg_estimators = rbind(WLS_NOint_mat_reg_estimators,
                                         data.frame(SACE, SACE_conditional, WLS_NOint_matching_reg_estimators) )
    
    WLS_YESint_mat_reg_estimators = rbind(WLS_YESint_mat_reg_estimators,
                                          data.frame(SACE, SACE_conditional, WLS_YESint_matching_reg_estimators) )
    
    OLS_NOint_mat_reg_estimators = rbind(OLS_NOint_mat_reg_estimators,
                                         data.frame(SACE, SACE_conditional, OLS_NOint_matching_reg_estimators) )
    
    OLS_YESint_mat_reg_estimators = rbind(OLS_YESint_mat_reg_estimators,
                                          data.frame(SACE, SACE_conditional, OLS_YESint_matching_reg_estimators) )
    
    # CI of matching plain estimators and regression matching estimators
    CI_mat = rbind( CI_mat, data.frame(SACE, SACE_conditional
                                       ,CI_naives_before_matching, CI_matching_estimators
                                       ,WLS_NOint_matching_reg_estimators_CI, WLS_YESint_matching_reg_estimators_CI
                                       ,OLS_NOint_matching_reg_estimators_CI, OLS_YESint_matching_reg_estimators_CI) )
    
    
    # list_means_by_subset
    list_means_by_subset[[i]] = means_by_subset_mat
    
    #list_means_by_g
    list_mean_by_g[[i]] = apply(as.matrix(arrange(mean_by_g, g)), 2, as.numeric)
    
    # check ties in BC caliper
    list_BCclpr[[i]] = BCclpr_untrt_surv_matched_untrt_matched_trt[[2]]
    
    # if the em process converges, add 1 to the index
    i = i + 1; real_iter_ind = real_iter_ind + 1
    end_time1 <- Sys.time()
    print(paste0("one iteration lasts ", difftime(end_time1, start_time1)))
    print(difftime(end_time1, start_time1))
  }
  # out of for loop for all the samples (all in all: param_n_sim samples)
  
  # only summary of EM coefficients
  if(only_EM_bool){
    return(list(list_beta_S0=list_beta_S0, list_coeff_ah=list_coeff_ah, list_coeff_pro=list_coeff_pro))
  }
  
  # summary of mat_param_estimators: mean and sd
  param_SACE = mean(mat_param_estimators$SACE)
  MSE_fun <- function (x) mean((x-param_SACE)^2) # param_SACE SACE mean(x)
  mat_param_estimators = rbind(mat_param_estimators,
                               mean = apply(mat_param_estimators, 2, mean), 
                               med = apply(mat_param_estimators, 2, median),
                               SD = apply(mat_param_estimators, 2, sd)
                               , MSE <- list.cbind(lapply(mat_param_estimators, FUN = MSE_fun))
  )
  
  rownames(mat_param_estimators) = 
    c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  
  #TODO if only_naive_bool==TRUE, i.e. we want only naive, stop here
  if(only_naive_bool==TRUE){
    CI_mat$SACE = mean(CI_mat$SACE)
    return(list(mat_param_estimators=mat_param_estimators, CI_mat=CI_mat))
  }
  
  # summary of mat_param_estimators: mean and sd
  list_reg_LS = list(WLS_NOint_mat_reg_estimators=WLS_NOint_mat_reg_estimators, WLS_YESint_mat_reg_estimators=WLS_YESint_mat_reg_estimators,
                     OLS_NOint_mat_reg_estimators=OLS_NOint_mat_reg_estimators, OLS_YESint_mat_reg_estimators=OLS_YESint_mat_reg_estimators)
  for(i in 1:length(list_reg_LS)){
    list_reg_LS[[i]] = rbind(list_reg_LS[[i]],
                             mean = apply(list_reg_LS[[i]], 2, mean), 
                             med = apply(list_reg_LS[[i]], 2, median),
                             SD = apply(list_reg_LS[[i]], 2, sd)
                             , MSE <- list.cbind(lapply(list_reg_LS[[i]], FUN = MSE_fun))
    )
    rownames(list_reg_LS[[i]]) = 
      c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  }
  
  apply(CI_mat[,c(1:2)], 2, as.numeric)
  CI_mat$SACE = mean(CI_mat$SACE); CI_mat$SACE_conditional  = mean(CI_mat$SACE_conditional )
  
  # calculating EM estimators and arrange in a data frame
  # TODO put this thing in a function
  coeffs_ah = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_ah, "[[", l)
  })
  coeffs_pro = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_pro, "[[", l)
  })
  coeffs = data.table( rbind( list.rbind(coeffs_ah), list.rbind(coeffs_pro) ) )
  #TODO genefilter package: coeffs[, `:=` (SD = rowSds(as.matrix(coeffs)), mean = rowMeans(coeffs))]
  coeffs = data.frame(coeffs, SD = apply(coeffs, 1, sd), mean = apply(coeffs, 1, mean))
  coeffs$parameter = c(gamma_ah, gamma_pro)
  coeffs$diff = coeffs$mean - coeffs$parameter 
  coeffs$perc = coeffs$diff / abs(coeffs$parameter)
  coeffs_df = data.frame(coeffs)
  rownames(coeffs_df) = colnames(mat_gamma)
  #rownames(coeffs_df) = c("coeff_ah_0", "coeff_ah_1","coeff_pro_0", "coeff_pro_1")
  
  # list_means_by_subset
  mean_list_means_by_subset = calculate_mean_repeated_as_and_pro(list_means_by_subset, FALSE)
  
  # list_means_by_subset
  #mean_list_by_g = calculate_mean_repeated_as_and_pro(list_mean_by_g, FALSE)
  mean_list_by_g = apply(simplify2array(list_mean_by_g), 2, rowMeans, na.rm = TRUE)
  #mean_list_by_g[,1] = mapvalues(as.numeric(mean_list_by_g[,1]), from = c(0:3), to = c("har", "as", "ns", "pro"))
  
  return(list(mat_param_estimators = mat_param_estimators, 
              WLS_NOint_mat_reg_estimators = list_reg_LS$WLS_NOint_mat_reg_estimators,
              WLS_YESint_mat_reg_estimators = list_reg_LS$WLS_YESint_mat_reg_estimators,
              OLS_NOint_mat_reg_estimators = list_reg_LS$OLS_NOint_mat_reg_estimators,
              OLS_YESint_mat_reg_estimators = list_reg_LS$OLS_YESint_mat_reg_estimators,
              CI_mat = CI_mat,
              coeffs_df = coeffs_df, 
              mean_list_means_by_subset = mean_list_means_by_subset,
              mean_list_by_g=mean_list_by_g,
              list_EM_not_conv = list_EM_not_conv,
              list_BCclpr = list_BCclpr
  ))
}


calculate_mean_repeated_as_and_pro = function(list_of_lists, mean_repeated_as_and_pro_boll=TRUE){
  temp = lapply(1:length(list_of_lists), function(i){
    list_of_lists[[i]] = as.matrix(list_of_lists[[i]])
    apply(list_of_lists[[i]], 2, as.numeric)
  })
  mean_repeated_as_and_pro = data.frame(apply(simplify2array(temp), 1:2, mean))
  rownames(mean_repeated_as_and_pro) = rownames(list_of_lists[[1]])
  if(mean_repeated_as_and_pro_boll == TRUE){
    as_rep = subset(mean_repeated_as_and_pro,
                    select = grep("as_", colnames(mean_repeated_as_and_pro)))
    pro_rep = subset(mean_repeated_as_and_pro,
                     select = grep("pro_", colnames(mean_repeated_as_and_pro)))
    #mean_as = as.matrix(as_rep) %*% c(1:ncol(as_rep)) / ifelse(as_rep>0, 1, 0) %*% c(1:ncol(as_rep))
    mean_as = as.matrix(as_rep) %*% c(1:ncol(as_rep)) / apply(as_rep, 1, sum)
    mean_pro = as.matrix(pro_rep) %*% c(1:ncol(pro_rep)) / apply(pro_rep, 1, sum)
    mean_repeated_as_and_pro = data.frame(mean_repeated_as_and_pro, mean_as = mean_as, mean_pro = mean_pro)
  }
  return(mean_repeated_as_and_pro)
}
  

