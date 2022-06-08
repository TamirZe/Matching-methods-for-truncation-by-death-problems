simulate_data_run_EM_and_match = function(only_EM_bool=FALSE, return_EM_PS=FALSE, index_set_of_params, gamma_ah, gamma_pro, gamma_ns, xi, xi_est, 
                                          two_log_models=TRUE, two_log_est_EM=FALSE,
                                          misspec_PS, misspec_outcome=0, transform_x=0, 
                                          funcform_factor_sqr, funcform_factor_log, 
                                          param_n, param_n_sim, iterations, epsilon_EM = 0.001,
                                          caliper, match_on = NULL, mu_x_fixed=FALSE, x_as, only_naive_bool=FALSE){
  
  X_sub_cols = paste0("X", c(1:(dim_x)))
  list_coeff_ah <- list_coeff_pro <- list_beta_S0 <- list_EM_not_conv <- 
    list_std_mean_diff <- list_means_by_subset <- list_BCclpr <- list_mean_by_g <- list()
  WLS_NOint_mat_reg_estimators <- WLS_YESint_mat_reg_estimators <-
    OLS_NOint_mat_reg_estimators <- OLS_YESint_mat_reg_estimators <-
    mat_param_estimators <- mat_excluded_included_matching <- CI_mat <- mat_std_mean_diff <-  NULL
  
  # run over param_n_sim different samples, each with param_n observations
  # i is the index. we dont consider iteration when EM does not converge
  i = 1; real_iter_ind = 1;  index_EM_not_conv = 0
  while (i <= param_n_sim) {
    #for (i in 1:param_n_sim)
    print(paste0("this is index_set_of_params ", index_set_of_params))
    print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match. ",
                 "index_EM_not_conv: ", index_EM_not_conv, ". real number of iterations: "  , real_iter_ind, "."))
    start_time1 <- Sys.time()
    # simulate data
    list_data_for_EM_and_X = simulate_data_function(seed_num=NULL, gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, 
                                                    xi=xi, xi_est=xi_est, two_log_models=two_log_models, param_n=param_n,
                                                    misspec_PS=misspec_PS, misspec_outcome=misspec_outcome, transform_x=transform_x,
                                                    funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log)
    
    data_for_EM = list_data_for_EM_and_X$dt
    mean_by_g = list_data_for_EM_and_X$mean_by_g
    x = list_data_for_EM_and_X$x_obs; x_PS = data.frame(list_data_for_EM_and_X$x_PS)
    x_outcome = data.frame(list_data_for_EM_and_X$x_outcome)
    OBS_table = list_data_for_EM_and_X$OBS_table
    pis = list_data_for_EM_and_X$pis; pis_est = list_data_for_EM_and_X$pis_est
    vec_OBS_table = t(c(OBS_table)); colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
    # "real" SACE parameter
    SACE = list_data_for_EM_and_X$true_SACE
    # naive estimators
    naive_sace_estimation = naive_sace_estimation_func(data_for_EM)
    # naive (composite)
    most_naive_est = naive_sace_estimation$most_naive_est
    most_naive_est_se = naive_sace_estimation$most_naive_est_se
    # survivors naive
    sur_naive_est = naive_sace_estimation$sur_naive_est
    sur_naive_est_se = naive_sace_estimation$sur_naive_est_se
    # CI of naive (composite) and survivors naive (S1)
    CI_naives_before_matching = naive_sace_estimation$CI_naives_before_matching
    
    if(only_naive_bool==TRUE){
      # all naive estimators together in the current row of mat_param_estimators
      mat_param_estimators = rbind( mat_param_estimators,
                                    data.frame(SACE, most_naive_est, most_naive_est_se, sur_naive_est, sur_naive_est_se,
                                               pis, t(pis_est) ))
      CI_mat = rbind( CI_mat, data.frame(SACE, CI_naives_before_matching) )
      i = i + 1
      next()
    }
    
    # EM and PS estiamtion
    # If beta_S0=NULL, employ two logistic regressions during the EM
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
    list_beta_S0[[i]] = beta_S0; list_coeff_ah[[i]] = coeff_ah; list_coeff_pro[[i]] = coeff_pro
    
    PS_est = est_ding_lst$ps.score
    data_with_PS = data.table(data_for_EM, PS_est)
    
    # Ding estimator DL plain estimator and model assisted estimator
    DING_est = est_ding_lst$AACE
    DING_model_assisted_est_ps = est_ding_lst$AACE.reg
    
    # if PS_est contains NAS, it probably implies that the EM process has not converged, so skip this iteration and go to the next
    if( sum(is.na(PS_est)) > 0 | (est_ding_lst$iter == iterations+1 & est_ding_lst$error>=epsilon_EM)){ # or if(est_ding_lst$iter == iterations+1 & est_ding_lst$error>=epsilon_EM)
      print("EM has not convereged")
      index_EM_not_conv = index_EM_not_conv + 1
      list_EM_not_conv$probs[[index_EM_not_conv]] = PS_est
      list_EM_not_conv$coeffs[[index_EM_not_conv]] = data.frame(rbind(coeff_ah=coeff_ah, coeff_pro=coeff_pro))
      colnames(list_EM_not_conv$coeffs[[index_EM_not_conv]]) = X_sub_cols
      list_EM_not_conv$probs_nas[[index_EM_not_conv]] = c(total_na = sum(is.na(PS_est)), 
                                                          prop_na = sum(is.na(PS_est)) / ( nrow(PS_est) * ncol(PS_est) ) ) %>% round(3)
      real_iter_ind = real_iter_ind + 1
      next()
    }
    
    # calculate weights: O11_prior_ratio, O11_posterior_ratio W_1_as, and W_1_as_true
    weights_lst = add_weights_func(data_with_PS=data_with_PS, pis=pis, pis_est=pis_est)
    O11_prior_ratio = weights_lst$O11_prior_ratio
    O11_prior_ratio_true = weights_lst$O11_prior_ratio_true
    data_with_PS = weights_lst$data_with_PS 
      
    # return only EM coefficients
    if(only_EM_bool){
      i = i + 1
      next()
    }
    
    # EM summary
    if(return_EM_PS == TRUE){
      #TODO function?
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
                  em_NOT_conv = est_ding_lst$iter == iterations+1 & est_ding_lst$error>=epsilon_EM, real_iter_ind=real_iter_ind))
    }
    
    # matching
    # perform matching on all the 3 dataset (full, wout A=0,S=0, and only S=1)
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
    
    
    # WLS ONLY FOR replace_vec == TRUE # rep_bool_false_true = 2 for WLS and 1 for OLS
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
    
    #TODO Coverage CI
    
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
    
    # add the current iteration (i.e. param_n_sim) to the big lists/matrices that contain info re all iterations
    mat_param_estimators = rbind( mat_param_estimators,
                                  data.frame(SACE
                                             ,DING_est, DING_model_assisted_est_ps, matching_estimators
                                             ,most_naive_est, most_naive_est_se, sur_naive_est, sur_naive_est_se
                                             ,pis, t(pis_est), vec_OBS_table
                                  ))
    
    # regression estimators
    WLS_NOint_mat_reg_estimators = rbind(WLS_NOint_mat_reg_estimators, data.frame(SACE, WLS_NOint_matching_reg_estimators))
    WLS_YESint_mat_reg_estimators = rbind(WLS_YESint_mat_reg_estimators, data.frame(SACE, WLS_YESint_matching_reg_estimators))
    OLS_NOint_mat_reg_estimators = rbind(OLS_NOint_mat_reg_estimators, data.frame(SACE, OLS_NOint_matching_reg_estimators))
    OLS_YESint_mat_reg_estimators = rbind(OLS_YESint_mat_reg_estimators, data.frame(SACE, OLS_YESint_matching_reg_estimators))
    
    # CI of matching crude, BC and regression estimators after matching
    CI_mat = rbind( CI_mat, data.frame(SACE
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
  } # out of for loop for all the samples (all in all: param_n_sim samples)
  
  # only summary of EM coefficients
  if(only_EM_bool){
    return(list(list_beta_S0=list_beta_S0, list_coeff_ah=list_coeff_ah, list_coeff_pro=list_coeff_pro))
  }
  
  # summary of mat_param_estimators: mean, med, empirical sd and MSE
  param_SACE = mean(mat_param_estimators$SACE)
  MSE_fun <- function (x) mean((x-param_SACE)^2) # param_SACE SACE mean(x)
  mat_param_estimators = rbind(mat_param_estimators,
         mean = apply(mat_param_estimators, 2, mean), med = apply(mat_param_estimators, 2, median),
         SD = apply(mat_param_estimators, 2, sd), MSE <- list.cbind(lapply(mat_param_estimators, FUN = MSE_fun))
  )
  rownames(mat_param_estimators) = c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  
  # summary of regression estimators: mean, med, empirical sd and MSE
  list_reg_LS = list(WLS_NOint_mat_reg_estimators=WLS_NOint_mat_reg_estimators, WLS_YESint_mat_reg_estimators=WLS_YESint_mat_reg_estimators,
                     OLS_NOint_mat_reg_estimators=OLS_NOint_mat_reg_estimators, OLS_YESint_mat_reg_estimators=OLS_YESint_mat_reg_estimators)
  for(i in 1:length(list_reg_LS)){
    list_reg_LS[[i]] = rbind(list_reg_LS[[i]],
                             mean = apply(list_reg_LS[[i]], 2, mean), med = apply(list_reg_LS[[i]], 2, median),
                             SD = apply(list_reg_LS[[i]], 2, sd), MSE <- list.cbind(lapply(list_reg_LS[[i]], FUN = MSE_fun))
    )
    rownames(list_reg_LS[[i]]) = c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  }
  
  CI_mat$SACE = mean(CI_mat$SACE)
  
  #TODO if only_naive_bool==TRUE, i.e. we want only naive, stop here
  if(only_naive_bool==TRUE){
    return(list(mat_param_estimators=mat_param_estimators, CI_mat=CI_mat))
  }
  
  # calculating EM coefficient estimators and arrange in a data frame
  EM_coeffs_df = arrange_EM_coeffs(list_coeff_ah=list_coeff_ah, list_coeff_pro=list_coeff_pro, dim_x=dim_x)
  
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
              coeffs_df = EM_coeffs_df, 
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

# calculating EM coefficient estimators and arrange in a data frame
arrange_EM_coeffs = function(list_coeff_ah, list_coeff_pro, dim_x){
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
  EM_coeffs_df = data.frame(coeffs)
  rownames(EM_coeffs_df) = colnames(mat_gamma)
  #rownames(EM_coeffs_df) = c("coeff_ah_0", "coeff_ah_1","coeff_pro_0", "coeff_pro_1")
  return(EM_coeffs_df)
}

