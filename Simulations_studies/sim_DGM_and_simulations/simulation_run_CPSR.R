gamma_pro = rep(0, dim_x)
gamma_as = as.numeric(mat_gamma[1, c(1:dim_x)])
gamma_ns =  as.numeric(mat_gamma[1, (dim_x+1): (2*dim_x)])

# misspec_PS: 0 <- NO mis, 2: add transformations to PS model, and remain original X's in ourcome model.
simulate_data_function = function(seed_num=NULL, gamma_as, gamma_ns, gamma_pro, param_n,
                                  misspec_PS, funcform_mis_out=FALSE,
                                  funcform_factor_sqr=0, funcform_factor_log=0, only_mean_x_bool=FALSE){
  if(is.null(seed_num)!=TRUE){set.seed(seed_num)}
  
  # draw covariate matrix
  x_obs <- matrix( c( rep(1,param_n), 
                      mvrnorm(param_n, mu=mean_x, Sigma = diag(var_x, cont_x))), nrow = param_n )
  # add categorial variable, if needed (needed if categ_x > 0)
  if(categ_x > 0){
    x_categorial = data.table(list.cbind(lapply(vec_p_categ, function(x) rbinom(n=param_n,prob=x,size=1))))
    x_obs = as.matrix(cbind(x_obs, x_categorial))
  }
  colnames(x_obs) = paste0("X", c(1:dim_x))
  
  # misspecation
  if(misspec_PS == 0){
    x = x_obs; x_PS = x_obs
    gamma_as_adj = gamma_as; gamma_ns_adj = gamma_ns; gamma_pro_adj = gamma_pro; betas_GPI_adj = betas_GPI 
  }
  
  # misspec2: replace 2 X's with x^2 and ~log(X), to PS model and possibly to outcome model (if funcform_mis_out == TRUE) 
  if(misspec_PS == 2){
    x_misspec = as.matrix(data.frame(x_obs[,(ncol(x_obs) - 1)]^2,
                                     log(x_obs[,ncol(x_obs)] - (min(x_obs[,ncol(x_obs)]) - 0.1)))) %>% as.data.frame
    colnames(x_misspec) = c("X_sqr", "X_log")
    # PS true model covariates
    x_PS = as.matrix( data.frame( x_obs[,-c((ncol(x_obs) - 1) ,ncol(x_obs))], x_misspec ) )
    # Y true model covariates:
    # if funcform_mis_out == FALSE (default), Y on original (obs) X 
    # if funcform_mis_out == TRUE, Y on the transformation of X 
    if(funcform_mis_out == FALSE){
      x = x_obs
    }else{x = x_PS}
    gamma_as_adj = c(gamma_as[-c((ncol(x_obs) - 1) ,ncol(x_obs))], funcform_factor_sqr*gamma_as[2], funcform_factor_log*gamma_as[2])
    gamma_ns_adj = c(gamma_ns[-c((ncol(x_obs) - 1) ,ncol(x_obs))], funcform_factor_sqr*gamma_ns[2], funcform_factor_log*gamma_ns[2])
    gamma_pro_adj = rep(0, length(gamma_as_adj))
    betas_GPI_adj = betas_GPI
    colnames(betas_GPI_adj) = rep("", ncol(betas_GPI_adj))
  }
  
  # vector of probabilities
  vProb = cbind(exp(x_PS%*%gamma_as_adj), exp(x_PS%*%gamma_pro_adj), exp(x_PS%*%gamma_ns_adj)) 
  prob = vProb / apply(vProb, 1, sum) 
  probs_mean = apply(vProb, 2, mean) / sum(apply(vProb, 2, mean))
  # multinomial draws
  mChoices = t(apply(prob, 1, rmultinom, n = 1, size = 1))
  g_vec_num = apply(mChoices, 1, function(z) which(z==1))
  g_vec = ifelse(g_vec_num == 1, "as", ifelse(g_vec_num == 2, "pro", 
                                              ifelse(g_vec_num == 3, "ns", "har")))
  # descriptive of the principal scores
  pi = table(g_vec) / param_n
  pi = t(c(pi)); colnames(pi) = paste0("pi_", colnames(pi))
  
  # generate data ####
  # data is going to be used in the EM first, in simulate_data_run_EM_and_match. Thus, data contains the "obs" X.
  data = data.frame(prob = prob, x_obs, g = g_vec, g_num = g_vec_num,
                    A = rbinom(param_n, 1, prob_A))
  data$S = ifelse(data$g == "as", 1, ifelse( data$g == "pro" & data$A == 1, 1,
                                             ifelse( data$g == "har" & data$A == 0, 1, 0 ) ))
  mean_by_g = data.table(data)[, lapply(.SD, mean), by="g"]
  mean_by_A_g = data.table(data)[, lapply(.SD, mean), by=c("A", "g")] %>% arrange(g,A)
  x_as = filter(mean_by_g, g=="as") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_pro = filter(mean_by_g, g=="pro") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  x_ns = filter(mean_by_g, g=="ns") %>% subset(select = grep("X|^A$", colnames(mean_by_g))) %>% as.matrix
  if(only_mean_x_bool==TRUE){
    return(list(x_as=x_as, x_pro=x_pro, x_ns=x_ns, pi=pi, mean_by_A_g=mean_by_A_g))
  }
  
  # models for Y(1) & y(0) ####
  PO_by_treatment_and_stratum = function(coeffs, sigma_square_param){
    return(rnorm(n, mean = x %*% matrix(coeffs,nrow = dim_x, ncol = 1), sd = sqrt(sigma_square_param)))
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
      PO_by_treatment_and_stratum(betas_GPI_adj[l,], sigma_square_ding[l])
    })
    two_PO = data.frame(list.cbind(two_PO))
  }
  
  colnames(two_PO) = c("Y1", "Y0")
  mean(two_PO$Y1, na.rm = T); mean(two_PO$Y0, na.rm = T)
  dt = data.frame(data, two_PO)
  
  # generate Y with SUTVA
  dt$Y = (dt$A * dt$Y1 + (1 - dt$A) * dt$Y0) * dt$S
  dt = data.frame(id = c(1:param_n), dt)
  dt = data.table(dt)
  # senity check
  as_1 = filter(dt, g=="as", A==1); mean(as_1$Y)
  pro_1 = filter(dt, g=="pro", A==1); mean(pro_1$Y)
  dt$OBS = paste0("O(", dt$A, ",", dt$S, ")")
  #OBS table
  OBS_values = data.frame(unique(cbind(dt$A, dt$S))); colnames(OBS_values) = c("S", "A")
  obs_table = table(dt$OBS)
  OBS_table = matrix(c(obs_table[1], obs_table[3], obs_table[2], obs_table[4]), nrow = 2, ncol = 2)
  OBS_table = OBS_table/ param_n
  rownames(OBS_table) = c("A=0", "A=1"); colnames(OBS_table) = c("S=0", "S=1")
  
  # estimation of strata proportions (under mono) 
  pi_as_est = mean(filter(dt, A==0)$S)
  pi_ns_est = 1 - mean(filter(dt, A==1)$S)
  pi_pro_est = mean(filter(dt, A==1)$S) - mean(filter(dt, A==0)$S)
  pis_est = c(pi_as_est_func = pi_as_est, pi_ns_est_func = pi_ns_est, pi_pro_est_func = pi_pro_est)
  
  return(list(dt=dt, x_obs=x_obs, x_PS=x_PS, x_outcome=x,
              OBS_table=OBS_table, pi=pi, pis_est=pis_est, probs_mean=probs_mean))
}


simulate_data_run_EM_and_match = function(return_EM_PS = FALSE, index_set_of_params, gamma_as, gamma_ns, gamma_pro,
                                          misspec_PS, funcform_mis_out=FALSE, funcform_factor_sqr=0, funcform_factor_log=0, 
                                          match_and_reg_watch_true_X=FALSE, param_n, param_n_sim, iterations, epsilon_EM = 0.001,
                                          caliper, match_on = NULL, mu_x_fixed=FALSE, x_as){
  
  X_sub_cols = paste0("X", c(1:(dim_x)))
  list_dat_EM <- list_coeff_as <- list_coeff_ns <- list()
  # run over param_n_sim different samples, each with param_n observations
  WLS_NOint_mat_reg_estimators <- WLS_YESint_mat_reg_estimators <-
    OLS_NOint_mat_reg_estimators <- OLS_YESint_mat_reg_estimators <-
    mat_param_estimators <- mat_excluded_included_matching <- 
    CI_mat <- mat_diff_distance_aspr_asas <- mat_std_mean_diff <-  NULL
  list_repeated_as_and_pro <- list_diff_distance_aspr_asas <- list_matched_units <- 
    list_std_mean_diff <- list_means_by_subset <- list_EM_not_conv <- list_BCclpr <- list()
  
  
  # run over param_n_sim different samples, each with param_n observations
  # i is the index. we dont consider iteration when EM does not converge
  i = 1; real_iter_ind = 1;  index_EM_not_conv = 0
  while (i <= param_n_sim) {
    #for (i in 1:param_n_sim)
    print(paste0("this is index_set_of_params ", index_set_of_params))
    print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match. ",
                 "index_EM_not_conv: ", index_EM_not_conv, ". real number of iterations: "  , real_iter_ind, "."))
    start_time1 <- Sys.time()
    list_data_for_EM_and_X = simulate_data_function(seed_num=NULL, gamma_as, gamma_ns, gamma_pro, param_n,
                                                    misspec_PS, funcform_mis_out, funcform_factor_sqr, funcform_factor_log)
    data_for_EM = list_data_for_EM_and_X$dt
    x = list_data_for_EM_and_X$x_obs; x_PS = data.frame(list_data_for_EM_and_X$x_PS)
    x_outcome = data.frame(list_data_for_EM_and_X$x_outcome)
    OBS_table = list_data_for_EM_and_X$OBS_table
    pis = list_data_for_EM_and_X$pi; pis_est_from_func = list_data_for_EM_and_X$pis_est
    vec_OBS_table = t(c(OBS_table))
    colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
    # pi estimation
    pi_as_est = mean(filter(data_for_EM, A==0)$S)
    pi_ns_est = 1 - mean(filter(data_for_EM, A==1)$S)
    pi_pro_est = mean(filter(data_for_EM, A==1)$S) - mean(filter(data_for_EM, A==0)$S)
    pis_est = c(pi_as_est = pi_as_est, pi_ns_est = pi_ns_est, pi_pro_est = pi_pro_est)
    
    # real parameter
    SACE = mean(data_for_EM[g_num==1 , Y1]) - mean(data_for_EM[g_num==1, Y0])
    SACE_conditional = mean(data_for_EM[A==1 & g_num==1 , Y]) - mean(data_for_EM[A==0 & g_num==1, Y])
    
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
      # TODO put all naive estimators together in the current row of mat_param_estimators
      mat_param_estimators = rbind( mat_param_estimators,
                                    data.frame(SACE, most_naive_est, most_naive_est_se, sur_naive_est, sur_naive_est_se,
                                               pis, t(pis_est), t(pis_est_from_func) ))
      CI_mat = rbind( CI_mat, data.frame(SACE, CI_naives_before_matching) )
      i = i + 1
      next()
    }
    
    # Ding estimator
    #TODO in ding the pis order id PROB[i,] = c(prob.c, prob.a, prob.n)/sum
    start_timeDing <- Sys.time()
    est_ding_lst = PSPS_M_weighting(Z=data_for_EM$A, D=data_for_EM$S,
                                    X=as.matrix(subset(data_for_EM, 
                                                       select = grep(paste(X_sub_cols[-1], collapse="|"), colnames(data_for_EM)))),  
                                    Y=data_for_EM$Y, trc = TRUE, ep1 = 1, ep0 = 1, beta.a = NULL, beta.n = NULL,
                                    iter.max = iterations , error0 = epsilon_EM) 
    end_timeDing <- Sys.time()
    print(paste0("Ding EM lasts ", difftime(end_timeDing, start_timeDing)))
    # adjust the cols the same order as in myEM: my order is: as, ns, pro. ding order: c(prob.c, prob.a, prob.n)
    PS_est = data.frame(est_ding_lst$PROB[,2], est_ding_lst$PROB[,3], est_ding_lst$PROB[,1])
    colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
    data_with_PS = data.table(data_for_EM, PS_est)
    
    # if PS_est contains NAS, it probably implies that the EM process diverged, skip this iteration and go to the next
    if( sum(is.na(PS_est)) > 0){
      index_EM_not_conv = index_EM_not_conv + 1
      list_EM_not_conv$probs[[index_EM_not_conv]] = PS_est
      coeff_as = est_ding_lst$beta.a ; coeff_ns = est_ding_lst$beta.n
      list_EM_not_conv$coeffs[[index_EM_not_conv]] = data.frame(rbind(coeff_as=coeff_as, coeff_ns=coeff_ns))
      colnames(list_EM_not_conv$coeffs[[index_EM_not_conv]]) = X_sub_cols
      list_EM_not_conv$probs_nas[[index_EM_not_conv]] = c(total_na = sum(is.na(PS_est)), 
                                                          prop_na = sum(is.na(PS_est)) / ( nrow(PS_est) * ncol(PS_est) ) ) %>% round(3)
      real_iter_ind = real_iter_ind + 1
      next()
    }
    
    DING_est = est_ding_lst$AACE
    DING_model_assisted_est_ps = est_ding_lst$AACE.reg
    coeff_as = est_ding_lst$beta.a ; coeff_ns = est_ding_lst$beta.n
    list_coeff_as[[i]] = coeff_as; list_coeff_ns[[i]] = coeff_ns
    
    # EM summary
    if(return_EM_PS == TRUE){
      PS_true_EM_compr = subset( data_with_PS, select = grep("^id$|^g$|prob.|EM",colnames(data_with_PS)) )
      PS_true_EM_compr = rapply(object = PS_true_EM_compr, f = round, classes = "numeric", how = "replace", digits = 3)
      colnames(PS_true_EM_compr) = mgsub(colnames(PS_true_EM_compr), c("prob.1", "prob.2", "prob.3"), c("prob_as", "prob_pro", "prob_ns"))
      PS_true_EM_compr = data.frame(id = PS_true_EM_compr$id, g = PS_true_EM_compr$g,
                                    prob_as = PS_true_EM_compr$prob_as, EMest_p_as=PS_true_EM_compr$EMest_p_as, diff = PS_true_EM_compr$prob_as - PS_true_EM_compr$EMest_p_as,
                                    prob_pro = PS_true_EM_compr$prob_pro, EMest_p_pro=PS_true_EM_compr$EMest_p_pro, prob_ns = PS_true_EM_compr$prob_ns, EMest_p_ns=PS_true_EM_compr$EMest_p_ns)
      return(list(PS_true_EM_compr=PS_true_EM_compr,OBS_table=OBS_table, pis=pis, EM_coeffs=EM_coeffs))
    }
    
    # calculate O11_posterior_ratio and O11_prior_ratio
    # pis order: as, ns, pro
    pis = data.frame(pis) 
    O11_prior_ratio = pis[which(names(pis)=="pi_as")] / (pis[which(names(pis)=="pi_as")] + pis[which(names(pis)=="pi_pro")]) 
    data_with_PS$O11_prior_ratio = O11_prior_ratio
    data_with_PS[, `:=` (O11_posterior_ratio = EMest_p_as / (EMest_p_as + EMest_p_pro),
                         W_1_as = ( EMest_p_as / (EMest_p_as + EMest_p_pro) ) / O11_prior_ratio)]
    
    # run for all options (3 options - full dataset, wout A=0,S=0, only S=1)
    data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 
    lst_matching_estimators_end_excluded_included = list()
    replace_vec = c(FALSE, TRUE)
    for(j in c(1:length(replace_vec))){
      lst_matching_estimators_end_excluded_included[[j]] =
        lapply(1:length(data_list), function(l){
          my_matching_func_multiple(match_on = match_on, X_sub_cols, data_list[[l]],
                                    weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2,
                                    min_PS = min_PS, min_diff_PS = min_diff_PS,
                                    caliper = caliper, OBS_table, change_id = TRUE, mu_x_fixed=mu_x_fixed, x_as=x_as, pass_tables_matched_units=FALSE)
        })
    }
    
    
    # TODO 1.
    matching_estimators = lapply(1:length(replace_vec), function(j){
      data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]], head, 8))))) # 8
    })
    
    # MULTIPLE MATCHING WITH BC
    matching_estimators = list.cbind(matching_estimators)
    colnames(matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(matching_estimators)/2)),
                                           "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), each=2, times=3*2*2),
                                           rep(c("_PS", "_maha", "","_HL","_BC","_BCclpr", "_BC_inter","_BCclpr_inter"), each=6, times=2),
                                           rep(c("_est","_SE"), times=length(matching_estimators)/2)) 
    
    CI_matching_estimators = lapply(1:length(replace_vec), function(j){
      as.vector(t(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
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
        lapply(lst_matching_estimators_end_excluded_included[[j]], "[[",
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
    
    
    # TODO 2.
    # excluded_included_matching from matching function
    # if we have replace FALSE or TRUE, use only the FALSE replace,
    # since when replacing this info is not informative
    excluded_included_matching = 
      as.numeric(sapply(lst_matching_estimators_end_excluded_included[[1]], 
                        "[[", "included_excluded_in_matching"))
    names(excluded_included_matching) = 
      paste0(rep(c("all", "wout_O_0_0", "S1" ), each = 5), "_", 
             rep(colnames(t(sapply(lst_matching_estimators_end_excluded_included[[1]],
                                   "[[", "included_excluded_in_matching")))))
    
    # TODO 3.
    # repeated_as_and_pro
    repeated_as_and_pro_lst = lapply(1:length(replace_vec), function(j){
      lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "repeated_as_and_pro")})
    repeated_as_and_pro_lst = unlist(repeated_as_and_pro_lst, recursive = FALSE)
    
    # TODO 6. tables_matched_units
    matched_units_lst = lapply(1:length(replace_vec), function(j){
      lapply(lst_matching_estimators_end_excluded_included[[j]],  
             "[[", "tables_matched_units")})
    matched_units_lst = unlist(matched_units_lst, recursive = FALSE)
    
    # TODO 4.
    # diff_distance_aspr_asas, 1 number per each part, so we have 6 numbers in total
    diff_distance = lapply(1:length(replace_vec), function(j){
      data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                            "[[", "diff_distance_aspr_asas")))))
    })
    diff_distance = list.cbind(diff_distance)
    colnames(diff_distance) = paste0("d_", paste0( "rep", rep(substr(replace_vec,1,1), each=3), "_",
                                                   paste0(rep(c("MATCH", "MATCH"), each = 3),
                                                          "_", c("all", "wout_O_0_0", "S1"))))
    
    # TODO 5.
    # standardized mean diff as2pro, as2as per each data set, with and wout replacements
    # each element in  the list is a matrix. with all covariates std diff per each part, 6 parts in total
    std_mean_diff_lst = lapply(1:length(replace_vec), function(j){
      lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "std_diff_2_cols")})
    std_mean_diff_lst = unlist(std_mean_diff_lst, recursive = FALSE)
    
    # TODO 7. means_by_subset
    means_by_subset_lst = lapply(1:length(replace_vec), function(j){
      lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "means_by_subset")})
    means_by_subset_mat = list.rbind(unlist(means_by_subset_lst, recursive = FALSE))
    rownames(means_by_subset_mat) = paste0(rep(c("reF_", "reT_"), each=18), 
                                           rep(c("all", "wout_0_0", "S1"), each=6), "_", rownames(means_by_subset_mat))
    
    
    # check ties in BC caliper
    BCclpr_untrt_surv_matched_untrt_matched_trt = lapply(1:length(replace_vec), function(j){
      list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                        "[[", "BCclpr_untrt_surv_matched_untrt_matched_trt"))
    })
    
    
    # TODO 1. put all results together in the current row of mat_param_estimators
    mat_param_estimators = rbind( mat_param_estimators,
                                  data.frame(SACE, SACE_conditional,
                                             DING_est, DING_model_assisted_est_ps
                                             ,matching_estimators
                                             , diff_distance,
                                             most_naive_est, most_naive_est_se, sur_naive_est, sur_naive_est_se,
                                             pis, t(pis_est), t(pis_est_from_func), vec_OBS_table
                                  ))
    
    # TODO regression estimators
    WLS_NOint_mat_reg_estimators = rbind(WLS_NOint_mat_reg_estimators,
                                         data.frame(SACE, SACE_conditional, WLS_NOint_matching_reg_estimators) )
    
    WLS_YESint_mat_reg_estimators = rbind(WLS_YESint_mat_reg_estimators,
                                          data.frame(SACE, SACE_conditional, WLS_YESint_matching_reg_estimators) )
    
    OLS_NOint_mat_reg_estimators = rbind(OLS_NOint_mat_reg_estimators,
                                         data.frame(SACE, SACE_conditional, OLS_NOint_matching_reg_estimators) )
    
    OLS_YESint_mat_reg_estimators = rbind(OLS_YESint_mat_reg_estimators,
                                          data.frame(SACE, SACE_conditional, OLS_YESint_matching_reg_estimators) )
    
    CI_mat = rbind( CI_mat, data.frame(SACE, SACE_conditional
                                       ,CI_naives_before_matching, CI_matching_estimators
                                       ,WLS_NOint_matching_reg_estimators_CI, WLS_YESint_matching_reg_estimators_CI
                                       ,OLS_NOint_matching_reg_estimators_CI, OLS_YESint_matching_reg_estimators_CI) )
    
    # TODO 2. put all results together in the current row of mat_excluded_included_matching
    mat_excluded_included_matching = 
      rbind(mat_excluded_included_matching, excluded_included_matching)
    
    # TODO 3. list of repeated_as_and_pro
    list_repeated_as_and_pro[[i]] = list.rbind(repeated_as_and_pro_lst)
    rownames(list_repeated_as_and_pro[[i]])=
      paste0("rep", rep(substr(replace_vec,1,1), each=3), "_", c("all", "wout_O_0_0", "S1"))
    
    # TODO 6. tables_matched_units
    list_matched_units[[i]] = rbindlist(matched_units_lst
                                        #, fill=TRUE
                                        #,use.names=FALSE
    )
    rownames(list_matched_units[[i]])=
      paste0("rep", rep(substr(replace_vec,1,1), each=3), "_", c("all", "wout_O_0_0", "S1"))
    
    # TODO 4. matrix of diff_abs_distance (number paer each part, so 6 in total)
    # this info is also in mat_param_estimators
    mat_diff_distance_aspr_asas = rbind(mat_diff_distance_aspr_asas, diff_distance)
    
    # TODO 5. list of std_diffs per all X's; size of list as the samples
    list_std_mean_diff[[i]] = std_mean_diff_lst
    
    # TODO 7. list_means_by_subset
    list_means_by_subset[[i]] = means_by_subset_mat
    
    # check ties in BC caliper
    list_BCclpr[[i]] = BCclpr_untrt_surv_matched_untrt_matched_trt[[2]]
    
    # if the em process converges, add 1 to the index
    i = i + 1; real_iter_ind = real_iter_ind + 1
    end_time1 <- Sys.time()
    print(paste0("one iteration lasts ", difftime(end_time1, start_time1)))
    print(difftime(end_time1, start_time1))
  }
  # out of for loop for all the samples (all in all: param_n_sim samples)
  
  
  # TODO 1.
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
  
  #TODO if only_naive_bool==TRUE, i.e. we ant only naive, stop here
  if(only_naive_bool==TRUE){
    CI_mat$SACE = mean(CI_mat$SACE)
    return(list(mat_param_estimators=mat_param_estimators, CI_mat=CI_mat))
  }
  
  # TODO 1.b
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
  coeffs_as = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_as, "[[", l)
  })
  coeffs_ns = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_ns, "[[", l)
  })
  coeffs = data.table( rbind( list.rbind(coeffs_as), list.rbind(coeffs_ns) ) )
  
  #TODO FIX THIS HERE
  print("mean and sd")
  coeffs = data.table(coeffs)
  #TODO genefilter package: coeffs[, `:=` (SD = rowSds(as.matrix(coeffs)), mean = rowMeans(coeffs))]
  coeffs = data.frame(coeffs, SD = apply(coeffs, 1, sd), mean = apply(coeffs, 1, mean))
  coeffs$parameter = c(gamma_as, gamma_ns)
  coeffs$diff = coeffs$mean - coeffs$parameter 
  coeffs$perc = coeffs$diff / abs(coeffs$parameter)
  coeffs_df = data.frame(coeffs)
  rownames(coeffs_df) = colnames(mat_gamma)
  #rownames(coeffs_df) = c("coeff_as_0", "coeff_as_1","coeff_ns_0", "coeff_ns_1")
  
  # TODO 2. mean and sd for mat_excluded_included_matching
  mat_excluded_included_matching = rbind(mat_excluded_included_matching, 
                                         mean = apply(mat_excluded_included_matching, 2, mean),
                                         sd = apply(mat_excluded_included_matching, 2, sd))
  # change order of mat_excluded_included_matching
  mat_excluded_included_matching = data.frame(mat_excluded_included_matching)
  ind1 = c(1,6,11)
  mat_excluded_included_matching = subset(mat_excluded_included_matching,
                                          select = c(ind1 , ind1+1, ind1+2, ind1+3, ind1+4))
  
  
  # TODO 3. list of repeated_as_and_pro
  mean_list_repeated_as_and_pro = 
    calculate_mean_repeated_as_and_pro(list_repeated_as_and_pro, TRUE)
  
  # TODO 6 THE SAME AS FOR list_repeated_as_and_pro to list_matched_units
  mean_list_matched_units = 
    calculate_mean_repeated_as_and_pro(list_matched_units, FALSE)
  
  
  # TODO 4. matrix of diff_abs_distance (number paer each part, so 6 in total)
  # I can calculate mean, but I think there is no need
  
  # TODO 5. list of std_diffs per all X's; size of list as the samples
  mean_list_std_mean_diff = calculate_mean_of_dfs_in_lists_of_lists(list_std_mean_diff)
  
  # TODO 7. list_means_by_subset
  mean_list_means_by_subset = calculate_mean_repeated_as_and_pro(list_means_by_subset, FALSE)
  
  
  return(list(mat_param_estimators = mat_param_estimators, 
              WLS_NOint_mat_reg_estimators = list_reg_LS$WLS_NOint_mat_reg_estimators,
              WLS_YESint_mat_reg_estimators = list_reg_LS$WLS_YESint_mat_reg_estimators,
              OLS_NOint_mat_reg_estimators = list_reg_LS$OLS_NOint_mat_reg_estimators,
              OLS_YESint_mat_reg_estimators = list_reg_LS$OLS_YESint_mat_reg_estimators,
              CI_mat = CI_mat,
              coeffs_df = coeffs_df, 
              mat_excluded_included_matching = mat_excluded_included_matching, 
              mean_list_repeated_as_and_pro = mean_list_repeated_as_and_pro,
              mat_diff_distance_aspr_asas = mat_diff_distance_aspr_asas,
              mean_list_std_mean_diff = mean_list_std_mean_diff,
              mean_list_matched_units = mean_list_matched_units,
              mean_list_means_by_subset = mean_list_means_by_subset,
              list_EM_not_conv = list_EM_not_conv,
              list_BCclpr = list_BCclpr
  ))
}


