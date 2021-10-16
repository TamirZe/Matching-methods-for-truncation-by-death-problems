# LIN REG LARGE DATASET
simulate_mtach_and_lin_reg = function(param_n, seed){
  #set.seed(100 + seed)
  list_data_for_EM_and_X = simulate_data_function(gamma_as, gamma_ns, param_n, epsilon_1_GPI = epsilon_1_GPI)
  
  data_for_EM = list_data_for_EM_and_X[[1]]; x = list_data_for_EM_and_X[[2]]
  OBS_table = list_data_for_EM_and_X[[3]]; pis = list_data_for_EM_and_X[[4]]
  vec_OBS_table = t(c(OBS_table))
  colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
  ###########################################################
  # real parameter
  # SACE = mean(data_for_EM[g_num==1 , Y1]) - mean(data_for_EM[g_num==1, Y0])
  # SACE_conditional = mean(data_for_EM[A==1 & g_num==1 , Y]) - mean(data_for_EM[A==0 & g_num==1, Y])
  # 
  # # naive estimators
  # most_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y])
  # most_naive_est_se = sqrt(( var(data_for_EM[A==1, Y]) + var(data_for_EM[A==0, Y]) ) / nrow(data_for_EM))
  # 
  # sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
  # sur_naive_est_se = sqrt(( var(data_for_EM[A==1 & S==1, Y]) + var(data_for_EM[A==0 & S==1, Y]) ) / nrow(data_for_EM))

  ###########################################################
  start_time2 <- Sys.time()
  EM_list = function_my_EM(data_for_EM, iterations, epsilon_EM)
  end_time2 <- Sys.time()
  print(paste0("function_my_EM lasts ", difftime(end_time2, start_time2)))
  dat_EM = EM_list[[1]]
  # after running EM, merge both data tables
  data_with_PS = data.table(merge(x = data_for_EM,
                                  y = subset(dat_EM, select = c(id, p_as, p_ns, p_pro, max_strata_per_subj)),
                                  by = "id", all.x = TRUE))
  
  # EM coeffs
  # coeff_as, coeff_ns, coeff_pro
  coeff_as = EM_list[[2]] ; coeff_ns = EM_list[[3]]
  #list_dat_EM[[i]] = dat_EM
  list_coeff_as[[i]] = coeff_as; list_coeff_ns[[i]] = coeff_ns
  
  #########################################################################################
  # calculating PS from the M step in the EM, not from the E step
  # the E step takes into account also the cells, and we don't want to do such thing here yet
  PS_est = cbind(exp(x%*%coeff_as), exp(x%*%coeff_ns), exp(x%*%gamma_pro))
  PS_est = PS_est / apply(PS_est, 1, sum)
  colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
  data_with_PS = data.table(data.frame(data_with_PS, PS_est))
  
  # DING estimator
  # O11_prior_ratio = pis[1] / (pis[1] + pis[3]) 
  # data_with_PS[, `:=` (O11_posterior_ratio = EMest_p_as / (EMest_p_as + EMest_p_pro),
  #                      O11_prior_ratio = O11_prior_ratio, 
  #                      W_1_as = ( EMest_p_as / (EMest_p_as + EMest_p_pro) ) / O11_prior_ratio)]
  # 
  # data_with_PS[, W_1_as_Y := W_1_as * Y]
  # DING_est = mean(data_with_PS[A==1 & S == 1, W_1_as_Y]) - mean(data_with_PS[A==0 & S == 1, Y])
  # 
  # ##### DING model assisted, 3 options, only 1 for now: 
  # DING_model_assisted_est_ps = DING_model_assisted_func(data_with_PS, x)
  
  X_sub_cols = paste0("X", c(1:(dim_x)))
  m_data = data_with_PS[S==1]

  # matching
  replace = TRUE; estimand = "ATC"; change_id = TRUE; mahal_match = 2; M=1; caliper = 0.05
  if(change_id == TRUE){
    print("change id")
    m_data$id = c(1:nrow(m_data))
  }
  vec_caliper = c(rep(1000, length(X_sub_cols[-1])), caliper)
  #w_mat = diag(length(X_sub_cols[-1]) + 1) / length(X_sub_cols[-1]) 
  w_mat = diag(length(X_sub_cols[-1]) + 1) / 
    length(X_sub_cols[-1])  * c(apply(subset(m_data, select = X_sub_cols[-1]), 2, var), 1)
  w_mat[nrow(w_mat), ncol(w_mat)]  = 0
  ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                         #, X=m_data[,"est_p_as"]
                         , X = subset(m_data, select = c(X_sub_cols[-1], "EMest_p_as"))
                         #,ties=FALSE
                         ,caliper = vec_caliper
                         ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                         ,Weight.matrix = w_mat
  )
  
  print(ATE_MATCH_PS$estimand)
  #ATE_MATCH_PS$est; ATE_MATCH_PS$se; ATE_MATCH_PS$nobs; ATE_MATCH_PS$index.dropped
  ncols  = ncol(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                       select = c("id", "p_as", "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1]))) + 1
  dt_match = data.table(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                               select = c("id", "p_as", "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1])),
                        ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control,
                        subset(m_data[ATE_MATCH_PS$index.control, ], 
                               select = c("id", "p_as", "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1])))
  colnames(dt_match)[(ncols + 1): (2 * ncols)] = 
    paste0("A0_", colnames(dt_match)[(ncols + 1): (2 * ncols)])
  colnames(dt_match)[c(ncols: (ncols+1))] = c("id_trt", "id_ctrl")
  unique(dt_match$A0_id) %>% length() == nrow(dt_match)
  unique(dt_match$id_trt) %>% length()
  identical(as.numeric(dt_match$id), dt_match$id_trt)
  sum(m_data$id %in% dt_match$id)
  # keep only S = 1
  dt_match_S1 = filter(dt_match, S == 1 & A0_S==1)
  
  ###### estimation
  #  diff of means
  SACE_matching_est = mean(dt_match_S1$Y) - mean(dt_match_S1$A0_Y)
  
  #  WLS
  WLS_NO =regression_adjusted_function(dt_match_S1, m_data, covariates = X_sub_cols[-1], 
                     reg_covariates = X_sub_cols[-1], interactions_bool = FALSE, LS="WLS")
  WLS_YES =regression_adjusted_function(dt_match_S1, m_data, covariates = X_sub_cols[-1], 
                     reg_covariates = X_sub_cols[-1], interactions_bool = TRUE, LS="WLS")
  
  #  OLS
  OLS_NO =regression_adjusted_function(dt_match_S1, m_data, covariates = X_sub_cols[-1], 
                     reg_covariates = X_sub_cols[-1], interactions_bool = FALSE, LS="OLS")
  OLS_YES = regression_adjusted_function(dt_match_S1, m_data, covariates = X_sub_cols[-1], 
                     reg_covariates = X_sub_cols[-1], interactions_bool = TRUE, LS="OLS")
  
  EST = data.frame(WLS_NO=WLS_NO[[1]][1], WLS_YES=WLS_YES[[1]][1], OLS_NO=OLS_NO[[1]][1], OLS_YES=OLS_YES[[1]][1])
  SE = data.frame(WLS_NO=WLS_NO[[1]][3], WLS_YES=WLS_YES[[1]][3], OLS_NO=OLS_NO[[1]][2], OLS_YES=OLS_YES[[1]][2])
  colnames(EST) <- colnames(SE) <- c("WLS_NO", "WLS_YES", "OLS_NO", "OLS_YES")
  return(list(EST=EST, SE=SE, SACE_matching_est=SACE_matching_est, 
              WLS_NO=WLS_NO[[1]], WLS_YES=WLS_YES[[1]], OLS_NO=OLS_NO[[1]], OLS_YES=OLS_YES[[1]]))

}

replications = 10
lst_results = list()
for(i in 1:replications){
  print(i)
  lst_results[[i]] = simulate_mtach_and_lin_reg(param_n, seed=i)
}

EST_results = data.frame(sapply(lst_results, "[[", "EST") %>% apply(2, as.numeric) %>% round(3))
EST_sum = data.frame(Mean = apply(EST_results, 1, mean), Emp_sd = apply(EST_results, 1, sd))
SE_results = data.frame(sapply(lst_results, "[[", "SE") %>% apply(2, as.numeric) %>% round(3))
SE_sum = data.frame(Est_sd = apply(SE_results, 1, mean), X = apply(SE_results, 1, sd))
EST_SE_sum = data.frame(EST_sum, SE_sum)
rownames(EST_SE_sum) = rownames(sapply(lst_results, "[[", "EST"))

