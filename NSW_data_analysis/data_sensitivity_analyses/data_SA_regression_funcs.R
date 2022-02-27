#########################################################################################
#func: regression on the matched dataset with original Y ####
regression_function_one_model = function(reg_data_matched_SA, data_reg, reg_after_match, repl = TRUE){
  
  #TODO regression wout A-X interactions
  f_wout_intercations = as.formula(paste0("Y ~ ", paste(c("A", reg_after_match), collapse = " + ")))
  model_wout_intercations = lm(f_wout_intercations, reg_data_matched_SA)  # reg_data_matched_SA
  coeffs_wout_intercations = model_wout_intercations$coefficients
  if(repl){
    vcov_wout_interactions = vcovCL(model_wout_intercations, cluster = ~ pair + id)
    se_wout_intercations = sqrt(diag(vcov_wout_interactions))
  }else{
    vcov_wout_interactions = vcovCL(model_wout_intercations, cluster = ~ pair)
    se_wout_intercations = sqrt(diag(vcov_wout_interactions))
  }
  
  TE_wout_intercations = coeffs_wout_intercations["A"]
  
  #TODO regression with A-X interactions
  f_with_intercations = as.formula(paste0("Y ~ ", paste(c("A", reg_after_match), collapse = " + "), " + ",
                                          paste(rep("A*",5), reg_after_match, collapse=" + "))) 
  model_with_intercations = lm(f_with_intercations, reg_data_matched_SA)  # reg_data_matched_SA
  mean_x = apply(subset(filter(reg_data_matched_SA, A==0), select = reg_after_match),2, mean) # reg_data_matched_SA
  coeffs_with_interactions = model_with_intercations$coefficients
  coeffs_interactions = coeffs_with_interactions[grep(":", names(coeffs_with_interactions))]
  TE_with_intercations = coeffs_with_interactions["A"] + coeffs_interactions %*% mean_x
  if(repl){
    se_with_intercations = sqrt(diag(vcovCL(model_with_intercations, cluster = ~ pair + id)))
  }else{
    se_with_intercations = sqrt(diag(vcovCL(model_with_intercations, cluster = ~ pair)))
  }
  
  
  return(list(coeffs_wout_intercations=coeffs_wout_intercations, se_wout_intercations=se_wout_intercations,
              vcov_wout_interactions=vcov_wout_interactions,
              coeffs_with_interactions=coeffs_with_interactions, se_with_intercations=se_with_intercations))
}
#########################################################################################

#########################################################################################
#func: two separate regression models (trt and ctr) on the matched dataset with original Y ####
regression_function_two_models = function(reg_data_matched_SA, data_reg, reg_after_match, repl = TRUE){
  f = as.formula(paste0("Y ~ ", paste(reg_after_match, collapse = " + ")))
  
  #TODO regression treated
  model_trt = lm(f, filter(reg_data_matched_SA, A==1)) 
  coeffs_trt = model_trt$coefficients; summ_trt = summary(model_trt)
  identical(vcovHC(model_trt, cluster = ~ id, type = "HC"), sandwich(model_trt))
  
  if(repl == TRUE){ # WLS
    vcov_beta_trt <- vcovCL(model_trt, cluster = ~ id) #  pair NO NEED, because these are only the treated
    #vcov_beta_trt2 = sandwich(model_trt)
  }else{ # OLS
    #vcov_beta_trt = vcov(model_trt) 
  }
  
  #TODO regression untreated
  model_untrt = lm(f, filter(reg_data_matched_SA, A==0)) 
  coeffs_untrt = model_untrt$coefficients; summ_untrt = summary(model_untrt)
  vcov_beta_untrt = vcov(model_untrt)
  
  return(list(coeffs_trt=coeffs_trt, vcov_beta_trt=vcov_beta_trt,
              coeffs_untrt=coeffs_untrt, vcov_beta_untrt=vcov_beta_untrt))
}
#########################################################################################

#########################################################################################
#func: PPI -            predictions for units from O(0,1), plugging A={0,1} ####
SACE_estimation_1LEARNER_PPI = function(matched_data, reg_after_match, eps_sensi_PPI=1, 
                                        coeffs_regression_one_model, coeffs_regression_two_models, two_models_bool = TRUE){
  # func: weighted SE with replacements
  func_est_wtd_var_matching_w_rep = function(matched_data, f_alpha){ # f_alpha=1 gives regular weighted var
    M = nrow(filter(matched_data, A==0))
    est_var_Y0_mean = (1/M) * var(filter(matched_data, A==0)$Y) 
    est_var_Y1_unq = var(filter(matched_data, A==1)$Y) 
    sum_weights_squared_alpha = sum( ( (filter(matched_data, A==1)$w)^2 ) / ( f_alpha^2 ) ) 
    est_var_Y1_mean = (1/M^2) * sum_weights_squared_alpha * est_var_Y1_unq
    #wtd.var(filter(matched_data, A==1)$Y, filter(matched_data, A==1)$w, na.rm = TRUE) / M
    est_var_after_match_w_rep = est_var_Y1_mean + est_var_Y0_mean 
    return(est_var_after_match_w_rep)
  }
  
  #TODO CRUDE
  # crude diff estimator, adjusted by eps_PPI
  A1_S1_data = filter(matched_data, A==1)
  A0_S1_data = filter(matched_data, A==0)
  #attach(A1_S1_data); #attach(A0_S1_data)
  f_alpha_A1 = ( A1_S1_data$e_1_as + ((1 - A1_S1_data$e_1_as) * eps_sensi_PPI) )
  f_alpha_A0 = ( A0_S1_data$e_1_as + ((1 - A0_S1_data$e_1_as) * eps_sensi_PPI) )
  crude_Y1_adj = A1_S1_data$Y / f_alpha_A1
  crude_est_adj = mean(crude_Y1_adj) - mean(A0_S1_data$Y)
  crude_est_adj_se = sqrt( func_est_wtd_var_matching_w_rep(matched_data, f_alpha_A1) )
  
  # wout interactions
  # regression coefficients
  coeffs_regression_wout_inter = coeffs_regression_one_model$coeffs_wout_intercations
  coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
  coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
  vcov_beta_wout_inter = coeffs_regression_one_model$vcov_wout_interactions
  # predictions wout interactions
  #print(apply(matched_data, 2, class))
  A0_S1_data$Y1_pred =  coeff_A_wout_inter + 
    as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match)) %>% apply(2, as.numeric)) %*% coeffs_wout_inter
  A0_S1_data$Y1_pred_adj =  A0_S1_data$Y1_pred / f_alpha_A0
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred_adj) - mean(A0_S1_data$Y0_pred)
  
  #TODO SE with f_alpha_A0
  g_alpha_A0 = 1 / f_alpha_A0 # g_alpha_A0 = (1 - f_alpha_A0) / f_alpha_A0 
  X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
  X_tilde = cbind(intercept = (1 - f_alpha_A0) * X[,"intercept"], A = 1, (1 - f_alpha_A0) * X[,-which(colnames(X)=="intercept")]) # X_tilde = cbind(intercept = X[,"intercept"], A = (1/(1-f_alpha_A0)), X[,-intercept])
  unit_diff_vec = g_alpha_A0 * (X_tilde %*% coeffs_regression_wout_inter)
  SACE_1LEARNER_adj_calc = mean(unit_diff_vec)
  print(SACE_1LEARNER_adj - SACE_1LEARNER_adj_calc)
  
  X_tilde_g_alpha = g_alpha_A0 * X_tilde
  sum_X_tilde_g_alpha = apply(X_tilde_g_alpha, 2, sum)
  var_sum_units_A0_S1_tilde_g_alpha = t(sum_X_tilde_g_alpha) %*%  vcov_beta_wout_inter %*% sum_X_tilde_g_alpha 
  SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_g_alpha )
  # if(eps_sensi_PPI != 1){
  #   SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_g_alpha )
  #   }else{
  #   SACE_1LEARNER_adj_se = as.numeric(coeffs_regression_one_model$se_wout_intercations["A"])
  # }
  
  # with interactions
  #1. if we calculated 1 regression model
  # TODO add SE. for now im not really using it, thus I didnt add SEs
  if(two_models_bool == FALSE){
    # regression prediction with interactions
    # with interactions
    coeffs_regression_with_inter = coeffs_regression_one_model$coeffs_with_interactions
    coeff_A_with_inter = coeffs_regression_with_inter["A"]; intercept_with_inter = coeffs_regression_with_inter["(Intercept)"]
    coeffs_with_inter_A1 = coeffs_regression_with_inter[-grep("A|(Intercept)", names(coeffs_regression_wout_inter))]
    coeffs_with_inter_A0 = coeffs_regression_with_inter[-grep(":|A", names(coeffs_regression_with_inter))]
    
    # predictions with interactions
    #mean_x = apply(subset(A0_S1_data, select = reg_after_match),2, mean)
    A0_S1_data$Y1_pred_inter = intercept_with_inter + coeff_A_with_inter +
      as.matrix(subset(A0_S1_data, select = reg_after_match)) %*% 
      (coeffs_with_inter_A1[-grep(":",names(coeffs_with_inter_A1))] + coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))])
    A0_S1_data$Y1_pred_inter_adj =  A0_S1_data$Y1_pred_inter / f_alpha_A0
    A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
    
    #SACE_1LEARNER_inter_chck = coeff_A_with_inter +
    #  sum( apply(subset(A0_S1_data, select = reg_after_match), 2, mean) * coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))] ) 
    SACE_1LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter)
    
    return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, crude_est_adj=crude_est_adj,
             SACE_1LEARNER=SACE_1LEARNER, SACE_1LEARNER_inter=SACE_1LEARNER_inter))
  }
  
  #2. if we calculated 2 regression models
  if(two_models_bool == TRUE){
    beta_trt = coeffs_regression_two_models$coeffs_trt; vcov_beta_trt = coeffs_regression_two_models$vcov_beta_trt
    beta_untrt = coeffs_regression_two_models$coeffs_untrt; vcov_beta_untrt = coeffs_regression_two_models$vcov_beta_untrt
    #f_alpha_A0 = ( A0_S1_data$e_1_as + ((1 - A0_S1_data$e_1_as) * eps_sensi_PPI) )
    X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
    X_adj = X / f_alpha_A0
    A0_S1_data = filter(matched_data, A==0 & S==1); A1_S1_data = filter(matched_data, A==1 & S==1)
    
    A0_S1_data$Y1_pred_inter = X %*% beta_trt
    A0_S1_data$Y1_pred_inter_adj = A0_S1_data$Y1_pred_inter / f_alpha_A0
    A0_S1_data$Y0_pred_inter = X %*% beta_untrt
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter)
    
    # SE, considering f_alpha_A0(X) as constant, i.e not considering gammas variance 
    sum_X = apply(X, 2, sum); sum_X_adj = apply(X_adj, 2, sum)
    var_sum_units_A0_S1 = t(sum_X_adj) %*%  vcov_beta_trt %*% sum_X_adj + t(sum_X) %*%  vcov_beta_trt %*% sum_X 
    SACE_1LEARNER_inter_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1 )
    
    return(list(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_adj_se=SACE_1LEARNER_adj_se,
                SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, SACE_1LEARNER_inter_adj_se=SACE_1LEARNER_inter_adj_se,
                crude_est_adj=crude_est_adj, crude_est_adj_se=crude_est_adj_se))
  }
}
#########################################################################################

#########################################################################################
#func: monotonoicity -  predictions for units from O(0,1), plugging A={0,1}
SACE_estimation_1LEARNER_mono = function(matched_data, reg_after_match, alpha0_mono=1, xi=0, 
                                         coeffs_regression_one_model, coeffs_regression_two_models, two_models_bool = TRUE){
  # func: weighted SE with replacements
  func_est_wtd_var_matching_w_rep = function(matched_data, f_alpha){ # f_alpha=1 gives regular weighted var
    M = nrow(filter(matched_data, A==0))
    est_var_Y0_mean = (1/M) * var(filter(matched_data, A==0)$Y) 
    est_var_Y1_unq = var(filter(matched_data, A==1)$Y) 
    sum_weights_squared = sum( (filter(matched_data, A==1)$w)^2 ) 
    est_var_Y1_mean = (1/M^2) * sum_weights_squared * est_var_Y1_unq
    #wtd.var(filter(matched_data, A==1)$Y, filter(matched_data, A==1)$w, na.rm = TRUE) / M
    est_var_after_match_w_rep = est_var_Y1_mean + (f_alpha^2 * est_var_Y0_mean) 
    return(est_var_after_match_w_rep)
  }
  
  #TODO CRUDE
  # crude diff estimator, adjusted by alpha0_mono
  A1_S1_data = filter(matched_data, A==1)
  A0_S1_data = filter(matched_data, A==0)
  f_alpha0 = (1 + xi) / (1 + (xi * alpha0_mono)) 
  
  crude_Y0_adj = A0_S1_data$Y * f_alpha0
  crude_est_adj = mean(A1_S1_data$Y) - mean(crude_Y0_adj)
  crude_est_adj_se = sqrt( func_est_wtd_var_matching_w_rep(matched_data, f_alpha0) )
  
  # wout interactions
  # regression coefficients
  coeffs_regression_wout_inter = coeffs_regression_one_model$coeffs_wout_intercations
  coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
  coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
  vcov_beta_wout_inter = coeffs_regression_one_model$vcov_wout_interactions
  # predictions wout interactions
  A0_S1_data$Y1_pred =  coeff_A_wout_inter + 
    as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match)) %>% apply(2, as.numeric)) %*% coeffs_wout_inter
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  A0_S1_data$Y0_pred_adj =  A0_S1_data$Y0_pred * f_alpha0
  
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred_adj)
  
  #TODO SE with f_alpha0
  g_alpha0 = 1 
  X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
  X_tilde = cbind(intercept = (1 - f_alpha0) * X[,"intercept"], A = 1, (1 - f_alpha0) * X[,-which(colnames(X)=="intercept")])
  unit_diff_vec = g_alpha0 * (X_tilde %*% coeffs_regression_wout_inter)
  SACE_1LEARNER_adj_calc = mean(unit_diff_vec)
  print(SACE_1LEARNER_adj - SACE_1LEARNER_adj_calc)
  
  X_tilde_g_alpha0 = g_alpha0 * X_tilde
  sum_X_tilde_g_alpha0 = apply(X_tilde_g_alpha0, 2, sum)
  var_sum_units_A0_S1_tilde_g_alpha0 = t(sum_X_tilde_g_alpha0) %*%  vcov_beta_wout_inter %*% sum_X_tilde_g_alpha0 
  SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_g_alpha0 )
  
  # with interactions
  #1. if we calculated 1 regression model
  # TODO add SE. for now im not really using it, thus I didnt add SEs
  if(two_models_bool == FALSE){
    # with interactions
    coeffs_regression_with_inter = coeffs_regression_one_model$coeffs_with_interactions
    coeff_A_with_inter = coeffs_regression_with_inter["A"]; intercept_with_inter = coeffs_regression_with_inter["(Intercept)"]
    coeffs_with_inter_A1 = coeffs_regression_with_inter[-grep("A|(Intercept)", names(coeffs_regression_wout_inter))]
    coeffs_with_inter_A0 = coeffs_regression_with_inter[-grep(":|A", names(coeffs_regression_with_inter))]
    
    # predictions with interactions
    A0_S1_data$Y1_pred_inter = intercept_with_inter + coeff_A_with_inter +
      as.matrix(subset(A0_S1_data, select = reg_after_match)) %*% 
      (coeffs_with_inter_A1[-grep(":",names(coeffs_with_inter_A1))] + coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))])
    A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
    A0_S1_data$Y0_pred_inter_adj =  A0_S1_data$Y0_pred_inter * f_alpha0 
    SACE_1LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter_adj)
    
    return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, crude_est_adj=crude_est_adj,
             SACE_1LEARNER=SACE_1LEARNER, SACE_1LEARNER_inter=SACE_1LEARNER_inter))
  }
  
  #2. if we calculated 2 regression models
  if(two_models_bool == TRUE){
    beta_trt = coeffs_regression_two_models$coeffs_trt; vcov_beta_trt = coeffs_regression_two_models$vcov_beta_trt
    beta_untrt = coeffs_regression_two_models$coeffs_untrt; vcov_beta_untrt = coeffs_regression_two_models$vcov_beta_untrt
    X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
    X_adj = X * f_alpha0
    
    A0_S1_data$Y1_pred_inter = X %*% beta_trt
    A0_S1_data$Y0_pred_inter = X %*% beta_untrt
    A0_S1_data$Y0_pred_inter_adj = A0_S1_data$Y0_pred_inter * f_alpha0
    SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter) -  mean(A0_S1_data$Y0_pred_inter_adj)
    
    # SE, considering f_alpha0(X) as constant, i.e not considering gammas variance 
    sum_X = apply(X, 2, sum); sum_X_adj = apply(X_adj, 2, sum)
    var_sum_units_A0_S1 = t(sum_X) %*%  vcov_beta_trt %*% sum_X + t(sum_X_adj) %*%  vcov_beta_trt %*% sum_X_adj 
    SACE_1LEARNER_inter_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1 )
    
    return(list(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_adj_se=SACE_1LEARNER_adj_se,
                SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, SACE_1LEARNER_inter_adj_se=SACE_1LEARNER_inter_adj_se,
                crude_est_adj=crude_est_adj, crude_est_adj_se=crude_est_adj_se))
  }
}
#########################################################################################
