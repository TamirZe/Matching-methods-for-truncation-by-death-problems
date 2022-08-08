#########################################################################################
#func: regression on the matched dataset with original Y ####
regression_function_one_model = function(reg_data_matched, reg_after_match, repl=TRUE){
  
  #TODO regression wout A-X interactions
  f_wout_intercations = as.formula(paste0("Y ~ ", paste(c("A", reg_after_match), collapse = " + ")))
  model_wout_intercations = lm(f_wout_intercations, reg_data_matched) 
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
  model_with_intercations = lm(f_with_intercations, reg_data_matched) 
  mean_x = apply(subset(filter(reg_data_matched, A==0), select = reg_after_match),2, mean) 
  coeffs_with_interactions = model_with_intercations$coefficients
  coeffs_interactions = coeffs_with_interactions[grep(":", names(coeffs_with_interactions))]
  TE_with_intercations = coeffs_with_interactions["A"] + coeffs_interactions %*% mean_x
  if(repl){
    vcov_with_interactions = vcovCL(model_with_intercations, cluster = ~ pair + id)
    se_with_intercations = sqrt(diag(vcov_with_interactions))
  }else{
    vcov_with_interactions = vcovCL(model_with_intercations, cluster = ~ pair)
    se_with_intercations = sqrt(diag(vcov_with_interactions))
  }
  
  
  return(list(coeffs_wout_intercations=coeffs_wout_intercations, vcov_wout_interactions=vcov_wout_interactions,
              se_wout_intercations=se_wout_intercations,
              coeffs_with_interactions=coeffs_with_interactions, vcov_with_interactions=vcov_with_interactions,
              se_with_intercations=se_with_intercations
              ))
}
#########################################################################################

#########################################################################################
#func: two separate regression models (trt and ctr separately) on the matched dataset with original Y ####
regression_function_two_models = function(reg_data_matched, reg_after_match, repl=TRUE){
  f = as.formula(paste0("Y ~ ", paste(reg_after_match, collapse = " + ")))
  
  #TODO regression treated
  model_trt = lm(f, filter(reg_data_matched, A==1)) 
  coeffs_trt = model_trt$coefficients; summ_trt = summary(model_trt)
  #identical(vcovHC(model_trt, cluster = ~ id, type = "HC"), sandwich(model_trt))
  
  if(repl == TRUE){ # WLS
    vcov_beta_trt <- vcovCL(model_trt, cluster = ~ id) #  NO NEED for pair, because these are only the treated
    #vcov_beta_trt2 = sandwich(model_trt)
  }else{ # OLS
    #vcov_beta_trt = vcov(model_trt) 
  }
  
  #TODO regression untreated
  model_untrt = lm(f, filter(reg_data_matched, A==0)) 
  coeffs_untrt = model_untrt$coefficients; summ_untrt = summary(model_untrt)
  vcov_beta_untrt = vcov(model_untrt)
  
  return(list(coeffs_trt=coeffs_trt, vcov_beta_trt=vcov_beta_trt,
              coeffs_untrt=coeffs_untrt, vcov_beta_untrt=vcov_beta_untrt))
}
#########################################################################################

#########################################################################################
# func: weighted SE with replacements
func_est_wtd_var_matching_w_rep = function(matched_data, alpha, xi=NULL, SA_bool){ # f_alpha=1 gives regular weighted var
  M = nrow(filter(matched_data, A==0)) #= sum(filter(matched_data, A==1)$w)
  est_var_Y0_mean = (1/M) * var(filter(matched_data, A==0)$Y) 
  est_var_Y1_unq = var(filter(matched_data, A==1)$Y) 
  if(SA_bool == "PPI"){
    f_alpha1 = 1 / ( filter(matched_data, A==1)$pi_tilde_as1 + ((1 - filter(matched_data, A==1)$pi_tilde_as1) * alpha) )
    sum_weights_squared_alpha1 = sum( ( f_alpha1^2 ) * ( (filter(matched_data, A==1)$w)^2 ) ) 
    est_var_Y1_mean = (1/M^2) * est_var_Y1_unq * sum_weights_squared_alpha1 
    #wtd.var(filter(matched_data, A==1)$Y, filter(matched_data, A==1)$w, na.rm = TRUE) / M
    est_var_after_match_w_rep = est_var_Y1_mean + est_var_Y0_mean  #TODO 06.08.2022 what about -2*COV 
  }
  if(SA_bool == "mono"){
    f_xi_alpha0 = (1 + xi) / (1 + (xi * alpha0))
    sum_weights_squared = sum( (filter(matched_data, A==1)$w)^2 ) 
    est_var_Y1_mean = (1/M^2) * est_var_Y1_unq * sum_weights_squared 
    est_var_after_match_w_rep = est_var_Y1_mean + (f_xi_alpha0^2 * est_var_Y0_mean) #TODO 06.08.2022 what about -2*COV 
  }
  return(est_var_after_match_w_rep)
}
#########################################################################################

#########################################################################################
#func: PPI -            predictions for units from O(0,1) (plugging A={0,1}) ####
SACE_estimation_LEARNER_PPI = function(reg_data_matched, reg_after_match, alpha1=1, 
            coeffs_regression_one_model, coeffs_regression_two_models, two_models_bool=TRUE){
  
  # CRUDE ####
  # crude diff estimator, adjusted by f_alpha_A1
  A1_S1_data = filter(reg_data_matched, A==1)
  A0_S1_data = filter(reg_data_matched, A==0)
  f_alpha1_A1 = 1 / ( A1_S1_data$pi_tilde_as1 + ((1 - A1_S1_data$pi_tilde_as1) * alpha1) )
  f_alpha1_A0 = 1 / ( A0_S1_data$pi_tilde_as1 + ((1 - A0_S1_data$pi_tilde_as1) * alpha1) )
  crude_Y1_adj = f_alpha1_A1 * A1_S1_data$Y 
  crude_est_adj = mean(crude_Y1_adj) - mean(A0_S1_data$Y)
  matched_data = reg_data_matched %>% distinct(id, .keep_all = TRUE) %>% subset(select = -pair)
  #weighted.mean(filter(matched_data, A==1)$Y, filter(matched_data, A==1)$w)
  crude_est_adj_se = sqrt( func_est_wtd_var_matching_w_rep(matched_data=matched_data, alpha=alpha1, SA_bool="PPI") )
  
  # wout interactions - one model  ####
  # regression coefficients from coeffs_regression_one_model - a model that was fitted using both treated and untreated
  coeffs_regression_wout_inter = coeffs_regression_one_model$coeffs_wout_intercations
  coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
  coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
  vcov_beta_wout_inter = coeffs_regression_one_model$vcov_wout_interactions
  # predictions wout interactions
  #print(apply(reg_data_matched, 2, class))
  A0_S1_data$Y1_pred =  coeff_A_wout_inter + 
    as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match)) %>% apply(2, as.numeric)) %*% coeffs_wout_inter
  A0_S1_data$Y1_pred_adj = f_alpha1_A0 * A0_S1_data$Y1_pred
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred_adj) - mean(A0_S1_data$Y0_pred)
  
  #TODO SE with f_alpha1_A0
  X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
  X_tilde = cbind(intercept = ( 1 - (1/f_alpha1_A0) ) * X[,"intercept"], A = 1, 
        ( 1 - (1/f_alpha1_A0) ) * X[,-which(colnames(X)=="intercept")]) 
  # checkings
  unit_diff_vec = f_alpha1_A0 * (X_tilde %*% coeffs_regression_wout_inter)
  SACE_1LEARNER_adj_calc = mean(unit_diff_vec)
  print(SACE_1LEARNER_adj - SACE_1LEARNER_adj_calc)
  
  X_tilde_f_alpha1 = f_alpha1_A0 * X_tilde
  sum_X_tilde_f_alpha1 = apply(X_tilde_f_alpha1, 2, sum)
  var_sum_units_A0_S1_tilde_f_alpha1 = t(sum_X_tilde_f_alpha1) %*%  vcov_beta_wout_inter %*% sum_X_tilde_f_alpha1 
  SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_f_alpha1 )
  # if(alpha1 != 1){ SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_f_alpha1 )
  #   }else{ SACE_1LEARNER_adj_se = as.numeric(coeffs_regression_one_model$se_wout_intercations["A"]) }
  
  
  # with interactions ####
  #1. if we calculated 1 regression model
  # TODO add SE. for now I'm not really using it, thus I didnt add SEs
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
    A0_S1_data$Y1_pred_inter_adj =  f_alpha1_A0 * A0_S1_data$Y1_pred_inter
    A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
    
    #SACE_1LEARNER_inter_chck = coeff_A_with_inter +
    #  sum( apply(subset(A0_S1_data, select = reg_after_match), 2, mean) * coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))] ) 
    SACE_LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
    SACE_LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter)
    
    return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_LEARNER_inter_adj=SACE_LEARNER_inter_adj, 
             crude_est_adj=crude_est_adj,
             SACE_1LEARNER=SACE_1LEARNER, SACE_LEARNER_inter=SACE_LEARNER_inter))
  }
  
  #2. if we calculated 2 regression models
  if(two_models_bool == TRUE){
    # regression coefficients from coeffs_regression_two_models - models that were fitted separately on treated and untreated
    beta_trt = coeffs_regression_two_models$coeffs_trt
    vcov_beta_trt = coeffs_regression_two_models$vcov_beta_trt
    beta_untrt = coeffs_regression_two_models$coeffs_untrt
    vcov_beta_untrt = coeffs_regression_two_models$vcov_beta_untrt
    X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
    X_f_alpha1 = f_alpha1_A0 * X 
    
    A0_S1_data$Y1_pred_inter = X %*% beta_trt
    A0_S1_data$Y1_pred_inter_adj = f_alpha1_A0 * A0_S1_data$Y1_pred_inter
    A0_S1_data$Y0_pred_inter = X %*% beta_untrt
    SACE_LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter)
    
    # SE, considering f_alpha1_A0(X) as constant, i.e not considering gammas variance 
    sum_X = apply(X, 2, sum)
    sum_X_f_alpha1 = apply(X_f_alpha1, 2, sum)
    
    var_sum_units_A0_S1 = ( t(sum_X_f_alpha1) %*% vcov_beta_trt %*% sum_X_f_alpha1 ) + 
      ( t(sum_X) %*% vcov_beta_untrt %*% sum_X ) # this is the formula if estimated coefficients from the two model are independent
    SACE_LEARNER_inter_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1 )
    
    return(list(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_adj_se=SACE_1LEARNER_adj_se,
                SACE_LEARNER_inter_adj=SACE_LEARNER_inter_adj, SACE_LEARNER_inter_adj_se=SACE_LEARNER_inter_adj_se,
                crude_est_adj=crude_est_adj, crude_est_adj_se=crude_est_adj_se))
  }
}
#########################################################################################

#########################################################################################
#func: monotonoicity -  predictions for units from O(0,1) (plugging A={0,1})
SACE_estimation_LEARNER_mono = function(reg_data_matched, reg_after_match, alpha0=1, xi=0, 
           coeffs_regression_one_model, coeffs_regression_two_models, two_models_bool=TRUE){
  
  #TODO CRUDE ####
  # crude diff estimator, adjusted by alpha0
  A1_S1_data = filter(reg_data_matched, A==1)
  A0_S1_data = filter(reg_data_matched, A==0)
  f_xi_alpha0 = (1 + xi) / (1 + (xi * alpha0)) 
  
  crude_Y0_adj = f_xi_alpha0 * A0_S1_data$Y
  crude_est_adj = mean(A1_S1_data$Y) - mean(crude_Y0_adj)
  matched_data = reg_data_matched %>% distinct(id, .keep_all = TRUE) %>% subset(select = -pair)
  crude_est_adj_se = sqrt( func_est_wtd_var_matching_w_rep(matched_data=matched_data, alpha=alpha0, xi=xi, SA_bool="mono") )
  
  # wout interactions - one model  ########
  # regression coefficients
  coeffs_regression_wout_inter = coeffs_regression_one_model$coeffs_wout_intercations
  coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
  coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
  vcov_beta_wout_inter = coeffs_regression_one_model$vcov_wout_interactions
  # predictions wout interactions
  A0_S1_data$Y1_pred =  coeff_A_wout_inter + 
    as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match)) %>% apply(2, as.numeric)) %*% coeffs_wout_inter
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  A0_S1_data$Y0_pred_adj = f_xi_alpha0 * A0_S1_data$Y0_pred 
  
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  coeff_A_wout_inter - SACE_1LEARNER
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred_adj)
  
  #TODO SE with f_xi_alpha0
  X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
  X_tilde = cbind(intercept = (1 - f_xi_alpha0) * X[,"intercept"], A = 1, 
                  (1 - f_xi_alpha0) * X[,-which(colnames(X)=="intercept")])
  # checkings
  unit_diff_vec = (X_tilde %*% coeffs_regression_wout_inter)
  SACE_1LEARNER_adj_calc = mean(unit_diff_vec)
  print(SACE_1LEARNER_adj - SACE_1LEARNER_adj_calc)
  
  #X_tilde_g_alpha0 = g_alpha0 * X_tilde
  sum_X_tilde_f_xi_alpha0 = apply(X_tilde, 2, sum) # apply(X_tilde_g_alpha0, 2, sum)
  var_sum_units_A0_S1_tilde_f_xi_alpha0 = t(sum_X_tilde_f_xi_alpha0) %*%  vcov_beta_wout_inter %*% sum_X_tilde_f_xi_alpha0 
  SACE_1LEARNER_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1_tilde_f_xi_alpha0 )
  
  # with interactions ####
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
    A0_S1_data$Y0_pred_inter_adj =  f_xi_alpha0 * A0_S1_data$Y0_pred_inter  
    SACE_LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
    SACE_LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter_adj)
    
    return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_LEARNER_inter_adj=SACE_LEARNER_inter_adj, crude_est_adj=crude_est_adj,
             SACE_1LEARNER=SACE_1LEARNER, SACE_LEARNER_inter=SACE_LEARNER_inter))
  }
  
  #2. if we calculated 2 regression models
  if(two_models_bool == TRUE){
    beta_trt = coeffs_regression_two_models$coeffs_trt
    vcov_beta_trt = coeffs_regression_two_models$vcov_beta_trt
    beta_untrt = coeffs_regression_two_models$coeffs_untrt
    vcov_beta_untrt = coeffs_regression_two_models$vcov_beta_untrt
    X = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) 
    X_f_xi_alpha0 = f_xi_alpha0 * X
    
    A0_S1_data$Y1_pred_inter = X %*% beta_trt
    A0_S1_data$Y0_pred_inter = X %*% beta_untrt
    A0_S1_data$Y0_pred_inter_adj = f_xi_alpha0 * A0_S1_data$Y0_pred_inter
    SACE_LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter_adj)
    
    # SE, considering f_xi_alpha0(X) as constant, i.e not considering gammas variance 
    sum_X = apply(X, 2, sum)
    sum_X_f_xi_alpha0 = apply(X_f_xi_alpha0, 2, sum)
    var_sum_units_A0_S1 = (t(sum_X) %*%  vcov_beta_trt %*% sum_X) + 
                          (t(sum_X_f_xi_alpha0) %*%  vcov_beta_untrt %*% sum_X_f_xi_alpha0) # = f_xi_alpha0^2 * (t(sum_X) %*%  vcov_beta_untrt %*% sum_X) 
    SACE_LEARNER_inter_adj_se = sqrt( (1 / nrow(A0_S1_data)^2) * var_sum_units_A0_S1 )
    
    return(list(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_adj_se=SACE_1LEARNER_adj_se,
                SACE_LEARNER_inter_adj=SACE_LEARNER_inter_adj, SACE_LEARNER_inter_adj_se=SACE_LEARNER_inter_adj_se,
                crude_est_adj=crude_est_adj, crude_est_adj_se=crude_est_adj_se))
  }
}
#########################################################################################


