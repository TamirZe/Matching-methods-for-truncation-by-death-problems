regression_adjusted_function = function(m_data, matched_data, matched_pairs,
                                        reg_covariates, interactions_bool, LS, mu_x_fixed=FALSE, x_as=NULL){
  
  estimation_with_interactions = function(lin_reg_fit_matched, coeffs_table, data_matched_reg,  
                                          lin_reg_fit_wls_matched=NULL, coeffs_table_wls=NULL, 
                                          LS_reg_name, reg_covariates, mu_x_fixed=FALSE, x_as=NULL){
    coeffs = coeffs_table[,1]
    names_coeffs_interactions = c("A", names(coeffs)[grep(":", names(coeffs))])
    coeffs_interactions = coeffs[names_coeffs_interactions]
    var_mu_x = apply(subset(filter(data_matched_reg, A==0), select = reg_covariates), 2, var) / 
      nrow(filter(data_matched_reg, A==0))
    
    x = subset(data_matched_reg, select = c("A", reg_covariates))
    # the 1 is for the add to intercerpt for the treated (which is the coeff of A)
    if(mu_x_fixed == TRUE){mu_x_vec=x_as}else{
      mu_x_vec = c(A = 1, apply(subset(filter(x, A==0), select = -A), 2, mean))
    }
    beta_est = mu_x_vec %*% coeffs_interactions
    
    if(LS_reg_name=="OLS"){
      cluster_coeffs = data.frame(coef_test(lin_reg_fit_matched, vcov = "CR1", cluster = data_matched_reg$pair))
      vcov_clstr_interactions = vcovCL(lin_reg_fit_matched, cluster = ~ pair)[names_coeffs_interactions, names_coeffs_interactions]
      vcov_clstr_interactions2 <- 
        vcovCR(lin_reg_fit_matched, cluster = data_matched_reg$pair, type = "CR1")[names_coeffs_interactions, names_coeffs_interactions] 
      se_beta_est_clstr <- 
        sqrt( matrix(mu_x_vec, nrow=1) %*% vcov_clstr_interactions %*% matrix(mu_x_vec, ncol=1) )
      
      # DELTA METHOD for OLS
      se_beta_est_clstr_DM = sqrt( se_beta_est_clstr^2 + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
      estimators = c(beta_est, se_beta_est_clstr_DM)
      }
    
    if(LS_reg_name=="WLS"){
      se = coeffs_table[,2]
      se_interactions = se[names_coeffs_interactions]
      vcov = vcov(lin_reg_fit_wls_matched) 
      vcov_interactions_naive = vcov[names_coeffs_interactions, names_coeffs_interactions]
      var_beta_est_naive = matrix(mu_x_vec, nrow=1) %*% vcov_interactions_naive %*% matrix(mu_x_vec, ncol=1) 
      se_beta_est_naive = sqrt(var_beta_est_naive)
      se_beta_est_naive_DM = sqrt( var_beta_est_naive + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
     
       # robust clustered SE estimator
      vcov_clstr_interactions = vcovCL(lin_reg_fit_matched, cluster = ~ pair + id)[names_coeffs_interactions, names_coeffs_interactions]
      se_beta_est_clstr <- 
        sqrt( matrix(mu_x_vec, nrow=1) %*% vcov_clstr_interactions %*% matrix(mu_x_vec, ncol=1) )
      vcov_sand_interactions = sandwich(lin_reg_fit_wls_matched)[names_coeffs_interactions, names_coeffs_interactions]
      se_beta_est_sand <- 
        sqrt( matrix(mu_x_vec, nrow=1) %*% vcov_sand_interactions %*% matrix(mu_x_vec, ncol=1) )
      
      # DELTA METHOD for WLS
      se_beta_est_clstr_DM = sqrt( se_beta_est_clstr^2 + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
      se_beta_est_sand_DM = sqrt( se_beta_est_sand^2 + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
    
      estimators = c(beta_est, se_beta_est_naive_DM, se_beta_est_sand_DM, se_beta_est_clstr_DM) 
  }
    return(estimators)
}
  
  #(comparedf(matched_data[,colnames(reg_data_matched)] %>% arrange(pair,id), reg_data_matched %>% arrange(pair,id))) #identical 
  reg_data_matched = matched_data %>% arrange(pair, A)
  
  # reg model's formula - wout and with A-X interaction
  if(interactions_bool == FALSE){
    print("NO interactions in model")
    f = as.formula(paste0("Y ~ ", paste(c("A", reg_covariates), collapse = " + ")))
  }else{
    print("interactions in model")
    # if only A in the model, there are actually no interactions
    if(is.null(reg_covariates)){
      print("only A in the model, there are actually no interactions")
      f = as.formula(paste0("Y ~ ", paste(c("A", reg_covariates), collapse = " + "))) 
    }else{
      f = as.formula(paste0("Y ~ ", paste(c("A", reg_covariates), collapse = " + "), " + ",
                            paste("A *", reg_covariates, collapse=" + "))) 
    }
  }
  
  # OLS
  if(LS=="OLS"){
    print("OLS")
    lin_reg_matched = lm(formula = f , data = reg_data_matched)
    OLS_sum = summary(lin_reg_matched)
    coeffs_table = OLS_sum$coefficients
    
    # WITHOUT INTERACTIONS
    if(interactions_bool == FALSE){
      beta_est = coeffs_table["A","Estimate"]
      se_beta_est = coeftest(lin_reg_matched, vcov. = vcovCL, cluster = ~ pair)["A", "Std. Error", drop = FALSE]
      #se_beta_est2 = data.frame(coef_test(lin_reg_matched, vcov = "CR1", cluster = reg_data_matched$pair))["A", "SE"]
      estimator_and_se = c(beta_est, se_beta_est)
    }
    # WITH INTERACTIONS
    if(interactions_bool == TRUE){
      estimator_and_se = estimation_with_interactions(
        lin_reg_fit_matched=lin_reg_matched, coeffs_table=coeffs_table, data_matched_reg=reg_data_matched, 
        LS_reg_name="OLS", reg_covariates=reg_covariates, mu_x_fixed=mu_x_fixed, x_as=x_as)
    }
    
    estimator_and_se = data.frame(t(estimator_and_se))
    colnames(estimator_and_se) = c("OLS_estimator", "OLS_se")
    CI_LS = round(estimator_and_se$OLS_estimator + c(-1,1) * 1.96 * estimator_and_se$OLS_se, 3)
    CI_LS = paste(CI_LS, sep = ' ', collapse = " , ")
  }
  
  # WLS
  if(LS=="WLS"){
    print("WLS")
    # adjust for replacements
    weights = data.frame(table(matched_pairs$id))
    colnames(weights) = c("id", "weight")
    weights = apply(weights, 2, as.numeric)
    reg_data_matched_weights = merge(weights, reg_data_matched, by = "id") %>% arrange(id)
    reg_data_matched_weights = subset(reg_data_matched_weights, !duplicated(reg_data_matched_weights$id))
    
    lin_reg_matched = lm(formula = f , data = reg_data_matched)
    lin_reg_matched_wls = lm(formula = f , data = reg_data_matched_weights, weights = reg_data_matched_weights$weight)
    
    # WLS regression treatment coefficients estimator
    # WOUT INTERACTIONS
    if(interactions_bool == FALSE){
      beta_est <- lin_reg_matched$coefficients[2]
      # WLS regression treatment SE coefficients estimator
      se_beta_est_naive <- summary(lin_reg_matched)$coefficients[2, 2]
      # robust sandwich SE estimator
      se_beta_est_sand <- sqrt(sandwich(lin_reg_matched_wls)[2,2]) # robust is actually sandwich
      # robust clustered SE estimator
      se_beta_est_clstr <- coeftest(lin_reg_matched, vcov. = vcovCL, cluster = ~ pair + id)["A","Std. Error",drop = FALSE]
      estimator_and_se <- c(beta_est, se_beta_est_naive, se_beta_est_sand, se_beta_est_clstr)
    }
    # WITH INTERACTIONS
    if(interactions_bool == TRUE){
      # summary of WLS
      sum_matched = summary(lin_reg_matched)
      coeffs_table = sum_matched$coefficients
      sum_matched_wls = summary(lin_reg_matched_wls)
      coeffs_table_wls = sum_matched_wls$coefficients
      
      estimator_and_se = estimation_with_interactions(
        lin_reg_fit_matched=lin_reg_matched, coeffs_table=coeffs_table, data_matched_reg=reg_data_matched, 
        lin_reg_fit_wls_matched=lin_reg_matched_wls, coeffs_table_wls=coeffs_table_wls,
        LS_reg_name="WLS", reg_covariates=reg_covariates, mu_x_fixed=mu_x_fixed, x_as=x_as) 
    }
    estimator_and_se = data.frame(t(estimator_and_se))
    colnames(estimator_and_se) = paste0("WLS", "_", c("estimator", "naive_se", "sandwi_se", "clstr_se"))
    CI_LS = round(estimator_and_se$WLS_estimator +
                    c(-1,1) * 1.96 * estimator_and_se$WLS_clstr_se, 3)
    CI_LS = paste(CI_LS, sep = ' ', collapse = " , ")
  }
  
  return(list(estimator_and_se=estimator_and_se, CI_LS=CI_LS, lm_obj=lin_reg_matched))
}



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

