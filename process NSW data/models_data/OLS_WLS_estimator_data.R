
# # simulations function
# LS_NOinter_sim_new = 
#   regression_adjusted_function(rep_bool=replace, dt_match_S1=dt_and_pairs_match_lst[[3]]$dt_match_S1,m_data=m_data,
#                                matched_pairs=dt_and_pairs_match_lst[[3]]$matched_pairs, covariates=reg_cov,
#                                interactions_bool=T, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
# # DATA function
# LS_NOinter_data =
#   regression_adjusted_function_data_old(dt_match_S1=dt_and_pairs_match_lst[[3]]$dt_match_S1, m_data=m_data,
#                                matched_pairs=dt_and_pairs_match_lst[[3]]$matched_pairs, covariates = reg_cov,
#                                interactions_bool=T, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)


regression_adjusted_function = function(rep_bool, dt_match_S1, m_data,
                        matched_pairs, covariates = reg_cov,
                        interactions_bool, LS, mu_x_fixed = FALSE, x_as){
  estimation_with_interactions = function(lin_reg_fit_matched, coeffs_table, data_matched_reg, 
                                          lin_reg_fit_wls_matched=NULL, coeffs_table_wls=NULL, # data_matched_reg_weights=NULL,
                                          LS_for_estimation, x_as){
    coeffs = coeffs_table[,1]
    #coeftest(lin_reg_fit_matched, vcov. = vcovHC, type = "HC1")
    # se is sqrt(vcov)
    names_coeffs_interactions = c("A", names(coeffs)[grep(":", names(coeffs))])
    coeffs_interactions = coeffs[names_coeffs_interactions]
    var_mu_x = apply(subset(filter(data_matched_reg, A==0), select = covariates), 2, var) / 
      nrow(filter(data_matched_reg, A==0))
    
    #TODO NEW X
    #x = expandRows(data_matched_reg, "weight") 
    x = subset(data_matched_reg, select = c("A", covariates))
    # the 1 is for the add to intercerpt for the treated (which is the coeff of A)
    x_check = subset(data_matched_reg, select = c("A", covariates))
    mu_x_check = c(A = 1, apply(subset(filter(x_check, A==0), select = -A), 2, mean))
    
    #mu_x_vec = c(A = 1, apply(subset(x, select = -A), 2, mean))
    if(mu_x_fixed == TRUE){mu_x_vec=x_as}else{
      mu_x_vec = c(A = 1, apply(subset(filter(x, A==0), select = -A), 2, mean))
    }
    beta_est = mu_x_vec %*% coeffs_interactions
    
    if(LS_for_estimation=="OLS"){
      cluster_coeffs = data.frame(coef_test(lin_reg_fit_matched, vcov = "CR1", cluster = data_matched_reg$pair))
      #TODO 02.08.2021
      vcov_clstr_interactions = vcovCL(lin_reg_fit_matched, cluster = ~ pair)[names_coeffs_interactions, names_coeffs_interactions]
      vcov_clstr_interactions2 <- 
        vcovCR(lin_reg_fit_matched, cluster = data_matched_reg$pair, type = "CR1")[names_coeffs_interactions, names_coeffs_interactions] 
      se_beta_est_clstr <- 
        sqrt( matrix(mu_x_vec, nrow=1) %*% vcov_clstr_interactions %*% matrix(mu_x_vec, ncol=1) )
      # ADD from DELTA METHOD for WLS
      se_beta_est_clstr_DM = sqrt( se_beta_est_clstr^2 + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
      estimators = c(beta_est, se_beta_est_clstr_DM)
    }
    
    if(LS_for_estimation=="WLS"){
      se = coeffs_table[,2]
      se_interactions = se[names_coeffs_interactions]
      vcov = vcov(lin_reg_fit_wls_matched) # vcov(lin_reg_fit_matched)
      vcov_interactions_naive = vcov[names_coeffs_interactions, names_coeffs_interactions]
      var_beta_est_naive = matrix(mu_x_vec, nrow=1) %*% vcov_interactions_naive %*% matrix(mu_x_vec, ncol=1) 
      se_beta_est_naive = sqrt(var_beta_est_naive)
      se_beta_est_naive_DM = sqrt( var_beta_est_naive + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
      
      # robust clustered SE estimator
      #identical(vcovHC(lin_reg_fit_matched, type = "HC"), sandwich(lin_reg_fit_matched))
      vcov_clstr_interactions = vcovCL(lin_reg_fit_matched, cluster = ~ pair + id)[names_coeffs_interactions, names_coeffs_interactions]
      se_beta_est_clstr <- 
        sqrt( matrix(mu_x_vec, nrow=1) %*% vcov_clstr_interactions %*% matrix(mu_x_vec, ncol=1) )
      #a = coeftest(lin_reg_fit_matched, vcov. = vcovCL, cluster = ~ pair + id)[,"Std. Error",drop = FALSE]
      #sqrt(diag(b)) == a
      vcov_sand_interactions = sandwich(lin_reg_fit_wls_matched)[names_coeffs_interactions, names_coeffs_interactions]
      se_beta_est_sand <- 
        sqrt( matrix(mu_x_vec, nrow=1) %*% vcov_sand_interactions %*% matrix(mu_x_vec, ncol=1) )
      #sqrt(matrix(mu_x_vec, nrow=1) %*% vcov(lin_reg_fit_matched)[names_coeffs_interactions, names_coeffs_interactions] %*% matrix(mu_x_vec, ncol=1) )
      
      # ADD from DELTA METHOD for WLS
      #var_mu_x = apply(subset(filter(x, A==0), select = -A), 2, sd) / nrow(filter(x, A==0))
      se_beta_est_clstr_DM = sqrt( se_beta_est_clstr^2 + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
      se_beta_est_sand_DM = sqrt( se_beta_est_sand^2 + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
      
      estimators = c(beta_est, se_beta_est_naive_DM, se_beta_est_sand_DM, se_beta_est_clstr_DM) 
      # estimators = data.frame(t(c(beta_est, se_beta_est_naive)))
      # colnames(estimators) = c("OLS_estimator", "OLS_se")
      # estimators = data.frame(t(c(estimators, se_beta_est_sand)))
      # colnames(estimators) = paste0(rep(c("appearance_over_unq"), each=ncol(estimators)),
      #                                    "_", c("estimator", "naive_se", "sandwi_se"))
    }
    return(estimators)
  }
  
  # run only OLS when matching is wout replacement, and WLS when matching is with replacement
  if( (rep_bool == T & LS == "OLS") | (rep_bool == F & LS == "WLS") ){return("not the regression model we need")}
  
  # reg covariares
  reg_data_matched_wout_pairs = filter(m_data, id %in% 
                                         c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) )
  #reg_data_matched2 = m_data[c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) ,]
  reg_data_matched_wout_pairs = subset(reg_data_matched_wout_pairs, select = c("id", "A", covariates, "Y"))
  reg_data_matched = merge(matched_pairs, reg_data_matched_wout_pairs, by="id") %>% arrange(pair, A)
  length(unique(dt_match_S1$id)); length(dt_match_S1$id); length(unique(dt_match_S1$A0_id)); length(dt_match_S1$A0_id); nrow(dt_match_S1)
  nrow(filter(m_data, A==0&S==1))
  
  if(interactions_bool == FALSE){
    print("NO interactions in model")
    f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + ")))
  }else{
    print("interactions in model")
    # TODO if only A in the model, there are actually no interactions
    if(is.null(covariates)){
      print("only A in the model, there are actually no interactions")
      f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + "))) 
    }else{
      f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + "), " + ",
                            paste("A *", covariates, collapse=" + "))) # rep("A*",length(covariates))
    }
  }
  
  # OLS
  if(LS=="OLS"){
    print("OLS")
    lin_reg_matched = lm(formula = f , data = reg_data_matched)
    OLS_sum = summary(lin_reg_matched)
    coeffs_table = OLS_sum$coefficients
    
    # TODO WITHOUT INTERACTIONS
    if(interactions_bool == FALSE){
      #beta_est = lin_reg_matched$coefficients[2]
      beta_est = coeffs_table["A","Estimate"]
      #TODO 02.08.2021
      se_beta_est = coeftest(lin_reg_matched, vcov. = vcovCL, cluster = ~ pair)["A","Std. Error",drop = FALSE]
      se_beta_est2 = data.frame(coef_test(lin_reg_matched, vcov = "CR1", cluster = reg_data_matched$pair))["A", "SE"]
      #se_beta_est_naive = coeffs_table["A","Std. Error"]
      estimator_and_se = c(beta_est, se_beta_est)
      
      # TODO WITH INTERACTIONS
    }else{
      estimator_and_se = estimation_with_interactions(lin_reg_fit_matched=lin_reg_matched, coeffs_table= coeffs_table,
                                                      LS_for_estimation="OLS", data_matched_reg=reg_data_matched, x_as=x_as)
    }
    estimator_and_se = data.frame(t(estimator_and_se))
    colnames(estimator_and_se) = c("OLS_estimator", "OLS_se")
    CI_LS = round(estimator_and_se$OLS_estimator + c(-1,1) * 1.96 * estimator_and_se$OLS_se, 3)
    CI_LS = paste(CI_LS, sep = ' ', collapse = " , ")
  }
  
  
  # WLS
  if(LS=="WLS"){
    print("WLS")
    # TODO adjust for replacements
    # TODO Regression adjusted matching_from_real_data on the matched set
    # TODO regression_adjusted_function(dt_match_S1, reg_covariates = X_sub_cols[-1])
    
    weights_trt = data.frame(table(dt_match_S1$id_trt))
    colnames(weights_trt) = c("id", "weight")
    weights_trt$id = as.numeric(as.character(weights_trt$id))
    weights_ctr = data.frame(id = dt_match_S1$id_ctrl , weight = rep(1, length(dt_match_S1$id_ctrl)))
    weights = rbind(weights_trt, weights_ctr)
    reg_data_matched_weights = merge(weights, reg_data_matched_wout_pairs, by= "id") #%>% arrange(pair, A)
    #weights$appearance_over_unique_id = weights$weight / nrow(weights)
    #weights$AbadieImbens785_overN0 = weights$weight / length(weights_ctr$weight) #  length(weights_ctr$weight) = sum(weights_ctr$weight)
    #weights$AbadieImbens785_over_matched_units = weights$AbadieImbens785_overN0 / 2
    
    # TODO? reg covariares instead of covariates
    lin_reg_matched = lm(formula = f , data = reg_data_matched)
    lin_reg_matched_wls = lm(formula = f , data = reg_data_matched_weights, weights = reg_data_matched_weights$weight)
    # summary of WLS
    sum_matched = summary(lin_reg_matched)
    sum_matched_wls = summary(lin_reg_matched_wls)
    coeffs_table = sum_matched$coefficients
    coeffs_table_wls = sum_matched_wls$coefficients
    
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
      
    # WITH INTERACTIONS
    }else{
      estimator_and_se = estimation_with_interactions(lin_reg_fit_matched=lin_reg_matched, coeffs_table=coeffs_table, data_matched_reg=reg_data_matched,
                                                      lin_reg_fit_wls_matched=lin_reg_matched_wls, coeffs_table_wls=coeffs_table_wls,
                                                      LS_for_estimation="WLS", x_as=x_as) 
    }
    estimator_and_se = data.frame(t(estimator_and_se))
    # colnames(estimator_and_se) = paste0(rep(c("appearance_over_unq"), each=ncol(estimator_and_se)), "_", c("estimator", "naive_se", "sandwi_se"))
    colnames(estimator_and_se) = paste0("WLS", "_", c("estimator", "naive_se", "sandwi_se", "clstr_se"))
    CI_LS = round(estimator_and_se$WLS_estimator +
                    c(-1,1) * 1.96 * estimator_and_se$WLS_clstr_se, 3)
    CI_LS = paste(CI_LS, sep = ' ', collapse = " , ")
  }
  
  return(list(estimator_and_se=estimator_and_se, CI_LS=CI_LS, coeffs_table=coeffs_table))
}


regression_adjusted_function_data_old = function(rep_bool=NULL, dt_match_S1, m_data,
                        matched_pairs, covariates = reg_cov,
                        interactions_bool, LS, mu_x_fixed = FALSE, x_as){
  estimation_with_interactions_new = function(lin_reg_matched, coeffs_table, LS_for_estimation,
                                              reg_data_matched_OLS_WLS, x_as){
    coeffs = coeffs_table[,1]
    # TODO library(lmtest)
    #coeftest(lin_reg_matched, vcov. = vcovHC, type = "HC1")
    
    # se is sqrt(vcov)
    se = coeffs_table[,2]
    vcov = vcov(lin_reg_matched)
    names_coeffs_interactions = c("A", names(coeffs)[grep(":", names(coeffs))])
    coeffs_interactions = coeffs[names_coeffs_interactions]
    se_interactions = se[names_coeffs_interactions]
    vcov_interactions = vcov[names_coeffs_interactions, names_coeffs_interactions]
    var_mu_x = apply(subset(filter(reg_data_matched_OLS_WLS, A==0), select = covariates), 2, var) / 
      nrow(filter(reg_data_matched_OLS_WLS, A==0))
    
    if(LS_for_estimation=="OLS"){
      x = subset(reg_data_matched_OLS_WLS, select = c("A", covariates))
      
      #mu_x_vec = c(A = 1, apply(subset(x, select = -A), 2, mean))
      if(mu_x_fixed == TRUE){mu_x_vec=x_as}else{
        mu_x_vec = c(A = 1, apply(subset(filter(x, A==0), select = -A), 2, mean))
      }
    }
    
    if(LS_for_estimation=="WLS"){
      # the 1 is for the add to intercerpt for the treated (which is the coeff of A)
      x_check = subset(reg_data_matched_OLS_WLS, select = c("A", covariates)) # reg_data_matched_weigts
      mu_x_check = c(A = 1, apply(subset(filter(x_check, A==0), select = -A), 2, mean))
      #TODO NEW X
      x = expandRows(reg_data_matched_OLS_WLS, "weight") # reg_data_matched_weigts
      x = subset(x, select = c("A", covariates))
      
      #mu_x_vec = c(A = 1, apply(subset(x, select = -A), 2, mean))
      if(mu_x_fixed == TRUE){
        mu_x_vec=x_as
      }else{
        mu_x_vec = c(A = 1, apply(subset(filter(x, A==0), select = -A), 2, mean))
      }
      
      #approx_mean_x == mu_x_vec
      #c(1,mean_as) - mu_x_vec
      
      # robust SE estimator
      #identical(vcovHC(lin_reg_matched, type = "HC"), sandwich(lin_reg_matched))
      v_cov_sandwich_interactions = sandwich(lin_reg_matched)[names_coeffs_interactions, names_coeffs_interactions]
      se_beta_est_robust <- 
        sqrt( matrix(mu_x_vec, nrow=1) %*% v_cov_sandwich_interactions %*% matrix(mu_x_vec, ncol=1) )
      #sqrt(matrix(mu_x_vec, nrow=1) %*% vcov(lin_reg_matched)[names_coeffs_interactions, names_coeffs_interactions] %*% matrix(mu_x_vec, ncol=1) )
      
      # ADD from DELTA METHOD for WLS
      #var_mu_x = apply(subset(filter(x, A==0), select = -A), 2, sd) / nrow(filter(x, A==0))
      se_beta_est_robust_DM = sqrt( se_beta_est_robust^2 + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
    }
    
    beta_est = mu_x_vec %*% coeffs_interactions 
    var_beta_est_naive = matrix(mu_x_vec, nrow=1)%*%vcov_interactions%*%matrix(mu_x_vec, ncol=1) 
    se_beta_est_naive = sqrt(var_beta_est_naive)
    se_beta_est_naive_DM = sqrt( var_beta_est_naive + sum( (coeffs_interactions[-1])^2 * var_mu_x ) )
    estimators = c(beta_est, se_beta_est_naive_DM)
    # estimators = data.frame(t(c(beta_est, se_beta_est_naive)))
    # colnames(estimators) = c("OLS_estimator", "OLS_se")
    
    # add robust SE estimator to WLS
    if(LS_for_estimation=="WLS"){
      estimators = c(estimators, se_beta_est_robust_DM)
      # estimators = data.frame(t(c(estimators, se_beta_est_robust)))
      # colnames(estimators) = paste0(rep(c("appearance_over_unq"), each=ncol(estimators)),
      #                                    "_", c("estimator", "naive_se", "sandwi_se"))
    }
    return(estimators)
  }
  
  # reg covariares
  reg_data_matched = filter(m_data, id %in% 
                              c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) )
  #reg_data_matched2 = m_data[c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) ,]
  reg_data_matched =  subset(reg_data_matched, select = c("id", "A", covariates, "Y"))
  
  if(interactions_bool == FALSE){
    print("NO interactions in model")
    f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + ")))
  }else{
    print("interactions in model")
    # TODO if only A in the model, there are actually no interactions
    if(is.null(covariates)){
      print("only A in the model, there are actually no interactions")
      f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + "))) 
    }else{
      f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + "), " + ",
                            paste(rep("A*",5), covariates, collapse=" + ")))
    }
  }
  
  # OLS
  if(LS=="OLS"){
    print("OLS")
    if(is.null(matched_pairs)==FALSE){reg_data_matched = merge(matched_pairs, reg_data_matched, by="id")}
    lin_reg_matched = lm(formula = f , data = reg_data_matched)
    OLS_sum = summary(lin_reg_matched)
    coeffs_table = OLS_sum$coefficients
    
    # TODO WITHOUT INTERACTIONS
    if(interactions_bool == FALSE){
      #beta_est = lin_reg_matched$coefficients[2]
      beta_est = coeffs_table["A","Estimate"]
      
      # if we have pairs, it means there are no replacements in matching (so we use OLS, and not WLS!)
      if(is.null(matched_pairs)==FALSE){
        # Imbens Spiess 2019: Clustered standard errors
        cluster_coeffs = data.frame(coef_test(lin_reg_matched, vcov = "CR1", cluster = reg_data_matched$pair))
        se_beta_est = cluster_coeffs["A", "SE"]
        # until I cahnge it, OLS reg runs inside the matching function even if we have replacements. I delete the redundent est in sim1.
        # if we dont have pairs, it means we have replacements in matching, and that we are going to delete this se est anyway. 
      }else{
        # TODO @@@ simple OLS wout inter SE
        se_beta_est = coeffs_table["A","Std. Error"]
        # TODO sandwich OLS wout inter SE
        # se_beta_est = sqrt(sandwich(lin_reg_matched)["A","A"]) 
      }
      OLS_est_se = c(beta_est, se_beta_est)
      
      # TODO WITH INTERACTIONS
    }else{
      OLS_est_se = estimation_with_interactions_new(lin_reg_matched, coeffs_table, LS_for_estimation="OLS", 
                                                reg_data_matched_OLS_WLS=reg_data_matched, x_as=x_as)
    }
    OLS_est_se = data.frame(t(OLS_est_se))
    colnames(OLS_est_se) = c("OLS_estimator", "OLS_se")
    OLS_CI_by_SE_and_Z_val = round(OLS_est_se$OLS_estimator + c(-1,1) * 1.96 * OLS_est_se$OLS_se, 3)
    CI_LS = paste(OLS_CI_by_SE_and_Z_val, sep = ' ', collapse = " , ")
    estimator_and_se_estimator = OLS_est_se
    
  }
  
  
  # WLS
  if(LS=="WLS"){
    print("WLS")
    # TODO adjust for replacements
    # TODO Regression adjusted matching_from_real_data on the matched set
    # TODO regression_adjusted_function(dt_match_S1, reg_covariates = X_sub_cols[-1])
    
    weights_trt = data.frame(table(dt_match_S1$id_trt))
    colnames(weights_trt) = c("id", "weight")
    weights_trt$id = as.numeric(as.character(weights_trt$id))
    weights_ctr = data.frame(id = dt_match_S1$id_ctrl , weight = rep(1, length(dt_match_S1$id_ctrl)))
    weights = rbind(weights_trt, weights_ctr)
    weights$appearance_over_unique_id = weights$weight / nrow(weights)
    weights$AbadieImbens785_overN0 = weights$weight / length(weights_ctr$weight) #  length(weights_ctr$weight) = sum(weights_ctr$weight)
    weights$AbadieImbens785_over_matched_units = weights$AbadieImbens785_overN0 / 2
    
    
    length(unique(dt_match_S1$id)); length(dt_match_S1$id); length(unique(dt_match_S1$A0_id)); nrow(dt_match_S1)
    
    reg_data_matched_weigts = merge(weights, reg_data_matched, by= "id")
    
    sum(weights$id %in% reg_data_matched$id) == nrow(reg_data_matched)
    sum(reg_data_matched$id %in% weights$id) == nrow(reg_data_matched)
    identical(as.numeric(reg_data_matched$id), weights$id[order(weights$id)])
    
    
    # TODO? reg covariares instead of covariates
    
    lin_reg_matched = lm(formula = f , data = reg_data_matched_weigts
                         ,weights = reg_data_matched_weigts$appearance_over_unique_id)
    sum(reg_data_matched$Y==0)
    # summary of WLS
    sum_matched = summary(lin_reg_matched)
    coeffs_table = sum_matched$coefficients
    # WLS regression treatment coefficients estimator
    # WOUT INTERACTIONS
    if(interactions_bool == FALSE){
      beta_est <- lin_reg_matched$coefficients[2]
      # WLS regression treatment SE coefficients estimator
      se_beta_est_naive <- summary(lin_reg_matched)$coefficients[2, 2]
      # robust SE estimator
      se_beta_est_robust <- sqrt(sandwich(lin_reg_matched)[2,2])
      weighted_est_se <- c(beta_est, se_beta_est_naive, se_beta_est_robust)
      # WITH INTERACTIONS
    }else{
      weighted_est_se = estimation_with_interactions_new(lin_reg_matched, coeffs_table, LS_for_estimation="WLS", 
                                                     reg_data_matched_OLS_WLS=reg_data_matched_weigts, x_as=x_as)
    }
    weighted_est_se = data.frame(t(weighted_est_se))
    # colnames(weighted_est_se) = paste0(rep(c("appearance_over_unq"), each=ncol(weighted_est_se)),
    #                  "_", c("estimator", "naive_se", "sandwi_se"))
    colnames(weighted_est_se) = paste0("WLS", "_", c("estimator", "naive_se", "sandwi_se"))
    
    WLS_CI_by_SE_and_Z_val = round(weighted_est_se$WLS_estimator +
                                     c(-1,1) * 1.96 * weighted_est_se$WLS_sandwi_se, 3)
    CI_LS = paste(WLS_CI_by_SE_and_Z_val, sep = ' ', collapse = " , ")
    
    estimator_and_se_estimator = weighted_est_se
    
  }
  
  return(list(estimator_and_se=estimator_and_se_estimator, CI_LS=CI_LS, coeffs_table=coeffs_table))
}






