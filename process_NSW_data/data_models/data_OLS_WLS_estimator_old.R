# TODO Regression adjusted matching_from_real_data on the matched set
# TODO adjust for replacements

# TODO for now I use covariates. I dont really use reg_covariates.
# TODO change it for the general case where the model is misspecified, so we dont use the same x's as the true model
regression_adjusted_function = function(dt_match_S1, m_data, matched_pairs=NULL,
                                        covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                        interactions_bool = TRUE, LS="OLS", mu_x_fixed = FALSE, x_as){
  estimation_with_interactions = function(lin_reg_matched, coeffs_table, LS_for_estimation,
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
      OLS_est_se = estimation_with_interactions(lin_reg_matched, coeffs_table, LS_for_estimation="OLS", 
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
      weighted_est_se = estimation_with_interactions(lin_reg_matched, coeffs_table, LS_for_estimation="WLS", 
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
  
  return(list(estimator_and_se_estimator=estimator_and_se_estimator, CI_LS=CI_LS, coeffs_table=coeffs_table))
}






regression_adjusted_measures_function = function(dt_match_S1, m_data,
             covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1], interactions_bool = TRUE){
  # TODO adjust for replacements
  # TODO Regression adjusted matching_from_real_data on the matched set
  # TODO regression_adjusted_function(dt_match_S1, reg_covariates = X_sub_cols[-1])
  
  weights_trt = data.frame(table(dt_match_S1$id_trt))
  colnames(weights_trt) = c("id", "weight");
  weights_trt$id = as.numeric(as.character(weights_trt$id))
  weights_ctr = data.frame(id = dt_match_S1$id_ctrl , weight = rep(1, length(dt_match_S1$id_ctrl)))
  weights = rbind(weights_trt, weights_ctr)
  weights$appearance_over_unique_id = weights$weight / nrow(weights)
  weights$AbadieImbens785_overN0 = weights$weight / length(weights_ctr$weight) #  length(weights_ctr$weight) = sum(weights_ctr$weight)
  weights$AbadieImbens785_over_matched_units = weights$AbadieImbens785_overN0 / 2
  
  
  length(unique(dt_match_S1$id)); length(dt_match_S1$id); length(unique(dt_match_S1$A0_id)); nrow(dt_match_S1)
  
  # reg covariares
  reg_data_matched = filter(m_data, id %in% 
                              c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) )
  #reg_data_matched2 = m_data[c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) ,]
  reg_data_matched =  subset(reg_data_matched, select = c("id", "A", covariates, "Y"))
  reg_data_matched_weigts = merge(weights, reg_data_matched, by= "id")
  
  #####################################################################
  # reg covariares instead of covariates
  if(interactions_bool == FALSE){
    print("NO interactions in WLS model")
    f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + ")))
  }else{
    print("interactions in WLS model")
    f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + "), " + ",
                          paste(rep("A*",5), covariates, collapse=" + ")))
  }
  lin_reg_matched = lm(formula = f , data = reg_data_matched_weigts
                       ,weights = reg_data_matched_weigts$weight)
  sum(reg_data_matched$Y==0)
  # summary of WLS
  sum_matched = summary(lin_reg_matched)
  
  # WLS estimator, model with interactions 
  x = subset(reg_data_matched_weigts, select = c("A", covariates))
  mu_x_vec = c(A = 1, apply(subset(x, select = -A), 2, mean))
  beta_est = mu_x_vec %*% coeffs_interactions 
  
  #TODO checking: mu_x_m_data_as is the real E[x|g=as]; betas_GPI are the real parameters coefficients
  #TODO check if E[x|g=as] %*% beta_interactions is unbiased to the SACE
  
  #@@@
  m_data_as = filter(m_data, g=="as")
  mu_x_m_data_as = c(A=1, apply(subset(m_data_as, select = covariates), 2, mean))
  # OF COURSE, UNBIASED 
  param_real_mu_as_real_coeff = mu_x_m_data_as %*% (betas_GPI[1,] - betas_GPI[2,]) 
  #@@@
  
  #@@@
  coeffs = sum_matched$coefficients[,1]
  names_coeffs_interactions = c("A", names(coeffs)[grep(":", names(coeffs))])
  coeffs_interactions = coeffs[names_coeffs_interactions]
  param_real_mu_as_est_coeff = mu_x_m_data_as %*% coeffs_interactions #unbiased
  #@@@
  
  weighted_est_se = c(beta_est, param_real_mu_as_real_coeff, param_real_mu_as_est_coeff)
  
  #weighted_est_se = data.frame(as.matrix(weighted_est_se, nrow=1, ncol=9))
  weighted_est_se = data.frame(t(weighted_est_se))
  colnames(weighted_est_se) = paste0(rep(c("appearance", "appearance_over_unq", "AbaImb_overN0", "AbaImb_over2N0"), each=3),
                                     "_", c("estimator", "naive_se", "sandwi_se"))
  return(weighted_est_se)
}




