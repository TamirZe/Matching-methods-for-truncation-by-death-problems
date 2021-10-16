# regression_adjusted_function(dt_match_S1, m_data, 
#                              covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
#                              interactions_bool = FALSE)

regression_adjusted_function_multiple_weights = function(dt_match_S1, m_data,
                                        covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                        interactions_bool = TRUE){
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
  
  # reg covariares
  reg_data_matched = filter(m_data, id %in% 
                              c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) )
  #reg_data_matched2 = m_data[c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) ,]
  reg_data_matched =  subset(reg_data_matched, select = c("id", "A", covariates, "Y"))
  reg_data_matched_weigts = merge(weights, reg_data_matched, by= "id")
  
  sum(weights$id %in% reg_data_matched$id) == nrow(reg_data_matched)
  sum(reg_data_matched$id %in% weights$id) == nrow(reg_data_matched)
  identical(as.numeric(reg_data_matched$id), weights$id[order(weights$id)])
  #####################################################################
  #TODO? reg covariares instead of covariates
  wts = c("weight", "appearance_over_unique_id", "AbadieImbens785_overN0", "AbadieImbens785_over_matched_units")
  weighted_est_se = c()
  for(i in 1:4){
    if(interactions_bool == FALSE){
      print("NO interactions in WLS model")
      f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + ")))
    }else{
      print("interactions in WLS model")
      f = as.formula(paste0("Y ~ ", paste(c("A", covariates), collapse = " + "), " + ",
                            paste(rep("A*",5), covariates, collapse=" + ")))
    }
    lin_reg_matched = lm(formula = f , data = reg_data_matched_weigts
                         ,weights = reg_data_matched_weigts[,wts[i]])
    sum(reg_data_matched$Y==0)
    # summary of WLS
    sum_matched = summary(lin_reg_matched)
    # WLS regression treatment coefficients estimator
    # WOUT INTERACTIONS
    if(interactions_bool == FALSE){
      beta_est <- lin_reg_matched$coefficients[2]
      # WLS regression treatment SE coefficients estimator
      se_beta_est_naive <- summary(lin_reg_matched)$coefficients[2, 2]
      # robust SE estimator
      se_beta_est_robust <- sqrt(sandwich(lin_reg_matched)[2,2])
    # WITH INTERACTIONS
    }else{
      coeffs_table = sum_matched$coefficients
      coeffs = coeffs_table[,1]
      # se is sqrt(vcov)
      se = coeffs_table[,2]
      vcov = vcov(lin_reg_matched)
      names_coeffs_interactions = c("A", names(coeffs)[grep(":", names(coeffs))])
      coeffs_interactions = coeffs[names_coeffs_interactions]
      se_interactions = se[names_coeffs_interactions]
      vcov_interactions = vcov[names_coeffs_interactions, names_coeffs_interactions]
      # the 1 is for the add to intercerpt for the treated (which is the coeff of A)
      x1 = subset(reg_data_matched_weigts, select = c("A", covariates))
      mu_x_vec1 = c(A = 1, apply(subset(x1, select = -A), 2, mean))
      # TODO NEW X
      x = expandRows(reg_data_matched_weigts, "weight")
      x = subset(x, select = c("A", covariates))
      mu_x_vec = c(A = 1, apply(subset(x, select = -A), 2, mean))
      #approx_mean_x == mu_x_vec
      #c(1,mean_as) - mu_x_vec
      beta_est = mu_x_vec %*% coeffs_interactions
      var_beta_est_naive = matrix(mu_x_vec, nrow=1)%*%vcov_interactions%*%matrix(mu_x_vec, ncol=1)
      se_beta_est_naive = sqrt(var_beta_est_naive)
      # robust SE estimator
      v_cov_sandwich_interactions = sandwich(lin_reg_matched)[names_coeffs_interactions, names_coeffs_interactions]
      se_beta_est_robust <-
        sqrt(matrix(mu_x_vec, nrow=1) %*% v_cov_sandwich_interactions %*% matrix(mu_x_vec, ncol=1) )
    }
    weighted_est_se = c(weighted_est_se, beta_est, se_beta_est_naive, se_beta_est_robust)
  }
  
  cols = c("estimator", "naive_se", "sandwi_se")
  weighted_est_se = data.frame(t(weighted_est_se))
  colnames(weighted_est_se) = paste0(rep(wts, each=length(cols)),
                                     "_", c("estimator", "naive_se", "sandwi_se"))

  return(weighted_est_se)
}


