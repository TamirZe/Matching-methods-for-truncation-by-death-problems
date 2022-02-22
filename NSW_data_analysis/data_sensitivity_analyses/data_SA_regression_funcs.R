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
