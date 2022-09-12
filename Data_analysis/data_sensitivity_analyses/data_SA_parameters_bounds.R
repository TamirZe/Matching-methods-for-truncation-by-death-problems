alpha_bounds = function(dataset_arm, reg_variables){
  Y_model = lm(as.formula(paste0("Y ~ ", paste(reg_variables, collapse = " + "))), 
               data = dataset_arm)
  Y_hat = predict(Y_model)
  Y_hat_sort = sort(Y_hat) 
  min_alpha = min(Y_hat_sort)
  max_alpha = max(Y_hat_sort)
  alpha_lower_bound = min_alpha / max_alpha
  alpha_upper_bound = max_alpha / min_alpha
  
  return(c(alpha_lower_bound=alpha_lower_bound, alpha_upper_bound=alpha_upper_bound))
}


