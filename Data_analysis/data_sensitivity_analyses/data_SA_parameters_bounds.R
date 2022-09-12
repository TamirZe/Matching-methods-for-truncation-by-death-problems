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


reg_variables = c("age", "education", "re75", "black", "hispanic", "married", "nodegree", "emp75")
reg_variables = reg_after_match[-1]

alpha1_bounds = alpha_bounds(dataset_arm = data_with_PS[data_with_PS$S == 1,] %>% filter(A==1), 
                             reg_variables = reg_after_match[-1])
alpha0_bounds = alpha_bounds(dataset_arm = data_with_PS[data_with_PS$S == 1,] %>% filter(A==0), 
                             reg_variables = reg_after_match[-1])
