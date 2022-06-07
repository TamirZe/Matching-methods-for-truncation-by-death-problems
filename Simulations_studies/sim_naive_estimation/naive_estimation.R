naive_sace_estimation_func = function(data_for_EM){
  # naive
  most_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y]) 
  most_naive_est_se = sqrt(  ( var(data_for_EM[A==1, Y])  / nrow(data_for_EM[A==1, ]) ) + 
                               ( var(data_for_EM[A==0, Y])  / nrow(data_for_EM[A==0, ]) )  )  
  CI_by_SE_and_Z_val_most_naive = round(most_naive_est + c(-1,1) * 1.96 * most_naive_est_se, 3)
  CI_by_SE_and_Z_val_most_naive = paste(CI_by_SE_and_Z_val_most_naive, sep = ' ', collapse = " , ")
  
  # survivors naive
  sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
  sur_naive_est_se = sqrt(  ( var(data_for_EM[A==1 & S==1, Y])  / nrow(data_for_EM[A==1 & S==1, ]) ) + 
                              ( var(data_for_EM[A==0 & S==1, Y])  / nrow(data_for_EM[A==0 & S==1, ]) )  )
  CI_by_SE_and_Z_val_sur_naive = round(sur_naive_est + c(-1,1) * 1.96 * sur_naive_est_se, 3)
  CI_by_SE_and_Z_val_sur_naive = paste(CI_by_SE_and_Z_val_sur_naive, sep = ' ', collapse = " , ")
  
  CI_naives_before_matching = data.frame(CI_by_SE_and_Z_val_most_naive, CI_by_SE_and_Z_val_sur_naive)
  colnames(CI_naives_before_matching) = c("naive_without_matching", "survivors_naive_without_matching")
  
  return(list(most_naive_est=most_naive_est, most_naive_est_se=most_naive_est_se, CI_by_SE_and_Z_val_most_naive=CI_by_SE_and_Z_val_most_naive,
              sur_naive_est=sur_naive_est, sur_naive_est_se=sur_naive_est_se, CI_by_SE_and_Z_val_sur_naive=CI_by_SE_and_Z_val_sur_naive,
              CI_naives_before_matching=CI_naives_before_matching))
  
}
