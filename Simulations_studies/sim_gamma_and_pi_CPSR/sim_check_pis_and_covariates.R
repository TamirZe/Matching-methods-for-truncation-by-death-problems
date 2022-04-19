#############################################################################################
# extract strata proportion and mean covariates ####
extract_pis_from_scenarios = function(nn=250000, xi=0, misspec_PS=0){
  big_lst = list(); mat_x_as <- mat_pis <- mat_x_by_g_A <- NULL
  for( k in c(1 : nrow(mat_gamma)) ){
    gamma_ah = as.numeric(mat_gamma[k, c(1:dim_x)])
    gamma_pro =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
    gamma_ns = gamma_ns
    lst_mean_x_and_pi = simulate_data_function(gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns, xi=xi, two_log_models=TRUE,
                                               param_n=nn, misspec_PS=misspec_PS, funcform_mis_out=FALSE, 
                                               funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, only_mean_x_bool=TRUE)
    big_lst[[k]] = lst_mean_x_and_pi
    mat_x_as = rbind(mat_x_as, lst_mean_x_and_pi$x_as)
    mat_pis = rbind(mat_pis, lst_mean_x_and_pi$pi)
    mat_x_by_g_A = rbind(mat_x_by_g_A, data.frame(Scenar = k, lst_mean_x_and_pi$mean_by_A_g))
  }
  #mat_pis = data.frame(pi_as=mat_pis[,1], pi_pro=mat_pis[,3], pi_ns=mat_pis[,2])
  round(mat_pis,3)
  return(list(mat_pis=mat_pis, mat_x_by_g_A=mat_x_by_g_A, big_lst=big_lst, mat_x_as=mat_x_as))
}
mat_gamma[,c(1,2,dim_x+1,dim_x+2)]
extract_pis_lst = extract_pis_from_scenarios(nn=1000000, xi=xi, misspec_PS=2); mat_pis_per_gamma = extract_pis_lst$mat_pis
mat_pis_per_gamma


big_lst = list(); big_mat_x_by_g_A=NULL
for(i in 1:100){
  print(i)
  big_lst[[i]] = extract_pis_from_scenarios(nn=2000, xi=xi, misspec_PS=0)
  big_mat_x_by_g_A = rbind(big_mat_x_by_g_A, big_lst[[i]]$mat_x_by_g_A)
}
big_mat_x_by_g_A = subset(big_mat_x_by_g_A, select = c(Scenar,A,g, grep("X", colnames(big_mat_x_by_g_A))))
big_mat_x_by_g_A = data.table(big_mat_x_by_g_A)[, lapply(.SD, mean), by=c("Scenar", "A", "g")] %>% arrange(Scenar, g, A)
#############################################################################################

#############################################################################################
# extract mean covariates of as ####
mat_x_as = NULL
for( k in c(1 : nrow(mat_gamma)) ){
  gamma_as = as.numeric(mat_gamma[k, c(1:dim_x)]); gamma_ns =  as.numeric(mat_gamma[k, (dim_x+1): (2*dim_x)])
  lst_mean_x_and_pi = simulate_data_function(gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, xi=xi,
                                             param_n=250000, misspec_PS=0, funcform_mis_out=FALSE, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, only_mean_x_bool=TRUE)
  mat_x_as = rbind(mat_x_as, lst_mean_x_and_pi$x_as)
}
#############################################################################################