##################################################################################################################
# check PS sd over N
mat_pis_per_gamma = mat_pis_per_gamma_mis_PS = NULL
sd_lst = mean_lst = list()
n_vec = c(2000, 20000, 100000)
for (j in 1:length(n_vec)){
  print(j)
  for (i in 1:100) {
    extract_pis_lst = extract_pis_from_scenarios(nn=n_vec[j], mat_gamma=mat_gamma, xi=xi, misspec_PS=0, two_log_models_DGM=T)
    mat_pis_per_gamma = rbind(mat_pis_per_gamma, 
                              cbind(extract_pis_lst$mat_pis[1,], extract_pis_lst$mat_pis[2,]))
    
    extract_pis_lst_mis_PS = extract_pis_from_scenarios(nn=n_vec[j], mat_gamma=mat_gamma, xi=xi, misspec_PS=2, two_log_models_DGM=T)
    mat_pis_per_gamma_mis_PS = rbind(mat_pis_per_gamma_mis_PS, 
                                     cbind(extract_pis_lst_mis_PS$mat_pis[1,], extract_pis_lst_mis_PS$mat_pis[2,]))
    
  }
  
  # sd and mean over all iterations, by g, mis/correct and N
  sd_pis = data.frame(g = rownames(mat_pis_per_gamma), mat_pis_per_gamma) %>% group_by(g) %>% summarise_each(sd) 
  sd_pis_mis_PS = data.frame(g = rownames(mat_pis_per_gamma_mis_PS), mat_pis_per_gamma_mis_PS) %>% group_by(g) %>% summarise_each(sd)
  sd_lst[[j]] = cbind(crc = sd_pis, mis = sd_pis_mis_PS)
  
  mean_pis = data.frame(g = rownames(mat_pis_per_gamma), mat_pis_per_gamma) %>% group_by(g) %>% summarise_each(mean) 
  mean_pis_mis_PS = data.frame(g = rownames(mat_pis_per_gamma_mis_PS), mat_pis_per_gamma_mis_PS) %>% group_by(g) %>% summarise_each(mean)
  mean_lst[[j]] = cbind(crc = mean_pis, mis = mean_pis_mis_PS)
  
  names(sd_lst)[[j]] = names(mean_lst)[[j]] = paste0("n = ", n_vec[j])
  
  # sd_over_n
  sd_over_n = list.rbind(mean_lst)[,-1] %>% group_by(mis.g) %>% summarise_each(sd) 
}
a=sd_lst; b=mean_lst
##################################################################################################################
