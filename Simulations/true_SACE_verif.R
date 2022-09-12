##################################################################################################################
# check PS sd over iterations and over N ####
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

##################################################################################################################


###############################################################################################
# one_large_simulation VS multiple small simulations ####

# true SACE parameter from one large simulation
one_large_simulation = simulate_data_func(
  gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns,
  xi=xi, two_log_models_DGM=two_log_models_DGM, param_n=200000, 
  misspec_PS=2, misspec_outcome=0, transform_x=transform_x,
  funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
  betas_GPI=betas_GPI, var_GPI=var_GPI, rho_GPI_PO=rho_GPI_PO)
true_SACE = one_large_simulation$true_SACE

# SACE parameter from mean of multiple simulations of sample size param_n

SACE_vec = Y1_vec = Y0_vec = 
  Y1_as_vec = Y0_as_vec = Y1_pro_vec = Y0_pro_vec = Y1_ns_vec = Y0_ns_vec = vector(length = 100)
prob_as_vec = prob_pro_vec = prob_ns_vec = prob_as_vec2 = prob_pro_vec2 = prob_ns_vec2 = c() 
pi_mat = matrix(nrow = length(SACE_vec), ncol = 3)
mean_x_mat = NULL
big_lst = mean_by_g_lst = list()
for (i in 1:length(SACE_vec)){
  print(i)
  one_small_simulation = simulate_data_func(
    gamma_ah=gamma_ah, gamma_pro=gamma_pro, gamma_ns=gamma_ns,
    xi=xi, two_log_models_DGM=two_log_models_DGM, param_n=param_n, 
    misspec_PS=2, misspec_outcome=0, transform_x=transform_x,
    funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
    betas_GPI=betas_GPI, var_GPI=var_GPI, rho_GPI_PO=rho_GPI_PO)
  
  big_lst[[i]] = one_small_simulation
  pi_mat[i,] = one_small_simulation$pis
  mean_by_g_lst[[i]] = apply(one_small_simulation$mean_by_g, 2, as.numeric)
  crnt_mean_x = apply(data.frame(one_small_simulation$dt[,-c("id","g","Y1","Y0","Y","OBS")], one_small_simulation$x_misspec), 2,mean)
  mean_x_mat = rbind(mean_x_mat, crnt_mean_x); colnames(mean_x_mat) = names(crnt_mean_x)
  SACE_vec[i] = one_small_simulation$true_SACE
  Y1_vec[i] = mean(one_small_simulation$dt[, Y1])
  Y0_vec[i] = mean(one_small_simulation$dt[, Y0])
  
  #Y1_as_vec[i] = mean(one_small_simulation$dt[g == "as", Y1])
  #Y0_as_vec[i] = mean(one_small_simulation$dt[g == "as", Y0])
  # Y1_pro_vec[i] = mean(one_small_simulation$dt[g == "pro", Y1])
  # Y0_pro_vec[i] = mean(one_small_simulation$dt[g == "pro", Y0])
  # Y1_ns_vec[i] = mean(one_small_simulation$dt[g == "ns", Y1])
  # Y0_ns_vec[i] = mean(one_small_simulation$dt[g == "ns", Y0])
  Y1_as_vec = c(Y1_as_vec, one_small_simulation$dt[g == "as", Y1])
  Y0_as_vec = c(Y0_as_vec, one_small_simulation$dt[g == "as", Y0])
  Y1_pro_vec = c(Y1_pro_vec, one_small_simulation$dt[g == "pro", Y1])
  Y0_pro_vec = c(Y0_pro_vec, one_small_simulation$dt[g == "pro", Y0])
  Y1_ns_vec = c(Y1_ns_vec, one_small_simulation$dt[g == "ns", Y1])
  Y0_ns_vec = c(Y0_ns_vec, one_small_simulation$dt[g == "ns", Y0])
  
  prob_as_vec = c(prob_as_vec, one_small_simulation$dt[, prob_as])
  prob_pro_vec = c(prob_pro_vec, one_small_simulation$dt[, prob_pro])
  prob_ns_vec = c(prob_ns_vec, one_small_simulation$dt[, prob_ns])
  prob_as_vec2 = mean(one_small_simulation$dt[, prob_as])
  prob_pro_vec2 = mean(one_small_simulation$dt[, prob_pro])
  prob_ns_vec2 = mean(one_small_simulation$dt[, prob_ns])
  
  
}

# true SACE
true_SACE
mean(SACE_vec); mean(Y1_as_vec) - mean(Y0_as_vec)

# pis
one_large_simulation$pis
apply(pi_mat, 2, mean)

# conditional probabilities
c(mean(one_large_simulation$dt$prob_as), mean(one_large_simulation$dt$prob_ns), mean(one_large_simulation$dt$prob_pro))
c(mean(prob_as_vec), mean(prob_ns_vec), mean(prob_pro_vec))
c(mean(prob_as_vec2), mean(prob_ns_vec2), mean(prob_pro_vec2))

# mean of X by G 
# X_ps includes all transformations x_misspec # x (outcome) includes c("X_sqr", "X_log")
#TODO MEAN OF ORIGINAL X'S ARE DIFFERENT IN THE LARGE AND SMALL SIM @@@@
one_large_simulation$mean_by_g
mean_by_g_sum = data.frame(apply(simplify2array(mean_by_g_lst), 2, rowMeans, na.rm = TRUE))
mean_by_g_sum[,1] = c(mapvalues(as.numeric(mean_by_g_sum[,1]), from = c(0:3), to = c("har", "as", "ns", "pro")))
all_dts = data.frame(do.call("rbind", lapply(big_lst, "[[", "dt")), do.call("rbind", lapply(big_lst, "[[", "x_misspec")))[,colnames(mean_by_g_sum)]
mean_by_g_all_dts = data.table(all_dts)[, lapply(.SD, mean), by="g"] %>% arrange(g)

# mean of X by G - marginal 
apply(data.frame(one_large_simulation$dt, one_large_simulation$x_misspec)[,colnames(mean_by_g_sum)][,-1], 2, mean)
apply(mean_x_mat, 2, mean)


#PO's - marginal and given G's

#Y1
mean(one_large_simulation$dt[,Y1])
mean(Y1_vec) 
mean(one_large_simulation$dt[g == "as", Y1])
mean(Y1_as_vec)
mean(one_large_simulation$dt[g == "pro", Y1])
mean(Y1_pro_vec)
mean(one_large_simulation$dt[g == "ns", Y1])
mean(Y1_ns_vec)
#Y0
mean(one_large_simulation$dt[, Y0])
mean(Y0_vec)
mean(one_large_simulation$dt[g == "as", Y0])
mean(Y0_as_vec)
mean(one_large_simulation$dt[g == "pro", Y0])
mean(Y0_pro_vec)
mean(one_large_simulation$dt[g == "ns", Y0])
mean(Y0_ns_vec)

# prob_as hist
p_as_dat = data.frame(c(large = one_large_simulation$dt$prob_as, small = prob_as_vec))
p_as_dat = data.frame(sim = substr(rownames(p_as_dat),1,5), p_as = p_as_dat[,1])
ggplot(p_as_dat, aes(p_as, fill = as.factor(sim))) + geom_density(alpha = 0.5) + labs(fill = "sim") + 
  ggtitle(paste0("prob as, mean: large-", round(mean(one_large_simulation$dt$prob_as),3), ". small-", round(mean(prob_as_vec),3)))
#ggplot(p_as_dat, aes(p_as, fill = as.factor(sim))) + geom_histogram(alpha = 0.5, aes(y = (..count..)/sum(..count..)), position = 'identity', bins=18) + labs(fill = "sim") # aes(y = ..density..) # (..count..)/sum(..count..)

# prob_pro hist
hist(one_large_simulation$dt$prob_pro, xlim=c(0,1.1), col="blue", cex.main=0.9, main = 
       paste0("prob pro, mean: large (blue: one large mean)-", round(mean(one_large_simulation$dt$prob_pro),3), ". small-", round(mean(prob_pro_vec),3)))
hist(prob_pro_vec, add=T, col=rgb(0, 1, 0, 0.5))

# prob_ns hist
hist(one_large_simulation$dt$prob_ns, xlim=c(0,1.1), col="blue", cex.main=0.9, main = 
       paste0("prob ns, mean: large (blue: one large mean)-", round(mean(one_large_simulation$dt$prob_ns),3), ". small-", round(mean(prob_ns_vec),3)))
hist(prob_ns_vec, add=T, col=rgb(0, 1, 0, 0.5))
###############################################################################################