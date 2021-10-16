# data_with_PS from sim1.R
# SACE ~ 5.687 (1M obs); 5.690121 (3M);
# SACE~ 5.69-5.7 (200K); 5.73 (250K)
# DING_model_assisted_est_ps ~ 5.716
# DING_est ~ 5.745242

# TODO for manual run of the large dataset estimators
# data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 
# SACE = mean(data_with_PS[g_num==1 , Y1]) - mean(data_with_PS[g_num==1, Y0])
# X_sub_cols = paste0("X", c(1:(dim_x)))
# replace_vec = c(FALSE, TRUE)
#change_id_vec = c(FALSE, TRUE)


###########################################################################################
# TODO run of the functions without running in lopp on seeds and gamma values
seed= 222; gamma_pro = rep(0, dim_x); gamma_as = as.numeric(mat_gamma[2, c(1:dim_x)]);
gamma_ns =  as.numeric(mat_gamma[2, (dim_x+1): (2*dim_x)]); param_n = 10000

preprocess_for_large =  fun_simulate_data(seed, gamma_as, gamma_ns, param_n)
data_list = preprocess_for_large[[1]]; X_sub_cols = preprocess_for_large[[2]]; SACE = preprocess_for_large[[3]] 
DING_est = preprocess_for_large[[4]]; DING_model_assisted_est_ps = preprocess_for_large[[5]]
# match_on = "O11_posterior_ratio"
large_dataset_estimators = fun_large_dataset_estimators(match_on = NULL, data_list, 
                                          X_sub_cols, SACE, DING_est, DING_model_assisted_est_ps)
match_list_change_id_before_matching = large_dataset_estimators[[1]]
large_matching_estimators = large_dataset_estimators[[2]]; large_matching_WLS_est = large_dataset_estimators[[3]]
mat_mean_by_subset = large_dataset_estimators[[4]]; replace_F_T_yes_interactions = large_dataset_estimators[[5]]

save(match_list_change_id_before_matching, file = "match_list_change_id_before_matching.RData")
save(large_matching_estimators, file = "large_matching_estimators.RData")
save(large_matching_WLS_est, file = "large_matching_WLS_est.RData")
save(mat_mean_by_subset , file = "mat_mean_by_subset.RData")
save(replace_F_T_yes_interactions , file = "replace_F_T_yes_interactions.RData")
###########################################################################################


###########################################################################################
mat_gamma = matrix(c(rep(0.1, dim_x),rep(1.25, dim_x)
                     ,rep(1.5, dim_x), rep(-1, dim_x)
                     ,rep(1, dim_x),rep(0.25, dim_x))
                   ,nrow = 3, byrow = T)
colnames(mat_gamma) = paste0( "gamma", paste(rep(c(0:(dim_x-1)), times = 2)), rep(c("as", "ns"), each = dim_x) )
seed_values = seq(101, 151, 50)
lst_of_lst_large_dataset_estimators = list()
for(i in 1 :nrow(mat_gamma)){
  print(paste0("gamma index is ", i))
  lst_gamma_val = list()
  for(j in 1 : length(seed_values)){
    print(paste0("seed index is ", j))
    preprocess_for_large =  fun_simulate_data(seed_values[j], gamma_as = as.numeric(mat_gamma[i, c(1:dim_x)]), 
                              gamma_ns = as.numeric(mat_gamma[i, (dim_x+1): (2*dim_x)]), param_n)
    data_list = preprocess_for_large[[1]]; X_sub_cols = preprocess_for_large[[2]]; SACE = preprocess_for_large[[3]] 
    DING_est = preprocess_for_large[[4]]; DING_model_assisted_est_ps = preprocess_for_large[[5]]
    lst_gamma_val[[j]] = fun_large_dataset_estimators(data_list, X_sub_cols, SACE, 
                                                            DING_est, DING_model_assisted_est_ps)
    }
  lst_of_lst_large_dataset_estimators[[i]] = lst_gamma_val
}
#save(lst_of_lst_large_dataset_estimators , file = "lst_of_lst_large_dataset_estimators.RData")

# extract only large_matching_estimators
lst_large_matching_estimators = list()
for(i in 1:length(lst_of_lst_large_dataset_estimators)){ # nrow(mat_gamma)
  temp = list()
  for(j in 1:length(lst_of_lst_large_dataset_estimators[[i]])){ # length(seed_values)
    temp[[j]] = lst_of_lst_large_dataset_estimators[[i]][[j]]$large_matching_estimators
  }
  mat_temp = rbindlist(temp)
  lst_large_matching_estimators[[i]] = mat_temp
  assign(paste0("gamma_", i, "_matching_estimators"),mat_temp)
  t(mat_temp) %>% xtable(digits=c(3))
}
t(mat_gamma) %>% xtable(digits=c(3)); t(seed_values) %>% xtable(digits=c(3))
t(lst_large_matching_estimators[[1]]) %>% xtable(digits=c(3))
t(lst_large_matching_estimators[[2]]) %>% xtable(digits=c(3))
t(lst_large_matching_estimators[[3]]) %>% xtable(digits=c(3))
###########################################################################################



###########################################################################################
fun_large_dataset_estimators(match_on = NULL, data_list, 
                             X_sub_cols, SACE, DING_est, DING_model_assisted_est_ps)



fun_large_dataset_estimators = function(match_on = NULL, data_list, X_sub_cols, SACE,
    DING_est = DING_est, DING_model_assisted_est_ps, replace_vec = c(FALSE, TRUE)){
  match_list_change_id_before_matching = list()
  for(j in c(1:length(replace_vec))){
    print(j)
    match_list_change_id_before_matching[[j]] =
      lapply(1:length(data_list), function(l){
        # what is the first argument?
        my_matching_func_basic(match_on = match_on, X_sub_cols, data_list[[l]],
         weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2,
         min_PS = min_PS, min_diff_PS = min_diff_PS,
         caliper = caliper, OBS_table)
      })
  }
  
  # DING estimator
  # prior_ratio = pis[1] / (pis[1] + pis[3]) 
  # data_with_PS[, `:=` (posterior_ratio = est_p_as / (est_p_as + est_p_pro),
  #                      prior_ratio = prior_ratio, 
  #                      W_1_as = ( est_p_as / (est_p_as + est_p_pro) ) / prior_ratio)]
  # 
  # data_with_PS[, W_1_as_Y := W_1_as * Y]
  # DING_est = mean(data_with_PS[A==1 & S == 1, W_1_as_Y]) - mean(data_with_PS[A==0 & S == 1, Y])
  # ##### DING model assisted
  # DING_model_assisted_est_ps = DING_model_assisted_func(data_with_PS, x)
  
  
  large_matching_estimators = lapply(1:length(replace_vec), function(j){
    data.frame(t(unlist(list.rbind(lapply(match_list_change_id_before_matching[[j]],
                                          head, 2)))))
  })
  large_matching_estimators = list.cbind(large_matching_estimators)
  colnames(large_matching_estimators) = paste0(rep(paste0( "rep", rep(substr(replace_vec,1,1), each=6), "_",
                     paste0(rep(c("MATCH", "MATCH"), each = 3),
                            "_", c("all", "wout_O_0_0", "S1")),
                     rep(c(rep("", 3), rep("_HL", 3)), 2)), each = length(c("EST","SE"))), "_", 
         rep(c("est","SE"), times=length(large_matching_estimators)/2))
  
  large_matching_estimators = data.frame(SACE, DING_est, DING_model_assisted_est_ps, 
                                         large_matching_estimators)
  
  # wout interactions
  NOint_matching_reg_estimators = lapply(2:2, function(j){
    lapply(match_list_change_id_before_matching[[j]], "[[",
           "WLS_NOinteractions_reg_adj_estimators_and_se")})
  NOint_matching_reg_estimators = NOint_matching_reg_estimators[[1]][[3]][[1]]
  colnas = paste0(rep(c("all", "wout_O_0_0", "S1"), each=length(NOint_matching_reg_estimators[[1]])), "_",
                  rep(colnames(NOint_matching_reg_estimators[[1]]), times=3))
  NOint_matching_reg_estimators = list.cbind(NOint_matching_reg_estimators)
  colnames(NOint_matching_reg_estimators) = colnas
  # TODO check it, probably dont need it
  NOint_matching_reg_estimators = subset(NOint_matching_reg_estimators, 
           select = grep("SACE|AbaImb_overN0", colnames(NOint_matching_reg_estimators)))
  
  # with interactions
  YESint_matching_reg_estimators = lapply(2:2, function(j){
    lapply(match_list_change_id_before_matching[[j]], "[[",
           "YESinteractions_reg_adj_estimators_and_se")})
  YESint_matching_reg_estimators = YESint_matching_reg_estimators[[1]][[3]][[1]]
  colnas = paste0(rep(c("all", "wout_O_0_0", "S1"), each=length(YESint_matching_reg_estimators[[1]])), "_",
                  rep(colnames(YESint_matching_reg_estimators[[1]]), times=3))
  YESint_matching_reg_estimators = list.cbind(YESint_matching_reg_estimators)
  colnames(YESint_matching_reg_estimators) = colnas
  # TODO check it, probably dont need it
  YESint_matching_reg_estimators = subset(YESint_matching_reg_estimators, 
             select = grep("SACE|AbaImb_overN0", colnames(YESint_matching_reg_estimators)))
  
  large_matching_WLS_est = rbind("NOint_matching_reg_estimators" = NOint_matching_reg_estimators, 
                                 "YESint_matching_reg_estimators" = YESint_matching_reg_estimators)
  #large_matching_WLS_est = large_matching_WLS_est[,c(1:3,13:15,25:27)]
  
  
  # means_by_subset
  mat_mean_by_subset = lapply(1:2, function(j){
    lapply(match_list_change_id_before_matching[[j]], "[[", "means_by_subset")})
  mat_mean_by_subset = rbind("rep_F" = list.rbind(mat_mean_by_subset[[1]]), 
                             "rep_T" = list.rbind(mat_mean_by_subset[[2]]))
  rownames(mat_mean_by_subset) = paste0(rep(c("reF_", "reT_"), each=18), 
                rep(c("all", "wout_0_0", "S1"), each=6), "_", rownames(mat_mean_by_subset))
  
  # WLS regression alos wout replacements in the matching
  replace_F_T = lapply(1:2, function(j){
    lapply(match_list_change_id_before_matching[[j]], "[[",
           "YESinteractions_reg_adj_estimators_and_se")})
  replace_F_T = rbind("rep_F" = list.cbind(replace_F_T[[1]]), "rep_T" = list.cbind(replace_F_T[[2]]))
  colnames(replace_F_T) = colnas
  replace_F_T_yes_interactions = replace_F_T[,c(1:3,13:15,25:27)]
  
  return(list(match_list_change_id_before_matching = match_list_change_id_before_matching,
    large_matching_estimators = large_matching_estimators, large_matching_WLS_est = large_matching_WLS_est,
    mat_mean_by_subset = mat_mean_by_subset, replace_F_T_yes_interactions = replace_F_T_yes_interactions))
}


fun_simulate_data = function(seed, gamma_as, gamma_ns, param_n,
                             iterations = 12, epsilon = 0.001){
  start_time1 <- Sys.time()
  set.seed(seed)
  list_data_for_EM_and_X = simulate_data_function(gamma_as, gamma_ns, param_n)
  data_for_EM = list_data_for_EM_and_X[[1]]; x = list_data_for_EM_and_X[[2]]
  OBS_table = list_data_for_EM_and_X[[3]]; pis = list_data_for_EM_and_X[[4]]
  vec_OBS_table = t(c(OBS_table))
  colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
  ###########################################################
  # real parameter
  SACE = mean(data_for_EM[g_num==1 , Y1]) - mean(data_for_EM[g_num==1, Y0])
  SACE_conditional = mean(data_for_EM[A==1 & g_num==1 , Y]) - mean(data_for_EM[A==0 & g_num==1, Y])
  
  # naive estimators
  
  most_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y])
  # first version
  #omit_Y0_naive_est = mean(data_for_EM[A==1 & Y != 0, Y]) - mean(data_for_EM[A==0 & Y != 0, Y])
  # second version
  omit_Y_less0_naive_est = mean(data_for_EM[A==1 & Y > 0, Y]) - mean(data_for_EM[A==0 & Y > 0, Y])
  
  # sur_naive_est naive is (almost) like restricted analysis, since we ignore subjects with Y = 0
  # only almost in this case, since some of the subjects hs negative outcome
  sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
  
  ###########################################################
  start_time2 <- Sys.time()
  EM_list = function_my_EM(data_for_EM, iterations, epsilon)
  end_time2 <- Sys.time()
  print(paste0("function_my_EM lasts ", difftime(end_time2, start_time2)))
  dat_EM = EM_list[[1]]
  # after running EM, merge both data tables
  data_with_PS = data.table(merge(x = data_for_EM,
                                  y = subset(dat_EM, select = c(id, p_as, p_ns, p_pro, max_strata_per_subj)),
                                  by = "id", all.x = TRUE))
  
  # EM coeffs
  # coeff_as, coeff_ns, coeff_pro
  coeff_as = EM_list[[2]] ; coeff_ns = EM_list[[3]]
  #list_dat_EM[[i]] = dat_EM
  #list_coeff_as[[i]] = coeff_as; list_coeff_ns[[i]] = coeff_ns
  
  #########################################################################################
  # calculating PS from the M step in the EM, not from the E step
  # the E step takes into account also the cells, and we don't want to do such thing here yet
  PS_est = cbind(exp(x%*%coeff_as), exp(x%*%coeff_ns), exp(x%*%gamma_pro))
  PS_est = PS_est / apply(PS_est, 1, sum)
  colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
  data_with_PS = data.table(data.frame(data_with_PS, PS_est))
  
  # DING estimator
  O11_prior_ratio = pis[1] / (pis[1] + pis[3]) 
  data_with_PS[, `:=` (O11_posterior_ratio = EMest_p_as / (EMest_p_as + EMest_p_pro),
                       O11_prior_ratio = O11_prior_ratio, 
                       W_1_as = ( EMest_p_as / (EMest_p_as + EMest_p_pro) ) / O11_prior_ratio)]
  
  data_with_PS[, W_1_as_Y := W_1_as * Y]
  DING_est = mean(data_with_PS[A==1 & S == 1, W_1_as_Y]) - mean(data_with_PS[A==0 & S == 1, Y])
  
  ##### DING model assisted, 3 options, only 1 for now: 
  DING_model_assisted_est_ps = DING_model_assisted_func(data_with_PS, x)
  
  
  #########################################################################################
  # MATCHING and estimation 
  
  # run for all options (3 options)
  # m_data = data_with_PS; data_with_PS[OBS != "O(0,0)"]; m_data = data_with_PS[S==1]
  data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 
  X_sub_cols = paste0("X", c(1:(dim_x)))
  return(list(data_list = data_list, X_sub_cols = X_sub_cols,  SACE = SACE, 
              DING_est = DING_est, DING_model_assisted_est_ps = DING_model_assisted_est_ps))
}
###########################################################################################




