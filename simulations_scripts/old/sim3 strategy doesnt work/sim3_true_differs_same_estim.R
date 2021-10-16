gamma_pro = rep(0, dim_x)
gamma_as = as.numeric(mat_gamma[1, c(1:dim_x)])
gamma_ns =  as.numeric(mat_gamma[1, (dim_x+1): (2*dim_x)])


# simulate data per each misspec_PS (0,1,2) and use simulate_data_function from sim2 
# then, run EM and estimation process on all 3 possibilities together (its the same obs data actually)
# TODO simulate_datasets_by_misspec_PS:
# function that generates 3 lists per each misspec_PS_index type.
# each list includes param_n_sim lists (for instabnce, 750 simulations)
# every list of this 750 lists is also a list, with the output of simulate_data_function from sim2
simulate_datasets_by_misspec_PS = function(seed_vec, gamma_as, gamma_ns, gamma_pro, 
                                           param_n_sim, first_misspec=0, last_misspec=2,
                                           U_factor, funcform_factor_sqr, funcform_factor_log){
multi_list_data_for_EM_and_X = list()
  for(misspec_PS_index in first_misspec:last_misspec){
    #set.seed(102)
    print(paste0("misspec_PS_index is ", misspec_PS_index))
    temp_misspec_PS_lts = list()
    for (i in 1:param_n_sim){
      seed_vec[i]
      print(paste0("param_n_sim is ", param_n_sim))
      # simulate_data_function from sim2
      temp_misspec_PS_lts[[i]] = simulate_data_function(seed_num=seed_vec[i],
        gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro, param_n=param_n,
        misspec_PS = misspec_PS_index, misspec_outcome_funcform,
        U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
        epsilon_1_GPI = epsilon_1_GPI)
    }
    multi_list_data_for_EM_and_X[[(misspec_PS_index+1)]] = temp_misspec_PS_lts 
  }
names(multi_list_data_for_EM_and_X) = paste0("misspec_PS", c(first_misspec:last_misspec))
return(multi_list_data_for_EM_and_X)
}

# seed_vec is a vector in length param_n_sim
set.seed(101)
seed_vec = sample(c(1:10000), param_n_sim, replace = FALSE)
multi_list_data_for_EM_and_X = simulate_datasets_by_misspec_PS(seed_vec, gamma_as, gamma_ns, gamma_pro, 
                 param_n_sim, first_misspec=0, last_misspec=2, U_factor, funcform_factor_sqr, funcform_factor_log)
all_nsim_dt_misspec0_lst = list.rbind(lapply(multi_list_data_for_EM_and_X[[1]], "[[", "dt"))
x0 =  subset(all_nsim_dt_misspec0_lst, select = grep("X", colnames(all_nsim_dt_misspec0_lst), ignore.case = FALSE))
all_nsim_dt_misspec1_lst = list.rbind(lapply(multi_list_data_for_EM_and_X[[2]], "[[", "dt"))
x1 =  subset(all_nsim_dt_misspec1_lst, select = grep("X", colnames(all_nsim_dt_misspec1_lst), ignore.case = FALSE))
all_nsim_dt_misspec2_lst = list.rbind(lapply(multi_list_data_for_EM_and_X[[3]], "[[", "dt"))
x2 =  subset(all_nsim_dt_misspec2_lst, select = grep("X", colnames(all_nsim_dt_misspec2_lst), ignore.case = FALSE))
# TODO IMPO@@@ check the x's are equal for all misspec types
View(x0); View(x1); View(x2)
identical(x0,x1); identical(x1,x2)


# implement EM: the same one for all 3 misspec types,
# since they observed the same X'S, if set.seed(102) worked well
implement_EM_onetime_for_all_misspec_type = function(param_n_sim, iterations, epsilon_EM, gamma_pro,
                                                     all_nsim_dt_misspec0_lst){
  print(paste0("length of all_nsim_dt_misspec0_lst is ", length(all_nsim_dt_misspec0_lst)))
  print(paste0("param_n_sim is ", param_n_sim))
  print("they 2 should be the same")
  dat_EM_list = list(); coeff_as_list = list(); coeff_ns_list = list()
  for (i in 1:param_n_sim) {
    data_for_EM = all_nsim_dt_misspec0_lst[[i]]
    start_time2 <- Sys.time()
    EM_list = function_my_EM(data_for_EM, iterations, epsilon_EM)
    end_time2 <- Sys.time()
    print(paste0("function_my_EM lasts ", difftime(end_time2, start_time2)))
    dat_EM = EM_list[[1]]
    
    # EM coeffs
    # coeff_as, coeff_ns, coeff_pro
    coeff_as = unlist(EM_list[[2]]) ; coeff_ns = unlist(EM_list[[3]])
    coeff_as_list[[i]] = coeff_as; coeff_ns_list[[i]] = coeff_ns
    
    # calculating PS from the M step in the EM, not from the E step
    # the E step takes into account also the cells, and we don't want to do such thing here yet
    x = as.matrix(subset(dat_EM, select = grep("X", colnames(dat_EM), ignore.case = FALSE)))
    PS_est = cbind(exp(x%*%coeff_as), exp(x%*%coeff_ns), exp(x%*%gamma_pro))
    PS_est = PS_est / apply(PS_est, 1, sum)
    colnames(PS_est) = c("EMest_p_as", "EMest_p_ns", "EMest_p_pro")
    data_with_PS = data.table(data.frame(id = dat_EM$id, PS_est))
    dat_EM_list[[i]] = data_with_PS
  }
  return(list(dat_EM_list=dat_EM_list, coeff_as_list=coeff_as_list, coeff_ns_list=coeff_ns_list))
}

# seed_vec is a vector in length param_n_sim
simulate_multi_misspec_by_calculate_EM_forall_together = 
  function(seed_vec, index_set_of_params, gamma_as, gamma_ns, gamma_pro, 
           param_n_sim,first_misspec=0,last_misspec=2,
           U_factor, funcform_factor_sqr, funcform_factor_log, match_on){
    multi_list_data_for_EM_and_X = 
      simulate_datasets_by_misspec_PS(seed_vec=seed_vec, 
                                      gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro,
                                      U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log,
                                      param_n_sim=param_n_sim, first_misspec=first_misspec, last_misspec=last_misspec)
    all_nsim_dt_misspec0_lst = lapply(multi_list_data_for_EM_and_X[[1]], "[[", "dt")
    EM_lst = implement_EM_onetime_for_all_misspec_type(param_n_sim, iterations, epsilon_EM, gamma_pro,
                                                       all_nsim_dt_misspec0_lst)
    # SAVE Pper mat gamma
    #assign(paste0("multi_list_data_for_EM_and_X_gamma_", index_set_of_params, ".RData"), multi_list_data_for_EM_and_X)
    #save(multi_list_data_for_EM_and_X, file = paste0("multi_list_data_for_EM_and_X_gamma_", index_set_of_params, ".RData"))
    #save(EM_lst, file = paste0("EM_lst_gamma_", index_set_of_params, ".RData"))
    #saveRDS(multi_list_data_for_EM_and_X, file = paste0("multi_list_data_for_EM_and_X_gamma_", index_set_of_params, ".rds"))
    
    all_misspec_big_sim_lst = list()
    for(misspec_PS_index in 1:length(multi_list_data_for_EM_and_X)){
      all_misspec_big_sim_lst[[misspec_PS_index]] = 
        simulate_data_run_EM_and_match_multi_misspec_PS(
          multi_list_data_for_EM_and_X[[misspec_PS_index]], EM_lst,
          return_EM_PS = FALSE, index_set_of_params=index_set_of_params, gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro,
          misspec_PS=misspec_PS_index, misspec_outcome_funcform=FALSE
          #,U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, 
          ,match_and_reg_watch_true_X=FALSE, param_n, param_n_sim, iterations = iterations, epsilon_EM = epsilon_EM,
          caliper, epsilon_1_GPI = epsilon_1_GPI, match_on = match_on, mu_x_fixed=FALSE, x_as)
    }
    return(all_misspec_big_sim_lst)
  }


# set.seed(101)
# seed_vec = sample(c(1:10000), param_n_sim, replace = FALSE)
# A = simulate_multi_misspec_by_calculate_EM_forall_together(seed_vec=seed_vec, index_set_of_params=index_set_of_params,
#    gamma_as=gamma_as, gamma_ns=gamma_ns, gamma_pro=gamma_pro,
#    param_n_sim,first_misspec=0,last_misspec=2,
#    U_factor=U_factor, funcform_factor_sqr=funcform_factor_sqr, funcform_factor_log=funcform_factor_log, match_on=match_on)



  
# TODO misspec_PS: 0 <- NO, 1:add 2 new U's to PS model and to outcome model
# 2: add transformations to PS model, and remain original X's in ourcome model.
# when misspec_PS = 1 & U_factor=0, there is no misspecification. U_factor=1 means the coeffs of U are the same as X, in the PS model

list_data_for_EM_and_X_all_sim_for_specific_misspec = multi_list_data_for_EM_and_X[[1]]


simulate_data_run_EM_and_match_multi_misspec_PS = function(
        list_data_for_EM_and_X_all_sim_for_specific_misspec, EM_lst,
        return_EM_PS = FALSE, index_set_of_params=-101, gamma_as, gamma_ns, gamma_pro,
        misspec_PS, misspec_outcome_funcform=FALSE
        #, U_factor=0, funcform_factor_sqr=0, funcform_factor_log=0
        , match_and_reg_watch_true_X=FALSE, param_n, param_n_sim, iterations = 12, epsilon_EM = 0.001,
        caliper, epsilon_1_GPI = 1, match_on = NULL, mu_x_fixed=FALSE, x_as){
  
  list_dat_EM <- list_coeff_as <- list_coeff_ns <- list()
  # run over param_n_sim different samples, each with param_n observations
  WLS_NOint_mat_reg_estimators <- WLS_YESint_mat_reg_estimators <-
    OLS_NOint_mat_reg_estimators <- OLS_YESint_mat_reg_estimators <-
    mat_param_estimators <- mat_excluded_included_matching <- 
    CI_mat <- mat_diff_distance_aspr_asas <- mat_std_mean_diff <-  NULL
  list_repeated_as_and_pro <- list_diff_distance_aspr_asas <- 
    list_matched_units <- list_std_mean_diff <- list_means_by_subset <- list()
  
  # run over param_n_sim different samples, each with param_n observations
  for (i in 1:param_n_sim){
    #set.seed(100 + i)
    print(paste0("this is index_set_of_params ", index_set_of_params))
    print(paste0("this is n_sim ", i, " in simulate_data_run_EM_and_match"))
    start_time1 <- Sys.time()
    #set.seed(258)
    list_data_for_EM_and_X = list_data_for_EM_and_X_all_sim_for_specific_misspec[[i]]
    data_for_EM = list_data_for_EM_and_X$dt
    x = list_data_for_EM_and_X$x_obs; x_PS = data.frame(list_data_for_EM_and_X$x_PS)
    x_outcome = data.frame(list_data_for_EM_and_X$x_outcome)
    OBS_table = list_data_for_EM_and_X$OBS_table; pis = list_data_for_EM_and_X$pi
    vec_OBS_table = t(c(OBS_table))
    colnames(vec_OBS_table) = c("A0_S0", "A1_S0", "A0_S1", "A1_S1")
    
    # EM elements from pre calcultaed EM_lst
    dat_EM = EM_lst$dat_EM_list[[i]]; coeff_as = EM_lst$coeff_as_list[[i]]; coeff_ns = EM_lst$coeff_ns_list[[i]]
    list_coeff_as[[i]] = coeff_as; list_coeff_ns[[i]] = coeff_ns
    ###########################################################
    # real parameter
    SACE = mean(data_for_EM[g_num==1 , Y1]) - mean(data_for_EM[g_num==1, Y0])
    SACE_conditional = mean(data_for_EM[A==1 & g_num==1 , Y]) - mean(data_for_EM[A==0 & g_num==1, Y])
    
    # naive estimators
    
    most_naive_est = mean(data_for_EM[A==1, Y]) - mean(data_for_EM[A==0, Y])
    most_naive_est_se = sqrt(( var(data_for_EM[A==1, Y]) + var(data_for_EM[A==0, Y]) ) / nrow(data_for_EM))
    CI_by_SE_and_Z_val_most_naive = round(most_naive_est + c(-1,1) * 1.96 * most_naive_est_se, 3)
    CI_by_SE_and_Z_val_most_naive = paste(CI_by_SE_and_Z_val_most_naive, sep = ' ', collapse = " , ")
    # first version
    #omit_Y0_naive_est = mean(data_for_EM[A==1 & Y != 0, Y]) - mean(data_for_EM[A==0 & Y != 0, Y])
    # second version
    omit_Y_less0_naive_est = mean(data_for_EM[A==1 & Y > 0, Y]) - mean(data_for_EM[A==0 & Y > 0, Y])
    
    # sur_naive_est naive is (almost) like restricted analysis, since we ignore subjects with Y = 0
    # only almost in this case, since some of the subjects has negative outcome
    sur_naive_est = mean(data_for_EM[A==1 & S == 1, Y]) - mean(data_for_EM[A==0 & S == 1, Y])
    sur_naive_est_se = sqrt(( var(data_for_EM[A==1 & S==1, Y]) + var(data_for_EM[A==0 & S==1, Y]) ) / nrow(data_for_EM))
    CI_by_SE_and_Z_val_sur_naive = round(sur_naive_est + c(-1,1) * 1.96 * sur_naive_est_se, 3)
    CI_by_SE_and_Z_val_sur_naive = paste(CI_by_SE_and_Z_val_sur_naive, sep = ' ', collapse = " , ")
    
    CI_naives_before_matching = data.frame(CI_by_SE_and_Z_val_most_naive, CI_by_SE_and_Z_val_sur_naive)
    colnames(CI_naives_before_matching) = c("naive_without_matching", "survivors_naive_without_matching")
    ###########################################################
    # after running EM, merge both data tables
    data_with_PS = data.table(merge(data_for_EM, dat_EM, by = "id", all.x = TRUE))
    
    if(return_EM_PS == TRUE){
      PS_true_EM_compr = subset( data_with_PS, select = grep("^id$|^g$|prob.|EM",colnames(data_with_PS)) )
      PS_true_EM_compr = rapply(object = PS_true_EM_compr, f = round, classes = "numeric", how = "replace", digits = 3)
      colnames(PS_true_EM_compr) = mgsub(colnames(PS_true_EM_compr), c("prob.1", "prob.2", "prob.3"),
                                         c("prob_as", "prob_pro", "prob_ns"))
      PS_true_EM_compr = data.frame(id = PS_true_EM_compr$id, g = PS_true_EM_compr$g,
                                    prob_as = PS_true_EM_compr$prob_as, EMest_p_as=PS_true_EM_compr$EMest_p_as, 
                                    diff = PS_true_EM_compr$prob_as - PS_true_EM_compr$EMest_p_as,
                                    prob_pro = PS_true_EM_compr$prob_pro, EMest_p_pro=PS_true_EM_compr$EMest_p_pro, 
                                    prob_ns = PS_true_EM_compr$prob_ns, EMest_p_ns=PS_true_EM_compr$EMest_p_ns)
      return(list(PS_true_EM_compr=PS_true_EM_compr,OBS_table=OBS_table, pis=pis))
    }
    
    
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
    X_sub_cols = paste0("X", c(1:(dim_x)))
    # TODO if I want to allow matching (mahalanobis, euclead) and rfegression to see the true x
    # TODO i.e. misspecification ONLY in the PS model!!!
    # TODO usually, it should be FALSE.
    # TODO @@@ CHECK THAT IT DOES WHAT I MEANT TO
    if(match_and_reg_watch_true_X == TRUE & misspec_PS==1){
      data_with_PS = data.table(subset(data_with_PS, 
                select = -grep(paste(X_sub_cols, collapse="|"), colnames(data_with_PS))), x_PS)
    }
    
    # run for all options (3 options)
    # m_data = data_with_PS; data_with_PS[OBS != "O(0,0)"]; m_data = data_with_PS[S==1]
    data_list = list(data_with_PS, data_with_PS[OBS != "O(0,0)"], data_with_PS[S==1]) 
    #save(data_list50, file = "data_list50.RData")
    
    if(!exists("one_type_replace")){
      print("NO one_type_replace")
      lst_matching_estimators_end_excluded_included = list()
      replace_vec = c(FALSE, TRUE)
      for(j in c(1:length(replace_vec))){
        lst_matching_estimators_end_excluded_included[[j]] =
          lapply(1:length(data_list), function(l){
            # my_matching_func_basic # my_matching_func_multiple
            my_matching_func_multiple(match_on = match_on, X_sub_cols, data_list[[l]],
                                      weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2,
                                      min_PS = min_PS, min_diff_PS = min_diff_PS,
                                      caliper = caliper, OBS_table, change_id = TRUE, mu_x_fixed=mu_x_fixed, x_as=x_as)
          })
      }
      
      
      # TODO 1.
      # matching estimators from all data sets, with and wout replacements
      
      # TODO @@@ WOUT bias correction
      # matching_estimators = lapply(1:length(replace_vec), function(j){
      #   data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]], head, 2)))))
      # })
      #matching_estimators = list.cbind(matching_estimators)
      #colnames(matching_estimators) = paste0(rep(paste0( "rep", rep(substr(replace_vec,1,1), each=6), "_",
      #               paste0(rep(c("MATCH", "MATCH"), each = 3),
      #               "_", c("all", "wout_O_0_0", "S1")),
      #               rep(c(rep("", 3), rep("_HL", 3)), 2)), each = length(c("EST","SE"))), "_", 
      #               rep(c("est","SE"), times=length(matching_estimators)/2))
      
      
      # with bias correction
      matching_estimators = lapply(1:length(replace_vec), function(j){
        data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]], head, 6)))))
      })
      
      # TODO BASIC MATCHING, with only 1 matching type and BC
      # matching_estimators = list.cbind(matching_estimators)
      # colnames(matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(matching_estimators)/2)),
      #             "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), each=2, times=3*2),
      #             rep(c("","_HL","_BC","_BCclpr"), each=6, times=2),
      #             rep(c("_est","_SE"), times=length(matching_estimators)/2))
      
      # MULTIPLE MATCHING WITH BC
      matching_estimators = list.cbind(matching_estimators)
      colnames(matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(matching_estimators)/2)),
                                             "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), each=2, times=3*2*2),
                                             rep(c("_PS", "_maha", "","_HL","_BC","_BCclpr"), each=6, times=2),
                                             rep(c("_est","_SE"), times=length(matching_estimators)/2)) 
      
      
      # CI naive and HL
      # CI_matching_estimators_lst = lapply(1:length(replace_vec), function(j){
      #   lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "CI_crude_HL_BC")})
      # 
      # CI_matching_estimators = lapply(1:length(replace_vec), function(j){
      #   data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
      #                                         "[[", "CI_crude_HL_BC")))))
      # })
      # 
      
      CI_matching_estimators = lapply(1:length(replace_vec), function(j){
        as.vector(t(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                               "[[", "CI_crude_HL_BC"))))))
      })
      CI_matching_estimators = data.frame((unlist(CI_matching_estimators))) %>% t
      rownames(CI_matching_estimators) = ""
      
      # TODO WOUT BC
      # colnames(CI_matching_estimators) = paste0( "rep", rep(substr(replace_vec,1,1), each=6), "_",
      #            paste0(rep(c("MATCH", "MATCH"), each = 3), "_", c("all", "wout_O_0_0", "S1")), rep(c(rep("", 3), rep("_HL", 3)), 2))
      # 
      
      # TODO BASIC MATCHING WITH BC
      # colnames(CI_matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(CI_matching_estimators)/2)),
      #   "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), times=3*2), rep(c("","_HL","_BC","_BCclpr"), each=3, times=2))
      
      # MULTIPLE MATCHING WITH BC
      colnames(CI_matching_estimators) = paste0(paste0("rep", rep(substr(replace_vec,1,1), each=length(CI_matching_estimators)/2)),
                                                "_", "MATCH_", rep(c("all", "wout_O_0_0", "S1"), times=3*2), rep(c("_PS", "_maha","","_HL","_BC","_BCclpr"), each=3, times=2))
      
      
      # rep_bool_false_true = 2 for WLS and 1 for OLS
      #CI_or_TABLE_EST_SE = "CI_LS"
      # rep_bool_false_true = 2; estimator_str = "WLS_NOinteractions_reg_adj_estimators_and_se"
      arrange_lin_reg_estimators = function(rep_bool_false_true, estimator_str,
                                            CI_or_TABLE_EST_SE="estimator_and_se_estimator", name=""){
        LS_lin_reg_estimators = lapply(rep_bool_false_true:rep_bool_false_true, function(j){
          lapply(lst_matching_estimators_end_excluded_included[[j]], "[[",
                 estimator_str)})
        LS_lin_reg_estimators = lapply(LS_lin_reg_estimators[[1]], "[[", CI_or_TABLE_EST_SE)
        if(CI_or_TABLE_EST_SE=="estimator_and_se_estimator"){
          colnas = paste0(rep(colnames(LS_lin_reg_estimators[[1]]), times=3), "_",
                          rep(c("all", "wout_O_0_0", "S1"), each=length(LS_lin_reg_estimators[[1]])))
          LS_lin_reg_estimators = list.cbind(LS_lin_reg_estimators)
          colnames(LS_lin_reg_estimators) = colnas
          
        }else{
          colnas = name
          LS_lin_reg_estimators = list.cbind(LS_lin_reg_estimators)
          colnames(LS_lin_reg_estimators) =  paste0(colnas, "_", c("all", "wout_O_0_0", "S1"))
        }
        
        return(LS_lin_reg_estimators)
      }
      
      # WLS ONLY FOR replace_vec == TRUE
      # WLS wout interactions
      WLS_NOint_matching_reg_estimators = arrange_lin_reg_estimators(2, "WLS_NOinteractions_reg_adj_estimators_and_se")
      WLS_NOint_matching_reg_estimators_CI = arrange_lin_reg_estimators(2, "WLS_NOinteractions_reg_adj_estimators_and_se", "CI_LS", name="WLS_NOint")
      # WLS with interactions
      WLS_YESint_matching_reg_estimators = arrange_lin_reg_estimators(2, "WLS_YESinteractions_reg_adj_estimators_and_se")
      WLS_YESint_matching_reg_estimators_CI = arrange_lin_reg_estimators(2, "WLS_YESinteractions_reg_adj_estimators_and_se", "CI_LS", name="WLS_YESint")
      # OLS ONLY FOR replace_vec == FALSE
      # OLS wout interactions
      OLS_NOint_matching_reg_estimators = arrange_lin_reg_estimators(1, "OLS_NOinteractions_reg_adj_estimators_and_se")
      OLS_NOint_matching_reg_estimators_CI = arrange_lin_reg_estimators(1, "OLS_NOinteractions_reg_adj_estimators_and_se", "CI_LS", name="OLS_NOint")
      # OLS with interactions
      OLS_YESint_matching_reg_estimators = arrange_lin_reg_estimators(1, "OLS_YESinteractions_reg_adj_estimators_and_se")
      OLS_YESint_matching_reg_estimators_CI = arrange_lin_reg_estimators(1, "OLS_YESinteractions_reg_adj_estimators_and_se", "CI_LS", name="OLS_YESint")
      
      
      # TODO Coverage CI
      
      
      # TODO 2.
      # excluded_included_matching from matching function
      # if we have replace FALSE or TRUE, use only the FALSE replace,
      # since when replacing this info is not informative
      excluded_included_matching = 
        as.numeric(sapply(lst_matching_estimators_end_excluded_included[[1]], 
                          "[[", "included_excluded_in_matching"))
      names(excluded_included_matching) = 
        paste0(rep(c("all", "wout_O_0_0", "S1" ), each = 5), "_", 
               rep(colnames(t(sapply(lst_matching_estimators_end_excluded_included[[1]],
                                     "[[", "included_excluded_in_matching")))))
      
      # TODO 3.
      # repeated_as_and_pro
      repeated_as_and_pro_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "repeated_as_and_pro")})
      repeated_as_and_pro_lst = unlist(repeated_as_and_pro_lst, recursive = FALSE)
      
      # TODO 6. tables_matched_units
      matched_units_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]],  
               "[[", "tables_matched_units")})
      matched_units_lst = unlist(matched_units_lst, recursive = FALSE)
      
      # TODO 4.
      # diff_distance_aspr_asas, 1 number per each part, so we have 6 numbers in total
      diff_distance = lapply(1:length(replace_vec), function(j){
        data.frame(t(unlist(list.rbind(lapply(lst_matching_estimators_end_excluded_included[[j]],
                                              "[[", "diff_distance_aspr_asas")))))
      })
      diff_distance = list.cbind(diff_distance)
      colnames(diff_distance) = paste0("d_", paste0( "rep", rep(substr(replace_vec,1,1), each=3), "_",
                                                     paste0(rep(c("MATCH", "MATCH"), each = 3),
                                                            "_", c("all", "wout_O_0_0", "S1"))))
      
      
      # TODO 5.
      # standardized mean diff as2pro, as2as per each data set, with and wout replacements
      # each element in  the list is a matrix
      # with all covariates std diff per each part, 6 parts in total
      std_mean_diff_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]], tail, 1)})
      std_mean_diff_lst = unlist(std_mean_diff_lst, recursive = FALSE)
      
      # TODO 7. means_by_subset
      means_by_subset_lst = lapply(1:length(replace_vec), function(j){
        lapply(lst_matching_estimators_end_excluded_included[[j]], "[[", "means_by_subset")})
      means_by_subset_mat = list.rbind(unlist(means_by_subset_lst, recursive = FALSE))
      rownames(means_by_subset_mat) = paste0(rep(c("reF_", "reT_"), each=18), 
                                             rep(c("all", "wout_0_0", "S1"), each=6), "_", rownames(means_by_subset_mat))
      
    }
    
    
    # TODO run my_matching_func_pairmatch similar process to the case 
    # where if(exists("one_type_replace") == TRUE)
    if(pairmatch_bool == TRUE){
      print("pairmatch")
      matching_estimators_end_excluded_included_pairmatch = 
        lapply(1:length(data_list), function(l){
          my_matching_func_pairmatch(X_sub_cols, data_list[[l]], 
                                     M=1, replace = FALSE,
                                     estimand = "ATC", mahal_match = 2, caliper = caliper, OBS_table)
        })
      
      matching_estimators_pairmatch = 
        list.rbind(lapply(matching_estimators_end_excluded_included_pairmatch,
                          head, 2))
      matching_estimators_pairmatch = 
        data.frame(t(unlist(matching_estimators_pairmatch)))
      colnames(matching_estimators_pairmatch) = 
        paste0(rep(c("pairMATCH_keepS0", "pairMATCH_removeS0"), each = 3),
               "_", c("all", "wout_O_0_0", "S1"))
      
      # excluded_included_matching from matching function
      excluded_included_matching_pairmatch = 
        as.numeric(sapply(matching_estimators_end_excluded_included_pairmatch,
                          "[[", 3))
      
      names(excluded_included_matching_pairmatch) = paste0(rep(c("all", "wout_O_0_0", "S1" ),each = 5), "_",
                                                           rep(colnames(t(sapply(matching_estimators_end_excluded_included_pairmatch,
                                                                                 "[[", 3)))))
      
      # TODO add std_mean_diff
      
    }
    
    
    # TODO 1. put all results together in the current row of mat_param_estimators
    mat_param_estimators = rbind( mat_param_estimators,
                                  data.frame(SACE, SACE_conditional, 
                                             DING_est, DING_model_assisted_est_ps
                                             #,matching_estimators_pairmatch
                                             ,matching_estimators
                                             , diff_distance,
                                             most_naive_est, most_naive_est_se, omit_Y_less0_naive_est, sur_naive_est, sur_naive_est_se,
                                             pis, vec_OBS_table
                                  ))
    
    # TODO 1. b
    
    WLS_NOint_mat_reg_estimators = rbind(WLS_NOint_mat_reg_estimators,
                                         data.frame(SACE, SACE_conditional, WLS_NOint_matching_reg_estimators) )
    
    WLS_YESint_mat_reg_estimators = rbind(WLS_YESint_mat_reg_estimators,
                                          data.frame(SACE, SACE_conditional, WLS_YESint_matching_reg_estimators) )
    
    OLS_NOint_mat_reg_estimators = rbind(OLS_NOint_mat_reg_estimators,
                                         data.frame(SACE, SACE_conditional, OLS_NOint_matching_reg_estimators) )
    
    OLS_YESint_mat_reg_estimators = rbind(OLS_YESint_mat_reg_estimators,
                                          data.frame(SACE, SACE_conditional, OLS_YESint_matching_reg_estimators) )
    
    CI_mat = rbind( CI_mat, data.frame(SACE, SACE_conditional
                                       ,CI_naives_before_matching, CI_matching_estimators
                                       ,WLS_NOint_matching_reg_estimators_CI, WLS_YESint_matching_reg_estimators_CI
                                       ,OLS_NOint_matching_reg_estimators_CI, OLS_YESint_matching_reg_estimators_CI) )
    
    # TODO 2. put all results together in the current row of mat_excluded_included_matching
    mat_excluded_included_matching = 
      rbind(mat_excluded_included_matching, excluded_included_matching)
    
    # TODO 3. list of repeated_as_and_pro
    list_repeated_as_and_pro[[i]] = list.rbind(repeated_as_and_pro_lst)
    rownames(list_repeated_as_and_pro[[i]])=
      paste0("rep", rep(substr(replace_vec,1,1), each=3), "_", c("all", "wout_O_0_0", "S1"))
    
    # TODO 6. tables_matched_units
    list_matched_units[[i]] = rbindlist(matched_units_lst, use.names=FALSE)
    rownames(list_matched_units[[i]])=
      paste0("rep", rep(substr(replace_vec,1,1), each=3), "_", c("all", "wout_O_0_0", "S1"))
    
    # TODO 4. matrix of diff_abs_distance (number paer each part, so 6 in total)
    # this info is also in mat_param_estimators
    mat_diff_distance_aspr_asas = rbind(mat_diff_distance_aspr_asas, diff_distance)
    
    # TODO 5. list of std_diffs per all X's; size of list as the samples
    list_std_mean_diff[[i]] = std_mean_diff_lst
    
    # TODO 7. list_means_by_subset
    list_means_by_subset[[i]] = means_by_subset_mat
    
    end_time1 <- Sys.time()
    print(paste0("one iteration lasts ", difftime(end_time1, start_time1)))
    print(difftime(end_time1, start_time1))
  }
  ################# out of for loop for all the samples (all in all: param_n_sim samples) #######
  
  
  # TODO 1.
  # summary of mat_param_estimators: mean and sd
  param_SACE = mean(mat_param_estimators$SACE)
  MSE_fun <- function (x) mean((x-param_SACE)^2) # param_SACE SACE mean(x)
  mat_param_estimators = rbind(mat_param_estimators,
                               mean = apply(mat_param_estimators, 2, mean), 
                               med = apply(mat_param_estimators, 2, median),
                               SD = apply(mat_param_estimators, 2, sd)
                               , MSE <- list.cbind(lapply(mat_param_estimators, FUN = MSE_fun))
  )
  
  rownames(mat_param_estimators) = 
    c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  # c(c(1:param_n_sim),c("mean","med","sd","MSE))
  #rownames(mat_param_estimators)[(nrow(mat_param_estimators) - 1): nrow(mat_param_estimators)] = c("mean", "sd")
  
  # TODO 1.b
  # summary of mat_param_estimators: mean and sd
  list_reg_LS = list(WLS_NOint_mat_reg_estimators=WLS_NOint_mat_reg_estimators, WLS_YESint_mat_reg_estimators=WLS_YESint_mat_reg_estimators,
                     OLS_NOint_mat_reg_estimators=OLS_NOint_mat_reg_estimators, OLS_YESint_mat_reg_estimators=OLS_YESint_mat_reg_estimators)
  for(i in 1:length(list_reg_LS)){
    list_reg_LS[[i]] = rbind(list_reg_LS[[i]],
                             mean = apply(list_reg_LS[[i]], 2, mean), 
                             med = apply(list_reg_LS[[i]], 2, median),
                             SD = apply(list_reg_LS[[i]], 2, sd)
                             , MSE <- list.cbind(lapply(list_reg_LS[[i]], FUN = MSE_fun))
    )
    rownames(list_reg_LS[[i]]) = 
      c(c(1:param_n_sim),c("mean","med","sd","MSE"))
  }
  
  apply(CI_mat[,c(1:2)], 2, as.numeric)
  #CI_mat = data.frame(rbind(CI_mat, apply(CI_mat, 2, function (x) if(is.numeric(x)){mean(x)}else{0})))
  #CI_mat[(param_n_sim + 1),] = c(mean(CI_mat$SACE), mean(CI_mat$SACE_conditional))
  CI_mat$SACE = mean(CI_mat$SACE); CI_mat$SACE_conditional  = mean(CI_mat$SACE_conditional )
  
  # calculating EM estimators and arrange in a data frame
  # TODO put this thing in a function
  coeffs_as = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_as, "[[", l)
  })
  coeffs_ns = lapply(c(1:dim_x), function(l){
    sapply(list_coeff_ns, "[[", l)
  })
  coeffs = data.table( rbind( list.rbind(coeffs_as), list.rbind(coeffs_ns) ) )
  
  # coeff_as_0 = sapply(list_coeff_as, "[[", 1); coeff_as_1 = sapply(list_coeff_as, "[[", 2)
  # coeff_as_2 = sapply(list_coeff_as, "[[", 3); coeff_as_3 = sapply(list_coeff_as, "[[", 4)
  # coeff_ns_0 = sapply(list_coeff_ns, "[[", 1); coeff_ns_1 = sapply(list_coeff_ns, "[[", 2)
  # coeffs = data.table(rbind(coeff_as_0, coeff_as_1, coeff_ns_0, coeff_ns_1))
  
  print("mean and sd")
  coeffs[, `:=` (SD = rowSds(as.matrix(coeffs)), mean = rowMeans(coeffs))]
  coeffs$parameter = c(gamma_as, gamma_ns)
  coeffs$diff = coeffs$mean - coeffs$parameter 
  coeffs$perc = coeffs$diff / abs(coeffs$parameter)
  coeffs_df = data.frame(coeffs)
  rownames(coeffs_df) = colnames(mat_gamma)
  #rownames(coeffs_df) = c("coeff_as_0", "coeff_as_1","coeff_ns_0", "coeff_ns_1")
  
  # TODO 2. mean and sd for mat_excluded_included_matching
  mat_excluded_included_matching = rbind(mat_excluded_included_matching, 
                                         mean = apply(mat_excluded_included_matching, 2, mean),
                                         sd = apply(mat_excluded_included_matching, 2, sd))
  # change order of mat_excluded_included_matching
  mat_excluded_included_matching = data.frame(mat_excluded_included_matching)
  ind1 = c(1,6,11)
  mat_excluded_included_matching = subset(mat_excluded_included_matching,
                                          select = c(ind1 , ind1+1, ind1+2, ind1+3, ind1+4))
  
  
  # TODO 3. list of repeated_as_and_pro
  # TODO could be nice to check this strait after the matching (before removing pairs with non-survivor) 
  # and the change before and after
  # TODO and also by proportion: if pi_as = 0.6, pi_pro = 0.3, so adjust to that, and then check
  # repeated_as_and_pro_lst
  mean_list_repeated_as_and_pro = 
    calculate_mean_repeated_as_and_pro(list_repeated_as_and_pro, TRUE)
  
  # TODO 6 THE SAME AS FOR list_repeated_as_and_pro to list_matched_units
  mean_list_matched_units = 
    calculate_mean_repeated_as_and_pro(list_matched_units, FALSE)
  
  
  # TODO 4. matrix of diff_abs_distance (number paer each part, so 6 in total)
  # I can calculate mean, but I think there is no need
  
  # TODO 5. list of std_diffs per all X's; size of list as the samples
  mean_list_std_mean_diff = calculate_mean_of_dfs_in_lists_of_lists(list_std_mean_diff)
  
  # TODO 7. list_means_by_subset
  mean_list_means_by_subset = calculate_mean_repeated_as_and_pro(list_means_by_subset, FALSE)
  
  
  return(list(mat_param_estimators = mat_param_estimators, 
              WLS_NOint_mat_reg_estimators = list_reg_LS$WLS_NOint_mat_reg_estimators,
              WLS_YESint_mat_reg_estimators = list_reg_LS$WLS_YESint_mat_reg_estimators,
              OLS_NOint_mat_reg_estimators = list_reg_LS$OLS_NOint_mat_reg_estimators,
              OLS_YESint_mat_reg_estimators = list_reg_LS$OLS_YESint_mat_reg_estimators,
              CI_mat = CI_mat,
              coeffs_df = coeffs_df, 
              mat_excluded_included_matching = mat_excluded_included_matching, 
              mean_list_repeated_as_and_pro = mean_list_repeated_as_and_pro,
              mat_diff_distance_aspr_asas = mat_diff_distance_aspr_asas,
              mean_list_std_mean_diff = mean_list_std_mean_diff,
              mean_list_matched_units = mean_list_matched_units,
              mean_list_means_by_subset = mean_list_means_by_subset
  ))
}










