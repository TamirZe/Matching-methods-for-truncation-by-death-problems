m_data = data_list[[1]]
# caliper is in sd
w_mat_bool = "NON-INFO"; replace = TRUE; estimand = "ATC"; change_id = TRUE; mahal_match = 2; M=1
reg_cov = reg_after_match # c("age", "education", "black", "hispanic", "married", "re74" ,"re75", "emp74", "emp75")
X_sub_cols = variables

matching_func_multiple_data = function(match_on = NULL,
           cont_cov_mahal = c("age", "education", "re74", "re75"),  reg_cov, X_sub_cols, 
           reg_BC = c("age", "education", "re74", "re75"), m_data, 
           w_mat_bool = "NON-INFO", M=1, replace, estimand = "ATC", mahal_match = 2, caliper = 0.25
           #,OBS_table
           ,change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units = FALSE, one_leraner_bool=FALSE, vertical_table = TRUE, rnd=1){
  if(change_id == TRUE){
    print("change id")
    m_data$id = c(1:nrow(m_data))
  }
  
  # check balance before 
  # cont and discrete
  disc_var = names(which(apply(subset(m_data, select = X_sub_cols[-1]), 2, function(x) { all(x %in% 0:1)})))
  cont_var = setdiff(X_sub_cols[-1], disc_var)
  # verically
  if(vertical_table == TRUE){
    balance_before_match = covarites_descriptive_table_cont_disc(dat = m_data, cov_descr = X_sub_cols, metric = "Surv", rnd=1)
  # horizontically
  } else{
     balance_before_match = subset(m_data, select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as",
                                                     X_sub_cols[-1]))[,lapply(.SD, mean), by="A"] 
     balance_before_match$N = m_data[, .N, by="A"]$N
     balance_before_match = rbind(balance_before_match[2,], balance_before_match[1,])
     balance_before_match = data.frame(Metric = "Surviors", subset(balance_before_match, select = c(A, N)),
                                      N_match = "", N_unq = "", subset(balance_before_match, select = -c(A, N)))
     est_var_x0 = apply(subset(filter(m_data, A==0), select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1])), 2, var)
     est_var_x1 = apply(subset(filter(m_data,A==0), select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1])), 2, var)
     sd_smd = sqrt(est_var_x0)
     #sd_smd = sqrt(0.5 * (est_var_x0 + est_var_x1))
     SMD = (balance_before_match[1,-c(1:5)] - balance_before_match[2,-c(1:5)]) / sd_smd[-1]
     balance_before_match[3,] = c(rep(" ",5), mutate_if(SMD, is.numeric, round, 3))
     balance_before_match = balance_before_match %>% mutate_if(is.numeric, round, 3)
  }

   
  print(paste0("replace is ", replace, " nrows is ", nrow(m_data)))
  #X_sub_cols = paste0("X", c(1:(dim_x)))
  # mahal_match for Weight = 2 for mahalanobis distance. 1 for inverse of variance
  vec_caliper = c(rep(1001, length(cont_cov_mahal)), caliper)
  # set weights matrix for mahalanobis for continous covariates
  
  if(w_mat_bool == "INVERSE_SD"){
    #TODO INVERSE SD
    w_mat = diag(length(cont_cov_mahal) + 1) /
      ( c(apply(subset(m_data, select = cont_cov_mahal), 2, sd), 1) )
    w_mat[nrow(w_mat), ncol(w_mat)]  = 0
    w_mat = w_mat * (1 / sum(w_mat))
  }else if(w_mat_bool == "NON-INFO"){
    #TODO NON-INFORMATIVE weights matrix
    w_mat = diag(length(cont_cov_mahal) + 1) / length(cont_cov_mahal)
    w_mat[nrow(w_mat), ncol(w_mat)]  = 0
  }else if(w_mat_bool == "RE74_75"){
    #TODO more weights to re74 re75 weights matrix
    w_mat = diag(length(cont_cov_mahal) + 1) 
    w_mat[nrow(w_mat), ncol(w_mat)]  = 0
    diag(w_mat)[grep("re74|re75", cont_cov_mahal)] = 2
    w_mat = w_mat / sum(diag(w_mat))
  }
  
  # TODO 2. MAHALANOBIS WITHOUT PS CALIPER
  print("MAHALANOBIS WITHOUT PS CALIPER")
  #set.seed(102)
  MATCH_MAHA_wout_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                               #, X=m_data[,"est_p_as"]
                               , X = subset(m_data, select = cont_cov_mahal)
                               ,ties=FALSE
                               #,caliper = vec_caliper ,Weight.matrix = w_mat
                               ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  )
  mala_wout_cal_lst = arrange_dataset_after_matching_DATA(m_data=m_data, match_obj=MATCH_MAHA_wout_PS,
                                                          replace_bool=replace, X_sub_cols=X_sub_cols)
  dt_match_S1_maha_wout_cal = mala_wout_cal_lst$dt_match_S1; matched_pairs_mala_wout_cal = mala_wout_cal_lst$matched_pairs
  diff_per_pair_maha_wout_cal = dt_match_S1_maha_wout_cal$Y - dt_match_S1_maha_wout_cal$A0_Y
  matched_pairs_maha_wout_cal = mala_wout_cal_lst$matched_pairs
  data_pairs_maha_wout_cal = merge(matched_pairs_maha_wout_cal, m_data, by="id", all.x = TRUE, all.y = FALSE) %>% arrange(pair, A)
  crude_inference_lst_mala_wout_cal = 
    crude_estimator_inference(match_obj=MATCH_MAHA_wout_PS, dt_match_S1_maha_wout_cal, diff_per_pair_maha_wout_cal, replace_bool=replace)
  est_crude_maha_wout_cal = crude_inference_lst_mala_wout_cal$SACE_matching_est
  se_crude_maha_wout_cal = crude_inference_lst_mala_wout_cal$SACE_matching_SE
  CI_by_SE_and_Z_val_naive_maha_wout_cal = crude_inference_lst_mala_wout_cal$CI_by_SE_and_Z_val_crude
  balance_maha_wout_cal = balance_after_matching(m_data=m_data, match_obj=MATCH_MAHA_wout_PS, dt_match=dt_match_S1_maha_wout_cal, 
                                                 X_sub_cols=X_sub_cols, metric="Mahal", replace=replace, smd_se="not_weighted", vertical_table=vertical_table)
  
  #TODO 1. MATCHING ONLY ONLY PS: EMest_p_as
  print(match_on)
  #set.seed(101)
  MATCH_PS_only  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                          , X = subset(m_data, select = match_on)
                          #, X=m_data[,"EMest_p_as"]         
                          #, X = subset(m_data, select = c(X_sub_cols[-1], "EMest_p_as"))
                          ,ties=FALSE
                          #, distance.tolerance = 1e-10
                          #,caliper = vec_caliper
                          ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                          #,Weight.matrix = w_mat
  )
  only_ps_lst = arrange_dataset_after_matching_DATA(m_data=m_data, match_obj=MATCH_PS_only, 
                          replace_bool=replace, X_sub_cols=X_sub_cols)
  dt_match_S1_only_ps = only_ps_lst$dt_match_S1; matched_pairs_only_ps = only_ps_lst$matched_pairs
  diff_per_pair_only_ps = dt_match_S1_only_ps$Y - dt_match_S1_only_ps$A0_Y 
  matched_pairs_only_ps = only_ps_lst$matched_pairs
  data_pairs_only_ps = merge(matched_pairs_only_ps, m_data, by="id", all.x = TRUE, all.y = FALSE) %>% arrange(pair, A)
  crude_inference_lst_only_ps = crude_estimator_inference(match_obj=MATCH_PS_only, dt_match_S1=dt_match_S1_only_ps,
                              diff_per_pair=diff_per_pair_only_ps, replace_bool=replace)
  est_crude_only_ps = crude_inference_lst_only_ps$SACE_matching_est
  se_crude_only_ps = crude_inference_lst_only_ps$SACE_matching_SE
  CI_by_SE_and_Z_val_naive_only_ps = crude_inference_lst_only_ps$CI_by_SE_and_Z_val_crude
  balance_only_ps = balance_after_matching(m_data=m_data, match_obj=MATCH_PS_only, dt_match=dt_match_S1_only_ps, 
     X_sub_cols=X_sub_cols, metric="PS", replace=replace, smd_se="not_weighted", vertical_table=vertical_table)
  
  
  # TODO 3. MAHALANOBIS WITH PS CALIPER
  print("MAHALANOBIS WITH PS CALIPER")
  #set.seed(103)
  ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                         #, X=m_data[,"est_p_as"]
                         , X = subset(m_data, select = c(cont_cov_mahal, match_on))
                         ,ties=FALSE
                         #, distance.tolerance = 1e-10
                         ,caliper = vec_caliper ,Weight.matrix = w_mat
                         ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  )
  
  weights = data.frame(table(c(ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control)))
  colnames(weights) = c("id", "w")
  matched_data_mahalPScal = expandRows(merge(weights, m_data, by="id"), "w") %>% arrange(id)
  ATE_MATCH_PS_lst = arrange_dataset_after_matching_DATA(m_data=m_data, match_obj=ATE_MATCH_PS, replace_bool=replace, X_sub_cols=X_sub_cols)
  dt_match_S1 = ATE_MATCH_PS_lst$dt_match_S1
  diff_per_pair = dt_match_S1$Y - dt_match_S1$A0_Y
  matched_pairs = ATE_MATCH_PS_lst$matched_pairs
  data_pairs_maha_cal_PS = merge(matched_pairs, m_data, by="id", all.x = TRUE, all.y = FALSE) %>% arrange(pair, A)
  crude_inference_lst = crude_estimator_inference(match_obj=ATE_MATCH_PS, dt_match_S1, diff_per_pair, replace_bool=replace)
  SACE_matching_est = crude_inference_lst$SACE_matching_est
  SACE_matching_SE = crude_inference_lst$SACE_matching_SE
  CI_by_SE_and_Z_val_naive = crude_inference_lst$CI_by_SE_and_Z_val_crude
  balance_maha_cal_PS = balance_after_matching(m_data=m_data, match_obj=ATE_MATCH_PS, dt_match=dt_match_S1,
     X_sub_cols=X_sub_cols, metric="Mahal cal", replace=replace, smd_se="not_weighted", vertical_table=vertical_table)
  # MATCH_obj = ATE_MATCH_PS
  
  extract_matched_set = function(MATCH_obj, m_data){
    amount_ctr = nrow(m_data) - sum(m_data$A)
    ind_ctr = MATCH_obj$index.control; ind_trt = MATCH_obj$index.treated
    ind_ctr_unq = unique(ind_ctr); ind_trt_unq = unique(ind_trt)
    print(paste0("len ctr_before match = ",amount_ctr, ", len ctr after match = ", length(ind_ctr), 
                 ", len ctr after match = ", length(ind_ctr_unq)))
    print(paste0("ind_ctr_unq and ind_ctr are identical: ", identical(ind_ctr_unq,ind_ctr)))
    weights = apply(data.frame(table(c(ind_ctr, ind_trt))), 2, as.numeric); 
    colnames(weights) = c("id", "w")
    matched_set = merge(weights, m_data, by = "id")
    print(paste0("cunique weights for controla are really = ", unique(filter(matched_set, A==0)$w)))
    return(matched_set=matched_set)
  }
  
  if(one_leraner_bool){
    # amount_ctr = nrow(m_data) - sum(m_data$A)
    # ind_ctr = ATE_MATCH_PS$index.control; ind_trt = ATE_MATCH_PS$index.treated
    # ind_ctr_unq = unique(ind_ctr); ind_trt_unq = unique(ind_trt)
    # print(paste0("len ctr_before match = ",amount_ctr, ", len ctr after match = ", length(ind_ctr), 
    #              ", len ctr after match = ", length(ind_ctr_unq)))
    # print(paste0("ind_ctr_unq and ind_ctr are identical: ", identical(ind_ctr_unq,ind_ctr)))
    # weights = apply(data.frame(table(c(ind_ctr, ind_trt))), 2, as.numeric); 
    # colnames(weights) = c("id", "w")
    # matched_set = merge(weights, m_data, by = "id")
    # print(paste0("unique weights for controla are really = ", unique(filter(matched_set, A==0)$w)))
    
    matched_set_PS_only = extract_matched_set(MATCH_PS_only, m_data)
    matched_set_MAHA_wout_PS = extract_matched_set(MATCH_MAHA_wout_PS, m_data)
    matched_set = extract_matched_set(ATE_MATCH_PS, m_data) # mahal with caliper
    matched_set_lst = list(PS=matched_set_PS_only, Mahal=matched_set_MAHA_wout_PS,
                           Mahal_PS_cal=matched_set)
    
    reg_data_matched_PS_only = merge(matched_pairs_only_ps, matched_set_PS_only, by="id") %>% arrange(pair, A)
    reg_data_matched_MAHA_wout_PS = merge(matched_pairs_maha_wout_cal, matched_set_MAHA_wout_PS, by="id") %>% arrange(pair, A)
    reg_data_matched = merge(matched_pairs, matched_set, by="id") %>% arrange(pair, A)
    reg_data_matched_lst = list(PS=reg_data_matched_PS_only, Mahal=reg_data_matched_MAHA_wout_PS,
                           Mahal_PS_cal=reg_data_matched)
    return(list(m_data=m_data, ATE_MATCH_PS_lst=ATE_MATCH_PS_lst,
                reg_data_matched_lst=reg_data_matched_lst, matched_set_lst=matched_set_lst))
  }
  
  # balance_maha_cal_PS = data.frame(
  # rbind(apply(subset(dt_match_S1, 
  #          select = paste0("A0_", c("A", "EMest_p_as", "EMest_p_pro", X_sub_cols[-1]))), 2, mean),
  #         apply(subset(dt_match_S1, 
  #          select = c("A", "EMest_p_as", "EMest_p_pro", X_sub_cols[-1])), 2, mean)) ) %>% round(3)
  # colnames(balance_maha_cal_PS) = substr(colnames(balance_maha_cal_PS), 4, 100)
  # balance_maha_cal_PS$N = c(nrow(filter(m_data, A==0)), nrow(filter(m_data, A==1)), "")
  # balance_maha_cal_PS$N_match = rep(nrow(dt_match_S1),3)
  # balance_maha_cal_PS$N_unq = c(length(unique(ATE_MATCH_PS$index.control)),length(unique(ATE_MATCH_PS$index.treated)))
  # balance_maha_cal_PS = data.frame(Metric = "Mahal cal", subset(balance_maha_cal_PS, select = c(A, N, N_match, N_unq)), 
  #                                  subset(balance_maha_cal_PS, select = -c(A, N, N_match, N_unq)))
  
  if(vertical_table==TRUE){
    balance_table = data.frame(Replacements = replace, cbind(filter(balance_before_match, Variable != "S"),
                                                             balance_only_ps$balance_match, balance_maha_wout_cal$balance_match, balance_maha_cal_PS$balance_match))
    balance_table = subset(balance_table, select = -grep("Variable.", colnames(balance_table)))
    balance_table_check = Reduce(function(x,y) merge(x = x, y = y, by = "Variable"), list(balance_before_match,
                                                                                          balance_only_ps$balance_match, balance_maha_wout_cal$balance_match, balance_maha_cal_PS$balance_match))
  }else{
    #rbind.fill
    balance_table = data.frame(Replacements = replace,
                               rbind(balance_before_match, balance_only_ps$balance_match, balance_maha_wout_cal$balance_match, balance_maha_cal_PS$balance_match))
    balance_table$Replacements = substr(balance_table$Replacements, 1, 1)
  }
  
  #summary(ATE_MATCH_PS)
  #ATE_MATCH_PS$est; ATE_MATCH_PS$se; ATE_MATCH_PS$nobs; ATE_MATCH_PS$index.dropped
  
  # TODO HL function
  func_wilcoxon_HL_est = function(boost_HL=FALSE, match_metric_lst, return_CI=FALSE){
    diff_per_pair = match_metric_lst$dt_match_S1$Y - match_metric_lst$dt_match_S1$A0_Y
    wilcoxon = wilcox.test(diff_per_pair, conf.int=T)
    SACE_matching_est_HL = as.numeric(wilcoxon$estimate) %>% round(3)
    SACE_matching_pval_HL = wilcoxon$p.value 
    # bootstrap for HL se
    if(boost_HL==TRUE){
      HL_boost_vec = c()
      for(i in 1:100){
        print(paste0("bootstrap ", i))
        d = match_metric_lst$dt_match_S1[sample(nrow(match_metric_lst$dt_match_S1), nrow(match_metric_lst$dt_match_S1), replace = T),]
        diff_per_pair_boost = d$Y - d$A0_Y
        HL_boost_vec[i] = wilcox.test(diff_per_pair_boost,conf.int=T)$estimate
      }
      SACE_matching_est_HL_bool = mean(HL_boost_vec)
      SACE_matching_se_HL = sd(HL_boost_vec)
    }else{
      # FAKE SE
      SACE_matching_se_HL = SACE_matching_pval_HL %>% round(3)
    }
    
    SACE_matching_CI_HL = c(as.character(as.numeric(round(wilcoxon$conf.int, 3))[1]), 
                            as.character(as.numeric(round(wilcoxon$conf.int, 3))[2]))
    SACE_matching_CI_HL = paste(SACE_matching_CI_HL, sep = ' ', collapse = " , ")
    
    if(return_CI==FALSE){
      return(c(HL_est=SACE_matching_est_HL, HL_se=SACE_matching_se_HL))
    }else{
      return(c(HL_est=SACE_matching_est_HL, HL_se=SACE_matching_se_HL, HL_CI=SACE_matching_CI_HL))
    }
  }
  
  # for each type of distance metric, run HL and OLS/WLS
  
  dt_and_pairs_match_lst = 
    list(only_ps=only_ps_lst, mala_wout_cal=mala_wout_cal_lst, mala_cal=ATE_MATCH_PS_lst)
  OLS_WLS_reg_lst <- coeffs_table <- HL_est_lst <- list()
  regression_est_se <-  regression_inter_est_se <- vector()
  
  for (i in 1:length(dt_and_pairs_match_lst)){
    print(names(dt_and_pairs_match_lst)[i])
    
    #TODO for now I use covariates. I dont really use reg_covariates.
    #TODO check what to do with matched_pairs
    if(replace == TRUE){
      # TODO WLS
      LS_NOinter =
        regression_adjusted_function(rep_bool=replace,
             dt_match_S1=dt_and_pairs_match_lst[[i]]$dt_match_S1, m_data=m_data,
             matched_pairs=dt_and_pairs_match_lst[[i]]$matched_pairs, covariates = reg_cov,
             interactions_bool = FALSE, LS="WLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
      LS_YESinter =
        regression_adjusted_function(rep_bool=replace,
             dt_match_S1=dt_and_pairs_match_lst[[i]]$dt_match_S1, m_data=m_data,
             matched_pairs=dt_and_pairs_match_lst[[i]]$matched_pairs, covariates = reg_cov,
             interactions_bool = TRUE, LS="WLS", mu_x_fixed=mu_x_fixed, x_as=x_as)

      # WLS_clstr_se # WLS_sandwi_se
      regression_est_se = c(regression_est_se, 
                            paste0(LS_NOinter$estimator_and_se$WLS_estimator %>% round(3), " (",
                            LS_NOinter$estimator_and_se$WLS_clstr_se %>% round(3), ")"))
      regression_inter_est_se = c(regression_inter_est_se, 
                            paste0(LS_YESinter$estimator_and_se$WLS_estimator %>% round(3),  " (",
                            LS_YESinter$estimator_and_se$WLS_clstr_se %>% round(3), ")"))
      
      # coeffs_NOinter = LS_NOinter$coeffs_table
      # coeffs_YESinter = LS_YESinter$coeffs_table
      
    }
    
    if(replace == FALSE){
      # TODO OLS
      LS_NOinter =
        regression_adjusted_function(rep_bool=replace,
               dt_match_S1=dt_and_pairs_match_lst[[i]]$dt_match_S1, m_data=m_data,
               matched_pairs=dt_and_pairs_match_lst[[i]]$matched_pairs, covariates = reg_cov, 
               interactions_bool = FALSE, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
      LS_YESinter =
        regression_adjusted_function(rep_bool=replace,
               dt_match_S1=dt_and_pairs_match_lst[[i]]$dt_match_S1, m_data=m_data,
               matched_pairs=dt_and_pairs_match_lst[[i]]$matched_pairs, covariates = reg_cov,
               interactions_bool = TRUE, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
      regression_est_se = c(regression_est_se, 
                            paste0(LS_NOinter$estimator_and_se$OLS_estimator %>% round(3), " (",
                            LS_NOinter$estimator_and_se$OLS_se %>% round(3) %>% round(3), ")"))
      regression_inter_est_se = c(regression_inter_est_se, 
                            paste0(LS_YESinter$estimator_and_se$OLS_estimator %>% round(3), " (",
                            LS_YESinter$estimator_and_se$OLS_se %>% round(3), ")"))
      # coeffs_NOinter = LS_NOinter$coeffs_table
      # coeffs_YESinter = LS_YESinter$coeffs_table
      }
  
    # OLS WLS estimators  
    OLS_WLS_reg_lst[[i]] = list(LS_NOinter=LS_NOinter, LS_YESinter=LS_YESinter)
    coeffs_table[[i]] = list(LS_NOinter=LS_NOinter$coeffs_table, LS_YESinter=LS_YESinter$coeffs_table)
    names(dt_and_pairs_match_lst)[i]                                               
    # HL estimator
    HL_est_lst[[i]] = 
      func_wilcoxon_HL_est(boost_HL=FALSE, dt_and_pairs_match_lst[[i]])
  }
  names(OLS_WLS_reg_lst) <- names(coeffs_table) <- names(HL_est_lst) <- names(dt_and_pairs_match_lst)
  print("finish lin reg!")
  crude_lst = list(crude_inference_lst_only_ps, crude_inference_lst_mala_wout_cal, crude_inference_lst)
  names(crude_lst) <- names(dt_and_pairs_match_lst)
  
  
  #TODO Bias Corrected estimator
  if(replace == TRUE){
    m_data_just_inter = m_data$A * subset(m_data, select = reg_BC)
    colnames(m_data_just_inter) = paste0("A_", colnames(m_data_just_inter))
    m_data_inter = data.table(m_data, m_data_just_inter)
    
    # PS only
    matchBC_PS <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = match_on),
                     Z = subset(m_data, select = reg_BC), BiasAdjust=TRUE
                     ,ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                     #, distance.tolerance = 1e-10
    )
    BCest_PS = round(as.numeric(matchBC_PS$est), 3)
    BCse_PS = round(matchBC_PS$se, 3)
    CI_by_SE_and_Z_val_BC_PS = round(BCest_PS + c(-1,1) * 1.96 * BCse_PS, 3)
    CI_by_SE_and_Z_val_BC_PS = paste(CI_by_SE_and_Z_val_BC_PS, sep = ' ', collapse = " , ")
    
    
    #set.seed(102)
    # TODO AI bias corrected, consider only when replace==TRUE
    matchBC <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = cont_cov_mahal),
                     Z = subset(m_data, select = reg_BC), BiasAdjust=TRUE
                     ,ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                     #, distance.tolerance = 1e-10
    )
    BCest = round(as.numeric(matchBC$est), 3)
    BCse = round(matchBC$se, 3)
    CI_by_SE_and_Z_val_BC = round(BCest + c(-1,1) * 1.96 * BCse, 3)
    CI_by_SE_and_Z_val_BC = paste(CI_by_SE_and_Z_val_BC, sep = ' ', collapse = " , ")
    
    #subset(m_data_inter, select = c(reg_BC, colnames(m_data_just_inter)))
    #set.seed(102)
    matchBC_inter <- Match(Y=m_data_inter[,Y], Tr=m_data_inter[,A], X = subset(m_data_inter, select = cont_cov_mahal), 
                           Z = subset(m_data_inter, select = c(reg_BC, colnames(m_data_just_inter))), BiasAdjust=TRUE,
                           ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
    )
    BCest_inter = round(as.numeric(matchBC_inter$est), 3)
    BCse_inter = round(matchBC_inter$se, 3)
    CI_by_SE_and_Z_val_BC_inter = round(BCest_inter + c(-1,1) * 1.96 * BCse_inter, 3)
    CI_by_SE_and_Z_val_BC_inter = paste(CI_by_SE_and_Z_val_BC_inter, sep = ' ', collapse = " , ")
    
    #set.seed(102)
    matchBC_clpr <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, 
                                                                  select = c(cont_cov_mahal, match_on)), 
                          Z = subset(m_data, select = reg_BC), BiasAdjust=TRUE
                          ,distance.tolerance = 1e-10
                          ,ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                          ,caliper = vec_caliper, Weight.matrix = w_mat
    )
    BCest_clpr = round(as.numeric(matchBC_clpr$est), 3)
    BCse_clpr = round(matchBC_clpr$se, 3)
    CI_by_SE_and_Z_val_BCclpr = round(BCest_clpr + c(-1,1) * 1.96 * BCse_clpr, 3)
    CI_by_SE_and_Z_val_BCclpr = paste(CI_by_SE_and_Z_val_BCclpr, sep = ' ', collapse = " , ")
    #summary.Match(matchBC, full=TRUE)
    
    #set.seed(102)
    matchBC_clpr_inter <- Match(Y=m_data_inter[,Y], Tr=m_data_inter[,A], X = subset(m_data_inter, 
                                                                  select = c(cont_cov_mahal, match_on)), 
                          Z = subset(m_data_inter, select = c(reg_BC, colnames(m_data_just_inter))), BiasAdjust=TRUE,
                          ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                          ,caliper = vec_caliper, Weight.matrix = w_mat
    )
    BCest_clpr_inter = round(as.numeric(matchBC_clpr_inter$est), 3)
    BCse_clpr_inter = round(matchBC_clpr_inter$se, 3)
    CI_by_SE_and_Z_val_BCclpr_inter = round(BCest_clpr_inter + c(-1,1) * 1.96 * BCse_clpr_inter, 3)
    CI_by_SE_and_Z_val_BCclpr_inter = paste(CI_by_SE_and_Z_val_BCclpr_inter, sep = ' ', collapse = " , ")
    
  }else{
    NO_BC_WOUT_REP = -101
    BCest_PS <- BCse_PS <- BCest <- BCse <- CI_by_SE_and_Z_val_BC <- BCest_inter <- BCse_inter <- CI_by_SE_and_Z_val_BC_inter <- 
      BCest_clpr <- BCse_clpr <- CI_by_SE_and_Z_val_BCclpr <- 
      BCest_clpr_inter <- BCse_clpr_inter <- CI_by_SE_and_Z_val_BCclpr_inter <- NO_BC_WOUT_REP
  }
  
  #TODO summary of matching estimators
  
  summary_table = data.frame( Replacements = replace, Metric = c("PS", "Mahal", "Mahal PS caliper"),                   
                             Crude_estimators = c(paste0(est_crude_only_ps %>% round(3), " (", se_crude_only_ps %>% round(3), ")"), 
                             paste0(est_crude_maha_wout_cal %>% round(3), " (", se_crude_maha_wout_cal %>% round(3), ")"),
                             paste0(SACE_matching_est %>% round(3), " (", SACE_matching_SE %>% round(3), ")")),
                             Regression = regression_est_se,
                             Regression_interactions = regression_inter_est_se,
                             HL = c(paste0(HL_est_lst$only_ps[1] %>% round(3), " (", HL_est_lst$only_ps[2] %>% round(3), ")"), 
                                    paste0(HL_est_lst$mala_wout_cal[1] %>% round(3), " (", HL_est_lst$mala_wout_cal[2] %>% round(3), ")"),
                                    paste0(HL_est_lst$mala_cal[1] %>% round(3), " (", HL_est_lst$mala_cal[2] %>% round(3), ")"))
                             )
                             
  summary_table = rbind( summary_table, c(replace, "BC PS",  paste0(BCest_PS, " (", BCse_PS, ")"), rep("",3)),
                         c(replace, "BC",  paste0(BCest, " (", BCse, ")"), rep("",3)),
                         c(replace, "BC inter",  paste0(BCest_inter, " (", BCse_inter, ")"), rep("",3)),
                         c(replace, "BC caliper",  paste0(BCest_clpr, " (", BCse_clpr, ")"), rep("",3)),
                         c(replace, "BC caliper inter",  paste0(BCest_clpr_inter, " (", BCse_clpr_inter, ")"), rep("",3)))
                          

  #TODO End of matching estimators!!!
    
  
  #TODO distribution of the x's; before matching and after matching
  # from now, see in the original function
  # descriprive before matching
  
  
  return(list(balance_table=balance_table, 
    balance_only_ps=balance_only_ps$balance_table1, balance_maha_wout_cal=balance_maha_wout_cal$balance_table1, balance_maha_cal_PS=balance_maha_cal_PS$balance_table1, 
    matched_data_mahalPScal=matched_data_mahalPScal, 
    data_pairs_lst = list(data_pairs_only_ps=data_pairs_only_ps, data_pairs_maha_wout_cal=data_pairs_maha_wout_cal, data_pairs_maha_cal_PS=data_pairs_maha_cal_PS),
    summary_table=summary_table, OLS_WLS_reg_lst=OLS_WLS_reg_lst, coeffs_table=coeffs_table))
}



balance_after_matching = function(m_data, match_obj, dt_match, X_sub_cols, metric, replace=FALSE,
                                  smd_se="weighted_simple_AI2011", vertical_table=TRUE, rnd=1){
  disc_var = names(which(apply(subset(m_data, select = X_sub_cols[-1]), 2, function(x) { all(x %in% 0:1)})))
  cont_var = setdiff(X_sub_cols[-1], disc_var)
  
  if(vertical_table == TRUE){
    m_data_trt_cont = subset(dt_match, select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", cont_var))
    m_data_untrt_cont = subset(dt_match, select = paste0("A0_", c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", cont_var)))
    trt = rbind(c( m_data[A==1, .N, by="A"]$N, 0), 
                data.frame(sapply(m_data_trt_cont, function(x) c(mean = mean(x), sd = sd(x))) %>% t))
    untrt = rbind(c( m_data[A==0, .N, by="A"]$N, 0), 
                  data.frame(sapply(m_data_untrt_cont, function(x) c(mean = mean(x), sd = sd(x))) %>% t))
    balance_match = data.frame(Variable = c("N", rownames(trt)[-1]),
                                      untrt = paste0(round(untrt$mean, rnd), ' (', round(untrt$sd, rnd), ')'),
                                      trt = paste0(round(trt$mean, rnd), ' (', round(trt$sd, rnd), ')'))
    balance_match$SMD = round( (untrt$mean - trt$mean) / untrt$sd, (rnd+1) )
    balance_match = filter(balance_match, !Variable %in% c("A", "age2")) 
    balance_match = rbind(balance_match[1,], c("N_match", nrow(dt_match), nrow(dt_match), ""),
          c("N_unq", length(unique(match_obj$index.control)), length(unique(match_obj$index.treated)), ""), balance_match[-1,])
    # discrete
    m_data_trt_disc = subset(dt_match, select = c("A", disc_var))
    m_data_untrt_disc = subset(dt_match, select = paste0("A0_", c("A", disc_var)))
    trt = data.frame(sapply(m_data_trt_disc, function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
    untrt = data.frame(sapply(m_data_untrt_disc, function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
    balance_match_disc = data.frame(Variable = rownames(trt),
                                           untrt = paste0(round(untrt$count, (rnd+1)), ' (', 100 * round(untrt$prop, (rnd+1)), '%)'),
                                           trt = paste0(round(trt$count, (rnd+1)), ' (', 100 * round(trt$prop, (rnd+1)), '%)'))
    balance_match_disc$SMD = round( (untrt$prop - trt$prop) / untrt$sd, (rnd+1) )
    balance_match = rbind(c("Metric", rep(metric, 3)), balance_match, balance_match_disc)
    balance_match = filter(balance_match, !Variable=="A")
  }
  else if(vertical_table == FALSE){
    balance_match = data.frame(
      rbind(apply(subset(dt_match, 
                       select = paste0("A0_", c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1]))), 2, mean),
            apply(subset(dt_match, 
                         select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1])), 2, mean)) ) 
    balance_match = mutate_if(balance_match, is.numeric, round, 2) 
    colnames(balance_match) = substr(colnames(balance_match), 4, 100)
    
    if(smd_se == "not_weighted"){ # se in denominator of smd is sd(X|A=0)
      sd_smd = apply(subset(dt_match, 
                select = paste0("A0_", c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1]))), 2, var)
      sd_smd = sqrt(sd_smd)
    }else if(smd_se == "weighted_simple_AI2011"){
      est_var_x0 = apply(subset(dt_match, 
                        select = paste0("A0_", c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1]))), 2, var)
      est_var_x1 = apply(subset(dt_match, 
                        select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1])), 2, var)
      sd_smd = sqrt(0.5 * (est_var_x0 + est_var_x1))
    }
    
    SMD = (balance_match[1,] - balance_match[2,]) / sd_smd
    balance_match[3,] = mutate_if(SMD, is.numeric, round, 3)
    balance_match$N = c(nrow(filter(m_data, A==0)), nrow(filter(m_data, A==1)), "")
    balance_match$N_match = c(rep(nrow(dt_match),2), "")
    balance_match$N_unq = c(length(unique(match_obj$index.control)),
                            length(unique(match_obj$index.treated)), "")
    balance_match = data.frame(Metric = metric, subset(balance_match, select = c(A, N, N_match, N_unq)), 
                                     subset(balance_match, select = -c(A, N, N_match, N_unq)))
    balance_match$A[3] = ""
    balance_match[c(1:2), c(6:ncol(balance_match))] = format(balance_match[c(1:2), c(6:ncol(balance_match))], nsmall = 2)
  }  
 
  # table1
  weights_trt = data.frame(table(dt_match$id_trt))
  colnames(weights_trt) = c("id", "weight")
  weights_trt$id = as.numeric(as.character(weights_trt$id))
  weights_ctr = data.frame(id = dt_match$id_ctrl , weight = rep(1, length(dt_match$id_ctrl)))
  weights = rbind(weights_trt, weights_ctr)
  data_matched_appear = merge(weights, filter(m_data, id %in% unique(c(dt_match$id_trt, dt_match$id_ctrl))), by= "id", all.x = TRUE)
  sapply(data_matched_appear, function(x) all(x) %in% c(0,1))
  catVars = c("black", "hispanic", "married", "nodegree", "emp74", "emp75")
  #catVars = which(colSums(apply(subset(data_matched_appear, select = -c(intercept, A, S)), 2, function(x) !x %in% c(0,1))) == 0)
  table1 <- CreateTableOne(vars = c("EMest_p_as", "EMest_p_pro", "e_1_as", X_sub_cols[-1]),
                           strata = "A", data = expandRows(data_matched_appear, count="weight"), 
                           smd = TRUE, factorVars = catVars) # factorVars = catVars
  balance_table1 = print(table1, smd = TRUE)
  return(list(balance_match=balance_match, balance_table1=balance_table1))
}


arrange_dataset_after_matching_DATA = function(m_data, match_obj, replace_bool, X_sub_cols){
  ncols  = ncol(subset(m_data[match_obj$index.treated, ], 
                       select = c("id"
                                  #, "p_as"
                                  , "EMest_p_as", "EMest_p_pro", "e_1_as", "Y", "A", "S", "g", X_sub_cols[-1]))) + 1
  dt_match = data.table(subset(m_data[match_obj$index.treated, ], 
                               select = c("id"
                                          #, "p_as"
                                          , "EMest_p_as", "EMest_p_pro", "e_1_as", "Y", "A", "S", "g", X_sub_cols[-1])),
                        match_obj$index.treated, match_obj$index.control,
                        subset(m_data[match_obj$index.control, ], 
                               select = c("id"
                                          #, "p_as"
                                          ,"EMest_p_as", "EMest_p_pro", "e_1_as", "Y", "A", "S", "g", X_sub_cols[-1])))
  colnames(dt_match)[(ncols + 1): (2 * ncols)] = 
    paste0("A0_", colnames(dt_match)[(ncols + 1): (2 * ncols)])
  colnames(dt_match)[c(ncols: (ncols+1))] = c("id_trt", "id_ctrl")
  unique(dt_match$A0_id) %>% length() == nrow(dt_match)
  unique(dt_match$id_trt) %>% length()
  identical(as.numeric(dt_match$id), dt_match$id_trt)
  sum(m_data$id %in% dt_match$id)
  # keep only S = 1
  dt_match_S1 = filter(dt_match, S == 1 & A0_S==1)
  
  # add pairs
  matched_pairs = data.frame(pair = c(1:length(match_obj$index.control)), 
                          ctr = match_obj$index.control, trt = match_obj$index.treated)
  matched_pairs = rbind(data.frame(id=matched_pairs$ctr, pair=matched_pairs$pair), 
                          data.frame(id=matched_pairs$trt, pair=matched_pairs$pair)) %>% arrange(pair) 
  #     
  # if(replace_bool == FALSE){
  #   matched_pairs = data.frame(pair = c(1:length(match_obj$index.control)), 
  #                              ctr = match_obj$index.control, trt = match_obj$index.treated)
  #   matched_pairs = rbind(data.frame(id=matched_pairs$ctr, pair=matched_pairs$pair), 
  #                         data.frame(id=matched_pairs$trt, pair=matched_pairs$pair)) %>% arrange(pair) 
  # }else{
  #   matched_pairs = NULL 
  # }
  return(list(dt_match_S1=dt_match_S1, matched_pairs=matched_pairs))
}

# crude_estimator_inference(match_obj=MATCH_PS_only, dt_match_S1=dt_match_S1_only_ps,
#                           diff_per_pair=diff_per_pair_only_ps, replace_bool=replace)
crude_estimator_inference = function(match_obj, dt_match_S1, diff_per_pair, replace_bool){
  # est
  SACE_matching_est = mean(dt_match_S1$Y) - mean(dt_match_S1$A0_Y)
  # SE 
  if(replace_bool==FALSE){
    SACE_matching_sd = sd(diff_per_pair)
    SACE_matching_SE = SACE_matching_sd / sqrt(nrow(dt_match_S1))
    # when there are rep, the smatched_sete est that match_obj, is only good for the matching with only S1, since otherwise we need to exclude subj's. 
  }else{
    #SACE_matching_SE = ifelse(replace_bool==TRUE, match_obj$se, match_obj$se.standard) 
    SACE_matching_SE = match_obj$se
    # need in matching on PS only, for some reason
    if(is.null(match_obj$se)){
      print(deparse(substitute(MATCH_PS_only)))
      SACE_matching_SE = match_obj$se.standard
    }
  }
  # CI
  CI_by_SE_and_Z_val_crude = round(SACE_matching_est + c(-1,1) * 1.96 * SACE_matching_SE, 3)
  CI_by_SE_and_Z_val_crude = paste(CI_by_SE_and_Z_val_crude, sep = ' ', collapse = " , ")
  return(list(SACE_matching_est=SACE_matching_est, SACE_matching_SE=SACE_matching_SE, 
              CI_by_SE_and_Z_val_crude=CI_by_SE_and_Z_val_crude))
}






