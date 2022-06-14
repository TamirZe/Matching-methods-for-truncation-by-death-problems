#m_data = data_with_PS
#m_data = data_with_PS[OBS != "O(0,0)"]
#m_data = data_with_PS[S==1]

#m_data = data_list[[3]]
# m_data = data_with_PS[1:3000,]
#m_data = data_with_PS[S==1,]

# caliper is in sd
#replace = T; estimand = "ATC"; change_id = TRUE; mahal_match = 2; M=1; caliper = 0.25


my_matching_func_multiple = function(match_on = NULL, X_sub_cols, m_data, weighting = FALSE,
                                     M=1, replace, estimand = "ATC", mahal_match = 2,
                                     min_PS = 0, min_diff_PS, caliper = 0.05, 
                                     OBS_table, change_id=TRUE, boost_HL=FALSE, mu_x_fixed, x_as,
                                     pass_tables_matched_units = FALSE){
  if(change_id == TRUE){
    print("change id")
    m_data$id = c(1:nrow(m_data))
  }
  print(paste0("replace is ", replace, " nrows is ", nrow(m_data)))
  # mahal_match for Weight = 2 for mahalanobis distance. 1 for inverse of variance
  vec_caliper = c(rep(1000, length(X_sub_cols[-1])), caliper)
  # w_mat = diag(length(X_sub_cols[-1]) + 1) / 
  #   length(X_sub_cols[-1])  * c(apply(subset(m_data, select = X_sub_cols[-1]), 2, var), 1)
  w_mat = diag(length(X_sub_cols[-1]) + 1) / length(X_sub_cols[-1]) 
  w_mat[nrow(w_mat), ncol(w_mat)]  = 0
  
  print(match_on)
  
  # TODO MATCHING ONLY ONLY PS: EMest_p_as
  print("MATCHING ON PS")
  MATCH_PS_only  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                          , X = subset(m_data, select = match_on)
                          ,ties=FALSE
                          ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  )
  only_ps_lst = arrange_dataset_after_matching(match_obj=MATCH_PS_only, m_data, replace_bool = replace, X_sub_cols)
  dt_match_S1_only_ps = only_ps_lst$dt_match_S1; matched_pairs_only_ps = only_ps_lst$matched_pairs
  diff_per_pair_only_ps = dt_match_S1_only_ps$Y - dt_match_S1_only_ps$A0_Y 
  crude_inference_lst_only_ps = 
    crude_estimator_inference(match_obj=MATCH_PS_only, dt_match_S1_only_ps, diff_per_pair_only_ps, replace_bool=replace)
  est_crude_only_ps = crude_inference_lst_only_ps$SACE_matching_est
  se_crude_only_ps = crude_inference_lst_only_ps$SACE_matching_SE
  CI_by_SE_and_Z_val_naive_only_ps = crude_inference_lst_only_ps$CI_by_SE_and_Z_val_crude
  
  # TODO MAHALANOBIS WITHOUT PS CALIPER
  print("MAHALANOBIS WITHOUT PS CALIPER")
  MATCH_MAHA_wout_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                               , X = subset(m_data, select = c(X_sub_cols[-1]))
                               ,ties=FALSE
                               ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  )
  mala_wout_cal_lst = arrange_dataset_after_matching(match_obj=MATCH_MAHA_wout_PS, m_data, replace_bool = replace, X_sub_cols)
  dt_match_S1_maha_wout_cal = mala_wout_cal_lst$dt_match_S1; matched_pairs_mala_wout_cal = mala_wout_cal_lst$matched_pairs
  diff_per_pair_maha_wout_cal = dt_match_S1_maha_wout_cal$Y - dt_match_S1_maha_wout_cal$A0_Y
  crude_inference_lst_mala_wout_cal = 
    crude_estimator_inference(match_obj=MATCH_MAHA_wout_PS, dt_match_S1_maha_wout_cal, diff_per_pair_maha_wout_cal, replace_bool=replace)
  est_crude_maha_wout_cal = crude_inference_lst_mala_wout_cal$SACE_matching_est
  se_crude_maha_wout_cal = crude_inference_lst_mala_wout_cal$SACE_matching_SE
  CI_by_SE_and_Z_val_naive_maha_wout_cal = crude_inference_lst_mala_wout_cal$CI_by_SE_and_Z_val_crude
  
  # TODO MAHALANOBIS WITH PS CALIPER
  print("MAHALANOBIS WITH PS CALIPER")
  ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                         , X = subset(m_data, select = c(X_sub_cols[-1], match_on))
                         ,ties=FALSE
                         ,caliper = vec_caliper ,Weight.matrix = w_mat
                         ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
  )
  ATE_MATCH_PS_lst = arrange_dataset_after_matching(match_obj=ATE_MATCH_PS, m_data, replace_bool = replace, X_sub_cols)
  dt_match_S1 = ATE_MATCH_PS_lst$dt_match_S1; ncols = ncol(dt_match_S1) / 2
  diff_per_pair = dt_match_S1$Y - dt_match_S1$A0_Y
  matched_pairs = ATE_MATCH_PS_lst$matched_pairs
  crude_inference_lst = crude_estimator_inference(match_obj=ATE_MATCH_PS, dt_match_S1, diff_per_pair, replace_bool=replace)
  SACE_matching_est = crude_inference_lst$SACE_matching_est
  SACE_matching_SE = crude_inference_lst$SACE_matching_SE
  CI_by_SE_and_Z_val_naive = crude_inference_lst$CI_by_SE_and_Z_val_crude
  dt_match_S1_pairs = merge(rbind(dt_match_S1[,1:(ncols-1)], dt_match_S1[,(ncols + 2) : (2 * ncols)], use.names=FALSE) %>% subset(select = c("id", "A", "Y")), 
                            matched_pairs, by="id") %>% arrange(pair)
  
  if(replace == TRUE){
    print("START BC")
    m_data_just_inter = m_data$A * subset(m_data, select = X_sub_cols[-1])
    colnames(m_data_just_inter) = paste0("A_", colnames(m_data_just_inter))
    m_data_inter = data.table(m_data, m_data_just_inter)
    
    # TODO AI bias corrected, consider only when replace==TRUE
    matchBC <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = c(X_sub_cols[-1])), 
                     Z = subset(m_data,select = X_sub_cols[-1]), BiasAdjust=TRUE,
                     ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
    )
    BCest = as.numeric(matchBC$est)
    BCse = ifelse(replace==TRUE, matchBC$se, matchBC$se.standard)
    CI_by_SE_and_Z_val_BC = round(BCest + c(-1,1) * 1.96 * BCse, 3)
    CI_by_SE_and_Z_val_BC = paste(CI_by_SE_and_Z_val_BC, sep = ' ', collapse = " , ")
    
    matchBC_inter <- Match(Y=m_data_inter[,Y], Tr=m_data_inter[,A], X = subset(m_data_inter, select = c(X_sub_cols[-1])), 
                           Z = subset(m_data_inter, select = c(X_sub_cols[-1], colnames(m_data_just_inter))), BiasAdjust=TRUE,
                           ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
    )
    BCest_inter = as.numeric(matchBC_inter$est)
    BCse_inter = matchBC_inter$se
    CI_by_SE_and_Z_val_BC_inter = round(BCest_inter + c(-1,1) * 1.96 * BCse_inter, 3)
    CI_by_SE_and_Z_val_BC_inter = paste(CI_by_SE_and_Z_val_BC_inter, sep = ' ', collapse = " , ")
    
    
    matchBC_clpr <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = c(X_sub_cols[-1], match_on)), 
                          Z = subset(m_data,select = X_sub_cols[-1]), BiasAdjust=TRUE,
                          ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                          ,caliper = vec_caliper, Weight.matrix = w_mat
    )
    BCclpr_untrt_surv_matched_untrt_matched_trt = data.frame( untrt_surv = sum(1-m_data$A),
                              matched_untrt = matchBC_clpr$index.control %>% length(),
                              unq_matched_untrt = unique(matchBC_clpr$index.control) %>% length(),
                              matched_trt = matchBC_clpr$index.treated %>% length(),
                              unq_matched_trt = unique(matchBC_clpr$index.treated) %>% length() )
    BCclpr_untrt_surv_matched_untrt_matched_trt$trt_added_by_ties = 
      BCclpr_untrt_surv_matched_untrt_matched_trt$matched_trt - BCclpr_untrt_surv_matched_untrt_matched_trt$matched_untrt 
    BCest_clpr = as.numeric(matchBC_clpr$est)
    BCse_clpr = matchBC_clpr$se
    CI_by_SE_and_Z_val_BCclpr = round(BCest_clpr + c(-1,1) * 1.96 * BCse_clpr, 3)
    CI_by_SE_and_Z_val_BCclpr = paste(CI_by_SE_and_Z_val_BCclpr, sep = ' ', collapse = " , ")
    #summary.Match(matchBC, full=TRUE)
    
    matchBC_clpr_inter <- Match(Y=m_data_inter[,Y], Tr=m_data_inter[,A], X = subset(m_data_inter, select = c(X_sub_cols[-1], match_on)), 
                                Z = subset(m_data_inter, select = c(X_sub_cols[-1], colnames(m_data_just_inter))), BiasAdjust=TRUE,
                                ties=TRUE ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                                ,caliper = vec_caliper, Weight.matrix = w_mat
    )
    BCest_clpr_inter = as.numeric(matchBC_clpr_inter$est)
    BCse_clpr_inter = matchBC_clpr_inter$se
    CI_by_SE_and_Z_val_BCclpr_inter = round(BCest_clpr_inter + c(-1,1) * 1.96 * BCse_clpr_inter, 3)
    CI_by_SE_and_Z_val_BCclpr_inter = paste(CI_by_SE_and_Z_val_BCclpr_inter, sep = ' ', collapse = " , ")
    
    print("END BC")
    
  }else{
    NO_BC_WOUT_REP = -101
    #BCest <- BCse <- CI_by_SE_and_Z_val_BC <- BCest_clpr <- BCse_clpr <- CI_by_SE_and_Z_val_BCclpr <- NO_BC_WOUT_REP
    BCest <- BCse <- CI_by_SE_and_Z_val_BC <- BCest_inter <- BCse_inter <- CI_by_SE_and_Z_val_BC_inter <- 
      BCest_clpr <- BCse_clpr <- CI_by_SE_and_Z_val_BCclpr <- BCest_clpr_inter <- BCse_clpr_inter <- 
      CI_by_SE_and_Z_val_BCclpr_inter <- BCclpr_untrt_surv_matched_untrt_matched_trt <- NO_BC_WOUT_REP
    
  }
  
  
  wilcoxon = wilcox.test(diff_per_pair,conf.int=T)
  SACE_matching_est_HL = wilcoxon$estimate
  SACE_matching_pval_HL = wilcoxon$p.value 
  # bootstrap for HL se
  if(boost_HL==TRUE){
    HL_boost_vec = c()
    for(i in 1:100){
      print(paste0("bootstrap ", i))
      d = dt_match_S1[sample(nrow(dt_match_S1), nrow(dt_match_S1), replace = T),]
      diff_per_pair_boost = d$Y - d$A0_Y
      HL_boost_vec[i] = wilcox.test(diff_per_pair_boost,conf.int=T)$estimate
    }
    SACE_matching_est_HL_bool = mean(HL_boost_vec)
    SACE_matching_se_HL = sd(HL_boost_vec)
  }else{
    SACE_matching_se_HL = SACE_matching_pval_HL
  }
  
  SACE_matching_CI_HL = c(as.character(as.numeric(round(wilcoxon$conf.int, 3))[1]), 
                          as.character(as.numeric(round(wilcoxon$conf.int, 3))[2]))
  SACE_matching_CI_HL = paste(SACE_matching_CI_HL, sep = ' ', collapse = " , ")
  
  #####################################################################
  # TODO adjust for replacements and with more than 1 to 1 matching
  # TODO change: if replace==TRUE: WLS, if replace==FALSE: OLS
  
  # TODO WLS
  WLS_NOinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(rep_bool = replace, dt_match_S1=dt_match_S1, m_data=m_data, matched_pairs=matched_pairs,
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = FALSE, LS="WLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  WLS_YESinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(rep_bool = replace, dt_match_S1=dt_match_S1, m_data=m_data, matched_pairs=matched_pairs,
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = TRUE, LS="WLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  
  # TODO OLS
  OLS_NOinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(rep_bool = replace, dt_match_S1=dt_match_S1, m_data=m_data, matched_pairs=matched_pairs, 
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = FALSE, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  OLS_YESinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(rep_bool = replace, dt_match_S1=dt_match_S1, m_data=m_data, matched_pairs=matched_pairs, # matched_pairs=NULL
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = TRUE, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  
  print("finish lin reg")
  
  #TODO distribution of the x's; before matching and after matching
  # descriprive before matching
  initial_data_x = subset(m_data, select = c("id", "A", "S", "g", X_sub_cols[-1]))
  initial_data_x_as = filter(initial_data_x, g=="as")
  initial_data_x_as_A0 = filter(initial_data_x, A==0, S==1) #only as here
  initial_data_x_as_A1 = filter(initial_data_x, A==1, S==1, g=="as")
  mean_as = apply(subset(initial_data_x_as, select = X_sub_cols[-1]), 2, mean)
  mean_A0_S1 = apply(subset(initial_data_x_as_A0, select = X_sub_cols[-1]), 2, mean)
  mean_A1_S1_as = apply(subset(initial_data_x_as_A1, select = X_sub_cols[-1]), 2, mean)
  
  # descriprive after matching
  dt_match_S1_A0 = subset(dt_match_S1, select = grep("A0|ctr", colnames(dt_match_S1)))
  mean_match_A0 = apply(subset(dt_match_S1_A0,select = paste0(rep("A0_"),X_sub_cols[-1])),2,mean)
  dt_match_S1_A1 = subset(dt_match_S1, select = -grep("A0|ctr", colnames(dt_match_S1)))
  mean_match_A1 = apply(subset(dt_match_S1_A1, select = X_sub_cols[-1]), 2, mean)
  approx_mean_x = (mean_match_A1 + mean_match_A0) / 2
  
  means_by_subset = 
    rbind(mean_as, mean_A0_S1, mean_A1_S1_as, mean_match_A0, mean_match_A1, approx_mean_x)
  
  # calculating the amount of as the matching process excluded
  OBS_table * param_n
  as_A0_matched = length(which(dt_match_S1$A0_g == "as"))
  # OBS_table[1,2] are the (A=0, S=1), thereare only as in this cell
  as_A0_unmatched = (OBS_table[1,2] * param_n) - as_A0_matched 
  nrow(filter(m_data, g=="as" & A==0)) - as_A0_matched # another options:
  setdiff(filter(m_data, g=="as" & A==0)$id, ATE_MATCH_PS$index.control) %>% length() # another options
  as_A1_matched = length(which(dt_match_S1$g == "as"))
  pro_A1_matched = length(which(dt_match_S1$g == "pro"))
  
 
  return(list(SACE_matching_est_SE_ps = c(est_crude_only_ps, se_crude_only_ps)
              ,SACE_matching_est_SE_maha = c(est_crude_maha_wout_cal,se_crude_maha_wout_cal)
              ,SACE_matching_est_SE_naive = c(SACE_matching_est, SACE_matching_SE)
              ,SACE_matching_est_HL = c(SACE_matching_est_HL, SACE_matching_se_HL)
              ,SACE_matching_est_SE_BC = c(BCest, BCse)
              ,SACE_matching_est_SE_BCclpr = c(BCest_clpr, BCse_clpr)
              ,SACE_matching_est_SE_BC_inter = c(BCest_inter, BCse_inter)
              ,SACE_matching_est_SE_BCclpr_inter = c(BCest_clpr_inter, BCse_clpr_inter)
              ,CI_crude_HL_BC = c(CI_by_SE_and_Z_val_naive_only_ps,CI_by_SE_and_Z_val_naive_maha_wout_cal,
                                  CI_by_SE_and_Z_val_naive, SACE_matching_CI_HL, CI_by_SE_and_Z_val_BC, CI_by_SE_and_Z_val_BCclpr,
                                  CI_by_SE_and_Z_val_BC_inter, CI_by_SE_and_Z_val_BCclpr_inter)
              ,WLS_NOinteractions_reg_adj_estimators_and_se = WLS_NOinteractions_reg_adj_estimators_and_se
              ,WLS_YESinteractions_reg_adj_estimators_and_se = WLS_YESinteractions_reg_adj_estimators_and_se
              ,OLS_NOinteractions_reg_adj_estimators_and_se = OLS_NOinteractions_reg_adj_estimators_and_se
              ,OLS_YESinteractions_reg_adj_estimators_and_se = OLS_YESinteractions_reg_adj_estimators_and_se
              ,means_by_subset = means_by_subset
              ,BCclpr_untrt_surv_matched_untrt_matched_trt=BCclpr_untrt_surv_matched_untrt_matched_trt
  ))
}



arrange_dataset_after_matching = function(match_obj, m_data, replace_bool, X_sub_cols){
  ncols  = ncol(subset(m_data[match_obj$index.treated, ], 
                       select = c("id"
                                  #, "p_as"
                                  , "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1]))) + 1
  dt_match = data.table(subset(m_data[match_obj$index.treated, ], 
                               select = c("id"
                                          #, "p_as"
                                          , "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1])),
                        match_obj$index.treated, match_obj$index.control,
                        subset(m_data[match_obj$index.control, ], 
                               select = c("id"
                                          #, "p_as"
                                          ,"EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1])))
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
  
  # return matched_pairs for OLS: when replace_bool == FALSE
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

crude_estimator_inference = function(match_obj, dt_match_S1, diff_per_pair, replace_bool){
  # est
  SACE_matching_est = mean(dt_match_S1$Y) - mean(dt_match_S1$A0_Y)
  # SE 
  if(replace_bool==FALSE){
    SACE_matching_sd = sd(diff_per_pair)
    SACE_matching_SE = SACE_matching_sd / sqrt(nrow(dt_match_S1))
    # when there are rep, the se est that match_obj, is only good for the matching with only S1, since otherwise we need to exclude subj's. 
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

tables_repeated = function(table_subjects){
  #table_treated_subjects = data.table(table(ATE_MATCH$index.treated))
  table_treated_subjects = data.table(apply(table_subjects, 2, as.numeric))
  colnames(table_treated_subjects)[1] = "id_trt"
  repeated_treated = table(table_treated_subjects$N)
  return(list(table_treated_subjects = table_treated_subjects, 
              repeated_treated = repeated_treated))
}

check_balance_function = function(data, cols){
  data = data.frame(data)
  diff_w = apply(data[ , cols], 2, mean) - 
    apply(data[ , paste0("A0", "_", cols)], 2, mean)
  std_diff_w = diff_w / 
    apply(data[ , paste0("A0", "_", cols)], 2, sd)  
  return(list(diff_w, std_diff_w))
}