post_matching_analysis_func = function(m_data, replace, all_measures_matched_lst){
  
 # crude estimator ####
  
  # MATCHING ONLY ONLY on the weights, O11_posterior_ratio
  crude_inference_only_ps_lst = 
    crude_estimator_inference(match_lst=all_measures_matched_lst$only_ps_lst, replace_bool=replace)
  
  # MAHALANOBIS WITHOUT PS CALIPER 
  crude_inference_maha_wout_cal_lst = 
    crude_estimator_inference(match_lst=all_measures_matched_lst$maha_wout_cal_lst, replace_bool=replace)
  
  # MAHALANOBIS WITH PS CALIPER
  crude_inference_maha_cal_lst = 
    crude_estimator_inference(match_lst=all_measures_matched_lst$maha_cal_lst, replace_bool=replace)
  dt_match_S1_maha_cal = all_measures_matched_lst$maha_cal_lst$dt_match_S1
  matched_data_maha_cal = all_measures_matched_lst$maha_cal_lst$matched_data
  matched_pairs_maha_cal = all_measures_matched_lst$maha_cal_lst$matched_pairs
  # dt_match_S1_pairs = merge(rbind(dt_match_S1[,1:(ncols-1)], dt_match_S1[,(ncols + 2) : (2 * ncols)], use.names=FALSE) %>% 
  #     subset(select = c("id", "A", "Y")), matched_pairs_maha_cal, by="id") %>% arrange(pair)
  
  
  # AI bias-corrected estimator, consider only when replace==TRUE ####
  #BC_matching_estimation_func()
  if(replace == TRUE){ 
    print("bias-corrected estimator")
    BC_inference_lst = BC_estimator_inference(match_obj=all_measures_matched_lst$matchBC_obj, m_data=m_data)
    BC_cal_inference_lst = BC_estimator_inference(match_obj=all_measures_matched_lst$matchBC_clpr_obj, m_data=m_data)
  }else{
    NO_BC_WOUT_REP = -101
    BC_inference_lst <- BC_cal_inference_lst <- 
      list(diff_per_pair=NO_BC_WOUT_REP, SACE_matching_est=NO_BC_WOUT_REP, SACE_matching_SE=NO_BC_WOUT_REP, 
           CI=NO_BC_WOUT_REP, BC_inference_lst=NO_BC_WOUT_REP)
  }
  
  # non-parametric tests ####
  # wilcoxon_and_HL only for MAHALANOBIS WITH PS CALIPER
  wilcoxon_and_HL = wilcoxon_and_HL_func(match_lst=all_measures_matched_lst$maha_cal_lst, boot_HL=FALSE)
  
  # linear regressions, only for MAHALANOBIS WITH PS CALIPER
  #TODO adjust for replacements and with more than 1 to 1 matching
  LS_bool = ifelse(replace == F, "OLS", "WLS")
  reg_wout_interactions_inference_lst =
    regression_adjusted_function(dt_match_S1=dt_match_S1_maha_cal, 
                                 m_data=m_data, matched_data=matched_data_maha_cal, matched_pairs=matched_pairs_maha_cal,
                                 covariates=X_sub_cols[-1], reg_covariates=X_sub_cols[-1],
                                 interactions_bool=FALSE, LS=LS_bool, mu_x_fixed=F, x_as=x_as)
  reg_with_interactions_inference_lst =
    regression_adjusted_function(dt_match_S1=dt_match_S1_maha_cal, 
                                 m_data=m_data, matched_data=matched_data_maha_cal, matched_pairs=matched_pairs_maha_cal,
                                 covariates =X_sub_cols[-1], reg_covariates=X_sub_cols[-1],
                                 interactions_bool=TRUE, LS=LS_bool, mu_x_fixed=F, x_as=x_as)
  
  
  return(list(crude_inference_only_ps_lst=crude_inference_only_ps_lst
              ,crude_inference_maha_wout_cal_lst=crude_inference_maha_wout_cal_lst
              ,crude_inference_maha_cal_lst=crude_inference_maha_cal_lst
              ,BC_inference_lst=BC_inference_lst
              ,BC_cal_inference_lst=BC_cal_inference_lst
              ,wilcoxon_and_HL=wilcoxon_and_HL
              ,reg_wout_interactions_inference_lst=reg_wout_interactions_inference_lst
              ,reg_with_interactions_inference_lst=reg_with_interactions_inference_lst
              ))
}

crude_estimator_inference = function(match_lst, replace_bool){
  match_obj = match_lst$match_obj; dt_match_S1 = match_lst$dt_match_S1
  
  diff_per_pair = dt_match_S1$Y - dt_match_S1$A0_Y
  # est
  # the SACE est of match_obj and SACE_matching_est are only equal for the matching with only S1, since otherwise we need to exclude subj's. 
  # SACE_matching_est_match_obj = match_obj$est
  SACE_matching_est = mean(dt_match_S1$Y) - mean(dt_match_S1$A0_Y)
  # SE 
  if(replace_bool==FALSE){
    # when there are rep, the se est of match_obj, is only good for the matching with only S1, since otherwise we need to exclude subj's. 
    SACE_matching_sd = sd(diff_per_pair)
    SACE_matching_SE = SACE_matching_sd / sqrt(nrow(dt_match_S1))
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
  CI = round(SACE_matching_est + c(-1,1) * 1.96 * SACE_matching_SE, 3)
  CI = paste(CI, sep = ' ', collapse = " , ")
  return(list(diff_per_pair=diff_per_pair
              ,SACE_matching_est=SACE_matching_est, SACE_matching_SE=SACE_matching_SE, CI=CI))
}

BC_estimator_inference = function(match_obj, m_data){
  diff_per_pair=-101 # there is no meaning for diff_per_pair in BC estimator, just for compatibility with the crude inference 
  # est
  SACE_matching_est = as.numeric(match_obj$est)
  # SE 
  SACE_matching_SE = match_obj$se
  # CI
  CI = round(SACE_matching_est + c(-1,1) * 1.96 * SACE_matching_SE, 3)
  CI = paste(CI, sep = ' ', collapse = " , ")
  
  # untrt_surv is TRUE untreated survivors only if m_data includes only S=1. Otherwise, use sum(1-m_data$A[m_data$S==1])
  BC_ties_multiple_treated = data.frame( untrt_surv = sum(1-m_data$A), 
                matched_untrt = match_obj$index.control %>% length(),
                unq_matched_untrt = unique(match_obj$index.control) %>% length(),
                matched_trt = match_obj$index.treated %>% length(),
                unq_matched_trt = unique(match_obj$index.treated) %>% length() )
  BC_ties_multiple_treated$trt_added_by_ties = BC_ties_multiple_treated$matched_trt - BC_ties_multiple_treated$matched_untrt
  
  return(list(diff_per_pair=diff_per_pair
              ,SACE_matching_est=SACE_matching_est, SACE_matching_SE=SACE_matching_SE, CI=CI
              ,BC_ties_multiple_treated=BC_ties_multiple_treated
              ))
}

wilcoxon_and_HL_func = function(match_lst, boot_HL=FALSE){
  dt_match_S1 = match_lst$dt_match_S1
  diff_per_pair = dt_match_S1$Y - dt_match_S1$A0_Y
  wilcoxon = wilcox.test(diff_per_pair,conf.int=T)
  SACE_matching_est = wilcoxon$estimate
  names(SACE_matching_est) = ""
  SACE_matching_pval_HL = wilcoxon$p.value 
  # bootstrap for HL se
  if(boot_HL==TRUE){
    HL_boot_vec = c()
    for(i in 1:100){
      print(paste0("bootstrap ", i))
      d = dt_match_S1[sample(nrow(dt_match_S1), nrow(dt_match_S1), replace = T),]
      diff_per_pair_boot = d$Y - d$A0_Y
      HL_boot_vec[i] = wilcox.test(diff_per_pair_boot,conf.int=T)$estimate
    }
    SACE_matching_est_HL_bool = mean(HL_boot_vec)
    SACE_matching_SE = sd(HL_boot_vec)
  }else{
    SACE_matching_SE = SACE_matching_pval_HL
  }
  
  CI = c(as.character(as.numeric(round(wilcoxon$conf.int, 3))[1]), 
                          as.character(as.numeric(round(wilcoxon$conf.int, 3))[2]))
  CI = paste(CI, sep = ' ', collapse = " , ")
  
  return(list(SACE_matching_est=SACE_matching_est, SACE_matching_SE=SACE_matching_SE, CI=CI
              ,SACE_matching_pval_HL=SACE_matching_pval_HL))
}

