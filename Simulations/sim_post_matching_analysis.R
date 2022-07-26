post_matching_analysis_func = function(m_data, replace, all_measures_matched_lst){
  
  # linear regressions, only for MAHALANOBIS WITH PS CALIPER
  LS_bool = ifelse(replace == F, "OLS", "WLS")
  
  # MATCHING ONLY ONLY on the weights, O11_posterior_ratio ####
  # when matching is on PS only (even with replacement), the SE estimator of Match object is not AI (2006) 
  crude_inference_ps_lst = 
    crude_estimator_inference(match_lst=all_measures_matched_lst$ps_lst, replace_bool=replace)
  matched_data_ps = all_measures_matched_lst$ps_lst$matched_data
  matched_pairs_ps = all_measures_matched_lst$ps_lst$matched_pairs
  
  reg_wout_interactions_inference_ps_lst =
    regression_adjusted_function(m_data=m_data, matched_data=matched_data_ps, matched_pairs=matched_pairs_ps,
                 reg_covariates=X_sub_cols[-1], interactions_bool=FALSE, LS=LS_bool, mu_x_fixed=F, x_as=NULL)
  reg_with_interactions_inference_ps_lst =
    regression_adjusted_function(m_data=m_data, matched_data=matched_data_ps, matched_pairs=matched_pairs_ps,
                 reg_covariates=X_sub_cols[-1], interactions_bool=TRUE, LS=LS_bool, mu_x_fixed=F, x_as=NULL)
  
  wilcoxon_and_HL_ps = wilcoxon_and_HL_func(match_lst=all_measures_matched_lst$ps_lst, boot_HL=FALSE)
  
  # MAHALANOBIS WITHOUT PS CALIPER ####
  crude_inference_mahal_lst = 
    crude_estimator_inference(match_lst=all_measures_matched_lst$mahal_lst, replace_bool=replace)
  matched_data_mahal = all_measures_matched_lst$mahal_lst$matched_data
  matched_pairs_mahal = all_measures_matched_lst$mahal_lst$matched_pairs
  
  reg_wout_interactions_inference_mahal_lst =
    regression_adjusted_function(m_data=m_data, matched_data=matched_data_mahal, matched_pairs=matched_pairs_mahal,
                reg_covariates=X_sub_cols[-1], interactions_bool=FALSE, LS=LS_bool, mu_x_fixed=F, x_as=NULL)
  reg_with_interactions_inference_mahal_lst =
    regression_adjusted_function(m_data=m_data, matched_data=matched_data_mahal, matched_pairs=matched_pairs_mahal,
                reg_covariates=X_sub_cols[-1], interactions_bool=TRUE, LS=LS_bool, mu_x_fixed=F, x_as=NULL)
  
  wilcoxon_and_HL_mahal = wilcoxon_and_HL_func(match_lst=all_measures_matched_lst$mahal_lst, boot_HL=FALSE)
  
  # MAHALANOBIS WITH PS CALIPER ####
  crude_inference_mahal_cal_lst = 
    crude_estimator_inference(match_lst=all_measures_matched_lst$mahal_cal_lst, replace_bool=replace)
  matched_data_mahal_cal = all_measures_matched_lst$mahal_cal_lst$matched_data
  matched_pairs_mahal_cal = all_measures_matched_lst$mahal_cal_lst$matched_pairs
  
  reg_wout_interactions_inference_mahal_cal_lst =
    regression_adjusted_function(m_data=m_data, matched_data=matched_data_mahal_cal, matched_pairs=matched_pairs_mahal_cal,
                 reg_covariates=X_sub_cols[-1], interactions_bool=FALSE, LS=LS_bool, mu_x_fixed=F, x_as=NULL)
  reg_with_interactions_inference_mahal_cal_lst =
    regression_adjusted_function(m_data=m_data, matched_data=matched_data_mahal_cal, matched_pairs=matched_pairs_mahal_cal,
                 reg_covariates=X_sub_cols[-1], interactions_bool=TRUE, LS=LS_bool, mu_x_fixed=F, x_as=NULL)
  
  wilcoxon_and_HL_mahal_cal = wilcoxon_and_HL_func(match_lst=all_measures_matched_lst$mahal_cal_lst, boot_HL=FALSE)
  
  # AI bias-corrected estimator, consider only when replace==TRUE ####
  if(replace == TRUE){ 
    print("bias-corrected estimator")
    BC_inference_ps_lst = BC_estimator_inference(match_obj=all_measures_matched_lst$matchBC_ps_obj, m_data=m_data)
    BC_inference_mahal_lst = BC_estimator_inference(match_obj=all_measures_matched_lst$matchBC_mahal_obj, m_data=m_data)
    BC_inference_mahal_cal_lst = BC_estimator_inference(match_obj=all_measures_matched_lst$matchBC_mahal_cal_obj, m_data=m_data)
  }else{
    NO_BC_WOUT_REP = -101
    BC_inference_ps_lst <- BC_inference_mahal_lst <- BC_inference_mahal_cal_lst <-
      list(diff_per_pair=NO_BC_WOUT_REP, SACE_matching_est=NO_BC_WOUT_REP, SACE_matching_SE=NO_BC_WOUT_REP, 
           CI=NO_BC_WOUT_REP, BC_inference_lst=NO_BC_WOUT_REP)
  }
  
  return(list(
               ps_estimators = list(
                crude_inference_lst=crude_inference_ps_lst
               ,BC_inference_lst=BC_inference_ps_lst
               ,wilcoxon_and_HL=wilcoxon_and_HL_ps
               ,reg_wout_interactions_inference_lst=reg_wout_interactions_inference_ps_lst
               ,reg_with_interactions_inference_lst=reg_with_interactions_inference_ps_lst)
               
               ,mahal_estimators = list(
               crude_inference_lst=crude_inference_mahal_lst
              ,BC_inference_lst=BC_inference_mahal_lst
              ,wilcoxon_and_HL=wilcoxon_and_HL_mahal
              ,reg_wout_interactions_inference_lst=reg_wout_interactions_inference_mahal_lst
              ,reg_with_interactions_inference_lst=reg_with_interactions_inference_mahal_lst)
              
              ,mahal_cal_estimators = list(
               crude_inference_lst=crude_inference_mahal_cal_lst
              ,BC_inference_lst=BC_inference_mahal_cal_lst
              ,wilcoxon_and_HL=wilcoxon_and_HL_mahal_cal
              ,reg_wout_interactions_inference_lst=reg_wout_interactions_inference_mahal_cal_lst
              ,reg_with_interactions_inference_lst=reg_with_interactions_inference_mahal_cal_lst)
              
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
    # need in matching on PS only, for some reason, meaning that for PS matching, there is no SE of AI (2006)
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
    wilc_na = -101
    SACE_matching_SE = wilc_na
  }
  
  CI = c(as.character(as.numeric(round(wilcoxon$conf.int, 3))[1]), 
                          as.character(as.numeric(round(wilcoxon$conf.int, 3))[2]))
  CI = paste(CI, sep = ' ', collapse = " , ")
  
  return(list(SACE_matching_est=SACE_matching_est, SACE_matching_SE=SACE_matching_SE, CI=CI
              ,SACE_matching_pval_HL=SACE_matching_pval_HL))
}

reg_estimator_per_measure = function(lst_one_measure, measure_name, replace){
  # est and SE
  reg_NOint_mat  = c(unlist(lst_one_measure$reg_wout_interactions_inference_lst$estimator_and_se))
  reg_YESint_mat  = c(unlist(lst_one_measure$reg_with_interactions_inference_lst$estimator_and_se))
  
  reg_matching_estimators = c(reg_NOint_mat[1], reg_YESint_mat[1])
  SE_colname = ifelse(replace==T, "WLS_clstr_se", "OLS_se")
  reg_matching_estimators_SE = c(reg_NOint_mat[SE_colname], reg_YESint_mat[SE_colname])
  
  # CI of regression estimators
  reg_matching_estimators_CI = c(
    unlist(lst_one_measure$reg_wout_interactions_inference_lst$CI_LS)
    ,unlist(lst_one_measure$reg_with_interactions_inference_lst$CI_LS))
  
  names(reg_matching_estimators) = names(reg_matching_estimators_SE) = names(reg_matching_estimators_CI) = 
    paste0(measure_name, rep(ifelse(replace==T, "_WLS", "_OLS"), 2), c("", "_int"))
  
  return(list(reg_matching_estimators=reg_matching_estimators, 
              reg_matching_estimators_SE=reg_matching_estimators_SE,
              reg_matching_estimators_CI=reg_matching_estimators_CI))
}

###########################################################################
#post_matching_analysis_lst[[1]][[1]]$reg_wout_interactions_inference_lst
# reg_estimator_per_measure2 = function(lst_one_measure_wout_rpl, lst_one_measure_with_rpl, measure_name){
#   # est and SE
#   OLS_NOint_mat  = c(unlist(lst_one_measure_wout_rpl$reg_wout_interactions_inference_lst$estimator_and_se))
#   OLS_YESint_mat  = c(unlist(lst_one_measure_wout_rpl$reg_with_interactions_inference_lst$estimator_and_se))
#   WLS_NOint_mat  = c(unlist(lst_one_measure_with_rpl$reg_wout_interactions_inference_lst$estimator_and_se))
#   WLS_YESint_mat  = c(unlist(lst_one_measure_with_rpl$reg_with_interactions_inference_lst$estimator_and_se))
#   
#   regression_matching_estimators = c(OLS_NOint_mat[1], OLS_YESint_mat[1], WLS_NOint_mat[1], WLS_YESint_mat[1])
#   regression_matching_estimators_SE = c(OLS_NOint_mat["OLS_se"], OLS_YESint_mat["OLS_se"], WLS_NOint_mat["WLS_clstr_se"], WLS_YESint_mat["WLS_clstr_se"])
#   names(regression_matching_estimators) = names(regression_matching_estimators_SE) = 
#     paste0(measure_name, c("_OLS", "_OLS_int", "_WLS", "_WLS_int"))
#   
#   # CI of regression estimators
#   regression_matching_estimators_CI = c(
#     OLS = unlist(lst_one_measure_wout_rpl$reg_wout_interactions_inference_lst$CI_LS)
#     ,OLS_int = unlist(lst_one_measure_wout_rpl$reg_with_interactions_inference_lst$CI_LS)
#     ,WLS = unlist(lst_one_measure_with_rpl$reg_wout_interactions_inference_lst$CI_LS)
#     ,WLS_int = unlist(lst_one_measure_with_rpl$reg_with_interactions_inference_lst$CI_LS))
#   
#   names(regression_matching_estimators) = names(regression_matching_estimators_SE) =  
#     names(regression_matching_estimators_CI) = paste0(measure_name, c("_OLS", "_OLS_int", "_WLS", "_WLS_int"))
#   
#   if(replace==TRUE){
#     
#   }
#   
#   return(list(regression_matching_estimators=regression_matching_estimators, 
#               regression_matching_estimators_SE=regression_matching_estimators_SE,
#               regression_matching_estimators_CI=regression_matching_estimators_CI))
# }
###########################################################################
