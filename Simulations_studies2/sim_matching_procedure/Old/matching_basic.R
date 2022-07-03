#m_data = data_with_PS
# m_data = data_with_PS[OBS != "O(0,0)"]
#m_data = data_with_PS[S==1]
#m_data = data_list[[3]]
# m_data = data_with_PS[1:3000,]
# m_data = data_with_PS[1:3000,][S==1]
# 
# 
# # caliper is in sd
#replace = TRUE; estimand = "ATC"; change_id = TRUE; mahal_match = 2; M=1; caliper = 0.05
# all.equal(dt_match_min_ps_w_scheme, dt_match_S1, check.attributes = FALSE)
# match_on = "O11_posterior_ratio"

# my_matching_func_basic(match_on = match_on, X_sub_cols, data_list[[l]],
#                        weighting = FALSE, M=1, replace = replace_vec[j], estimand = "ATC", mahal_match = 2,
#                        min_PS = min_PS, min_diff_PS = min_diff_PS,
#                        caliper = caliper, OBS_table, change_id = TRUE)

my_matching_func_basic = function(match_on = NULL, X_sub_cols, m_data, weighting = FALSE,
                M=1, replace, estimand = "ATC", mahal_match = 2,
                min_PS = 0, min_diff_PS, caliper = 0.05, 
                OBS_table, change_id=TRUE, boost_HL=FALSE, mu_x_fixed, x_as){
  if(change_id == TRUE){
    print("change id")
    m_data$id = c(1:nrow(m_data))
  }
  print(paste0("replace is ", replace, " nrows is ", nrow(m_data)))
  #X_sub_cols = paste0("X", c(1:(dim_x)))
  # mahal_match for Weight = 2 for mahalanobis distance. 1 for inverse of variance
  # TODO find a way to use caliper on the PS in the matching function
  vec_caliper = c(rep(1000, length(X_sub_cols[-1])), caliper)
  #w_mat = diag(length(X_sub_cols[-1]) + 1) / length(X_sub_cols[-1]) 
  w_mat = diag(length(X_sub_cols[-1]) + 1) / 
    length(X_sub_cols[-1])  * c(apply(subset(m_data, select = X_sub_cols[-1]), 2, var), 1)
  w_mat[nrow(w_mat), ncol(w_mat)]  = 0
  
  if(is.null(match_on) == TRUE){
    # TODO matching with caliper; if(caliper > 0){
    #set.seed(11) # , "EMest_p_as"
    print("match wout match_on")
    ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                           #, X=m_data[,"est_p_as"]
                           , X = subset(m_data, select = c(X_sub_cols[-1], "EMest_p_as"))
                           #,ties=FALSE
                           ,caliper = vec_caliper ,Weight.matrix = w_mat
                           ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
    )
  }else{
    print(match_on)
    ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                           , X = subset(m_data, select = match_on)
                           #, X=m_data[,"EMest_p_as"]         
                           #, X = subset(m_data, select = c(X_sub_cols[-1], "EMest_p_as"))
                           ,ties=FALSE
                           #,caliper = vec_caliper
                           ,M=M, replace = replace, estimand = estimand, Weight = mahal_match
                           #,Weight.matrix = w_mat
    )
  }
  
  # TODO AI bias corrected, consider only when replace==TRUE
  matchBC <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = c(X_sub_cols[-1])), 
                   Z = subset(m_data,select = X_sub_cols[-1]), BiasAdjust=TRUE,
                   ties=TRUE ,M=M, replace = replace, estimand = "ATC", Weight = mahal_match
    )
  BCest = as.numeric(matchBC$est); BCse = ifelse(replace==TRUE, matchBC$se, matchBC$se.standard)
  
  matchBC_clpr <- Match(Y=m_data[,Y], Tr=m_data[,A], X = subset(m_data, select = c(X_sub_cols[-1], "EMest_p_as")), 
                   Z = subset(m_data,select = X_sub_cols[-1]), BiasAdjust=TRUE,
                   ties=TRUE ,M=M, replace = replace, estimand = "ATC", Weight = mahal_match
                   ,caliper = vec_caliper, Weight.matrix = w_mat
  )
  BCest_clpr = as.numeric(matchBC_clpr$est); BCse_clpr = ifelse(replace==TRUE, matchBC_clpr$se, matchBC_clpr$se.standard)
  
  #summary.Match(matchBC, full=TRUE)
  #summary(ATE_MATCH_PS)
  
    
  #ATE_MATCH_PS$est; ATE_MATCH_PS$se; ATE_MATCH_PS$nobs; ATE_MATCH_PS$index.dropped
  ncols  = ncol(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                       select = c("id", "p_as", "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1]))) + 1
  dt_match = data.table(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                               select = c("id", "p_as", "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1])),
                        ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control,
                        subset(m_data[ATE_MATCH_PS$index.control, ], 
                               select = c("id", "p_as", "EMest_p_as", "Y", "A", "S", "g", X_sub_cols[-1])))
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
  if(replace == FALSE){
    matched_pairs = data.frame(pair = c(1:length(ATE_MATCH_PS$index.control)), 
                               ctr = ATE_MATCH_PS$index.control, trt = ATE_MATCH_PS$index.treated)
    matched_pairs = rbind(data.frame(id=matched_pairs$ctr, pair=matched_pairs$pair), 
                          data.frame(id=matched_pairs$trt, pair=matched_pairs$pair)) %>% arrange(pair) 
  }else{
    matched_pairs = NULL 
  }
  
  
  #TODO distribution of the x's; before matching and after matching
  # descriprive before matching
  initial_data_x = subset(m_data, select = c("id", "A", "S", "g", X_sub_cols[-1]))
  initial_data_x_as = filter(initial_data_x, g=="as")
  initial_data_x_as_A0 = filter(initial_data_x, A==0, S==1) #only as here
  initial_data_x_as_A1 = filter(initial_data_x, A==1, S==1, g=="as")
  mean_as = apply(subset(initial_data_x_as, select = X_sub_cols[-1]), 2, mean)
  mean_A0_S1 = apply(subset(initial_data_x_as_A0, select = X_sub_cols[-1]), 2, mean)
  mean_A1_S1_as = apply(subset(initial_data_x_as_A1, select = X_sub_cols[-1]), 2, mean)
  
  # histograms
  #apply(subset(initial_data_x_as, select = X_sub_cols[-1]), 2, function(i) hist(x, main = ))
  
  # lapply(c(1:length(X_sub_cols[-1])), function(i){
  #   hist(subset(initial_data_x_as_A1, select = X_sub_cols[-1])[,i], xlab = "covariate",
  #        main = paste0("initial_data_as_A1 " , X_sub_cols[-1][i]))
  # })
  
  # descriprive after matching
  dt_match_S1_A0 = subset(dt_match_S1, select = grep("A0|ctr", colnames(dt_match_S1)))
  mean_match_A0 = apply(subset(dt_match_S1_A0,select = paste0(rep("A0_"),X_sub_cols[-1])),2,mean)
  dt_match_S1_A1 = subset(dt_match_S1, select = -grep("A0|ctr", colnames(dt_match_S1)))
  mean_match_A1 = apply(subset(dt_match_S1_A1, select = X_sub_cols[-1]), 2, mean)
  approx_mean_x = (mean_match_A1 + mean_match_A0) / 2
  
  means_by_subset = 
    rbind(mean_as, mean_A0_S1, mean_A1_S1_as, mean_match_A0, mean_match_A1, approx_mean_x)
  
  ###### estimation
  # TODO diff of means
  SACE_matching_est = mean(dt_match_S1$Y) - mean(dt_match_S1$A0_Y)
  #ATE_MATCH_PS$est; ATE_MATCH_PS$est.noadj
  # TODO mean of diffs 
  diff_per_pair = dt_match_S1$Y - dt_match_S1$A0_Y
  #SACE_matching_est2 = mean(diff_per_pair)
  
  # TODO est se wout replacements
  if(replace==FALSE){ 
    SACE_matching_sd = sd(diff_per_pair)
    # SE 
    SACE_matching_SE = SACE_matching_sd / sqrt(nrow(dt_match_S1))
    #print("AI SE instead of naive SE est")
    #SACE_matching_SE = ATE_MATCH_PS$se.standard
    
    # TODO est se with replacements
  }else{
    weights_trt = data.frame(table(dt_match_S1$id_trt))
    colnames(weights_trt) = c("id", "weight")
    id_w_and_Y_trt = merge(weights_trt, subset(dt_match_S1, select = c(id, Y)), by="id") 
    id_w_and_Y_trt = id_w_and_Y_trt[!duplicated(id_w_and_Y_trt$id), ]
    
    func_est_var_withweights_matching_w_rep = function(dt_match_S1, id_w_and_Y_trt){
      M = nrow(dt_match_S1); N_1m = length(unique(dt_match_S1$id_trt))
      print(nrow(id_w_and_Y_trt) == N_1m)
      est_var_Y0 = var(dt_match_S1$A0_Y) # (sd(dt_match_S1$A0_Y))^2
      est_var_Y1 = var(id_w_and_Y_trt$Y) # (sd(id_w_and_Y_trt$Y))^2
      est_var_mean_Y1_with_weights = est_var_Y1 * sum((id_w_and_Y_trt$weight)^2)
      est_var_after_match_w_rep = ( (1/M^2) * est_var_mean_Y1_with_weights ) + ( (1/M) * est_var_Y0 )
      return(est_var_after_match_w_rep)
    }
    #SACE_matching_SE = sqrt( func_est_var_withweights_matching_w_rep(dt_match_S1, id_w_and_Y_trt) ) # SACE_matching_SE_my
    SACE_matching_SE = ATE_MATCH_PS$se
    
    # MISC
    # M = nrow(dt_match_S1)
    # est_var_Y0 = var(dt_match_S1$A0_Y)
    # #wtd.mean(id_w_and_Y_trt$Y, id_w_and_Y_trt$weight); mean(dt_match_S1$Y)
    # w_var <- wtd.var(id_w_and_Y_trt$Y, id_w_and_Y_trt$weight) /  M
    #SACE_matching_SE = sqrt(w_var + ( (1/M) * est_var_Y0 ))
    
  }
  
  
  
  # nrow(dt_match_S1) - 1 ; # https://www.youtube.com/watch?v=zD3VIBkwc-0
  t_val = qt(0.975, nrow(dt_match_S1) - 2, lower.tail = TRUE, log.p = FALSE)
  CI_by_SE_and_Z_val_naive = round(SACE_matching_est + c(-1,1) * 1.96 * SACE_matching_SE, 3)
  CI_by_SE_and_Z_val_naive = paste(CI_by_SE_and_Z_val_naive, sep = ' ', collapse = " , ")
  t_test = t.test(dt_match_S1$Y, dt_match_S1$A0_Y, 
                  var.equal = FALSE, paired = TRUE
                  #,alternative = c("two.sided")
  ) 
  CI_by_SE_and_Z_val_BC = round(BCest + c(-1,1) * 1.96 * BCse, 3)
  CI_by_SE_and_Z_val_BC = paste(CI_by_SE_and_Z_val_BC, sep = ' ', collapse = " , ")
  CI_by_SE_and_Z_val_BCclpr = round(BCest_clpr + c(-1,1) * 1.96 * BCse_clpr, 3)
  CI_by_SE_and_Z_val_BCclpr = paste(CI_by_SE_and_Z_val_BCclpr, sep = ' ', collapse = " , ")
  
  # TODO HL
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
  # TODO Regression adjusted matching_from_real_data on the matched set
  # TODO regression_adjusted_function(dt_match_S1, reg_covariates = X_sub_cols[-1])
  # TODO change: if replace==TRUE: WLS, if replace==TRUE: OLS
  
  # TODO WLS
  WLS_NOinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(dt_match_S1, m_data, matched_pairs=NULL,
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = FALSE, LS="WLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  WLS_YESinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(dt_match_S1, m_data, matched_pairs=NULL,
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = TRUE, LS="WLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  
  # TODO OLS
  OLS_NOinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(dt_match_S1, m_data, matched_pairs=matched_pairs, 
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = FALSE, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  OLS_YESinteractions_reg_adj_estimators_and_se =
    regression_adjusted_function(dt_match_S1, m_data, matched_pairs=NULL,
                                 covariates = X_sub_cols[-1], reg_covariates = X_sub_cols[-1],
                                 interactions_bool = TRUE, LS="OLS", mu_x_fixed=mu_x_fixed, x_as=x_as)
  
  
  ######## calculating the amount of as the matching process excluded
  OBS_table * param_n
  as_A0_matched = length(which(dt_match_S1$A0_g == "as"))
  # OBS_table[1,2] are the (A=0, S=1), thereare only as in this cell
  as_A0_unmatched = (OBS_table[1,2] * param_n) - as_A0_matched 
  as_A1_matched = length(which(dt_match_S1$g == "as"))
  pro_A1_matched = length(which(dt_match_S1$g == "pro"))
  
  # check the compisition in A1_S1: as and protected
  as_in_A1_S1 = m_data[g =="as" & A == 1,]
  pro_in_A1_S1 = m_data[g =="pro" & A == 1,]
  # senity check
  nrow(as_in_A1_S1) + nrow(pro_in_A1_S1) == OBS_table[2,2] * param_n
  as_A1_unmatched = nrow(as_in_A1_S1) - as_A1_matched
  included_excluded_in_matching = 
    data.frame(as_A0_matched, as_A0_unmatched, as_A1_matched, as_A1_unmatched, pro_A1_matched)
  colnames(included_excluded_in_matching) = 
    paste0("rep", substr(replace, 1,1), "_", colnames(included_excluded_in_matching))
  
  # checking covariates balance between treated and untreated
  dt_match_min_ps_A0 = dt_match_S1[ , c("A0_A", paste0("A0", "_", X_sub_cols[-1]))]
  dt_match_min_ps_A1 = dt_match_S1[ , c("A", X_sub_cols[-1])]
  
  
  # TODO repeated summary
  # TODO %%%%%%%%% NEXT LINE IS REALLY VERY IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  m_data$id_n = c(1:nrow(m_data))
  m_data = data.frame(id_n = m_data$id_n, subset(m_data, select = -id_n))
  befS1_table_treated_subjects = data.table(table(ATE_MATCH_PS$index.treated))
  befS1_list_table_treated = tables_repeated(befS1_table_treated_subjects)
  befS1_table_treated_subjects = befS1_list_table_treated[[1]]
  befS1_list_table_treated_repeated = befS1_list_table_treated[[2]]
  colnames(befS1_table_treated_subjects)[1] = "id_n"
  befS1_trt_matched_unq = merge(befS1_table_treated_subjects, subset(m_data, select = c(id_n, g)),
                                all.x = TRUE, all.y = FALSE, by = "id_n")
  befS1_trt_matched_total_by_g = ddply(befS1_trt_matched_unq, .(g), summarize, total_re = sum(N))
  # matched_treated_S1 = subset(dt_match_S1, select = c(id_trt, g))
  # matched_treated_S1 = arrange(matched_treated_S1[!duplicated(matched_treated_S1),], by=id_trt)
  # matched_treated_S1_repeats = merge(matched_treated_S1, table_treated_subjects,
  #                                   all.x = TRUE, all.y = FALSE, by = "id_trt")
  
  
  table_treated_subjects = data.table(table(dt_match_S1$id_trt))
  list_table_treated = tables_repeated(table_treated_subjects)
  table_treated_subjects = list_table_treated[[1]]
  repeated_treated = list_table_treated[[2]]
  
  matched_treated_S1 = subset(dt_match_S1, select = c(id_trt, g))
  matched_treated_S1 = arrange(matched_treated_S1[!duplicated(matched_treated_S1),], by=id_trt)
  matched_treated_S1_repeats = merge(matched_treated_S1, table_treated_subjects, 
                                     all.x = TRUE, all.y = FALSE, by = "id_trt")
  identical(matched_treated_S1_repeats[,-3], matched_treated_S1)
  
  S1_matched_distinguish_by_g =  
    ddply(matched_treated_S1_repeats, .(g), summarize, amount_of_subjects = length(id_trt))
  #m = m_data[unique(ATE_MATCH_PS$index.treated), c("id_n", "id", "g")]
  #S1_matched_distinguish_by_g_2 = ddply(m, .(g), summarize, total_re = length(id_n))
  
  #print("before repeated3")
  S1_matched_as = filter(matched_treated_S1_repeats, g=="as")
  list_table_treated_as = tables_repeated(subset(S1_matched_as, select = -g))
  repeated_treated_as = data.frame(list_table_treated_as[[2]])
  
  S1_matched_pro = filter(matched_treated_S1_repeats, g=="pro") 
  list_table_treated_pro = tables_repeated(subset(S1_matched_pro, select = -g))
  repeated_treated_pro = data.frame(list_table_treated_pro[[2]])
  
  as = repeated_treated_as; pro = repeated_treated_pro
  as = data.frame(Var1 = c(1:14), Freq = 0)
  REPS_as = as.numeric(as.character(repeated_treated_as$Var1))
  as$Freq[REPS_as[REPS_as<=14]] = repeated_treated_as$Freq
  as$Var1 = paste0("as_", as$Var1)
  rownames(as) = as$Var1
  
  pro = data.frame(Var1 = c(1:14), Freq = 0)
  REPS_pro = as.numeric(as.character(repeated_treated_pro$Var1))
  pro$Freq[REPS_pro[REPS_pro<=14]] = repeated_treated_pro$Freq
  pro$Var1 = paste0("pro_", pro$Var1)
  rownames(pro) = pro$Var1
  
  repeated_as_and_pro = data.frame(t(as), t(pro))
  repeated_as_and_pro = repeated_as_and_pro[-1,]
  
  #histogram 
  # hist(rep(as.numeric(repeated_treated_as$Var1), times = repeated_treated_as$Freq), 
  #      col='skyblue', border=F, xlab = "X2", main = colnames(d1)[i])
  # hist(rep(as.numeric(repeated_treated_pro$Var1), times = repeated_treated_pro$Freq),
  #      add=T, col=scales::alpha('green',.5), border=F, breaks = 10)
  #  legend('topright',c('as','pro'),
  #         fill = c('skyblue', 'green'), bty = 'n',
  #         border = NA)
  
  # }
  
  
  # checking covariates balance between treated and untreated 
  # without manu caliper matching
  balance_as_to_as = check_balance_function(
    filter(dt_match_S1, g == "as" & A0_g == "as"), X_sub_cols[-1])
  balance_as_to_pro = check_balance_function(
    filter(dt_match_S1, g == "pro" & A0_g == "as"), X_sub_cols[-1])
  std_diff_as_to_as = balance_as_to_as[[2]]
  std_diff_as_to_pro = balance_as_to_pro[[2]]
  
  std_diff_2_cols = data.frame(std_diff_as_to_pro, std_diff_as_to_as)
  diff_distance_aspr_asas = 
    mean(abs(std_diff_as_to_pro) - abs(std_diff_as_to_as))
  
  
  # TODO number of matchd units:
  # summary immediately after matching
  # m_data= only S1
  n_trt = nrow(filter(m_data, A==1)); n_ctr = nrow(filter(m_data, A==0))
  unq_trt_in_match = length(unique(ATE_MATCH_PS$index.treated))
  unq_ctr_in_match = length(unique(ATE_MATCH_PS$index.control))
  
  m_data_trt = m_data[unique(ATE_MATCH_PS$index.treated), ]
  m_data_ctr =m_data[unique(ATE_MATCH_PS$index.control), ]
  #m_data_trt = m_data[id %in% unique(ATE_MATCH_PS$index.treated) , ]
  #m_data_ctr = m_data[id %in% unique(ATE_MATCH_PS$index.control) , ]
  
  trt_tab_after_match_bef_S1 = c(table(m_data_trt$g))
  ctr_tab_after_match_bef_S1 = c(table(m_data_ctr$g))
  
  
  # after excluding pairs with nonsurvivor
  unq_after_S1 = filter(m_data, id %in%
                          c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) )
  table(unq_after_S1$A, unq_after_S1$g)
  trt_tab_after_S1 = t(table(filter(unq_after_S1, A==1)$g))
  ctr_tab_after_S1 = t(table(filter(unq_after_S1, A==0)$g))
  
  # including replacements, the total weights of each stratum
  # in the histogram (repeated as_pro) its the amount of re (1,2...) times how many (per stratum)
  # replaced this amounts (sum over all the amount of re (1,2...))
  d = data.table(befS1_trt_matched_total_by_g)
  d[g=="as",total_re]
  trt_total_as_Be = d[g=="as",total_re]; trt_total_pro_Be = d[g=="pro",total_re]
  trt_tab_after_S1 = c(as=trt_tab_after_S1[1], ns=0, pro=trt_tab_after_S1[2])
  ctr_tab_after_S1 = c(as=ctr_tab_after_S1[1], ns=0, pro=0)
  # if there is no ns in mathced treated before excluding non-survivors, it implies we in S1 MATCH
  if("ns" %in% as.character(d$g)){
    trt_total_ns_Be = d[g=="ns",total_re]
  }else{
    trt_total_ns_Be=0
    trt_tab_after_match_bef_S1 = c(as=trt_tab_after_match_bef_S1[1], ns=0, pro=trt_tab_after_match_bef_S1[2])
    ctr_tab_after_match_bef_S1 = c(as=ctr_tab_after_match_bef_S1[1], ns=0, pro=0)
    # pro=ifelse(is.na(ctr_tab_after_match_bef_S1[2]), 0, ctr_tab_after_match_bef_S1[2])
  }
  # FOR WOUT O(0,0)
  if(names(ctr_tab_after_match_bef_S1)[1] == 'as'){
    ctr_tab_after_match_bef_S1 = c(ctr_tab_after_match_bef_S1[1], ns=0, pro=0)
  }
  
  trt_total_as_Af = sum(as$Freq*c(1:length(as$Freq)))
  trt_total_pro_Af = sum(pro$Freq*c(1:length(pro$Freq)))
  
  tables_matched_units = data.frame(cbind(t(trt_tab_after_match_bef_S1), t(ctr_tab_after_match_bef_S1),
                                          t(trt_tab_after_S1), t(ctr_tab_after_S1), 
            trt_total_as_Be, trt_total_ns_Be, trt_total_pro_Be, trt_total_as_Af, trt_total_pro_Af))
  
  colnames(tables_matched_units) = c(paste0(rep(c("trt", "ctr"), each = 3), "_",
                        rep(c("as", "ns", "pro"), times=4), rep(c("BeS1", "AfS1"), each=6)),
                       "trt_total_as_Be", "trt_total_ns_Be", "trt_total_pro_Be", "trt_total_as_Af", "trt_total_pro_Af")

  
  # TODO histogram with ggplot
  # dt_match_min_ps_A0_n = dt_match_min_ps_A0
  # colnames(dt_match_min_ps_A0_n) = colnames(dt_match_min_ps_A1)
  # dt_match_compare = rbind(dt_match_min_ps_A0_n, dt_match_min_ps_A1)
  
  # ggplot(subset(dt_match_compare, select = c(A, X2)), aes(length, fill = A)) +
  #   geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
  
  # install.packages("stddif")
  # library(stddif)
  # stddiff.numeric(dt_match_compare, gcol = "A")
  
  # par(mfrow = c(1, length(X_sub_cols[-1])))
  # d0 = subset(dt_match_min_ps_A0, select = -A0_A);  d1 = subset(dt_match_min_ps_A1, select = -A)
  # for(i in c(1:length(X_sub_cols[-1]))){
  # hist(d0[,i], col='skyblue', border=F, xlab = "X2", main = colnames(d1)[i])
  # hist(d1[,i], add=T, col=scales::alpha('green',.5), border=F  )
  # legend('topleft',c('control','treatment'),
  #        fill = c('skyblue', 'green'), bty = 'n',
  #        border = NA)
  # }
  
  # SACE_matching_est_HL instead of SACE_matching_est_CI_HL
  return(list(SACE_matching_est_SE_naive = c(SACE_matching_est, SACE_matching_SE)
              ,SACE_matching_est_HL = c(SACE_matching_est_HL, SACE_matching_se_HL)
              ,SACE_matching_est_SE_BC = c(BCest, BCse)
              ,SACE_matching_est_SE_BCclpr = c(BCest_clpr, BCse_clpr)
              ,CI_crude_HL_BC = c(CI_by_SE_and_Z_val_naive,SACE_matching_CI_HL,CI_by_SE_and_Z_val_BC,CI_by_SE_and_Z_val_BCclpr)
              #,CI_by_SE_and_Z_val_naive = CI_by_SE_and_Z_val_naive
              ,WLS_NOinteractions_reg_adj_estimators_and_se = WLS_NOinteractions_reg_adj_estimators_and_se
              ,WLS_YESinteractions_reg_adj_estimators_and_se = WLS_YESinteractions_reg_adj_estimators_and_se
              ,OLS_NOinteractions_reg_adj_estimators_and_se = OLS_NOinteractions_reg_adj_estimators_and_se
              ,OLS_YESinteractions_reg_adj_estimators_and_se = OLS_YESinteractions_reg_adj_estimators_and_se
              ,included_excluded_in_matching = included_excluded_in_matching
              ,repeated_as_and_pro = repeated_as_and_pro
              ,diff_distance_aspr_asas = diff_distance_aspr_asas
              ,tables_matched_units = tables_matched_units
              ,means_by_subset = means_by_subset
              ,std_diff_2_cols = std_diff_2_cols
              ))
}

tables_repeated = function(table_subjects){
  #table_treated_subjects = data.table(table(ATE_MATCH$index.treated))
  table_treated_subjects = data.table(apply(table_subjects, 2, as.numeric))
  colnames(table_treated_subjects)[1] = "id_trt"
  repeated_treated = table(table_treated_subjects$N)
  return(list(table_treated_subjects = table_treated_subjects, 
              repeated_treated = repeated_treated))
}


#data = filter(dt_match_S1, g == "as" & A0_g == "as");  cols = X_sub_cols[-1]
#cols=X_sub_cols[-1]; data=dt_match_min_ps_w_scheme
check_balance_function = function(data, cols){
  data = data.frame(data)
  diff_w = apply(data[ , cols], 2, mean) - 
    apply(data[ , paste0("A0", "_", cols)], 2, mean)
  std_diff_w = diff_w / 
    apply(data[ , paste0("A0", "_", cols)], 2, sd)  
  return(list(diff_w, std_diff_w))
}


repeated_histogram = function(mat_all_repeated_as_and_pro){
  mat_all_repeated_as_and_pro_t = t(mat_all_repeated_as_and_pro)
  reT_mat_all_repeated_as_and_pro = subset(mat_all_repeated_as_and_pro_t,
                                           select = -grep("repF_", colnames(mat_all_repeated_as_and_pro_t)))
  
  reT_mat_all_repeated_hist = reT_mat_all_repeated_as_and_pro[-grep("mean", 
                                                                    rownames(reT_mat_all_repeated_as_and_pro)),] 
  len_each_strat = nrow(reT_mat_all_repeated_hist) / 2
  title_vec = colnames(reT_mat_all_repeated_hist)
  #par(mfrow = c(2,3))
  pdf(file = "repeated_by_str_hist.pdf")
  for(i in 1 : length(title_vec) ){
    #temp = t(data.frame(reT_mat_all_repeated_hist[,i]))
    #as = data.frame(re = c(1 : len_each_strat), count = temp[,grep("as_", colnames(temp))])
    #pro = data.frame(re = c(1 : len_each_strat), count = temp[,grep("pro", colnames(temp))])
    #ggplot(data = as) + geom_bar(aes(x = re, y = count), stat = "identity")
    
    temp = data.frame(reT_mat_all_repeated_hist[,i])
    temp = data.frame( count = temp[,1], re = rep(c(1 : len_each_strat), times=2),
                       stratum = rep(c("as","pro"), each=len_each_strat) )
    
    print( ggplot(data = temp, aes(x = re, fill = stratum)) +
             stat_identity(data = temp, aes(x = re, y = count), geom = "bar", alpha = 1) + 
             ggtitle(paste0("after matching and excluding pairs with non-surv ", 
                            title_vec[i])) )
    
    # print( ggplot(data = temp, aes(x = re, fill = stratum)) + 
    #          geom_bar(aes(x = re, y = count), stat = "identity") + geom_density(alpha = 0.5) + 
    #          ggtitle(paste0("after matching and excluding pairs with non-surv ", 
    #                         title_vec[i])) )
    # 
    
    
  } 
  dev.off()
}



