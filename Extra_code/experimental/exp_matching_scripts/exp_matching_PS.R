m_data = data_with_PS
#m_data = data_with_PS[OBS != "O(0,0)"]
#m_data = data_with_PS[S==1]
# m_data = d
#min_PS_weighted_match = 0.5
#min_diff_PS = 0.15

# caliper is in sd

my_matching_func_basic = function(replace_bool, j, X_sub_cols, m_data, weighting = FALSE,
                                  M=1, replace, estimand = "ATC", mahal_match = 2,
                                  min_PS = 0, min_diff_PS, caliper = 0.2, min_PS_weighted_match, 
                                  OBS_table){
  #X_sub_cols = paste0("X", c(1:(dim_x)))
  # mahal_match for Weight = 2 for mahalanobis distance. 1 for inverse of variance
  # TODO find a way to use caliper on the PS in the matching function
  ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
                         , X = subset(m_data, select = X_sub_cols[-1]), ties=FALSE
                         #, X=m_data[,est_p_as]
                         #,caliper = caliper
                         ,M=M, replace = replace, estimand = estimand, Weight = mahal_match)
  
  print(ATE_MATCH_PS$estimand)
  #ATE_MATCH_PS$est; ATE_MATCH_PS$se; ATE_MATCH_PS$nobs; ATE_MATCH_PS$index.dropped
  ncols  = ncol(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                       select = c("id", "p_as", "est_p_as", "Y", "A", "S", "g", X_sub_cols[-1]))) + 1
  dt_match = data.table(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                               select = c("id", "p_as", "est_p_as", "Y", "A", "S", "g", X_sub_cols[-1])),
                        ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control,
                        subset(m_data[ATE_MATCH_PS$index.control, ], 
                               select = c("id", "p_as", "est_p_as", "Y", "A", "S", "g", X_sub_cols[-1])))
  colnames(dt_match)[(ncols + 1): (2 * ncols)] = 
    paste0("A0_", colnames(dt_match)[(ncols + 1): (2 * ncols)])
  colnames(dt_match)[c(ncols: (ncols+1))] = c("id_trt", "id_ctrl")
  #dt_match_as_A0 = dt_match[A0_S==1, 6:10]
  
  ##### min_diff_PS (almost caliper) ####
  # TODO until I'll find how to fo it in the matching function, use caliper in simple way
  # which just remove pairs which not close enough
  dt_match_min_ps = dt_match[abs(est_p_as - A0_est_p_as) <= min_diff_PS , ]
  
  ##### min_PS  ##### 
  # use min_ps as min ps for the treated
  dt_match_min_ps = filter(dt_match_min_ps, est_p_as > min_PS)
  # use min_ps as min ps for the treated, not sure about it, since under mono, their all as
  #dt_match_min_ps = filter(dt_match, A0_est_p_as > min_PS)
  
  # keep only S = 1
  
  dt_match_min_ps = filter(dt_match_min_ps, S == 1 & A0_S==1)
  # estimation
  
  # TODO diff of means
  #SACE_matching_est = mean(dt_match_min_ps$Y) - mean(dt_match_min_ps$A0_Y) 
  
  
  # TODO mean of diffs 
  diff_per_pair = dt_match_min_ps$Y - dt_match_min_ps$A0_Y
  #SACE_matching_est = mean(diff_per_pair)
  
  # TODO HL
  wilcoxon = wilcox.test(diff_per_pair,conf.int=T)
  SACE_matching_est_HL = wilcoxon$estimate
  
  ######## calculating the amount of as the matching process excluded
  OBS_table * param_n
  as_A0_matched = length(which(dt_match_min_ps$A0_g == "as"))
  # OBS_table[1,2] are the (A=0, S=1), thereare only as in this cell
  as_A0_unmatched = (OBS_table[1,2] * param_n) - as_A0_matched 
  as_A1_matched = length(which(dt_match_min_ps$g == "as"))
  pro_A1_matched = length(which(dt_match_min_ps$g == "pro"))
  
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
  dt_match_min_ps_A0 = dt_match_min_ps[ , c("A0_A", paste0("A0", "_", X_sub_cols[-1]))]
  dt_match_min_ps_A1 = dt_match_min_ps[ , c("A", X_sub_cols[-1])]
  
  # TODO adjust for the case of repacements, so if we have 1000 treated,
  # some of them are the same subjects, see table_treated_subjects
  if(replace == TRUE){
    length(unique(ATE_MATCH_PS$index.treated))
    table_treated_subjects = data.table(table(ATE_MATCH_PS$index.treated))
    table_treated_subjects = data.table(apply(table_treated_subjects, 2, as.numeric))
    colnames(table_treated_subjects)[1] = "id_trt"
    repeated = table(table_treated_subjects$N)
    temp = subset(dt_match_min_ps, select = c(id_trt, g))
    matched_S1_unique = unique(merge(table_treated_subjects, temp, all.x = TRUE, all.y = FALSE, by = "id_trt"))
    S1_matched_distinguish_by_g =  
      ddply(matched_S1_unique, .(g), summarize, amount_of_subjects = length(id_trt))
    # repeated = ddply(table_treated_subjects, .(N), summarize, amount_of_subjects = length(id_trt))
  }
  
  
  
  #### wout manual caliper
  dt_match_min_ps_w_scheme = data.table(filter(dt_match, est_p_as > min_PS))
  ### keep only S = 1
  dt_match_min_ps_w_scheme = filter(dt_match_min_ps_w_scheme, S == 1 & A0_S==1)
  ### estimation
  SACE_matching_PS_weighted_match_est = 
    mean(dt_match_min_ps_w_scheme$Y) - mean(dt_match_min_ps_w_scheme$A0_Y)
  
  # checking covariates balance between treated and untreated 
  
  # with manu caliper matching
  # balance_as_to_as_manu_cal = check_balance_function(
  #   filter(dt_match_min_ps, g == "as" & A0_g == "as"), X_sub_cols[-1])
  # balance_as_to_pro_manu_cal = check_balance_function(
  #   filter(dt_match_min_ps, g == "pro" & A0_g == "as"), X_sub_cols[-1])
  # std_diff_w_as_to_as_manu_cal = balance_as_to_as_manu_cal[[2]]
  # std_diff_w_as_to_pro_manu_cal = balance_as_to_pro_manu_cal[[2]]
  
  # without manu caliper matching
  balance_as_to_as = check_balance_function(
    filter(dt_match_min_ps_w_scheme, g == "as" & A0_g == "as"), X_sub_cols[-1])
  balance_as_to_pro = check_balance_function(
    filter(dt_match_min_ps_w_scheme, g == "pro" & A0_g == "as"), X_sub_cols[-1])
  std_diff_w_as_to_as = balance_as_to_as[[2]]; std_diff_w_as_to_pro = balance_as_to_pro[[2]]
  diff_distance_aspr_asas = mean(abs(std_diff_w_as_to_pro) - abs(std_diff_w_as_to_as))
  
  #histogram with ggplot
  
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
  
  return(list(SACE_matching_est_HL, SACE_matching_PS_weighted_match_est, 
              included_excluded_in_matching, diff_distance_aspr_asas, 
              std_diff_w_as_to_as, std_diff_w_as_to_pro))
}

#cols=X_sub_cols[-1]; data=dt_match_min_ps_w_scheme
check_balance_function = function(data, cols){
  diff_w = apply(data[ , cols], 2, mean) - 
    apply(data[ , paste0("A0", "_", cols)], 2, mean)
  std_diff_w = diff_w / 
    apply(data[ , paste0("A0", "_", cols)], 2, sd)  
  return(list(diff_w, std_diff_w))
}





# 
# #@@@@@@weighting_scheme
# ###@@@@@@ min_PS FOR WEIGHTS
# dt_match_min_ps_w_scheme = data.table(filter(dt_match, est_p_as > min_PS_weighted_match))
# ###################################   weights by p_as  ######################################  
# dt_match_min_ps_w_scheme$diff = dt_match_min_ps_w_scheme[ ,Y] - dt_match_min_ps_w_scheme[ , A0_Y]
# ###@@@@@@ here were assign weights to the diffs and than doing calipper (1 - 0.6 = 0.4)
# ###@@@@@@ we can try the other way- first calipper and than calculate weights
# dt_match_min_ps_w_scheme$w = dt_match_min_ps_w_scheme[ , est_p_as] /
#   mean(dt_match_min_ps_w_scheme[ , est_p_as])
# 
# ###@@@@@@ sum(dt_match_w_scheme$w) == as.numeric(nrow(dt_match_w_scheme))
# dt_match_min_ps_w_scheme[ , weighted_diff := w * diff]
# ###@@@@@@ keep only S = 1
# dt_match_min_ps_w_scheme = filter(dt_match_min_ps_w_scheme, S == 1 & A0_S==1)
# ###@@@@@@ estimation
# SACE_matching_PS_weighted_match_est = mean(dt_match_min_ps_w_scheme$weighted_diff)
