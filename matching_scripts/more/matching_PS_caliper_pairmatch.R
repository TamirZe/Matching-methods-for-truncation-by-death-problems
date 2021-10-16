# examples in DOS cran and in rosenbaum book: design of obs studies
#caliper = 0.1
#m_data = data_with_PS
#m_data = data_with_PS[OBS != "O(0,0)"]
# m_data = data_with_PS[S==1]
# m_data = d
#min_PS_weighted_match = 0.5
#min_diff_PS = 0.15

my_matching_func_pairmatch = function(X_sub_cols, m_data,
                            M=1, replace = FALSE,
                            estimand = "ATC", mahal_match = 2, caliper, OBS_table){
  
  x <- subset(m_data, select = X_sub_cols[-1])
  # TODO adjust to the case of categorial (more than 1 ption) factor covariate X
  f = as.formula(paste0("A ~ ", paste(X_sub_cols[-1], collapse = " + ")
                        #, " + scores(est_p_as)"
  ))
  principal_score = m_data$est_p_as
  names(principal_score) <- rownames(m_data)
  # in this function, caliper is probably in absolute value, and not in sd(principal_score)
  mhd <- match_on(f, data = m_data) +
    caliper(match_on(principal_score, z = m_data$A), caliper)
   pm2 <- pairmatch(mhd, data = m_data
                   # controls = M,  
                    ) 
  summary(pm2); pm2 = na.omit(pm2)
  matched_data = data.table(m_data); matched_data$id = c(1:nrow(matched_data))
  matched_data = matched_data[as.numeric(names(pm2)) , ]
  matched_data$matched_pair = as.numeric(pm2)
  matched_data = data.table(arrange(matched_data, matched_pair, desc(A)))
  
  diff_each_pair = matched_data[, list(diff = -diff(Y)), by = matched_pair]
  SACE_matching_est_with_SO = mean(diff_each_pair$diff)
  
  # diff_each_pair2 = ddply(matched_data, .(matched_pair), summarize,
  #                        diff = Y[A==1] - Y[A==0])
  # diff_each_pair3 = ddply(matched_data, .(matched_pair), summarize,
  #            diff = -diff(Y))
  # all.equal(diff_each_pair3, diff_each_pair, check.attributes=FALSE)
  
  
  # remove pairs when one of the subjects is wout S=1
  matched_data = matched_data[, surv_pair := mean(S), by = matched_pair]
  matched_data_S1 = data.table(filter(matched_data, surv_pair==1))
  diff_each_pair_S1 = matched_data_S1[, list(diff = -diff(Y)), by = matched_pair]
  SACE_matching_est = mean(diff_each_pair_S1$diff)
    
  
  
  ######################################################################
  ######## calculating the amount of as the matching process excluded
  #OBS_table
  as_A0_matched = length(which(matched_data_S1$g == "as" & matched_data_S1$A==0))
  as_A0_unmatched = OBS_table[1,2] * param_n - as_A0_matched
  as_A1_matched = length(which(matched_data_S1$g == "as" & matched_data_S1$A==1))
  pro_A1_matched = length(which(matched_data_S1$g == "pro" & matched_data_S1$A==1))
  
  # ns senity check
  ns_A1_matched = length(which(matched_data$g == "ns" & matched_data$A==1))
  ns_A0_matched = length(which(matched_data$g == "ns" & matched_data$A==0))
  ns = ns_A1_matched + ns_A0_matched
 
  # check the compisition in A1_S1: as and protected
  as_in_A1_S1 = matched_data_S1[g =="as" & A == 1,]
  pro_in_A1_S1 = matched_data_S1[g =="pro" & A == 1,]
  
  as_A1_unmatched = 
  length(which(m_data$g == "as" & m_data$A==1)) - as_A1_matched
  
  included_excluded_in_matching = 
    data.frame(as_A0_matched, as_A0_unmatched, as_A1_matched, 
               as_A1_unmatched, pro_A1_matched)

  # TODO check for covariates balance
  
  return(list(SACE_matching_est_with_SO, SACE_matching_est,
              included_excluded_in_matching))
}





