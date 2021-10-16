# matching on Principal score
library("Matching")
#m_data = data_with_PS
#m_data = data_with_PS[OBS != "O(0,0)"]
#m_data = data_with_PS[S==1]
# m_data = d
#min_PS_weighted_match = 0.5
#min_diff_PS = 0.15

# TODO caliper is in sd
#replace = T; j = 1; reg_covariates = c("Y", "A", "re74", "age", "married")
my_matching_func_real_data = function(j, X_sub_cols, reg_covariates, 
          m_data, weighting = FALSE, M=1, replace, 
          estimand = "ATC", mahal_match = 2, caliper = 0.2, OBS_table){
  #X_sub_cols = paste0("X", c(1:(dim_x)))
  # mahal_match for Weight = 2 for mahalanobis distance. 1 for inverse of variance
  # TODO find a way to use caliper on the PS in the matching function
  # X_sub_cols[-1] est_p_as
  ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A]
               , X = subset(m_data, select = X_sub_cols[-1]), ties=FALSE
               #, X=m_data[,est_p_as]
               #,caliper = caliper
               ,M=M, replace = replace, estimand = estimand, Weight = mahal_match)
  
  print(ATE_MATCH_PS$estimand)
  # mb  <- MatchBalance(treat~age + I(age^2) + educ + I(educ^2) + black +
  #                       hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
  #                       u74 + u75, data=lalonde, match.out=rr, nboots=10)
  
  ncols  = ncol(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                       select = c("id", "p_as", "est_p_as", "Y", "A", "S", X_sub_cols[-1]))) + 1
  dt_match = data.table(subset(m_data[ATE_MATCH_PS$index.treated, ], 
                               select = c("id", "p_as", "est_p_as", "Y", "A", "S", X_sub_cols[-1])),
                        ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control,
                        subset(m_data[ATE_MATCH_PS$index.control, ], 
                               select = c("id", "p_as", "est_p_as", "Y", "A", "S", X_sub_cols[-1])))
  colnames(dt_match)[(ncols + 1): (2 * ncols)] = 
    paste0("A0_", colnames(dt_match)[(ncols + 1): (2 * ncols)])
  colnames(dt_match)[c(ncols: (ncols+1))] = c("id_trt", "id_ctrl")
  
  # keep only S = 1
  dt_match_S1 = filter(dt_match, S == 1 & A0_S==1)
  #dt_match_S1 = dt_match
  
  # estimation
  # TODO diff of means
  # SACE_matching_est_log = mean(dt_match_S1$Y) - mean(dt_match_S1$A0_Y)
  # sd_log_diff_per_pair =  sd(dt_match_S1$Y - dt_match_S1$A0_Y)
  
  # TODO mean of diffs 
  diff_per_pair = dt_match_S1$Y - dt_match_S1$A0_Y
  SACE_matching_est = mean(diff_per_pair)
  SACE_matching_sd = sd(diff_per_pair)
  # SACE_matching_SE =  sqrt((sd(dt_match_S1$Y))^2/ (nrow(dt_match_S1)) +
  #   (sd(dt_match_S1$A0_Y))^2 / (nrow(dt_match_S1)))
  SACE_matching_SE = sd(diff_per_pair) / sqrt(nrow(dt_match_S1))
  # nrow(dt_match_S1) - 1 ; # https://www.youtube.com/watch?v=zD3VIBkwc-0
  t_val = qt(0.975, nrow(dt_match_S1) - 2, lower.tail = TRUE, log.p = FALSE)
  CI_by_SE_and_t_val = SACE_matching_est + c(-1,1) * t_val * SACE_matching_SE 
  
  t_test = t.test(dt_match_S1$Y, dt_match_S1$A0_Y, var.equal = FALSE, paired = TRUE
                  #,alternative = c("two.sided")
  ) 
  # Regression adjusted
  # matched set
  reg_data_matched = m_data[c(unique(dt_match_S1$id), unique(dt_match_S1$A0_id)) ,]
  # reg_covariates
  reg_data_matched =  subset(reg_data_matched, select = reg_covariates)
  lin_reg_matched = lm(Y~., reg_data_matched)
  sum(reg_data_matched$Y==0)
  sum_matched = summary(lin_reg_matched)
  
  # Regression adjusted
  # all set
  reg_data =  subset(m_data, select = reg_covariates)
  lin_reg = lm(Y~., reg_data)
  sum(reg_data$Y==0)
  sum_all_set = summary(lin_reg)
  
  # HL
  wilcoxon = wilcox.test(diff_per_pair,conf.int=T)
  SACE_matching_est_HL = wilcoxon$estimate
  
  
  ######## calculating the amount of as the matching process excluded
  OBS_table_numbers = OBS_table * nrow(m_data)
  # checking covariates balance between treated and untreated
  dt_match_A0 = dt_match_S1[ , c("A0_A", paste0("A0", "_", X_sub_cols[-1]), "A0_p_as")]
  dt_match_A1 = dt_match_S1[ , c("A", X_sub_cols[-1], "est_p_as")]
  
  # TODO adjust for the case of repacements, so if we have 1000 treated,
  # some of them are the same subjects, see table_treated_subjects
  
  # if(replace == TRUE){
  #   length(unique(ATE_MATCH_PS$index.treated))
  #   length(unique(ATE_MATCH_PS$index.control))
  #   table_treated_subjects = data.table(table(ATE_MATCH_PS$index.treated))
  #   table_treated_subjects = data.table(apply(table_treated_subjects, 2, as.numeric))
  #   colnames(table_treated_subjects)[1] = "id_trt"
  #   repeated = table(table_treated_subjects$N)
  #   temp = subset(dt_match_min_ps, select = c(id_trt, g))
  #   matched_S1_unique = unique(merge(table_treated_subjects, temp, all.x = TRUE, all.y = FALSE, by = "id_trt"))
  #   S1_matched_distinguish_by_g =
  #     ddply(matched_S1_unique, .(g), summarize, amount_of_subjects = length(id_trt))
  #   # repeated = ddply(table_treated_subjects, .(N), summarize, amount_of_subjects = length(id_trt))
  # }
  
  
  # without manu caliper matching
  balance = check_balance_function(dt_match_S1, X_sub_cols[-1])
  mean_diff = balance[[1]]; std_dif = balance[[2]]
  data.frame(rbind(std_dif))
  #histogram with ggplot
  
  # dt_match_A0_n = dt_match_A0
  # colnames(dt_match_A0_n) = colnames(dt_match_A1)
  # dt_match_compare = rbind(dt_match_A0_n, dt_match_A1)
  # 
  # ggplot(subset(dt_match_compare, select = c(A, X2)), aes(length, fill = A)) +
  #  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')
  # 
  # install.packages("stddif")
  # library(stddif)
  # stddiff.numeric(dt_match_compare, gcol = "A")
  
  par(mfrow = c(1, length(X_sub_cols[-1])))
  par(mfrow = c(2, 2))
  d0 = subset(dt_match_A0, select = -A0_A);  d1 = subset(dt_match_A1, select = -A)
  for(i in c(1:length(X_sub_cols[-1]))){
    hist(d0[,i], col='skyblue', border=F, xlab = colnames(d1)[i], breaks=12, main = colnames(d1)[i])
    hist(d1[,i], add=T, col=scales::alpha('green',.5), border=F, breaks=12  )
    legend('topright',c('control','treatment'),
           fill = c('skyblue', 'green'), bty = 'n',
           border = NA)
  }
  
  return(list(SACE_matching_est = SACE_matching_est, SACE_matching_sd= SACE_matching_sd,
              SACE_matching_SE = SACE_matching_SE, CI_by_SE_and_t_val = CI_by_SE_and_t_val,
              SACE_matching_est_log = SACE_matching_est_log, sd_log_diff_per_pair = sd_log_diff_per_pair, 
              sum_matched = sum_matched, sum_all_set = sum_all_set,
              t_test = t_test, SACE_matching_est_HL = SACE_matching_est_HL, wilcoxon = wilcoxon,
              std_dif = std_dif))
}

#cols=X_sub_cols[-1]; data=dt_match_min_ps_w_scheme
check_balance_function = function(data, cols){
  diff_w = apply(data[ , cols], 2, mean) - 
    apply(data[ , paste0("A0", "_", cols)], 2, mean)
  std_diff_w = diff_w / 
    apply(data[ , paste0("A0", "_", cols)], 2, sd)  
  return(list(diff_w, std_diff_w))
}





