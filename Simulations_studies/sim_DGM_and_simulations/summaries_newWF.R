#############################################################################################
calculate_mean_repeated_as_and_pro = function(list_of_lists, mean_repeated_as_and_pro_boll=TRUE){
  temp = lapply(1:length(list_of_lists), function(i){
    list_of_lists[[i]] = as.matrix(list_of_lists[[i]])
    apply(list_of_lists[[i]], 2, as.numeric)
  })
  mean_repeated_as_and_pro = data.frame(apply(simplify2array(temp), 1:2, mean))
  rownames(mean_repeated_as_and_pro) = rownames(list_of_lists[[1]])
  if(mean_repeated_as_and_pro_boll == TRUE){
    as_rep = subset(mean_repeated_as_and_pro,
                    select = grep("as_", colnames(mean_repeated_as_and_pro)))
    pro_rep = subset(mean_repeated_as_and_pro,
                     select = grep("pro_", colnames(mean_repeated_as_and_pro)))
    #mean_as = as.matrix(as_rep) %*% c(1:ncol(as_rep)) / ifelse(as_rep>0, 1, 0) %*% c(1:ncol(as_rep))
    mean_as = as.matrix(as_rep) %*% c(1:ncol(as_rep)) / apply(as_rep, 1, sum)
    mean_pro = as.matrix(pro_rep) %*% c(1:ncol(pro_rep)) / apply(pro_rep, 1, sum)
    mean_repeated_as_and_pro = data.frame(mean_repeated_as_and_pro, mean_as = mean_as, mean_pro = mean_pro)
  }
  return(mean_repeated_as_and_pro)
}
#############################################################################################

#############################################################################################
calculate_coverage = function(CI_mat){
  bool_CI_contains_SACE = function(x, true_SACE){
    lower_bound =  as.numeric(sub(",.*", "", x)) #as.numeric(sapply(strsplit(as.character(x), ","), "[", 1)) 
    upper_bound = as.numeric(sub(".*,", "", x))
    bool = lower_bound <= true_SACE & upper_bound >= true_SACE   
    bool = ifelse(bool == T, 1, 0)
    return(bool)
  }
  bool_mat = data.frame(apply(CI_mat[,-1], 2, bool_CI_contains_SACE, true_SACE))
  bool_mat = data.frame(rbind(bool_mat, apply(bool_mat, 2, mean)))
  rownames(bool_mat) = c(1:nrow(CI_mat), "Coverage")
  return(bool_mat)
}
#############################################################################################

#############################################################################################
# summaries over all simulation iterations

index_EM_not_conv / param_n_sim # check number of iterations where the EM algo has not converged
i_EM_not_conv

# summary of matching_estimators_sum: mean, med, empirical sd and MSE ####
param_SACE = mean(matching_estimators_mat[,"SACE"])
MSE_fun <- function (x) mean((x-true_SACE)^2) # mean((x-param_SACE)^2) # mean((x-true_SACE)^2)
matching_estimators_sum = data.frame( rbind(na.omit(matching_estimators_mat), 
     mean = apply(na.omit(matching_estimators_mat), 2, mean), med = apply(na.omit(matching_estimators_mat), 2, median),
     SD = apply(na.omit(matching_estimators_mat), 2, sd), MSE = apply(na.omit(matching_estimators_mat), 2, FUN = MSE_fun)) )
rownames(matching_estimators_sum) = c(c(1:nrow(na.omit(matching_estimators_mat))), c("mean", "med", "sd", "MSE"))
matching_estimators_sum = matching_estimators_sum[c("mean", "sd", "MSE"),]

# summary of matching_estimators_SE_mat: mean ####
matching_estimators_SE_sum = data.frame( rbind(na.omit(matching_estimators_SE_mat), 
                                               SE = apply(na.omit(matching_estimators_SE_mat), 2, mean)) )
rownames(matching_estimators_SE_sum) = c(c(1:nrow(na.omit(matching_estimators_SE_mat))), "SE")
matching_estimators_SE_sum = matching_estimators_SE_sum[c("SE"),]

# summary of pis_pis_est_obs_mat: mean ####
pis_pis_est_obs_sum = data.frame( rbind(pis_pis_est_obs_mat, mean = apply(pis_pis_est_obs_mat, 2, mean)) )
rownames(pis_pis_est_obs_sum) = c(c(1:param_n_sim), "mean")

# summary of CI_mat ####
coverage_sum = calculate_coverage(CI_mat=na.omit(CI_mat))
coverage = data.frame(SACE = true_SACE, coverage_sum["Coverage",])

# check existence of ties when using the BC estimator ####
BC_ties_multiple_treated_sum = apply(na.omit(BC_ties_multiple_treated_mat), 2, mean)

# summary of EM coefficient estimators  ####
beta_S0_sum = apply(na.omit(beta_S0_mat), 2, mean)     # na.omit(beta_S0_mat) # coeff_ah_mat[-beta_S0_mat,]
coeff_ah_sum = apply(na.omit(coeff_ah_mat), 2, mean)   # na.omit(coeff_ah_mat) # coeff_ah_mat[-i_EM_not_conv,] 
coeff_pro_sum = apply(na.omit(coeff_pro_mat), 2, mean) # na.omit(coeff_pro_mat) # coeff_pro_mat[-i_EM_not_conv,] 

# summaries of x_obs (original covariates) ####

# summary of mean_by_g before matching
# mapvalues(mean_by_g$g, from = c("har", "as", "ns", "pro"), to = c(0:3))
mean_list_by_g_sum = apply(simplify2array(list_mean_by_g), 2, rowMeans, na.rm = TRUE)
mean_list_by_g_sum[,1] = c(mapvalues(as.numeric(mean_list_by_g_sum[,1]), from = c(0:3), to = c("har", "as", "ns", "pro")))
mean_list_by_g_sum = data.frame(mean_list_by_g_sum)

# summary of balance - mean (over all simulation iterations) of the X's means
# Filter(Negate... ) removes NULL lists (when EM has not converged) from the list_of_lists
balance_wout_rep_sum = calculate_mean_repeated_as_and_pro(
  list_of_lists=Filter(Negate(function(x) is.null(unlist(x))), balance_wout_rep_lst), FALSE)
balance_with_rep_sum = calculate_mean_repeated_as_and_pro(
  list_of_lists=Filter(Negate(function(x) is.null(unlist(x))), balance_with_rep_lst), FALSE)
#############################################################################################

#############################################################################################
# summary tables for LaTeX of the etimators
results_table = data.frame(t(rbind(matching_estimators_sum, matching_estimators_SE_sum, coverage))) %>% round(3)
results_table = results_table[!(row.names(results_table) %in% c("BC_rep_FALSE", "BC_cal_rep_FALSE")),] 
true_SACE
# Reduce(function(x, y) merge(x, y, by=0),
#        list(t(matching_estimators_sum), t(matching_estimators_SE_sum), t(coverage)))
#############################################################################################
