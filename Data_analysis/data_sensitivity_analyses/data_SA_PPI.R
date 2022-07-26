# seed because otherwise, when xi=0 in SA for monotonicity results will not be similar to results when alpha_1=1 in SA for PPI (for mahalanobis measure)
set.seed(101) 
matching_lst = matching_func_multiple_data(match_on = match_on,
       cont_cov_mahal=cont_cov_mahal,  reg_cov=reg_after_match, X_sub_cols=variables, 
       reg_BC=reg_BC, m_data=data_with_PS[S==1], # data_with_PS from data_main (the same data set used for analysis)
       w_mat_bool="NON-INFO", M=1, replace=TRUE, estimand="ATC", mahal_match=2, caliper=caliper,
       boost_HL=FALSE, one_leraner_bool=TRUE)

#TODO matched_set_lst for all distance metrics
matched_data_lst = matching_lst$matched_set_lst
reg_matched_lst = matching_lst$reg_data_matched_lst

eps_sensi_PPI_vec = seq(0.5,2,0.25)
eps_sensi_PPI_names = paste0("eps_PPI_",eps_sensi_PPI_vec)
reg_sensi_PPI <- NULL

#TODO run on all distance mesaures
for(ind_matched_set in c(1:length(reg_matched_lst))){ 
  #m_dat = matching_lst$m_data; dt_match = matching_lst$ATE_MATCH_PS_lst$dt_match_S1; matching_lst$m_data 
  matched_data = matched_data_lst[[ind_matched_set]]
  reg_data_matched_SA = reg_matched_lst[[ind_matched_set]]
  print(paste0("unique weights for control are really = ", unique(filter(matched_data, A==0)$w)))
  
  #regression on the matched dataset with original Y
  coeffs_regression_one_model = regression_function_one_model(reg_data_matched_SA=reg_data_matched_SA,
                                                              data_reg=matched_data, reg_after_match=reg_after_match) 
  coeffs_regression_two_models = regression_function_two_models(reg_data_matched_SA=reg_data_matched_SA,
                                                                data_reg=matched_data, reg_after_match=reg_after_match) 
  #########################################################################################
  
  #TODO calculate estimators for several values of sensitivity parameters for PPI (eps_sensi_PPI)
  #########################################################################################
  for (i in 1:length(eps_sensi_PPI_vec)) {
    print(eps_sensi_PPI_names[i])
    
    # ONE-LEARNER regression approach
    #predictions for units from O(0,1), plugging A={0,1} + 4. sensitivity adjustments
    reg_sensi_PPI = 
      rbind(reg_sensi_PPI, c(measure = names(matched_data_lst)[ind_matched_set], eps_PPI = eps_sensi_PPI_vec[i], 
       unlist(SACE_estimation_1LEARNER_PPI(matched_data=matched_data, reg_after_match=reg_after_match, eps_sensi_PPI=eps_sensi_PPI_vec[i],
       coeffs_regression_one_model=coeffs_regression_one_model, coeffs_regression_two_models=coeffs_regression_two_models,
       two_models_bool=TRUE))) )
  }
}

# process before plotting
reg_sensi_PPI = data.frame(reg_sensi_PPI)
reg_sensi_PPI[,-1] = apply(reg_sensi_PPI[,-1] , 2, as.numeric) %>% data.frame
reg_sensi_PPI[,-c(1,2)] = round(reg_sensi_PPI[,-c(1,2)]) 

reg_sensi_PPI_est = reg_sensi_PPI[,c(1:2,3,5,7)] %>% gather("Estimator", "Estimate", c(3:5)) %>% arrange(measure, eps_PPI)
reg_sensi_PPI_se = reg_sensi_PPI[,c(1:2,4,6,8)] %>% gather("Estimator", "SE", c(3:5)) %>% arrange(measure, eps_PPI)
reg_sensi_PPI_se$Estimator = mgsub(reg_sensi_PPI_se$Estimator, c("\\_se$"), "")
reg_sensi_PPI = merge(reg_sensi_PPI_est, reg_sensi_PPI_se, by=c("measure", "eps_PPI", "Estimator"))
reg_sensi_PPI$lower_CI = reg_sensi_PPI$Estimate - 1.96 * reg_sensi_PPI$SE
reg_sensi_PPI$upper_CI = reg_sensi_PPI$Estimate + 1.96 * reg_sensi_PPI$SE

legend_levels = c("Crude", "WLS", "WLS inter"); data_bool="LL"
reg_sensi_PPI$Estimator = mgsub(reg_sensi_PPI$Estimator,
          c("crude_est_adj", "SACE_1LEARNER_adj", "SACE_1LEARNER_inter_adj"), legend_levels)
reg_sensi_PPI$Estimator = factor(reg_sensi_PPI$Estimator, levels = legend_levels)
reg_sensi_PPI$set = data_bool
#########################################################################################

#########################################################################################
#plot ####
plot_SA_PPI = reg_sensi_PPI %>% filter(measure == "Mahal_PS_cal" & Estimator %in% c("WLS")) %>%
ggplot(aes(x=eps_PPI, y=Estimate)) + 
  geom_point(aes(col = Estimator, size = 7), size = 4) + theme_bw() + 
  scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + 
  geom_line(aes(col = Estimator, size = 2.5), size=2.5) + 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") +
  #scale_linetype_manual(values = c("dotted", "dotted")) +
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
  ) + 
  ylim(-2600,4000) + 
  ylab(label="Estimate") +
  xlab(label = bquote(alpha[1])) + # epsilon[PPI]
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5)), size=FALSE) + 
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    ,legend.position="none" # remove legend
    ) + 
  geom_hline(yintercept = 0)
#########################################################################################

