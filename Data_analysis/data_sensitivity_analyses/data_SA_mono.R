set.seed(101) 
# EM parameters
two_log_est_EM = FALSE
iterations_EM = 400; epsilon_EM = 1e-06
#########################################################################################
# sensitivity parameters ####
#bounds for xi
p1 = mean(data[A==1,S]); p0 = mean(data[A==0,S])
up_bound_xi = (1 - p1) / (p0 - (1 - p1))
xi_sensi_mono_vec = seq(0, 0.5, 0.1)
xi_sensi_mono_vec[length(xi_sensi_mono_vec)] = 0.48 #round(up_bound_xi,2)
xi_sensi_mono_names = paste0("xi_mono_", round(xi_sensi_mono_vec, 2))
# alpha0_mono
alpha0_mono_vec = seq(0.5, 2, 0.25) 
alpha0_mono_vec_names = paste0("alpha0_mono_", alpha0_mono_vec)  

# SA for each combinations of xi and alpha_0  and for all distance metrics ####
reg_sensi_mono <- NULL
for (j in 1:length(xi_sensi_mono_names)) {
  xi = xi_sensi_mono_vec[j]
  #########################################################################################
  tmp = data 
  # EM
  # calculate PS for the current xi
  
  if(two_log_est_EM == FALSE){
    #S(0)=1: Logistic regression S(0)=1 on X, using S|A=0
    fit_S0_in_A0 =
      glm(as.formula(paste0("S ~ ",paste(covariates_PS, collapse="+"))), data=filter(data, A==0), family="binomial")
    beta_S0 = fit_S0_in_A0$coefficients
  }else{beta_S0=NULL}
  est_ding_lst_SA_mono = xi_2log_PSPS_M_weighting(Z=tmp$A, D=tmp$S,
                           X=as.matrix(subset(tmp, select = covariates_PS)), Y=tmp$Y,
                           xi_est=xi, beta.S0=beta_S0, beta.ah=NULL, beta.c=NULL,
                           iter.max=iterations_EM, error0=epsilon_EM)
  # est_ding_lst_SA_mono = xi_PSPS_M_weighting_SA(Z=tmp$A, D=tmp$S, X=as.matrix(subset(tmp, select = covariates_PS)),  
  #                       Y=tmp$Y, eta=xi, beta.c = NULL, beta.n = NULL)
  DING_model_assisted_sensi = est_ding_lst_SA_mono$AACE.reg
  # Ding order: c(prob.c, prob.d, prob.a, prob.n)
  PS_est_sensi_mono = est_ding_lst_SA_mono$ps.score
  # adjust the cols the same order as in myEM: my order is: as, ns, pro, har.
  PS_est = data.frame(EMest_p_as = PS_est_sensi_mono[,3], EMest_p_ns = PS_est_sensi_mono[,4],
                      EMest_p_pro = PS_est_sensi_mono[,1], EMest_p_har = PS_est_sensi_mono[,2])
  tmp = data.table(tmp, PS_est)
  which(is.na(tmp)==TRUE)
  tmp$e_1_as = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_pro)
  tmp$e_0_as = tmp$EMest_p_as / (tmp$EMest_p_as + tmp$EMest_p_har)
  tmp = data.table(subset(tmp, select = -c(e_1_as, e_0_as)),  subset(tmp, select = c(e_1_as, e_0_as)))
  ##########################################################################

  # matching for the current xi ####
  # set.seed because otherwise the matching procedure will not yield the same results when alpha_0=1 under all values of xi, during the SA for monotonicity (for mahalanobis measure) 
  # also, otherwise, when xi=0 in SA for monotonicity results will not be similar to results when alpha_1=1 in SA for PPI (for mahalanobis measure)
  #TODO 22.02.22 even when we set.seed here matching on malanobis alone yields different results. Need to set.seed before every matching procedure for the same results
  set.seed(101) 
  matching_lst = matching_func_multiple_data(match_on = match_on
           ,cont_cov_mahal = cont_cov_mahal,  reg_cov = reg_after_match, X_sub_cols = variables
           ,reg_BC = reg_BC, m_data = tmp[S==1]
           ,w_mat_bool = "NON-INFO", M=1, replace=TRUE, estimand = "ATC", mahal_match = 2, caliper = caliper
           ,change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units=FALSE, one_leraner_bool=TRUE)
 
  # matched_set_lst for all distance metrics
  matched_data_lst = matching_lst$matched_set_lst
  reg_matched_lst = matching_lst$reg_data_matched_lst
  # run on all distance metrics
  for(ind_matched_set in c(1:length(reg_matched_lst))){  
    matched_data = matched_data_lst[[ind_matched_set]]
    reg_data_matched_SA = reg_matched_lst[[ind_matched_set]]
    print(paste0("unique weights for control are really = ", unique(filter(matched_data, A==0)$w)))
    
    #regression on the matched dataset with original Y
    coeffs_regression_one_model = regression_function_one_model(reg_data_matched_SA=reg_data_matched_SA,
      data_reg=matched_data, reg_after_match=reg_after_match) 
    coeffs_regression_two_models = regression_function_two_models(reg_data_matched_SA=reg_data_matched_SA,
      data_reg=matched_data, reg_after_match=reg_after_match) 
    
    for (i in 1:length(alpha0_mono_vec_names)) {
      print(alpha0_mono_vec_names[i])
      alpha0 = alpha0_mono_vec[i]; alpha1 = 1
      
      # ONE-LEARNER regression approach
      # predictions for units from O(0,1), plugging A={0,1} + 4. sensitivity adjustments
      reg_sensi_mono =
        rbind( reg_sensi_mono, c(measure = names(matched_data_lst)[ind_matched_set], xi_mono = xi, alpha0_mono = alpha0, 
        unlist(SACE_estimation_1LEARNER_mono(matched_data=matched_data, reg_after_match=reg_after_match, alpha0_mono=alpha0, xi=xi, 
        coeffs_regression_one_model=coeffs_regression_one_model, coeffs_regression_two_models=coeffs_regression_two_models, two_models_bool=TRUE))) )
      
      print(paste0(xi_sensi_mono_names[j], " ", alpha0_mono_vec_names[i]))
    }
  }
}
#########################################################################################

#########################################################################################\
# process before plotting ####
reg_sensi_mono = data.frame(reg_sensi_mono)
reg_sensi_mono[,-1] = apply(reg_sensi_mono[,-1] , 2, as.numeric)
reg_sensi_mono[,-c(1,2,3)] = round(reg_sensi_mono[,-c(1,2,3)])
reg_sensi_mono_est = reg_sensi_mono[,c(1:3,4,6,8)] %>% gather("Estimator", "Estimate", c(4:6)) %>% arrange(measure, xi_mono, alpha0_mono)
reg_sensi_mono_se = reg_sensi_mono[,c(1:3,5,7,9)] %>% gather("Estimator", "SE", c(4:6)) %>% arrange(measure, xi_mono, alpha0_mono)
reg_sensi_mono_se$Estimator = mgsub(reg_sensi_mono_se$Estimator, c("\\_se$"), "")
reg_sensi_mono = merge(reg_sensi_mono_est, reg_sensi_mono_se, by=c("measure", "xi_mono", "alpha0_mono", "Estimator"))
reg_sensi_mono$lower_CI = reg_sensi_mono$Estimate - 1.96 * reg_sensi_mono$SE
reg_sensi_mono$upper_CI = reg_sensi_mono$Estimate + 1.96 * reg_sensi_mono$SE

legend_levels = c("Crude", "WLS", "WLS inter")
reg_sensi_mono$Estimator = mgsub(reg_sensi_mono$Estimator,
                                 c("crude_est_adj", "SACE_1LEARNER_adj", "SACE_1LEARNER_inter_adj"), legend_levels)
reg_sensi_mono$Estimator = factor(reg_sensi_mono$Estimator, levels = legend_levels)
reg_sensi_mono$set = data_bool
#########################################################################################

#########################################################################################
# plot SA for monotonicity under PPI, as a function of xi and alpha_0 ####
plot_SA_mono <- ggplot(filter(reg_sensi_mono, measure == "Mahal_PS_cal" & Estimator %in% c("WLS")), 
                          aes(x=alpha0_mono, y=Estimate)) + 
  geom_point(aes(col = Estimator, size = 7), size = 3) + theme_bw() + 
  scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + 
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.1, position = "dodge", linetype="solid", color = "gray53") + 
  labs(colour = "Estimator"
       , size = 1) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(-2600,4000) +
  ylab(label="Estimate") + xlab(label = bquote(alpha[0])) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)
    ,legend.position="none" # remove legend
    ) + 
  geom_hline(yintercept = 0)

plot_SA_mono = plot_SA_mono +
  facet_grid(~ glue('xi*" = {xi_mono}"'), labeller = label_parsed) +
  theme(
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14), # X axis title
    axis.title.y=element_text(size=14), # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)   # Y axis text
  ) 
#########################################################################################

#########################################################################################
# plot SA for monotonicity under SPPI (alpha_0=1) under several xi values ####
reg_sensi_mono$measure = mgsub(as.character(reg_sensi_mono$measure), "_", " ")
reg_sensi_mono$measure = factor(reg_sensi_mono$measure, levels = c("Mahal", "Mahal PS cal", "PS"))

plot_SA_by_metric <- reg_sensi_mono %>% filter(alpha0_mono==1 & !Estimator=="WLS inter") %>% ggplot(aes(x=xi_mono, y=Estimate)) +
  geom_point(aes(col = Estimator, size = 7), size = 2) + 
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
  ) + 
  ylab(label="Estimate") +
  xlab(label = bquote(xi)) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  geom_hline(yintercept = 0 )

plot_SA_by_metric = plot_SA_by_metric + 
  scale_color_manual(name="Estimator", 
                     labels = legend_levels, 
                     values = c("Crude" = "orangered2", "WLS" = "green3", "WLS inter" = "cornflowerblue"))  +
  facet_grid(~ measure) + # set ~ measure # Metric ~ Set
  scale_x_continuous(breaks=xi_sensi_mono_vec) + 
  theme(
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=8, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  ) 

# plot_SA_mono_byMetric = plot_SA_mono + facet_wrap(~ Metric, ncol=3)

# EXTRACT LEGEND
library(cowplot); library(ggpubr)
lgnd_plt <- get_legend(plot_SA_by_metric) 
# Convert to a ggplot and print
as_ggplot(lgnd_plt)
plot_SA_woutLGND = plot_SA_by_metric + theme(legend.position = 'none') 
#########################################################################################


