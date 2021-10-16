#########################################################################################
#TODO 1-LEARNER regression approach:
#0. if SA for PPI: continue to 1 (set xi=0 to invoke monotonicity) 
#0. if SA for mono under PPI: set a value of xi. for this value:
#1. matching for O(0,1)
#2. regression on the matched dataset with original Y
#3. predictions for units from O(0,1), plugging A={0,1}
#4. sensitivity adjustments- Y1 in SA for PPI; Y0 in SA for mono under PPI
#5. if SA for PPI: weight with 1/N
#5. if SA for mono under PPI: 
#5.a weight with f_X(*|\xi)
#5.b ITERATE OVER SEVERAL VALUS OF xi
#########################################################################################

#2. func: regression on the matched dataset with original Y
#########################################################################################
regression_function = function(data_reg, reg_after_match){
  
  #TODO regression wout A-X interactions
  f_wout_intercations = as.formula(paste0("Y ~ ", paste(c("A", reg_after_match), collapse = " + ")))
  model_wout_intercations = lm(f_wout_intercations, data_reg, weights = w) 
  coeffs_wout_intercations = model_wout_intercations$coefficients
  TE_wout_intercations = coeffs_wout_intercations["A"]
  
  #TODO regression with A-X interactions
  f_with_intercations = as.formula(paste0("Y ~ ", paste(c("A", reg_after_match), collapse = " + "), " + ",
                                          paste(rep("A*",5), reg_after_match, collapse=" + "))) 
  model_with_intercations = lm(f_with_intercations, data_reg, weights = w) 
  mean_x = apply(subset(filter(data_reg, A==0), select = reg_after_match),2, mean)
  coeffs_with_interactions = model_with_intercations$coefficients
  coeffs_interactions = coeffs_with_interactions[grep(":", names(coeffs_with_interactions))]
  TE_with_intercations = coeffs_with_interactions["A"] + coeffs_interactions %*% mean_x
  
  return(list(coeffs_wout_intercations=coeffs_wout_intercations, coeffs_with_interactions=coeffs_with_interactions))
}
#########################################################################################


#3. func: predictions for units from O(0,1), plugging A={0,1}
#########################################################################################
SACE_estimation_1LEARNER = function(matched_data, coeffs_regression, reg_after_match, 
                                    eps_sensi_PPI=1, eps0_sensi_mono=1, eta=0, identific_xi=TRUE){
  # coefficients
  # wout interactions
  coeffs_regression_wout_inter = coeffs_regression$coeffs_wout_intercations
  coeff_A_wout_inter = coeffs_regression_wout_inter["A"]
  coeffs_wout_inter = coeffs_regression_wout_inter[-grep("A", names(coeffs_regression_wout_inter))]
  # with interactions
  coeffs_regression_with_inter = coeffs_regression$coeffs_with_interactions
  coeff_A_with_inter = coeffs_regression_with_inter["A"]; intercept_with_inter = coeffs_regression_with_inter["(Intercept)"]
  coeffs_with_inter_A1 = coeffs_regression_with_inter[-grep("A|(Intercept)", names(coeffs_regression_wout_inter))]
  coeffs_with_inter_A0 = coeffs_regression_with_inter[-grep(":|A", names(coeffs_regression_with_inter))]
  
  # predictions
  A0_S1_data = filter(matched_data, A==0 & S==1)
  attach(A0_S1_data)
  
  # wout interactions
  A0_S1_data$Y1_pred = coeff_A_wout_inter + 
    as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
  A0_S1_data$Y1_pred_adj =  A0_S1_data$Y1_pred / 
    ( ((1 - eps_sensi_PPI) * A0_S1_data$e_1_as) +  eps_sensi_PPI )
  
  A0_S1_data$Y0_pred = as.matrix(subset(A0_S1_data, select = c("intercept", reg_after_match))) %*% coeffs_wout_inter
 
  # with interactions
  #mean_x = apply(subset(A0_S1_data, select = reg_after_match),2, mean)
  A0_S1_data$Y1_pred_inter = coeff_A_with_inter + intercept_with_inter +
    as.matrix(subset(A0_S1_data, select = reg_after_match)) %*% 
    (coeffs_with_inter_A1[-grep(":",names(coeffs_with_inter_A1))] + coeffs_with_inter_A1[grep(":",names(coeffs_with_inter_A1))])
  A0_S1_data$Y1_pred_inter_adj =  A0_S1_data$Y1_pred_inter / 
    ( ((1 - eps_sensi_PPI) * A0_S1_data$e_1_as) +  eps_sensi_PPI )
  
  A0_S1_data$Y0_pred_inter = as.matrix(subset(A0_S1_data, select=c("intercept", reg_after_match))) %*% coeffs_with_inter_A0
  
  # Y0_adj
  if(identific_xi==TRUE){ # FALSE and TRUE must obtain the same result, F is the same term as T, when expressing eta(or xi) with eps0_sensi_mono
    A0_S1_data$Y0_pred_adj = (1 + eta) * A0_S1_data$Y0_pred /
      ( (1 - eps0_sensi_mono) +  ( (1 + eta) * eps0_sensi_mono ) )
    A0_S1_data$Y0_pred_inter_adj = (1 + eta) * A0_S1_data$Y0_pred_inter /
      ( (1 - eps0_sensi_mono) +  ( (1 + eta) * eps0_sensi_mono ) )
  }else if(identific_xi==FALSE){
    A0_S1_data$Y0_pred_adj =  A0_S1_data$Y0_pred /
      ( ((1 - eps0_sensi_mono) * A0_S1_data$e_0_as) +  eps0_sensi_mono )
    A0_S1_data$Y0_pred_inter_adj =  A0_S1_data$Y0_pred_inter /
      ( ((1 - eps0_sensi_mono) * A0_S1_data$e_0_as) +  eps0_sensi_mono )
  }
  
  # wout inter
  SACE_1LEARNER = mean(A0_S1_data$Y1_pred) - mean(A0_S1_data$Y0_pred) 
  SACE_1LEARNER_adj = mean(A0_S1_data$Y1_pred_adj) - mean(A0_S1_data$Y0_pred_adj) 
  # with inter
  SACE_1LEARNER_inter = mean(A0_S1_data$Y1_pred_inter) - mean(A0_S1_data$Y0_pred_inter) 
  SACE_1LEARNER_inter_adj = mean(A0_S1_data$Y1_pred_inter_adj) - mean(A0_S1_data$Y0_pred_inter_adj)
  
  # crude diff estimator, adjusted by eps_PPI
  crude_Y1_adj = matched_data$Y[matched_data$A==1] / 
    ( ((1 - eps_sensi_PPI) * matched_data$e_1_as[matched_data$A==1]) +  eps_sensi_PPI )
  crude_Y0_adj = (1 + eta) * matched_data$Y[matched_data$A==0] /
    ( (1 - eps0_sensi_mono) +  ( (1 + eta) * eps0_sensi_mono ) )
  crude_est_adj = mean(crude_Y1_adj) - mean(crude_Y0_adj)
  
  return(c(SACE_1LEARNER_adj=SACE_1LEARNER_adj, SACE_1LEARNER_inter_adj=SACE_1LEARNER_inter_adj, crude_est_adj=crude_est_adj,
           SACE_1LEARNER=SACE_1LEARNER, SACE_1LEARNER_inter=SACE_1LEARNER_inter))
}
#########################################################################################

# 1. matching for O(0,1), using parameters (reg_after_match for instance) from main
# mahalanobis with PS caliper
#########################################################################################
#TODO \XI
#TODO MY upper bound for xi/eta
p1 = mean(data[A==1,S]); p0 = mean(data[A==0,S])
up_bound_xi = (1 - p1) / (p1 - 1 + p0)
len = 5; win_len = up_bound_xi / len
xi_sensi_mono_vec = seq(0, up_bound_xi, win_len) %>% round(2)
# arrange
xi_sensi_mono_vec = c(seq(0, 0.4, 0.1), 0.48)
xi_sensi_mono_vec[length(xi_sensi_mono_vec)] = round(up_bound_xi,2)
# check
last(xi_sensi_mono_vec) == up_bound_xi
# names
xi_sensi_mono_names = paste0("xi_mono_", round(xi_sensi_mono_vec, 2))

# eps0
eps0_sensi_mono_vec = seq(0.5, 2, 0.25) 
# SPPI: eps0_sensi_mono_vec = 1
# names
eps0_sensi_mono_names = paste0("eps0_PPI_",eps0_sensi_mono_vec)  
# use xi in the identification from proposition 3, which does not requires PSs estimation, or use e_{0,as} which does require
identific_xi = TRUE

reg_sensi_mono <- NULL
for (j in 1:length(xi_sensi_mono_names)) {
  xi_tmp = xi_sensi_mono_vec[j]
  ##########################################################################
  tmp = data # data = adjust_data(data, 1000, data_bool=data_bool)
  attach(tmp)
  
  # calculate PS for the current xi
  tmp$g = ifelse( A==0 & S==1, "as", ifelse( A==1 & S==0, "har", ifelse(A==1 & S==1, "pro", "pro") )  )
  est_ding_lst_SA_mono = xi_PSPS_M_weighting_SA(Z=tmp$A, D=tmp$S,
                                             X=as.matrix(subset(tmp, select = covariates_PS)),  
                                             Y=tmp$Y, eta=xi_tmp, # eta = 0 implies monotonicity
                                             beta.c = NULL, beta.n = NULL)
  DING_model_assisted_sensi = est_ding_lst_SA_mono$AACE.reg
  # c(prob.c, prob.d, prob.a, prob.n)
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
  
  # matching for the current xi
  #####################################
  matching_lst = matching_func_multiple_data(match_on = match_on,
             cont_cov_mahal = cont_cov_mahal,  reg_cov = reg_after_match, X_sub_cols = variables, 
             reg_BC = reg_BC, m_data = tmp[S==1], 
             w_mat_bool = "NON-INFO", M=1, replace=TRUE, estimand = "ATC", mahal_match = 2, caliper = caliper
             #, OBS_table = descrip_all_data_OBS$OBS_table
             , change_id=TRUE, boost_HL=FALSE, pass_tables_matched_units=FALSE, one_leraner_bool=TRUE)
  #m_dat = matching_lst$m_data; dt_match = matching_lst$ATE_MATCH_PS_lst
  
  #matched_data = matching_lst$matched_set
  #TODO matched_set_lst for all distance metrics
  matched_data_lst = matching_lst$matched_set_lst
  #TODO run on all distance mterics
  for(ind_matched_set in c(1:length(matched_data_lst))){  
    matched_data = matched_data_lst[[ind_matched_set]]
    print(paste0("unique weights for control are really = ", unique(filter(matched_data, A==0)$w)))
  
    #2. regression on the matched dataset with original Y
    coeffs_regression = regression_function(data_reg=matched_data, reg_after_match=reg_after_match)
    #########################################################################################
  
    for (i in 1:length(eps0_sensi_mono_names)) {
      print(eps0_sensi_mono_names[i])
      eps0 = eps0_sensi_mono_vec[i]; eps = 1
      
      # DING estimator
      # est_ding_lst = PSPS_M_weighting(Z=data_with_PS$A, D=data_with_PS$S, 
      #   X=as.matrix(subset(data_with_PS, select = covariates_PS)), Y=data_with_PS$Y,
      #   trc = TRUE, ep1 = eps0_sensi_mono_vec[i], ep0 = 1, beta.a = NULL, beta.n = NULL, iter.max = iterations , error0 = epsilon_EM)
      # #DING_est_sensi_PPI = est_ding_lst$AACE
      # DING_model_assisted_sensi_PPI = est_ding_lst$AACE.reg
      
      
      # ONE-LEARNER regression approach
      #3. predictions for units from O(0,1), plugging A={0,1} + 4. sensitivity adjustments
      reg_sensi_mono =
        rbind(reg_sensi_mono, c(measure = names(matched_data_lst)[ind_matched_set], xi_mono = xi_tmp, eps0_mono = eps0, 
      SACE_estimation_1LEARNER(matched_data=matched_data, coeffs_regression=coeffs_regression, reg_after_match=reg_after_match, 
                            eps_sensi_PPI=eps, eps0_sensi_mono=eps0, eta=xi_tmp, identific_xi=identific_xi)) )
      print(paste0(xi_sensi_mono_names[j], " ", eps0_sensi_mono_names[i]))
    }
  }
}

reg_sensi_mono = data.frame(reg_sensi_mono)
assign(paste0(data_bool, "_TEMP") ,reg_sensi_mono)
get(paste0(data_bool, "_TEMP"))
reg_sensi_mono[,-1] = apply(reg_sensi_mono[,-1] , 2, as.numeric)
# process for ggplot
reg_sensi_mono = subset(reg_sensi_mono, select = -c(SACE_1LEARNER, SACE_1LEARNER_inter)) %>% data.frame()
reg_sensi_mono[,-c(1,2,3)] = round(reg_sensi_mono[,-c(1,2,3)])
reg_sensi_mono$measure = factor(reg_sensi_mono$measure , levels = c("Mahal_PS_cal", "Mahal", "PS"))
reg_sensi_mono = reg_sensi_mono %>% gather("Estimator", "Estimate", 4:6) %>% arrange(measure, xi_mono, eps0_mono)
legend_levels = c("Crude", "WLS", "WLS inter"); 
reg_sensi_mono$Estimator = mgsub(reg_sensi_mono$Estimator,
                                c("crude_est_adj", "SACE_1LEARNER_adj", "SACE_1LEARNER_inter_adj"), legend_levels)
reg_sensi_mono$Estimator = factor(reg_sensi_mono$Estimator, levels = legend_levels)
reg_sensi_mono$set = data_bool

# combine DW LL
reg_sensi_mono = rbind(DW_sensi_mono_under_PPI, LL_sensi_mono_under_PPI)
#reg_sensi_mono = reg_sensi_mono %>% filter(Estimator != "WLS inter")
reg_sensi_mono$Estimator = factor(reg_sensi_mono$Estimator, levels = c("Crude", "WLS","WLS inter")) # "WLS inter"
save(reg_sensi_mono, file = "reg_sensi_mono.RData")
#reg_sensi_mono = filter(reg_sensi_mono, measure == "Mahal_PS_cal")

# only mahal with caliper
reg_sensi_mono = reg_sensi_mono[as.character(reg_sensi_mono$measure) == "Mahal_PS_cal",]
# exclude WLS inter
reg_sensi_mono = filter(reg_sensi_mono, !Estimator=="WLS inter")
#########################################################################################

#gsub("^0\\.", "\\.", gsub("^-0\\.", "-\\.", reg_sensi_mono$eps0_mono))
#TODO plot monotonicity under SPPI or PPI
#########################################################################################
plot_sensi_PPI <- ggplot(reg_sensi_mono, aes(x=eps0_mono, y=Estimate)) + # eps0_mono # xi_mono
  geom_point(aes(col = Estimator, size = 7), size = 2) + 
  scale_color_manual(values = c("Crude" = "green3", "WLS" = "orangered2", "WLS inter" = "cornflowerblue")) + # WLS inter = "blue" # "WLS inter" = "cornflowerblue"
  geom_line(aes(col = Estimator, size = 1.5), size=1.5) + 
  #xlim("0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75") +
  labs(colour = "Estimator"
       , size = 1
       #, title=data_bool
  ) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab(label="Estimate") + xlab(label = bquote(epsilon[0])) + # epsilon[0] # epsilon[PPI] # xi
  #labs( y="Estimate", x=glue('esp[PPI]*" : {protected}"')) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=5))
         , size=FALSE
  ) + 
  scale_x_continuous(breaks=c(.5, 1, 1.5, 2), labels = c(".5", "1", "1.5", "2")) +  # xi_sensi_mono_vec # breaks=c(.5, 1, 2)
  theme(
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16)) + 
  geom_hline(yintercept = 0 )


# plot_sensi_PPI_DW_LL # plot_sensi_PPI_DW # plot_sensi_PPI_LL
#  # facet_wrap(~ set, ncol=2) # glue('xi*" : {xi_mono}"') # xi_mono
plot_sensi_PPI_DW_LL = plot_sensi_PPI +
  # facet_wrap(~ set, ncol=2) +
    facet_grid(~ glue('xi*" = {xi_mono}"'), labeller = label_parsed) + # set ~ glue('xi*" = {xi_mono}"') # ~ glue('xi*" = {xi_mono}"') # facet_wrap(~ set, ncol=2) +
  
  # theme(strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=10, face="bold"),
  #       strip.background = element_rect(colour="black", fill="white")) +
  theme(
    #legend.direction = "horizontal",
    strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="black", fill="white"), 
    axis.title.x=element_text(size=14),  # X axis title
    axis.title.y=element_text(size=14),  # Y axis title
    axis.text.x=element_text(size=10),  # X axis text
    axis.text.y=element_text(size=10)
  ) 

# EXTRACT LEGEND
library(cowplot); library(ggpubr)
lgnd_plt <- get_legend(plot_sensi_PPI)
# Convert to a ggplot and print
as_ggplot(lgnd_plt)
plot_sensi_PPI_woutLGND = lgnd_plt + theme(legend.position = 'none') 

# ggarrange
library(ggpubr)
ggarrange(plot_sensi_PPI_DW, plot_sensi_PPI_LL)
#########################################################################################


#########################################################################################
reg_sensi_mono$measure = mgsub(as.character(reg_sensi_mono$measure), "_", " ")
reg_sensi_mono$measure = factor(reg_sensi_mono$measure, levels = c("Mahal", "Mahal PS cal", "PS"))
reg_sensi_mono_sppi = filter(reg_sensi_mono, eps0_mono==1)
reg_sensi_mono_sppi = filter(reg_sensi_mono_sppi, !Estimator=="WLS inter")

plot_sensi <- ggplot(reg_sensi_mono_sppi, aes(x=xi_mono, y=Estimate)) +
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

plot_sens_byMetric = plot_sensi + 
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

# plot_sensi_PPI_byMetric = plot_sensi_PPI + facet_wrap(~ Metric, ncol=3)

# EXTRACT LEGEND
library(cowplot); library(ggpubr)
lgnd_plt <- get_legend(plot_sens_byMetric) # plot_sensi_DW_LL_line # plot_sens_byMetric
# Convert to a ggplot and print
as_ggplot(lgnd_plt)
plot_sensi_woutLGND = plot_sens_byMetric + theme(legend.position = 'none') 
#########################################################################################