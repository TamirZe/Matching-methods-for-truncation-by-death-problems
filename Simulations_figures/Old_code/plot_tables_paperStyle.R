library(cowplot); library(ggpubr); library(arsenal); library(tidyverse)

param_n=2000
mis_xi = 2
AX_interactions = T
misspec_outcome = 0 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_PS = 0 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
xi_values = c(0, 0.05, 0.1, 0.2) # sort(unique(full_results_table$xi)) # c(0, 0.05, 0.1, 0.2)

estimators_vec = c("maha_cal_rep_FALSE", "OLS_int", "maha_cal_rep_TRUE", "WLS_int", "BC_cal_rep_TRUE", "DL_MA_est")
legend_levels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA")

full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
      AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, 
      estimators_vec=estimators_vec, legend_levels=legend_levels)
small_large_pro = full_results_table %>% filter(xi==0.2)

#legend_levels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA")  # WLS or WLS inter
small_large_pro$shape = mgsub(small_large_pro$Estimator, c("OLS|WLS|BC", " inter| caliper"), c("Model-based", ""))
small_large_pro$SE_rel_bias = (small_large_pro$SE - small_large_pro$sd) / small_large_pro$sd

# TODO plot # ggpubr::show_point_shapes()
# Bias  ####
plot_tables_func_by_dimx_paperStyle(small_large_pro, param_n, mis_xi, AX_interactions, misspec_outcome, misspec_PS, 
                                  estimators_vec, legend_levels=NULL)



# WITH and WOUT repl with SE ####
small_large_pro %>%
  ggplot(aes(x = l, y = Bias)) + theme_bw() +
  geom_point(alpha = 0.65, aes(col = EstCombi, shape = EstCombi, size = sd)) + # ,size = #sd #MSE # Bias SE_rel_bias, # shape = as.character(shape)
  xlim("3", "5", "10") +
  xlab("Number of Covariates") + 
  labs(colour = "Estimator", shape = as.character("shape")) + 
  scale_colour_manual(name="",
                      breaks = c("DingLu MA", "Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      labels = c("DingLu MA", "Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      values = c("black", "forestgreen","dodgerblue3", "yellow", "blanchedalmond", "red4")) + # dodgerblue3 yellow  #c("forestgreen", "dodgerblue3",  "firebrick3","palevioletred3")
  scale_shape_manual(name="",
                     breaks = c("DingLu MA", "Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     labels = c("DingLu MA", "Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     values = c(20, 15, 16, 17, 18, 19)) +
  guides(col=guide_legend(nrow=1,byrow=TRUE), size=F) + 
  #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
  geom_hline(yintercept = 0) + 
  ylim(c(-0.75, 0.75)) + 
  facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {pi_as}"'), 
             labeller = label_parsed) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.575, 'cm'),
        axis.text = element_text(size  = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="white"))       
  
# SE relative Bias  ####
small_large_pro %>% filter(EstCombi != "DingLu MA") %>%
  ggplot(aes(x = l, y = SE_rel_bias)) + theme_bw() +
  geom_point(alpha = 0.65, size = 5, aes(col = EstCombi, shape = EstCombi,size = Bias)) + # ,size = Bias SE_rel_bias# , shape = as.character(shape)
  xlim("3", "5", "10") +
  xlab("Number of Covariates") + 
  labs(colour = "Estimator", shape = as.character("shape")) + 
  scale_colour_manual(name="",
                      breaks = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      labels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                      values = c("forestgreen","dodgerblue3", "yellow", "palevioletred3", "red4")) + # dodgerblue3 yellow  #c("forestgreen", "dodgerblue3",  "firebrick3","palevioletred3")
  scale_shape_manual(name="",
                     breaks = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     labels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With"),
                     values = c(15, 16, 17, 18, 19)) +
  guides(col=guide_legend(nrow=2,byrow=TRUE)) + 
  #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
  geom_hline(yintercept = 0) + 
  ylim(c(-0.75, 0.75)) + 
  facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {pi_as}"'), 
             labeller = label_parsed) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.575, 'cm'),
        axis.text = element_text(size  = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="white")) 















