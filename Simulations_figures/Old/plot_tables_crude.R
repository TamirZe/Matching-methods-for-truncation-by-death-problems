##################################################################################################
# estimators_vec: vec of estimators name, repl and estimator together. If NULL, keep all estimators
# legend_levels: the factor level argument for the legend of the ggplot

param_n=2000
mis_xi = 2
AX_interactions = T
misspec_outcome = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_PS = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
xi_values = c(0, 0.05, 0.1, 0.2) 

estimators_vec = c("PS_rep_FALSE", "maha_rep_FALSE", "maha_cal_rep_FALSE", "PS_rep_TRUE", "maha_rep_TRUE", "maha_cal_rep_TRUE", "DL_MA_est") # "DL_MA_est"
legend_levels = c("PS Crude With", "Mahal Crude With", "Crude With", "DingLu MA")   # "DingLu MA"
#legend_levels = c("PS Crude Wout", "Mahal Crude Wout", "Crude Wout", "PS Crude With", "Mahal Crude With", "Crude With", "DingLu MA")   # "DingLu MA"
legend_LABELS = c("PS", "Mahalanobis", "Mahalanobis with caliper", "Weighting")
colors_arg = c("palevioletred3","dodgerblue3", "red4", "forestgreen")
shapes_arg = c(15, 16, 17, 18)

full_results_table = combine_small_large_pro_func(param_n=param_n, xi_values=xi_values, mis_xi=mis_xi,
                                                  AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, 
                                                  estimators_vec=estimators_vec, legend_levels=legend_levels)
small_large_pro = full_results_table %>% filter(xi==0)

##################################################################################################
# TODO plot ####
# estimator by color and shape
plot_PS <- 
  small_large_pro %>% filter(EstCombi %in% legend_levels) %>%
  ggplot(aes(x=l, y=Bias)) + theme_bw() + 
  geom_point(alpha = 0.65, size = 5, aes(col = EstCombi, shape = EstCombi)) +
  xlim("3", "5", "10") +
  xlab("Number of Covariates") +
  expand_limits(y=c(-1.3,0.5)) + 
  scale_colour_manual(name="", breaks = legend_levels, labels = legend_LABELS, values = colors_arg) + # dodgerblue3 yellow  #c("forestgreen", "dodgerblue3",  "firebrick3","palevioletred3")
  scale_shape_manual(name="", breaks = legend_levels, labels = legend_LABELS, values = shapes_arg) +
  guides(col=guide_legend(nrow=1,byrow=TRUE), size=F)  +
  geom_hline(yintercept = 0 )


plot_PS = plot_PS + 
  facet_grid(glue('pi[pro]*" : {protected}"') ~ glue('pi[as]*" = {pi_as}"'), labeller = label_parsed) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.575, 'cm'),
        axis.text = element_text(size  = 12),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="white")) +
  theme(plot.title=element_text(size=14, color="black", hjust=0.5, lineheight=1.2
                                #,face="bold"
  ),  # title
  plot.subtitle=element_text(size=15, 
                             family="American Typewriter",
                             face="bold",
                             hjust=0.5),  # subtitle
  plot.caption=element_text(size=10),  # caption
  axis.title.x=element_text(size=22),  # X axis title
  axis.title.y=element_text(size=15),  # Y axis title
  axis.text.x=element_text(size=15, 
                           angle = 30,
                           vjust=.5),  # X axis text
  axis.text.y=element_text(size=15)) # Y axis text

# Extract the legend. Returns a gtable
lgnd_PS <- get_legend(plot_PS)
# Convert to a ggplot and print
as_ggplot(lgnd_PS)
plot_PS = plot_PS + theme(legend.position = 'none') 
##################################################################################################