##################################################################################################
# CRUDE matching estimators ####
# read data: final_tables_crude ####
tables_small_pro3 = get(load(paste0(small_pro_path3, "final_tables_crude.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro3 = get(load(paste0(large_pro_path3, "final_tables_crude.RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro5 = get(load(paste0(small_pro_path5, "final_tables_crude.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro5 = get(load(paste0(large_pro_path5, "final_tables_crude.RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro10 = get(load(paste0(small_pro_path10, "final_tables_crude.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro10 = get(load(paste0(large_pro_path10, "final_tables_crude.RData"))) #tables_large_pro = get(tables_large_pro)


legend_levels = c("Crude Wout", "PS Crude Wout", "Mahal Crude Wout", 
                  "Crude With", "PS Crude With", "Mahal Crude With")   # "DingLu MA"
estimators_vec = c("Crude Wout", "PS Crude Wout", "Mahal Crude Wout", "Crude With", "PS Crude With", "Mahal Crude With")
#c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA")

# data.frame(tables_large_pro3$`_S1`%>%select(-k), k = 1)
small_pro = rbind(
  add_bias_tables(list( data.frame(tables_small_pro3$`_S1`, k = 1) ), 
                  estimators_vec = estimators_vec, 
                  N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro5$`_S1`, k = 2) ), 
                  estimators_vec = estimators_vec,
                  N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro10$`_S1`, k = 3) ), 
                  estimators_vec = estimators_vec,
                  N_obs = param_n, num_of_x = 10, legend_levels=legend_levels)
)
small_pro$protected = "Low"

large_pro = rbind(
  add_bias_tables(list( data.frame(tables_large_pro3$`_S1`, k = 1) ), 
                  estimators_vec = estimators_vec,
                  N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro5$`_S1`,  k = 2) ), 
                  estimators_vec = estimators_vec,
                  N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro10$`_S1`, k = 3) ), 
                  estimators_vec = estimators_vec,
                  N_obs = param_n, num_of_x = 10, legend_levels=legend_levels)
)
large_pro$protected = "High"

small_large_pro = rbind(small_pro, large_pro)
small_large_pro$Pi_as = mgsub(small_large_pro$Scenario, c("A", "B"), c("0.5", "0.75"))
small_large_pro = data.table(arrange(small_large_pro, protected, Pi_as, k))

small_large_pro[, label := paste0(unique(SACE), collapse=","), by = c("protected", "Pi_as")]
temp = apply(data.frame(list.rbind(strsplit(small_large_pro$label, ","))), 2 , as.numeric) %>% round(2)
small_large_pro$label = paste0("SACE: ", apply(temp, 1, function(x) paste(x, collapse=", ")))


# TODO plot ####
# estimator by color and shape
plot_PS <- 
  small_large_pro %>% filter(EstCombi %in% c("Crude With", "PS Crude With", "Mahal Crude With")) %>%
  ggplot(aes(x=k, y=Bias)) + theme_bw() + 
  geom_point(alpha = 0.65, size = 5, aes(col = EstCombi, shape = EstCombi)) +
  xlim("3", "5", "10") +
  xlab("Number of Covariates") +
  expand_limits(y=c(-1.3,0.5)) + 
  scale_colour_manual(name="", 
                      breaks = c("Crude With", "PS Crude With", "Mahal Crude With"),
                      labels = c("Mahalanobis with caliper", "PS", "Mahalanobis"),
                      values = c("palevioletred3","dodgerblue3", "forestgreen")) + # dodgerblue3 yellow  #c("forestgreen", "dodgerblue3",  "firebrick3","palevioletred3")
  scale_shape_manual(name="", 
                     breaks = c("Crude With", "PS Crude With", "Mahal Crude With"),
                     labels = c("Mahalanobis with caliper", "PS", "Mahalanobis"),
                     values = c(15, 16, 17)) +
  guides(col=guide_legend(nrow=1,byrow=TRUE), size=F)  +
  geom_hline(yintercept = 0 )


plot_PS_final = plot_PS + 
  facet_grid(glue('pi[pro]*" : {protected}"') ~ glue('pi[as]*" = {Pi_as}"'), labeller = label_parsed) +
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