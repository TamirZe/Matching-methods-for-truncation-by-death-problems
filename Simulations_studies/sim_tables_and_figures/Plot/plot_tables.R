library(cowplot); library(ggpubr)
##################################################################################################
# TODO A func that add bias for a FIXED number of N ####
# lst:  final_tables_general or final_tables_crude, CAN BE ONLY 1 data frame (e.g. S1)
# estimators_vec: vec of estimators name, repl and estimator together. If NULL, keep all estimators
# legend_levels: the factor level argument for the legend of the ggplot

func_add_AND_remove_COLS = function(fin_tab, N, num_of_x){
  fin_tab$Estimator = mgsub(fin_tab$Estimator, "mahal", "Mahal")
  fin_tab$Replacements = mgsub(fin_tab$Replacements, c("Yes", "No"), c("With", "Wout"))
  if(!"N" %in% colnames(fin_tab)) { fin_tab$N = N }
  if(!"dim_x" %in% colnames(fin_tab)) { fin_tab$dim_x = num_of_x }
  #if(!"k" %in% colnames(fin_tab)) { fin_tab = fin_tab%>%select(-k) }
  return(fin_tab)
}

add_bias_tables = function(lst, estimators_vec=NULL, N_obs, num_of_x,
                           legend_levels = c("Crude Wout", "OLS", "Crude With", "WLS", "BC With","DingLu MA")){
  new_lst = list()
  for(i in 1:length(lst)){
    tmp = lst[[i]]
    tmp = func_add_AND_remove_COLS(tmp, N = N_obs, num_of_x = num_of_x)
    colnames(tmp)[grep("parameter", colnames(tmp))] = "Set_and_parameter"
    tmp = data.frame("Set_and_parameter" = tmp[,1], N = tmp$N, dim_x = tmp$dim_x, 
                     subset(tmp, select = !colnames(tmp) %in% c("Set_and_parameter", "N", "dim_x")) )
    # tmp[,1] IS tmp$`Set & parameter`
    tmp$Bias = tmp$Mean - as.numeric(substr(tmp[,1],3,10))  
    tmp = data.frame(Scenario = substr(tmp[,1],1,1), SACE = as.numeric(substr(tmp[,1],3,10)) %>% round(10), tmp)
    tmp$EstCombi = paste0(tmp$Estimator, " ", tmp$Replacements)
    tmp$EstCombi[tmp$EstCombi == "DingLu MA "] = "DingLu MA"
    if(length(unique(tmp$Scenario)) == 3){ # unique(tmp$Scenario) == c("A", "B", "C")
      tmp = tmp[tmp$Scenario!="A" ,]
      tmp$Scenario = mgsub(tmp$Scenario, c("B", "C"), c("A", "B"))
    }
    if(!is.null(estimators_vec)){
      tmp = filter(tmp, EstCombi %in% estimators_vec)
      # TODO make this levels more general, so we could consider all estimators,
      # probably the level should be an argument
      tmp$EstCombi = mgsub(tmp$EstCombi, c("OLS Wout", "OLS inter Wout", "WLS With", "WLS inter With", "BC caliper With"),
                           c("OLS", "OLS inter", "WLS", "WLS inter", "BC With"))
      tmp$EstCombi = factor(tmp$EstCombi, levels = legend_levels)
      #TODO order the cols as we want them
    }
    new_lst[[i]] = tmp
  }
  names(new_lst) = names(lst)
  if(length(new_lst)==1) new_lst = new_lst[[1]]
  return(new_lst)
}
##################################################################################################

##################################################################################################
# GENERAL matching estimators ####
# read data: matching on mahalanobis with PS caliper - several estimators (final_tables_general)
tables_small_pro3 = get(load(paste0(small_pro_path3, "final_tables_general_", ind,".RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro3 = get(load(paste0(large_pro_path3, "final_tables_general_", ind,".RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro5 = get(load(paste0(small_pro_path5, "final_tables_general_", ind,".RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro5 = get(load(paste0(large_pro_path5, "final_tables_general_", ind,".RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro10 = get(load(paste0(small_pro_path10, "final_tables_general_", ind,".RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro10 = get(load(paste0(large_pro_path10, "final_tables_general_", ind,".RData"))) #tables_large_pro = get(tables_large_pro)

legend_levels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA")  # WLS or WLS inter
estimators_vec = c("Crude Wout", "OLS inter Wout", "Crude With", "WLS inter With", "BC caliper With","DingLu MA")

# data.frame(tables_small_pro3$`_S1`%>%select(-k) , k = 1)
small_pro = rbind(
  add_bias_tables(list( data.frame(tables_small_pro3$`_S1`, k = 1) ), 
    estimators_vec = estimators_vec, N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro5$`_S1`, k = 2) ), 
    estimators_vec = estimators_vec, N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro10$`_S1`, k = 3) ), 
    estimators_vec = estimators_vec, N_obs = param_n, num_of_x = 10, legend_levels=legend_levels)
)
small_pro$protected = "Low"

large_pro = rbind(
  add_bias_tables(list( data.frame(tables_large_pro3$`_S1`, k = 1) ), 
    estimators_vec = estimators_vec, N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro5$`_S1`,  k = 2) ), 
    estimators_vec = estimators_vec, N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro10$`_S1`, k = 3) ), 
    estimators_vec = estimators_vec, N_obs = param_n, num_of_x = 10, legend_levels=legend_levels)
)
large_pro$protected = "High"

small_large_pro = rbind(small_pro, large_pro)
small_large_pro$Pi_as = mgsub(small_large_pro$Scenario, c("A", "B"), c("0.5", "0.75"))
small_large_pro = data.table(arrange(small_large_pro, protected, Pi_as, k))

# add label for the SACE, for each facet plot
small_large_pro[, label := paste0(unique(SACE), collapse=","), by = c("protected", "Pi_as")]
temp = apply(data.frame(list.rbind(strsplit(small_large_pro$label, ","))), 2 , as.numeric) %>% round(2)
small_large_pro$label = paste0("SACE: ", apply(temp, 1, function(x) paste(x, collapse=", ")))

# plot
# http://www.cookbook-r.com/Graphs/Facets_(ggplot2)/ # https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels # http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html

# estimator by color, emp sd by size
plot_general <- ggplot(small_large_pro, aes(x=k, y=Bias)) +
  geom_point(aes(col = EstCombi, size = 7)) + # shape = EstCombi
  xlim("3", "5", "10") +
  labs(colour = "Estimator"
       #, size = "Emp SD"
  ) + 
  guides(colour = guide_legend(order = 1, override.aes = list(size=7))
         , size=FALSE
  ) + 
  geom_hline(yintercept = 0 )

plot_general = plot_general + scale_color_manual(name="Estimator", 
   labels = legend_levels, 
    values = c("Crude Wout" = "forestgreen", "OLS inter" = "dodgerblue3",
        "Crude With" = "yellow1", "WLS inter" = "firebrick3", "BC With" = "palevioletred3","DingLu MA" = "black"))  +
  facet_grid(glue('pi[pro]*" : {protected}"') ~ glue('pi[as]*" = {Pi_as}"'), labeller = label_parsed) +
  theme(
    strip.text.x = element_text(size=16, face="bold"),
    strip.text.y = element_text(size=16, face="bold"),
    strip.background = element_rect(colour="black", fill="white")
  ) + 
  labs(y="Bias", x="k") +
  theme(plot.title=element_text(size=15, color="black", hjust=0.5,
                                lineheight=1.2),  # title
  plot.subtitle=element_text(size=15, 
                             family="American Typewriter",
                             face="bold",
                             hjust=0.5),  # subtitle
  plot.caption=element_text(size=10),  # caption
  axis.title.x=element_text(
    size=22),  # X axis title
  axis.title.y=element_text(size=22),  # Y axis title
  axis.text.x=element_text(size=15, 
                           angle = 30,
                           vjust=.5),  # X axis text
  axis.text.y=element_text(size=15)) # Y axis text

# Extract the legend. Returns a gtable
lgnd_general <- get_legend(plot_general)
# Convert to a ggplot and print
as_ggplot(lgnd_general)
plot_general = plot_general + theme(legend.position = 'none') 
##################################################################################################

##################################################################################################
# crude matching estimators ####
# read data: final_tables_crude
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


# TODO plot
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

