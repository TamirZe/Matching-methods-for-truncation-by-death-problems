# TODO A func that add bias for a FIXED number of N
# lst:  final_tables_general or final_tables_crude, CAN BE ONLY 1 DF (e.g. S1)
# estimators_vec: vec of estimators name, repl and estimator together. If NULL, keep all estimators
# legend_levels: the factor level argument for the legend of the ggplot

lst = list( data.frame(tables_large_pro5$`_S1`%>%select(-k),  k = 2) )
lst = llist( data.frame(tables_small_pro3$`_S1`, k = 1) )
func_add_AND_remove_COLS = function(fin_tab, N, num_of_x){
  if(!"N" %in% colnames(fin_tab)) { fin_tab$N = N }
  if(!"dim_x" %in% colnames(fin_tab)) { fin_tab$dim_x = num_of_x }
  #if(!"k" %in% colnames(fin_tab)) { fin_tab = fin_tab%>%select(-k) }
  return(fin_tab)
}
add_bias_tables = function(lst, estimators_vec=NULL, N_obs, num_of_x,
                           legend_levels = c("Crude No", "OLS", "Crude Yes", "WLS", "BC Yes","DingLu MA")){
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
      tmp$EstCombi = mgsub(tmp$EstCombi, c("OLS No", "OLS inter No", "WLS Yes", "WLS inter Yes"),
                           c("OLS", "OLS inter", "WLS", "WLS inter"))
      tmp$EstCombi = factor(tmp$EstCombi, levels = legend_levels)
      #TODO order the cols as we want them
    }
    new_lst[[i]] = tmp
  }
  names(new_lst) = names(lst)
  if(length(new_lst)==1) new_lst = new_lst[[1]]
  return(new_lst)
}

misspec_title = "ff misspec" # no misspec" "" # "U- S&Y misspec" # "ff misspec"
caliper = 0.25
##################################################################################################
# caliper 0.25 caliper 0.05
#TODO path to final tables,  
#TODO no misspec
# 3X
small_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/no misspec/small pro/3X/"
large_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/no misspec/large pro/3X/"
# 5X
small_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/no misspec/small pro/5X/"
large_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/no misspec/large pro/5X/"
# 10X
small_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/no misspec/small pro/10X/"
large_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/no misspec/large pro/10X/"

#TODO U- S&Y

#TODO Func Form
# 3X

small_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/Func Form/mis2 ff small pro -3 3/3X/"
large_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/Func Form/mis2 ff large pro -3 3/3X/"
# 5X
small_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/Func Form/mis2 ff small pro -3 3/5X/"
large_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/Func Form/mis2 ff large pro -3 3/5X/"
# 10X
small_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/Func Form/mis2 ff small pro -3 3/10X/"
large_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/caliper 0.25/model with interaction/Func Form/mis2 ff large pro -3 3/10X/"

##################################################################################################

##################################################################################
# TODO GENERAL matching estimators
# read data: final_tables_general
tables_small_pro3 = get(load(paste0(small_pro_path3, "final_tables_general.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro3 = get(load(paste0(large_pro_path3, "final_tables_general.RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro5 = get(load(paste0(small_pro_path5, "final_tables_general.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro5 = get(load(paste0(large_pro_path5, "final_tables_general.RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro10 = get(load(paste0(small_pro_path10, "final_tables_general.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro10 = get(load(paste0(large_pro_path10, "final_tables_general.RData"))) #tables_large_pro = get(tables_large_pro)


# TODO datasets with different k
# need to adjust k, the covariates, because of the ggplot. need to check that each k=1,2,3 corresponds to th coresponding real k.
# k = (dim_x-1)

legend_levels = c("Crude No", "OLS inter", "Crude Yes", "WLS inter", "BC Yes", "DingLu MA")  # WLS or WLS inter

# data.frame(tables_small_pro3$`_S1`%>%select(-k) , k = 1)
small_pro = rbind(
  add_bias_tables(list( data.frame(tables_small_pro3$`_S1`, k = 1) ), 
    estimators_vec = c("Crude No", "OLS inter No", "Crude Yes", "WLS inter Yes", "BC Yes","DingLu MA"), 
    N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro5$`_S1`, k = 2) ), 
    estimators_vec = c("Crude No", "OLS inter No", "Crude Yes", "WLS inter Yes", "BC Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro10$`_S1`, k = 3) ), 
    estimators_vec = c("Crude No", "OLS inter No", "Crude Yes", "WLS inter Yes", "BC Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 10, legend_levels=legend_levels)
)
#'mu[2]*" = 0.1"' 
#small_pro$protected = '"Low "pi[pro]*'
small_pro$protected = "Low"

large_pro = rbind(
  add_bias_tables(list( data.frame(tables_large_pro3$`_S1`, k = 1) ), 
    estimators_vec = c("Crude No", "OLS inter No", "Crude Yes", "WLS inter Yes", "BC Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro5$`_S1`,  k = 2) ), 
    estimators_vec = c("Crude No", "OLS inter No", "Crude Yes", "WLS inter Yes", "BC Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro10$`_S1`, k = 3) ), 
    estimators_vec = c("Crude No", "OLS inter No", "Crude Yes", "WLS inter Yes", "BC Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 10, legend_levels=legend_levels)
)
#large_pro$protected = '"High "pi[pro]*'
large_pro$protected = "High"

small_large_pro = rbind(small_pro, large_pro)
small_large_pro$Pi_as = mgsub(small_large_pro$Scenario, c("A", "B"), c("0.5", "0.75"))
# mgsub(small_large_pro$Scenario, c("A", "B", "C"), c("pi as = 0.25", "pi as = 0.5", "pi as = 0.75"))
#small_large_pro$Pi_as = mgsub(small_large_pro$Scenario, c("A", "B"), c('pi[as]*" = 0.5"', 'pi[as]*" = 0.75"'))
small_large_pro = data.table(arrange(small_large_pro, protected, Pi_as, k))

# add label for the SACE, for each facet plot
small_large_pro[, label := paste0(unique(SACE), collapse=","), by = c("protected", "Pi_as")]
#small_large_pro[, lapply(.SD, toString), by = c("protected", "Pi_as")]

# small_large_pro$label = paste0("SACE: ", as.numeric(sub(",.*", "", small_large_pro$label)) %>% round(2), ", ",
#                                as.numeric(sub(".*,", "", small_large_pro$label)) %>% round(2))
temp = apply(data.frame(list.rbind(strsplit(small_large_pro$label, ","))), 2 , as.numeric) %>% round(2)
small_large_pro$label = paste0("SACE: ", apply(temp, 1, function(x) paste(x, collapse=", ")))


# TODO plot
# http://www.cookbook-r.com/Graphs/Facets_(ggplot2)/
# http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
# https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels/

# sp <- ggplot(small_large_pro, aes(x=k, y=Bias, colour = factor(EstCombi))) + geom_point(shape=1)

# estimator by color, emp sd by size
plot_crude <- ggplot(small_large_pro, aes(x=k, y=Bias)) +
  geom_point(aes(col = EstCombi, size = Emp.SD)) + 
  #geom_point(aes(col = EstCombi, size = Emp.SD)) + # size by SD
  #scale_size_continuous(name = "Emp SD",range = c(3, 8)) +
  xlim("3", "5", "10") +
  labs(title="Bias of different estimators, N = 2000", y="Bias", x="k", caption=paste0("caliper = ", caliper)) +
  labs(colour = "Estimator"
       , size = "Emp SD"
       ) + 
  scale_color_discrete(name="Estimator")  +
  scale_size_continuous(name = "Emp SD",range = c(3, 8)) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=7))
         #, size=FALSE
         ) + 
  geom_hline(yintercept = 0 )
#geom_hline(yintercept = as.numeric(substr(tmp$Set...parameter,3,10)) )
#theme(legend.text=element_text(size=20))

# WLS or WLS inter

#facet_grid(glue('pi[pro]*" = {protected}"') ~ glue('pi[as]*" = {Pi_as}"'), labeller = label_parsed)
#facet_grid(protected ~ Pi_as, labeller = label_parsed)
plot_crude = plot_crude + scale_color_manual(name="Estimator", 
       labels = legend_levels, 
       values = c("Crude No" = "forestgreen", "OLS inter" = "dodgerblue3",
                  "Crude Yes" = "yellow1", "WLS inter" = "firebrick3",
                  "BC Yes" = "palevioletred3","DingLu MA" = "black"))  +
  facet_grid(glue('pi[pro]*" : {protected}"') ~ glue('pi[as]*" = {Pi_as}"'), labeller = label_parsed) +
  theme(strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="white")) + # fill="#CCCCFF"
  labs(title=paste0("Mahalanobis and PS caliper, N = 2000 ", misspec_title), y="Bias", x="k", caption=paste0("caliper = ", caliper)) +
  theme(plot.title=element_text(size=15, color="black", hjust=0.5,
                                lineheight=1.2
                                #, face="bold"
  ),  # title
  plot.subtitle=element_text(size=15, 
                             family="American Typewriter",
                             face="bold",
                             hjust=0.5),  # subtitle
  plot.caption=element_text(size=10),  # caption
  axis.title.x=element_text(
                            #vjust=10,  
                            size=22),  # X axis title
  axis.title.y=element_text(size=22),  # Y axis title
  axis.text.x=element_text(size=15, 
                           angle = 30,
                           vjust=.5),  # X axis text
  axis.text.y=element_text(size=15))  + # Y axis text
  geom_label(data = small_large_pro, aes(label=label), label.size = 0.7 # add SACE values per each facet cell
           ,x = Inf, y = -Inf
           ,hjust=1, vjust=0
           #,hjust="bottom", vjust="middle"
           ,inherit.aes = FALSE)



# TODO estimator by shape, emp sd by color

# cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# + scale_colour_manual(values=cbp1)

plt_crude_sd_col = ggplot(data = small_large_pro, aes(x=k, y=Bias))+
  geom_point(aes(size = 5, color = Emp.SD, shape = factor(EstCombi)))+
  guides(colour = guide_colourbar(order = 2),
         alpha = guide_legend(order = 1),
         size = guide_legend(order = 3)) +
  guides(size = FALSE) + 
  labs(shape="Estimator", colour="Emp SD") + 
  xlim("3", "5", "10") +
  labs(title=paste0("Mahalanobis and PS caliper, N = 2000, ", misspec_title), y="Bias", x="k", caption=paste0("caliper = ", caliper)) +
  #scale_color_discrete(name="Estimator") +
  guides(shape = guide_legend(order = 1, override.aes = list(size=7))) + 
  scale_colour_gradient(low = "red", high = "yellow", na.value = NA) + 
  #scale_color_gradientn(colours = rainbow(3)) + 
  #scale_color_brewer(palette = "Dark2") + 
  geom_hline(yintercept = 0 )

plt_crude_sd_col = plt_crude_sd_col + 
  #scale_shape_manual(values=c(19, 15, 13, 7, 17, 8)) + 
  scale_shape_manual(values=c(19, 17, 15, 25, 18, 8)) +
  facet_grid(protected ~ Pi_as) +
  theme(strip.text.x = element_text(size=16, face="bold"
                                    #, angle=75
  ),
  strip.text.y = element_text(size=12, face="bold"),
  strip.background = element_rect(colour="black", fill="white"))

##################################################################################




##################################################################################
# TODO crude matching estimators
# read data: final_tables_crude
tables_small_pro3 = get(load(paste0(small_pro_path3, "final_tables_crude.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro3 = get(load(paste0(large_pro_path3, "final_tables_crude.RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro5 = get(load(paste0(small_pro_path5, "final_tables_crude.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro5 = get(load(paste0(large_pro_path5, "final_tables_crude.RData"))) #tables_large_pro = get(tables_large_pro)
tables_small_pro10 = get(load(paste0(small_pro_path10, "final_tables_crude.RData"))) #tables_small_pro = get(tables_small_pro)
tables_large_pro10 = get(load(paste0(large_pro_path10, "final_tables_crude.RData"))) #tables_large_pro = get(tables_large_pro)

legend_levels = c("PS Crude No", "OLS", "PS Crude Yes", "WLS","DingLu MA")  
#c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA")

# data.frame(tables_large_pro3$`_S1`%>%select(-k), k = 1)
small_pro = rbind(
  add_bias_tables(list( data.frame(tables_small_pro3$`_S1`%>%select(-k), k = 1) ), 
    estimators_vec =c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA"), 
    N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro5$`_S1`, k = 2) ), 
    estimators_vec = c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_small_pro10$`_S1`, k = 3) ), 
    estimators_vec = c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 10, legend_levels=legend_levels)
)
small_pro$protected = "Low"

large_pro = rbind(
  add_bias_tables(list( data.frame(tables_large_pro3$`_S1`%>%select(-k), k = 1) ), 
    estimators_vec = c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 3, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro5$`_S1`,  k = 2) ), 
    estimators_vec = c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA"),
    N_obs = param_n, num_of_x = 5, legend_levels=legend_levels),
  add_bias_tables(list( data.frame(tables_large_pro10$`_S1`, k = 3) ), 
    estimators_vec = c("PS Crude No", "OLS No", "PS Crude Yes", "WLS Yes","DingLu MA"),
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
# estimator by color, emp sd by size
plot_PS <- ggplot(small_large_pro, aes(x=k, y=Bias)) +
  geom_point(aes(col = EstCombi, size = Emp.SD)) + 
  xlim("3", "5", "10") +
  labs(title="PS crude, N = 2000", y="Bias", x="k", caption=paste0("caliper = ", caliper)) +
  labs(colour="Estimator", size = "Emp SD") + 
  scale_color_discrete(name="Estimator")  +
  scale_size_continuous(name = "Emp SD",range = c(3, 8)) +
  guides(colour = guide_legend(order = 1, override.aes = list(size=7))) + 
  geom_hline(yintercept = 0 )
#geom_hline(yintercept = as.numeric(substr(tmp$Set...parameter,3,10)) )
#theme(legend.text=element_text(size=20))

plot_PS = plot_PS + scale_color_manual(name="Estimator", 
                                       labels = legend_levels, 
                                       values = c("PS Crude No" = "forestgreen", "OLS" = "dodgerblue3",
                                                  "PS Crude Yes" = "yellow1", "WLS" = "firebrick3",
                                                  "DingLu MA" = "black"))  +
  facet_grid(protected ~ Pi_as) +
  theme(strip.text.x = element_text(size=16, face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black", fill="white")) +
  labs(title=paste0("PS crude, N = 2000, ", misspec_title), y="Bias", x="k", caption=paste0("caliper = ", caliper)) +
  theme(plot.title=element_text(size=14, color="black", hjust=0.5, lineheight=1.2
                                #,face="bold"
                                ),  # title
        plot.subtitle=element_text(size=15, 
                                   family="American Typewriter",
                                   face="bold",
                                   hjust=0.5),  # subtitle
        plot.caption=element_text(size=10),  # caption
        axis.title.x=element_text(size=15),  # X axis title
        axis.title.y=element_text(size=15),  # Y axis title
        axis.text.x=element_text(size=15, 
                                 angle = 30,
                                 vjust=.5),  # X axis text
        axis.text.y=element_text(size=15))  # Y axis text


# TODO estimator by shape, emp sd by color

# cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# + scale_colour_manual(values=cbp1)

plt_PS_sd_col = ggplot(data = small_large_pro, aes(x=k, y=Bias))+
  geom_point(aes(size = 5, color = Emp.SD, shape = factor(EstCombi)))+
  guides(colour = guide_colourbar(order = 2),
         alpha = guide_legend(order = 1),
         size = guide_legend(order = 3)) +
  guides(size = FALSE) + 
  labs(shape="Estimator", colour="Emp SD") + 
  xlim("3", "5", "10") +
  labs(title=paste0("PS crude, N = 2000, ", misspec_title), y="Bias", x="k", caption="N = 2000") +
  #scale_color_discrete(name="Estimator") +
  guides(shape = guide_legend(order = 1, override.aes = list(size=7))) + 
  scale_colour_gradient(low = "red", high = "yellow", na.value = NA) + 
  #scale_color_gradientn(colours = rainbow(3)) + 
  #scale_color_brewer(palette = "Dark2") + 
  geom_hline(yintercept = 0 )

plt_PS_sd_col = plt_PS_sd_col + 
  scale_shape_manual(values=c(19, 17, 15, 25, 8)) +
  facet_grid(protected ~ Pi_as) +
  theme(strip.text.x = element_text(size=16, face="bold"
                                    #, angle=75
  ),
  strip.text.y = element_text(size=12, face="bold"),
  strip.background = element_rect(colour="black", fill="white"))
##################################################################################

grid.arrange(plot_crude, plt_crude_sd_col, plot_PS, plt_PS_sd_col, nrow=2)

# change names of misspec plots, to plot both no misspec and ff misspec
plot_crude_no_mis = plot_crude; plot_PS_no_mis = plot_PS
plot_crude_mis = plot_crude;  plt_crude_sd_col_mis = plt_crude_sd_col
plot_PS_mis = plot_PS; plt_PS_sd_col_mis = plt_PS_sd_col

# estimator by color
grid.arrange(plot_crude_no_mis, plot_crude_mis, plot_PS_no_mis, plot_PS_mis, nrow=2)
# estimator by shape
grid.arrange(plt_crude_sd_col, plt_crude_sd_col_mis, plt_PS_sd_col, plt_PS_sd_col_mis, nrow=2)



