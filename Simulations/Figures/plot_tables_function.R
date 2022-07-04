library(cowplot); library(ggpubr); library(textclean)

##################################################################################################
# path to results_tables ####
path_func = function(param_n, xi, AX_interactions, misspec_outcome, misspec_PS){
# Cluster/GH
#main_path = "/a/home/cc/stud_math/tamirzehavi/MatchingSACE/Simulation_studies/"
#path_data = paste0(main_path, "Data/")

# local
main_path = "~/A matching framework for truncation by death problems/"
path_data = paste0(main_path, "Data_cluster/")

path = paste0(path_data, "N=", param_n, 
              ifelse(AX_interactions==T, "/True_outcome_with_interactions/", "/True_outcome_wout_interactions/"),
              ifelse(misspec_outcome==0, "Correct_spec_outcome/", "Mis_spec_outcome/"),
              ifelse(misspec_PS==0, "Correct_spec_PS/", "Mis_spec_PS/"))
all_path_lst = list(
  small_pro_small_as_path3 = paste0(path, "3X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/"),
  small_pro_large_as_path3 = paste0(path, "3X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/"),
  large_pro_small_as_path3 = paste0(path, "3X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/"),
  large_pro_large_as_path3 = paste0(path, "3X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/"),
  small_pro_small_as_path5 = paste0(path, "5X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/"),
  small_pro_large_as_path5 = paste0(path, "5X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/"),
  large_pro_small_as_path5 = paste0(path, "5X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/"),
  large_pro_large_as_path5 = paste0(path, "5X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/"),
  small_pro_small_as_path10 = paste0(path, "10X/Low_pi_pro/pi_as_0.5/", "xi = ", xi, "/"),
  small_pro_large_as_path10 = paste0(path, "10X/Low_pi_pro/pi_as_0.75/", "xi = ", xi, "/"),
  large_pro_small_as_path10 = paste0(path, "10X/Large_pi_pro/pi_as_0.5/", "xi = ", xi, "/"),
  large_pro_large_as_path10 = paste0(path, "10X/Large_pi_pro/pi_as_0.75/", "xi = ", xi, "/")
)
  return(all_path_lst=all_path_lst)
  
}
##################################################################################################

##################################################################################################
# estimators_vec: vec of estimators name
# legend_levels: the factor level argument for the legend of the ggplot
add_bias_tables = function(res_tab, 
                           estimators_vec=NULL, legend_levels = c("Crude Wout", "OLS", "Crude With", "WLS", "BC With","DingLu MA"),
                           N, num_of_x){
  res_tab = data.frame(Estimator = rownames(res_tab), res_tab)
  #res_tab$Estimator = mgsub(res_tab$Estimator, "maha", "Mahal")
  if(!"N" %in% colnames(res_tab)) { res_tab$N = N }
  if(!"dim_x" %in% colnames(res_tab)) { res_tab$dim_x = num_of_x }
  res_tab$true_SACE = res_tab$mean[res_tab$Estimator == "true_SACE"] 
  res_tab = data.frame("true_SACE" = res_tab$true_SACE, N = res_tab$N, dim_x = res_tab$dim_x,  l = res_tab$l, pi_as = res_tab$pi_as,
                       subset(res_tab, select = !colnames(res_tab) %in% c("true_SACE", "N", "l", "dim_x", "pi_as")) )
  res_tab$Bias = res_tab$mean - as.numeric(res_tab$true_SACE)  
  res_tab$Rel_bias = res_tab$Bias / as.numeric(res_tab$true_SACE) 
  res_tab$true_SACE = round(res_tab$true_SACE, 4)
  
  # tmp$EstCombi = paste0(tmp$Estimator, " ", tmp$Replacements)
  # tmp$EstCombi[tmp$EstCombi == "DingLu MA "] = "DingLu MA"
  if(!is.null(estimators_vec)){
    res_tab = filter(res_tab, Estimator %in% estimators_vec)
    res_tab$EstCombi = mgsub(res_tab$Estimator, estimators_vec, legend_levels)
    res_tab$EstCombi = factor(res_tab$EstCombi, levels = legend_levels)
  }
  
  return(res_tab)
}
##################################################################################################


plot_tables_func = function(param_n, xi, ind, AX_interactions, misspec_outcome, misspec_PS, 
    estimators_vec = c("maha_cal_rep_FALSE", "OLS_int", "maha_cal_rep_TRUE", "WLS_int", "BC_cal_rep_TRUE", "DL_MA_est"), 
    legend_levels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA"),
    extract_legend_bool=FALSE){

  all_path_lst = path_func(param_n=param_n, xi=xi, AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS)

  # ind from plot_path, the index of xi in c(0,0.05,0.1,0.2) 
  ##################################################################################################
  # GENERAL matching estimators ####
  # read data: matching on mahalanobis with PS caliper - several estimators (final_tables_general)
  
  # 3X
  table_small_pro_small_as_3 = get(load(paste0(all_path_lst$small_pro_small_as_path3, "results_table_", ind,".RData"))) 
  table_small_pro_large_as_3 = get(load(paste0(all_path_lst$small_pro_large_as_path3, "results_table_", ind,".RData"))) 
  table_large_pro_small_as_3 = get(load(paste0(all_path_lst$large_pro_small_as_path3, "results_table_", ind,".RData"))) 
  table_large_pro_large_as_3 = get(load(paste0(all_path_lst$large_pro_large_as_path3, "results_table_", ind,".RData"))) 
  
  # 5X
  table_small_pro_small_as_5 = get(load(paste0(all_path_lst$small_pro_small_as_path5, "results_table_", ind,".RData"))) 
  table_small_pro_large_as_5 = get(load(paste0(all_path_lst$small_pro_large_as_path5, "results_table_", ind,".RData"))) 
  table_large_pro_small_as_5 = get(load(paste0(all_path_lst$large_pro_small_as_path5, "results_table_", ind,".RData"))) 
  table_large_pro_large_as_5 = get(load(paste0(all_path_lst$large_pro_large_as_path5, "results_table_", ind,".RData")))
  
  # 10X
  table_small_pro_small_as_10 = get(load(paste0(all_path_lst$small_pro_small_as_path10, "results_table_", ind,".RData"))) 
  table_small_pro_large_as_10 = get(load(paste0(all_path_lst$small_pro_large_as_path10, "results_table_", ind,".RData"))) 
  table_large_pro_small_as_10 = get(load(paste0(all_path_lst$large_pro_small_as_path10, "results_table_", ind,".RData"))) 
  table_large_pro_large_as_10 = get(load(paste0(all_path_lst$large_pro_large_as_path10, "results_table_", ind,".RData")))
  ##################################################################################################
  
  ##################################################################################################
  # create a combined table for the figure ####
  
  small_pro = rbind(
    add_bias_tables(res_tab = data.frame(table_small_pro_small_as_3, pi_as=0.5, l = 1), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 3),
    add_bias_tables(res_tab = data.frame(table_small_pro_large_as_3, pi_as=0.75, l = 1), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 3),
    
    add_bias_tables(res_tab = data.frame(table_small_pro_small_as_5, pi_as=0.5, l = 2), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 5),
    add_bias_tables(res_tab = data.frame(table_small_pro_large_as_5, pi_as=0.75, l = 2), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 5),
    
    add_bias_tables(res_tab = data.frame(table_small_pro_small_as_10, pi_as=0.5, l = 3), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 10),
    add_bias_tables(res_tab = data.frame(table_small_pro_large_as_10, pi_as=0.75, l = 3), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 10)
  )
  small_pro$protected = "Low"
  
  large_pro = rbind(
    add_bias_tables(res_tab = data.frame(table_large_pro_small_as_3, pi_as=0.5, l = 1), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 3),
    add_bias_tables(res_tab = data.frame(table_large_pro_large_as_3, pi_as=0.75, l = 1), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 3),
    
    add_bias_tables(res_tab = data.frame(table_large_pro_small_as_5, pi_as=0.5, l = 2), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 5),
    add_bias_tables(res_tab = data.frame(table_large_pro_large_as_5, pi_as=0.75, l = 2), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 5),
    
    add_bias_tables(res_tab = data.frame(table_large_pro_small_as_10, pi_as=0.5, l = 3), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 10),
    add_bias_tables(res_tab = data.frame(table_large_pro_large_as_10, pi_as=0.75, l = 3), 
                    estimators_vec = estimators_vec, legend_levels=legend_levels, N = param_n, num_of_x = 10)
  )
  large_pro$protected = "High"
  
  # combine small_pro and large_pro
  small_large_pro = rbind(small_pro, large_pro)
  small_large_pro = data.frame(subset(small_large_pro, select = c("true_SACE", "N", "l", "dim_x", "pi_as", "protected", "Estimator", "EstCombi")),
                               subset(small_large_pro, select = !colnames(res_tab) %in% c("true_SACE", "N", "l", "dim_x", "pi_as", "protected", "Estimator", "EstCombi")) )
  small_large_pro = data.table(arrange(small_large_pro, protected, pi_as, l))
  
  # add label for the SACE, for each facet plot
  small_large_pro[, label := paste0(unique(true_SACE), collapse=","), by = c("protected", "pi_as")]
  temp = apply(data.frame(list.rbind(strsplit(small_large_pro$label, ","))), 2 , as.numeric) %>% round(2)
  small_large_pro$label = paste0("SACE: ", apply(temp, 1, function(x) paste(x, collapse=", ")))
  ##################################################################################################
  
  ##################################################################################################
  # plot ####
  # http://www.cookbook-r.com/Graphs/Facets_(ggplot2)/ # https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels # http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
  
  #small_large_pro$Rel_bias[1:6] = 0
  plot_general <- ggplot(small_large_pro, aes(x=l, y=Bias)) + # Bias # Rel_bias 
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
                                                              "Crude With" = "yellow1", "WLS inter" = "firebrick3", "BC With" = "palevioletred3", "DingLu MA" = "black"))  +
    facet_grid(glue('pi[pro]*" : {protected}"') ~ glue('pi[as]*" = {pi_as}"'), labeller = label_parsed) +
    theme(
      strip.text.x = element_text(size=16, face="bold"),
      strip.text.y = element_text(size=16, face="bold"),
      strip.background = element_rect(colour="black", fill="white")
    ) + 
    labs(y="Bias", x="k") + 
    # param_n=param_n, xi=xi, AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
    #ggtitle(glue('xi" : 0"')) + 
    ggtitle(paste0(bquote(xi), "=", xi, ", interactions=", 
            AX_interactions, ", mis_Y=", misspec_outcome, ", mis_PS=", misspec_PS, ", N=", param_n)) + 
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
  
  if(extract_legend_bool == TRUE){
    # Extract the legend. Returns a gtable
    lgnd_general <- get_legend(plot_general)
    # Convert to a ggplot and print
    as_ggplot(lgnd_general)
    plot_general = plot_general + theme(legend.position = 'none') 
  }
  
  print(plot_general)
  return(list(plot_general=plot_general, small_large_pro=small_large_pro))
}

##############################################################
param_n=2000
xi_values = c(0, 0.05, 0.1, 0.2)
AX_interactions = T
misspec_outcome = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)
misspec_PS = 2 # no misspec (misspec_PS=0) # Func Form misspec (misspec_PS=2)

figures_xi_values_lst = list()
results_tables_xi_values_lst = list()
for (i in 1 : length(xi_values)) {
  xi=xi_values[i]
  ind=i-1
  figures_and_results_table = plot_tables_func(param_n=param_n, xi=xi, ind=ind,
                              AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS, 
                              estimators_vec = c("maha_cal_rep_FALSE", "OLS_int", "maha_cal_rep_TRUE", "WLS_int", "BC_cal_rep_TRUE", "DL_MA_est"), 
                              legend_levels = c("Crude Wout", "OLS inter", "Crude With", "WLS inter", "BC With", "DingLu MA"),
                              extract_legend_bool=FALSE) 
  figures_xi_values_lst[[i]] = figures_and_results_table$plot_general
  results_tables_xi_values_lst[[i]] = figures_and_results_table$small_large_pro
    
}
names(figures_xi_values_lst) = names(results_tables_xi_values_lst) = xi_values
##############################################################
  
  


