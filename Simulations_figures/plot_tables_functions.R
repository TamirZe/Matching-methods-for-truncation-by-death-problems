##################################################################################################
# path to results_tables ####
path_func_old = function(param_n, xi, AX_interactions, misspec_outcome, misspec_PS){
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

path_func = function(param_n, xi, mis_xi, AX_interactions, misspec_outcome, misspec_PS){
  # Cluster/GH
  #main_path = "/a/home/cc/stud_math/tamirzehavi/MatchingSACE/Simulation_studies/"
  #path_data = paste0(main_path, "Data/")
  
  # local
  main_path = "~/A matching framework for truncation by death problems/"
  path_data = paste0(main_path, "Data_cluster/")
  
  path = paste0(path_data, "Data_DGMseq/N=", param_n, "/",
                ifelse(mis_xi==0, "Data_crct_xi/", "Data_Mis_xi/"),
                ifelse(AX_interactions==T, "True_outcome_with_interactions/", "True_outcome_wout_interactions/"),
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
add_bias_tables = function(res_tab, estimators_vec=NULL, N, num_of_x){
  res_tab = data.frame(Estimator = rownames(res_tab), res_tab)
  if(!"N" %in% colnames(res_tab)) { res_tab$N = N }
  if(!"dim_x" %in% colnames(res_tab)) { res_tab$dim_x = num_of_x }
  res_tab$true_SACE = res_tab$mean[res_tab$Estimator == "true_SACE"] 
  res_tab = data.frame("true_SACE" = res_tab$true_SACE, N = res_tab$N, dim_x = res_tab$dim_x,  l = res_tab$l, pi_as = res_tab$pi_as,
                       subset(res_tab, select = !colnames(res_tab) %in% c("true_SACE", "N", "l", "dim_x", "pi_as")) )
  res_tab$Bias = res_tab$mean - as.numeric(res_tab$true_SACE)  
  res_tab$rel_bias = res_tab$Bias / as.numeric(res_tab$true_SACE) 
  res_tab$true_SACE = round(res_tab$true_SACE, 4)
  if(!is.null(estimators_vec)){
    res_tab = filter(res_tab, Estimator %in% estimators_vec)
    #res_tab$EstCombi = mgsub(res_tab$Estimator, estimators_vec, legend_levels)
    #res_tab$EstCombi = factor(res_tab$EstCombi, levels = legend_levels)
  }
  
  return(res_tab)
}
##################################################################################################

# estimators_vec: vec of estimators name
##################################################################################################
combine_small_large_pro_func = function(param_n, xi_values, mis_xi, AX_interactions, misspec_outcome, misspec_PS, estimators_vec){
  
  small_large_pro = NULL
  for(i in 1:length(xi_values)){
    xi = xi_values[i]
    ind=i-1
    all_path_lst = path_func(param_n=param_n, xi=xi, mis_xi=mis_xi, 
                             AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS)
    
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
    # l is actuallt k, but I used k for the row in mat_gamma
    small_pro = rbind(
      add_bias_tables(res_tab = data.frame(table_small_pro_small_as_3, pi_as=0.5, l = 1),  
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 3),
      add_bias_tables(res_tab = data.frame(table_small_pro_large_as_3, pi_as=0.75, l = 1), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 3),
      
      add_bias_tables(res_tab = data.frame(table_small_pro_small_as_5, pi_as=0.5, l = 2), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 5),
      add_bias_tables(res_tab = data.frame(table_small_pro_large_as_5, pi_as=0.75, l = 2), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 5),
      
      add_bias_tables(res_tab = data.frame(table_small_pro_small_as_10, pi_as=0.5, l = 3), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 10),
      add_bias_tables(res_tab = data.frame(table_small_pro_large_as_10, pi_as=0.75, l = 3), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 10)
    )
    small_pro$protected = "Low"
    
    large_pro = rbind(
      add_bias_tables(res_tab = data.frame(table_large_pro_small_as_3, pi_as=0.5, l = 1), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 3),
      add_bias_tables(res_tab = data.frame(table_large_pro_large_as_3, pi_as=0.75, l = 1), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 3),
      
      add_bias_tables(res_tab = data.frame(table_large_pro_small_as_5, pi_as=0.5, l = 2), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 5),
      add_bias_tables(res_tab = data.frame(table_large_pro_large_as_5, pi_as=0.75, l = 2), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 5),
      
      add_bias_tables(res_tab = data.frame(table_large_pro_small_as_10, pi_as=0.5, l = 3), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 10),
      add_bias_tables(res_tab = data.frame(table_large_pro_large_as_10, pi_as=0.75, l = 3), 
                      estimators_vec = estimators_vec, N = param_n, num_of_x = 10)
    )
    large_pro$protected = "High"
    
    # combine small_pro and large_pro
    small_large_pro_xi = rbind(small_pro, large_pro)
    small_large_pro_xi = data.frame(subset(small_large_pro_xi, select = c("true_SACE", "N", "l", "dim_x", "pi_as", "protected", "Estimator")), #  "EstCombi"
                                    subset(small_large_pro_xi, select = !colnames(small_large_pro_xi) %in% c("true_SACE", "N", "l", "dim_x", "pi_as", "protected", "Estimator")) ) #  "EstCombi"
    small_large_pro_xi = data.table(arrange(small_large_pro_xi, protected, pi_as, l))
    small_large_pro_xi$SE_rel_bias = (small_large_pro_xi$SE - small_large_pro_xi$sd) / small_large_pro_xi$sd
    small_large_pro_xi$SE_rel_bias[small_large_pro_xi$Estimator == "DL_MA_est"] = 0 # since SE were nor calculated for DL
    # add label for the SACE, for each facet plot, per scenario
    small_large_pro_xi[, label := paste0(unique(true_SACE), collapse=","), by = c("protected", "pi_as")]
    temp = apply(data.frame(list.rbind(strsplit(small_large_pro_xi$label, ","))), 2 , as.numeric) %>% round(2)
    small_large_pro_xi$label = paste0("SACE: ", apply(temp, 1, function(x) paste(x, collapse=", ")))
    small_large_pro = rbind(small_large_pro, data.frame(xi = xi, small_large_pro_xi))
  }
  ##################################################################################################
  return(small_large_pro)
}
##################################################################################################  

##################################################################################################
plot_tables_func_by_dimx = function(small_large_pro, param_n, mis_xi, AX_interactions, misspec_outcome, misspec_PS, 
                                    legend_levels, extract_legend_bool=FALSE){
  
  # plot ####
  # http://www.cookbook-r.com/Graphs/Facets_(ggplot2)/ # https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels # http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
  
  #small_large_pro$rel_bias[1:6] = 0
  plot_general <- ggplot(small_large_pro, aes(x=l, y=Bias)) + # Bias # rel_bias 
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
    # param_n=param_n, xi=unique(small_large_pro$xi), AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
    #ggtitle(glue('xi" : 0"')) + 
    ggtitle(paste0(bquote(xi), "=", unique(small_large_pro$xi), "(", ifelse(mis_xi, "Mis xi", "crct xi") , "), interactions=", 
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
  return(plot_general=plot_general)
}
##################################################################################################

##################################################################################################
plot_tables_func_by_xi = function(small_large_pro, param_n, mis_xi, AX_interactions, misspec_outcome, misspec_PS, 
                                  legend_levels, extract_legend_bool){
  # plot
  # http://www.cookbook-r.com/Graphs/Facets_(ggplot2)/ # https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels # http://r-statistics.co/Complete-Ggplot2-Tutorial-Part2-Customizing-Theme-With-R-Code.html
  
  #small_large_pro$rel_bias[1:6] = 0
  plot_general <- ggplot(small_large_pro, aes(x=as.factor(xi), y=Bias)) + # Bias # rel_bias 
    geom_point(aes(col = EstCombi, size = 7)) + # shape = EstCombi
    #xlim("0", ".05", "0.1", ".2") +
    labs(colour = "Estimator"
         #, size = "Emp SD"
    ) + 
    guides(colour = guide_legend(order = 1, override.aes = list(size=7))
           , size=FALSE
    ) + 
    geom_hline(yintercept = 0)
  
  plot_general = plot_general + scale_color_manual(name="Estimator", 
     labels = legend_levels, 
     values = c("Crude Wout" = "forestgreen", "OLS inter" = "dodgerblue3",
                "Crude With" = "yellow1", "WLS inter" = "firebrick3", "BC With" = "palevioletred3", "DingLu MA" = "black")) +
    facet_grid(glue('pi[pro]*" : {protected}"') ~ glue('pi[as]*" = {pi_as}"'), labeller = label_parsed) +
    theme(
      strip.text.x = element_text(size=16, face="bold"),
      strip.text.y = element_text(size=16, face="bold"),
      strip.background = element_rect(colour="black", fill="white")
    ) + 
    labs(y="Bias", x=bquote(xi)) + 
    # param_n=param_n, xi=xi, AX_interactions=AX_interactions, misspec_outcome=misspec_outcome, misspec_PS=misspec_PS
    #ggtitle(glue('xi" : 0"')) + 
    ggtitle(paste0(unique(small_large_pro$dim_x) ,"X, ", ifelse(mis_xi, "Mis xi", "Crct xi"), ", interactions=", 
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
  return(plot_general=plot_general)
}
##################################################################################################


##################################################################################################
plot_tables_func_by_dimx_paperStyle = function(small_large_pro, param_n, mis_xi, AX_interactions, misspec_outcome, misspec_PS,
                                               estimators_vec, legend_labels, colors_arg, shapes_arg, measure_to_plot="Bias"){
  
  #small_large_pro$shape = mgsub(small_large_pro$Estimator, c("OLS|WLS|BC", " inter| caliper"), c("Model-based", ""))
  
  # TODO plot # ggpubr::show_point_shapes()
  plot_general = small_large_pro %>% filter(Estimator %in% estimators_vec) %>% 
    ggplot(aes(x = l, y = Bias)) + theme_bw() +
    geom_point(alpha = 0.65, size = 5, aes(col = Estimator, shape = Estimator)) + # , shape = as.character(shape)
    xlim("3", "5", "10") +
    xlab("Number of Covariates") + 
    #labs(colour = "Estimator", shape = as.character("shape")) + 
    scale_colour_manual(name="", breaks = estimators_vec, labels = legend_labels, values = colors_arg) + 
    scale_shape_manual(name="", breaks = estimators_vec, labels = legend_labels, values = shapes_arg) +
    guides(col=guide_legend(nrow=2,byrow=TRUE)) + 
    #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
    geom_hline(yintercept = 0) + 
    #ylim(c(-0.75, 0.75)) + 
    facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {pi_as}"'), 
               labeller = label_parsed) +
    #ggtitle(paste0(measure_to_plot, ", ", 
    # bquote(xi), "=", unique(small_large_pro$xi), "(", ifelse(mis_xi, "Mis xi", "crct xi") , "), interactions=", 
    #               AX_interactions, ", mis_Y=", misspec_outcome, ", mis_PS=", misspec_PS, ", N=", param_n)) + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.575, 'cm'),
          axis.text = element_text(size  = 12),
          axis.title = element_text(size = 16),
          strip.text.x = element_text(size=16, face="bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_rect(colour="black", fill="white"))
  
  plot(plot_general)
  return(plot_general)
}
##################################################################################################

##################################################################################################
plot_tables_func_by_xi_paperStyle = function(small_large_pro, param_n, mis_xi, AX_interactions, misspec_outcome, misspec_PS, 
                                             estimators_vec, legend_labels, colors_arg, shapes_arg, measure_to_plot="Bias"){
  
  #small_large_pro$shape = mgsub(small_large_pro$Estimator, c("OLS|WLS|BC", " inter| caliper"), c("Model-based", ""))
  
  # TODO plot # ggpubr::show_point_shapes()
  plot_general = small_large_pro %>% filter(Estimator %in% estimators_vec) %>% 
    ggplot(aes(x = as.factor(xi), y = Bias)) + theme_bw() +
    geom_point(alpha = 0.65, size = 5, aes(col = Estimator, shape = Estimator)) + # , shape = as.character(shape)
    xlab(bquote(xi)) + 
    #labs(colour = "Estimator", shape = as.character("shape")) + 
    scale_colour_manual(name="", breaks = estimators_vec, labels = legend_labels, values = colors_arg) + 
    scale_shape_manual(name="", breaks = estimators_vec, labels = legend_labels, values = shapes_arg) +
    guides(col=guide_legend(nrow=2,byrow=TRUE)) + 
    #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
    geom_hline(yintercept = 0) + 
    #ylim(c(-0.75, 0.75)) + 
    facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {pi_as}"'), 
               labeller = label_parsed) +
    #ggtitle(paste0(measure_to_plot, ", ", 
    # unique(small_large_pro$dim_x) ,"X, ", ifelse(mis_xi, "Mis xi", "Crct xi"), ", interactions=", 
    #               AX_interactions, ", mis_Y=", misspec_outcome, ", mis_PS=", misspec_PS, ", N=", param_n)) + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.575, 'cm'),
          axis.text = element_text(size  = 12),
          axis.title = element_text(size = 16),
          strip.text.x = element_text(size=16, face="bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_rect(colour="black", fill="white"))
  
  plot(plot_general)
  return(plot_general)
}
##################################################################################################

##################################################################################################
plot_tables_func_by_dimx_paperStyle_flex = function(small_large_pro, param_n, mis_xi, AX_interactions, misspec_outcome, misspec_PS,
                                               estimators_vec, legend_labels, colors_arg, shapes_arg, measure_to_plot){
  
  #small_large_pro$shape = mgsub(small_large_pro$Estimator, c("OLS|WLS|BC", " inter| caliper"), c("Model-based", ""))
  
  # TODO plot # ggpubr::show_point_shapes()
  plot_general = small_large_pro %>% filter(Estimator %in% estimators_vec) %>% 
    ggplot(aes_string(x = "l", y = measure_to_plot)) + theme_bw() + # SE_rel_bias
    geom_point(alpha = 0.65, aes(col = Estimator, shape = Estimator)) + # , size = abs(SE_rel_bias)
    xlim("3", "5", "10") +
    xlab("Number of Covariates") + 
    #labs(colour = "Estimator", size = "SE_rel_bias") + 
    #labs(colour = "Estimator", shape = as.character("shape")) + 
    scale_colour_manual(name="", breaks = estimators_vec, labels = legend_labels, values = colors_arg) + 
    scale_shape_manual(name="", breaks = estimators_vec, labels = legend_labels, values = shapes_arg) +
    guides(col=guide_legend(nrow=2,byrow=TRUE)) + 
    #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = 0.95) +
    #ylim(c(-0.75, 0.75)) + 
    facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {pi_as}"'), 
               labeller = label_parsed) +
    ggtitle(paste0(measure_to_plot, ", ", 
                   bquote(xi), "=", unique(small_large_pro$xi), "(", ifelse(mis_xi, "Mis xi", "crct xi") , "), interactions=", 
                   AX_interactions, ", mis_Y=", misspec_outcome, ", mis_PS=", misspec_PS, ", N=", param_n)) + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.575, 'cm'),
          axis.text = element_text(size  = 12),
          axis.title = element_text(size = 16),
          strip.text.x = element_text(size=16, face="bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_rect(colour="black", fill="white"))
  
  plot(plot_general)
  return(plot_general)
}
##################################################################################################

##################################################################################################
plot_tables_func_by_xi_paperStyle_flex = function(small_large_pro, param_n, mis_xi, AX_interactions, misspec_outcome, misspec_PS, 
                                             estimators_vec, legend_labels, colors_arg, shapes_arg, measure_to_plot){
  
  #small_large_pro$shape = mgsub(small_large_pro$Estimator, c("OLS|WLS|BC", " inter| caliper"), c("Model-based", ""))
  
  # TODO plot # ggpubr::show_point_shapes()
  small_large_pro$xi = as.factor(small_large_pro$xi)
  plot_general = small_large_pro %>% filter(Estimator %in% estimators_vec) %>% 
    ggplot(aes_string(x = "xi", y = measure_to_plot)) + theme_bw() + # SE_rel_bias
    geom_point(alpha = 0.65, aes(col = Estimator, shape = Estimator)) + #, size = abs(SE_rel_bias)
    xlab(bquote(xi)) + 
    #labs(colour = "Estimator", shape = as.character("shape")) + 
    #labs(colour = "Estimator", size = "SE_rel_bias") + 
    scale_colour_manual(name="", breaks = estimators_vec, labels = legend_labels, values = colors_arg) + 
    scale_shape_manual(name="", breaks = estimators_vec, labels = legend_labels, values = shapes_arg) +
    guides(col=guide_legend(nrow=2,byrow=TRUE)) + 
    #guides(colour = guide_legend(order = 1, override.aes = list(size=7)), shape = guide_legend(title = "F"), size = FALSE) +
    geom_hline(yintercept = 0) + 
    geom_hline(yintercept = 0.95) + 
    #ylim(c(-0.75, 0.75)) + 
    facet_grid(glue::glue('pi[pro]*" : {protected}"') ~ glue::glue('pi[as]*" = {pi_as}"'), 
               labeller = label_parsed) +
    ggtitle(paste0(measure_to_plot, ", ", 
                   unique(small_large_pro$dim_x) ,"X, ", ifelse(mis_xi, "Mis xi", "Crct xi"), ", interactions=", 
                   AX_interactions, ", mis_Y=", misspec_outcome, ", mis_PS=", misspec_PS, ", N=", param_n)) + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.575, 'cm'),
          axis.text = element_text(size  = 12),
          axis.title = element_text(size = 16),
          strip.text.x = element_text(size=16, face="bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_rect(colour="black", fill="white"))
  
  plot(plot_general)
  return(plot_general)
}
##################################################################################################