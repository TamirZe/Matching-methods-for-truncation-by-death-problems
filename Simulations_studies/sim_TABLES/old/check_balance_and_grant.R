means_by_subset_sum = mat_all_means_by_subset[grep("S1|mean_as", 
                       rownames(mat_all_means_by_subset), ignore.case = F) , ] %>% round(4)
ctr_as_matched_after = sum(mat_all_repeated_as_and_pro["s1repT_S1",-c(29,30)] * c(c(1:14), c(1:14)))
S1_repF_T_mat_all_matched_units = data.frame(t(mat_all_matched_units[grep("F_S1|T_S1", rownames(mat_all_matched_units)),])) %>% round(3)
ctr_as_matched_after - S1_repF_T_mat_all_matched_units["ctr_asAfS1","s1repT_S1"]


#TODO path to final tables,  
# 3X
############### no misspec ############### 
# model with interaction # model wout interaction

small_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/no misspec/small pro/3X/"
large_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/no misspec/large pro/3X/"
# 5X
small_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/no misspec/small pro/5X/"
large_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/no misspec/large pro/5X/"
# 10X
small_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/no misspec/small pro/5X/"
large_pro_path10 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/no misspec/large pro/5X/"

#TODO Func Form
# 3X
small_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/Func Form/mis2 ff small pro -3 3/3X/"
large_pro_path3 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/Func Form/mis2 ff large pro -3 3/3X/"
# 5X
small_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/Func Form/mis2 ff small pro -3 3/5X/"
large_pro_path5 = "C:/Users/tamir/Documents/R projects/AA Thesis Causal inference/Simulations principal score estimation/R data/ding EM code/model wout interaction/Func Form/mis2 ff large pro -3 3/5X/"


#TODO read data: final_tables_general
# sum_mat_name = "mat_all_matched_units"  # mat_all_means_by_subset # deparse(substitute(sum_mat))
read_get_sum_mat = function(sum_mat_name){
  sum_mat_small_pro3 = get(load(paste0(small_pro_path3, sum_mat_name, ".RData")))
  sum_mat_large_pro3 = get(load(paste0(large_pro_path3, sum_mat_name, ".RData")))
  sum_mat_small_pro5 = get(load(paste0(small_pro_path5, sum_mat_name, ".RData"))) 
  sum_mat_large_pro5 = get(load(paste0(large_pro_path5, sum_mat_name, ".RData"))) 
  sum_mat_small_pro10 = get(load(paste0(small_pro_path10, sum_mat_name, ".RData"))) 
  sum_mat_large_pro10 = get(load(paste0(large_pro_path10, sum_mat_name, ".RData"))) 
  return(list(sum_mat_small_pro3=sum_mat_small_pro3, sum_mat_large_pro3=sum_mat_large_pro3,
              sum_mat_small_pro5=sum_mat_small_pro5, sum_mat_large_pro5=sum_mat_large_pro5,
              sum_mat_small_pro10=sum_mat_small_pro10, sum_mat_large_pro10=sum_mat_large_pro10))
}

#TODO arrange sum mat
# S1_repF_T_mat_all_matched_units
sum_mat_name = "mat_all_matched_units" 
S1_repF_T_mat_all_matched_units_small_pro3 = 
  data.frame(t(sum_mat_small_pro3[grep("F_S1|T_S1", rownames(sum_mat_small_pro3)),]))
S1_repF_T_mat_all_matched_units_large_pro3 = 
  data.frame(t(sum_mat_large_pro3[grep("F_S1|T_S1", rownames(sum_mat_large_pro3)),]))

# mat_all_means_by_subset
sum_mat_name = "mat_all_means_by_subset" 

arrange_means_by_subset_sum = function(sum_mat){
  means_by_subset = sum_mat[grep("S1|mean_as", 
     rownames(sum_mat), ignore.case = F) , ] %>% round(3)
  means_by_subset = means_by_subset[
    - intersect( grep("reF|approx_mean_x", rownames(means_by_subset), ignore.case = F),
      grep("match|approx_mean_x", rownames(means_by_subset)) ), ]
  return(means_by_subset)
}

means_by_subset_small_pro3 = arrange_means_by_subset_sum(sum_mat_small_pro3)
means_by_subset_large_pro3 = arrange_means_by_subset_sum(sum_mat_large_pro3)
means_by_subset_small_pro5 = arrange_means_by_subset_sum(sum_mat_small_pro5)
means_by_subset_large_pro5 = arrange_means_by_subset_sum(sum_mat_large_pro5)
 
 
#TODO print
# all_matched_units
print(S1_repF_T_mat_all_matched_units_small_pro3 %>% xtable(digits=c(2),
    caption = "Small pro 3X"), size="\\fontsize{11pt}{11pt}\\selectfont"
    #, include.rownames=F
    )
print(S1_repF_T_mat_all_matched_units_large_pro3 %>% xtable(digits=c(2),
    caption = "large pro 3X"), size="\\fontsize{11pt}{11pt}\\selectfont"
    #, include.rownames=F
    )

# means_by_subset
print(means_by_subset_small_pro3 %>% xtable(digits=c(3),
    caption = "Small pro 3X"), size="\\fontsize{11pt}{11pt}\\selectfont")
print(means_by_subset_large_pro3 %>% xtable(digits=c(3),
    caption = "large pro 3X"), size="\\fontsize{11pt}{11pt}\\selectfont")
print(means_by_subset_small_pro5 %>% xtable(digits=c(3),
    caption = "Small pro 5X"), size="\\fontsize{11pt}{11pt}\\selectfont")
print(means_by_subset_large_pro5 %>% xtable(digits=c(3),
    caption = "large pro 5X"), size="\\fontsize{11pt}{11pt}\\selectfont")




########################################################################################
# TODO GRANT
tables_small_pro3 = get(load(paste0(small_pro_path3, "final_tables_general.RData"))) 
tables_large_pro3 = get(load(paste0(large_pro_path3, "final_tables_general.RData"))) 
tables_small_pro5 = get(load(paste0(small_pro_path5, "final_tables_general.RData"))) 
tables_large_pro5 = get(load(paste0(large_pro_path5, "final_tables_general.RData")))

#tab = tables_small_pro5$`_S1`; tab$k = 5
arrange_tab_short = function(tab, 
   estimators_vec=c("Crude Yes", "OLS inter", "WLS inter", "Naive survivors", "DingLu MA"),
   N_obs, num_of_x, pi_pro_vec){
    tmp = tab
    tmp = tmp[!substring(tmp$`Set & parameter`,1,1)=="A" ,]
    tmp[,1] = mgsub(tmp[,1], c("B", "C"), c("A", "B"))
    tmp = func_add_AND_remove_COLS(tmp, N = N_obs, num_of_x = num_of_x)
    colnames(tmp)[grep("parameter", colnames(tmp))] = "Set_and_parameter"
    tmp = data.frame("Set_and_parameter" = tmp[,1], N = tmp$N, dim_x = tmp$dim_x, 
                     subset(tmp, select = !colnames(tmp) %in% c("Set_and_parameter", "N", "dim_x")) )
    # tmp[,1] IS tmp$`Set & parameter`
    tmp$Bias = tmp$Mean - as.numeric(substr(tmp[,1],3,10))  
    tmp = data.frame(Scenario = substr(tmp[,1],1,1), tmp)
    tmp$EstCombi = paste0(tmp$Estimator, " ", tmp$Replacements)
    tmp$EstCombi[tmp$EstCombi == "DingLu MA "] = "DingLu MA"
    tmp$EstCombi[tmp$EstCombi == "Naive survivors "] = "Naive survivors"
    tmp$pro[tmp$Scenario=="A"] = paste0("pi pro = ", pi_pro_vec[1])
    tmp$pro[tmp$Scenario=="B"] = paste0("pi pro = ", pi_pro_vec[2])
    if(!is.null(estimators_vec)){
      tmp$EstCombi = mgsub(tmp$EstCombi, c("OLS No", "OLS inter No", "WLS Yes", "WLS inter Yes"),
                           c("OLS", "OLS inter", "WLS", "WLS inter"))
      tmp = filter(tmp, EstCombi %in% estimators_vec)
    }
    tmp = data.frame(subset(tmp, select = c(Set_and_parameter, N, k, Scenario, pro, 
                                            EstCombi, Bias, Emp.SD)))
    tmp[,1] = paste0(substr(tmp[,1],1,1), " ", round(as.numeric(substr(tmp[,1],3,10)), 2) )
    colnames(tmp) = mgsub::mgsub(colnames(tmp), c("Set_and_parameter", "EstCombi", "Emp.SD"),
                 c("Set & parameter", "Estimator", "Emp SD"))
    tmp$Estimator = factor(tmp$Estimator, levels = estimators_vec)
    tmp$Bias = round(tmp$Bias, 3); tmp$`Emp SD` = round(tmp$`Emp SD`, 3)
    #tmp = tmp %>% mutate_if(is.numeric, round, digits=2)
  return(tmp)
}

estimators_vec=c("Crude No", "Crude Yes", "OLS", "WLS", "Naive survivors", "DingLu MA")

tables_small_pro3$`_S1`$k = 3; tables_small_pro5$`_S1`$k = 5
tables_large_pro3$`_S1`$k = 3; tables_large_pro5$`_S1`$k = 5

small_pro = rbind( arrange_tab_short(tables_small_pro3$`_S1`, estimators_vec=estimators_vec,
                      N_obs = param_n, num_of_x = 3, pi_pro_vec = c(0.1, 0.1)),
                    arrange_tab_short(tables_small_pro5$`_S1`, estimators_vec=estimators_vec,
                    N_obs = param_n, num_of_x = 5, pi_pro_vec = c(0.1, 0.1)) )

large_pro = rbind( arrange_tab_short(tables_large_pro3$`_S1`, estimators_vec=estimators_vec, 
                                      N_obs = param_n, num_of_x = 3, pi_pro_vec = c(0.35, 0.14)),
                    arrange_tab_short(tables_large_pro5$`_S1`, estimators_vec=estimators_vec, 
                                     N_obs = param_n, num_of_x = 5, pi_pro_vec = c(0.35, 0.14)) )

small_large_pro = rbind(small_pro, large_pro)
small_large_pro$Scenario = mgsub(small_large_pro$Scenario, c("A", "B"), c("pi as = 0.5", "pi as = 0.75"))
colnames(small_large_pro)[which(colnames(small_large_pro)=="Scenario")] = "as"
small_large_pro = arrange(small_large_pro, k)
small_large_pro$`Set & parameter` =
  paste0(rep(LETTERS[1 : (nrow(small_large_pro) / length(estimators_vec))], each =length(estimators_vec)), 
  " ",  substring(small_large_pro$`Set & parameter`, 3,10))
table_without_mis = small_large_pro
table_mis_funcform = small_large_pro


save(table_without_mis, file = "table_without_mis.RData")
save(table_mis_funcform, file = "table_mis_funcform.RData")

print(table_without_mis %>% xtable(digits=c(3), caption = "True model without intercations.
        Results when the PS mdoel is correctly specified,
        N = 2000, K= 3,5"), size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)
print(table_mis_funcform %>% xtable(digits=c(3), caption = "True model without intercations.
        Results when the PS mdoel is misspecified, 
        N = 2000, K= 3,5"), size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)
########################################################################################

