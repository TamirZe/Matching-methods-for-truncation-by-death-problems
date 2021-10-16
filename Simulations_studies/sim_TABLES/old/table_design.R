# mat_all_estimators; WLS_NOint_mat_regression_estimators; WLS_YESint_mat_regression_estimators
# OLS_NOint_mat_regression_estimators; OLS_YESint_mat_regression_estimators

################################################################
separate_to_sets = function(d){
  lst_sets = list()
  LETTERS_set = unique(substring(rownames(d),1,1))
  for(i in 1:length(LETTERS_set)){
    lst_sets[[i]] = d[grep(LETTERS_set[i], rownames(d), ignore.case=FALSE),]
  }
  return(lst_sets)
}
################################################################

################################################################
one_set = lst_separate_sets[[1]]; HL_and_naive="WLS"; estimators_names=c("WLS without interactions","OLS with interactions")


arrange_mean_and_sd_per_set = function(one_set, HL_and_naive="HL_and_naive", estimators_names=c("Naive", "HL")){
  if(HL_and_naive=="HL_and_naive"){
    rep_F_T = as.character(unique(one_set$Replacements))
    one_set_results = data.frame(t(subset(one_set, select = -Replacements)))
    # TODO better with grep
    empirical = one_set_results[c(1,3),]
    se = one_set_results[c(2,4), ]; se = data.frame(se[,1])
    colnames(empirical) = c("Mean", "Emp_SD", "MSE"); colnames(se) = c("Est_SE")
    set_name = unique(substr(colnames(one_set_results),1,1))
    one_set_rep_F_T = data.frame(Set = set_name, Replacements = rep_F_T,
         Matching_estimator = estimators_names,
         subset(empirical, select = -MSE), se, subset(empirical, select = MSE)) 
  }else{
    if (HL_and_naive=="OLS"){
      rep_F_T = "No"
      }
    else if (HL_and_naive=="WLS"){
      rep_F_T = "Yes"
      }
    one_set_results = data.frame(t(one_set))
    one_set_results = data.frame(rbind(one_set_results[c(3,4),], one_set_results[c(1,2),]))
    empirical = one_set_results[c(1,3),]
    se = one_set_results[c(2,4), ]; se = data.frame(se[,1])
    colnames(empirical) = c("Mean", "Emp_SD", "MSE"); colnames(se) = c("Est_SE")
    set_name = unique(substr(colnames(one_set_results),1,1))
    one_set_rep_F_T = data.frame(Set = set_name, Replacements = rep_F_T,
                                 Matching_estimator = estimators_names, empirical, se)
  }
  
  return(one_set_rep_F_T) 
}
################################################################

################################################################
# arrange_all_sets_est_se_sd(rep_FALSE_or_TRUE_data=tab_est, HL_and_naive="OLS", 
#            estimators_names=c("OLS without interactions", "OLS with interactions"))
tab_est = lst_rep_F_T[[1]]

#arrange_all_sets_est_se_sd(tab_est, HL_and_naive=HL_and_naive, estimators_names=estimators_names)

arrange_all_sets_est_se_sd = function(rep_FALSE_or_TRUE_data, HL_and_naive="HL_and_naive", 
                                      estimators_names=c("Naive", "HL")){
  lst_separate_sets = separate_to_sets(rep_FALSE_or_TRUE_data)
  lst_separate_sets_in_one_row = lapply(1:length(lst_separate_sets), function(l){
    arrange_mean_and_sd_per_set(lst_separate_sets[[l]], HL_and_naive=HL_and_naive, estimators_names=estimators_names)
  })
  df_separate_sets_in_one_row = rbindlist(lst_separate_sets_in_one_row)
  return(df_separate_sets_in_one_row)
}
################################################################


################################################################
#design_tables_one_dataset(tab_est, data_set, HL_and_naive="HL_and_naive", estimators_names=c("Naive", "HL"))
design_tables_one_dataset(lst_tab_est_LS[[l]], data_set, HL_and_naive=LS_str[l], 
                          estimators_names=c(paste0(LS_str[l], " without interactions"), paste0(LS_str[l], " with interactions")))
design_tables_one_dataset(tab_est = lst_tab_est_LS[[l]], data_set, HL_and_naive=LS_str[l], 
  estimators_names=c(paste0(LS_str[l], " without interactions"), paste0(LS_str[l], " with interactions")))
design_tables_one_dataset = function(tab_est, data_set, HL_and_naive="HL_and_naive", estimators_names=c("Naive", "HL")){
  tab_est = tab_est[-grep("med", rownames(tab_est)),]
  SACE = data.frame(Set  = unique(substring(rownames(tab_est),1,1)), SACE = tab_est[grep("mean", rownames(tab_est)) , "SACE"])
  if(HL_and_naive=="HL_and_naive"){
    print(HL_and_naive)
    tab_est = subset(tab_est, select = grep(paste(c("gamma","^SACE$", data_set), collapse = "|"), colnames(tab_est)))
    param_gammas = unique(subset(tab_est, select = grep("gamma", colnames(tab_est))))
    repF = data.frame(Replacements = "No", tab_est[, grep("repF_", colnames(tab_est))])
    repT = data.frame(Replacements = "Yes", tab_est[, grep("repT_", colnames(tab_est))])
    lst_rep_F_T = list(repF, repT)
    lst_final_tables_naive_HL = list()
    for (i in 1:length(lst_rep_F_T)) {
      lst_final_tables_naive_HL[[i]] = arrange_all_sets_est_se_sd(lst_rep_F_T[[i]],HL_and_naive=HL_and_naive,estimators_names=estimators_names)
    }
    df_final_tables_one_grp_of_estimators = rbindlist(lst_final_tables_naive_HL)
  }else{
    print(HL_and_naive)
    tab_est = subset(tab_est, select = grep(paste(c("^SACE$", data_set), collapse = "|"), colnames(tab_est)))
    tab_est = subset(tab_est, select = -SACE)
    param_gammas="LS"
    df_final_tables_one_grp_of_estimators = 
      arrange_all_sets_est_se_sd(tab_est, HL_and_naive=HL_and_naive, estimators_names=estimators_names)
  }
  df_final_tables_one_grp_of_estimators = arrange(merge(SACE, df_final_tables_one_grp_of_estimators, by="Set"), 
                                                  Set, Replacements)
  return(list(df_final_tables_one_grp_of_estimators = df_final_tables_one_grp_of_estimators
              ,param_gammas = param_gammas
              ))
}
################################################################


################################################################
arrange_LS_before_design_tables_one_dataset(LS_YESint=WLS_YESint_mat_regression_estimators,
                                            LS_NOint=WLS_NOint_mat_regression_estimators)
arrange_LS_before_design_tables_one_dataset = function(LS_YESint, LS_NOint){
  LS_YESint = subset(LS_YESint, select = -grep("SACE_conditional", colnames(LS_YESint))) 
  LS_NOint = subset(LS_NOint, select = -grep("SACE_conditional", colnames(LS_NOint))) 
  colnames(LS_YESint)[-1] = paste0("intYes_", colnames(LS_YESint)[-1]) 
  colnames(LS_NOint)[-1] = paste0("intNo_", colnames(LS_NOint)[-1])
  rows = rownames(LS_NOint)
  #LS_YESint$rows = rownames(LS_YESint); LS_NOint$rows = rownames(LS_NOint)
  tab_est = merge(LS_YESint, LS_NOint, by = 0) %>% arrange(rows)
  
  rownames(tab_est) = tab_est$Row.names 
  tab_est = subset(tab_est, select = -Row.names)
  if(identical(tab_est$SACE.x, tab_est$SACE.y)){
    tab_est = subset(tab_est, select = -grep("SACE.y", colnames(tab_est)))
    colnames(tab_est)[grep("SACE.x", colnames(tab_est))] = "SACE"
  }
  return(tab_est)
}
################################################################

################################################################
# design_tables_one_dataset_ALL_estimators(data_set = data_set_names_vec[1], mat_all_estimators,
#                                          OLS_YESint_mat_regression_estimators, OLS_NOint_mat_regression_estimators, 
#                                          WLS_YESint_mat_regression_estimators, WLS_NOint_mat_regression_estimators)

design_tables_one_dataset_ALL_estimators = function(data_set, mat_all_estimators,
  OLS_YESint_mat_regression_estimators, OLS_NOint_mat_regression_estimators, WLS_YESint_mat_regression_estimators, WLS_NOint_mat_regression_estimators){
  
  # mat_all_estimators: HL_and_naive
  tab_est = mat_all_estimators
  final_tables_naive_HL_ONE_dataset = 
    design_tables_one_dataset(tab_est, data_set, HL_and_naive="HL_and_naive", estimators_names=c("Crude", "HL"))
  df_final_tables_naive_HL_ONE_dataset = final_tables_naive_HL_ONE_dataset$df_final_tables_one_grp_of_estimators
  gamma = final_tables_naive_HL_ONE_dataset$param_gammas
  
  # LS
  #  OLS and WLS
  tab_est_OLS = arrange_LS_before_design_tables_one_dataset(OLS_YESint_mat_regression_estimators, OLS_NOint_mat_regression_estimators)
  tab_est_WLS = arrange_LS_before_design_tables_one_dataset(subset(WLS_YESint_mat_regression_estimators, select = -grep("naive", colnames(WLS_YESint_mat_regression_estimators))), 
                                                            subset(WLS_NOint_mat_regression_estimators, select = -grep("naive", colnames(WLS_NOint_mat_regression_estimators))))
  lst_tab_est_LS = list(tab_est_OLS, tab_est_WLS)
  LS_str = c("OLS", "WLS")
  lst_final_tables_OLS_WLS_dataset = lapply(1:length(lst_tab_est_LS), function(l){
    design_tables_one_dataset(lst_tab_est_LS[[l]], data_set, HL_and_naive=LS_str[l], 
      estimators_names=c(paste0(LS_str[l], " without interactions"), 
                         paste0(LS_str[l], " with interactions")))
    })
  df_final_tables_OLS_WLS_dataset = rbindlist(lapply(lst_final_tables_OLS_WLS_dataset, "[[", "df_final_tables_one_grp_of_estimators")) %>%
                arrange(Set, Replacements)
  
  # merge naive, OLS, WLS together
  df_final_tables_dataset = rbind(df_final_tables_naive_HL_ONE_dataset, df_final_tables_OLS_WLS_dataset) %>%
              arrange(Set, Replacements)
  return(list(df_final_tables_dataset = df_final_tables_dataset, gamma = gamma))
}
################################################################


# TODO @@@ run on all data_set
#paste(c("gamma","^SACE$", data_set), collapse = "|")
#data_set_names_vec = c("all_", "wout_O_0_0_", "S1_")[-2]
#data_set_names_vec = c("all", "wout_O_0_0", "S1")[-2]
data_set_names_vec = c("_all", "_wout_O_0_0", "_S1")[-2]

#data_set = data_set_names_vec[1] 
# a = design_tables_one_dataset_ALL_estimators(data_set, mat_all_estimators,
#     OLS_YESint_mat_regression_estimators, OLS_NOint_mat_regression_estimators,
#     WLS_YESint_mat_regression_estimators, WLS_NOint_mat_regression_estimators)

lst_final_tables_ALL_est_ALL_dataset_with_gamma = lapply(1:length(data_set_names_vec), function(l){
  design_tables_one_dataset_ALL_estimators(data_set_names_vec[l], mat_all_estimators,
      OLS_YESint_mat_regression_estimators, OLS_NOint_mat_regression_estimators, 
      WLS_YESint_mat_regression_estimators, WLS_NOint_mat_regression_estimators)
})
data_set_names_vec
lst_final_tables_ALL_est_ALL_dataset = lapply(lst_final_tables_ALL_est_ALL_dataset_with_gamma, "[[", "df_final_tables_dataset")
lst_final_tables_gammas = lapply(lst_final_tables_ALL_est_ALL_dataset_with_gamma, "[[", "gamma")

# TODO estimated MSE
#df = lst_final_tables_ALL_est_ALL_dataset[[1]]
calculate_MSE_per_row = function(df){
  df$my_MSE = (df$SACE - df$Mean)^2 + (df$Emp_SD)^2
  return(df)
}
for (i in 1:length(lst_final_tables_ALL_est_ALL_dataset)){
  lst_final_tables_ALL_est_ALL_dataset[[i]] = calculate_MSE_per_row(lst_final_tables_ALL_est_ALL_dataset[[i]])
  #gsub(lst_final_tables_ALL_est_ALL_dataset[[i]]$Matching_estimator, "Crude", "Naive")
}
names(lst_final_tables_ALL_est_ALL_dataset) = data_set_names_vec

# TODO duplicate cols for 2 param_n
  

#############################################################################
# COVERAGE RATE after calculating CI in sim1
sets_of_param = LETTERS[1:length(list_all_WLS_NOint_regression_estimators)] # nrow(mat_gamma)
list_all_CI_temp = list_all_CI # list_all_CI_temp from main
for (i in 1:length(list_all_CI_temp)) {
  rownames(list_all_CI_temp[[i]]) = paste0(sets_of_param[i], "_", rownames(list_all_CI_temp[[i]]))
}
#mat_all_CI_all_param = list.rbind(list_all_CI_temp)


calculate_coverage = function(mat_all_CI_one_param){
  rows = rownames(mat_all_CI_one_param)
  # coverage
  x = mat_all_CI_one_param$repF_MATCH_all # mat_all_CI_all_param$SACE
  bool_CI_contain_SACE = function(x){
    lower_bound =  as.numeric(sub(",.*", "", x)) #as.numeric(sapply(strsplit(as.character(x), ","), "[", 1)) 
    upper_bound = as.numeric(sub(".*,", "", x))
    bool = mat_all_CI_one_param$SACE >= lower_bound & mat_all_CI_one_param$SACE <= upper_bound 
    bool = ifelse(bool==T,1,0)
    return(bool)
  }
  
  bool_mat = data.frame(apply(mat_all_CI_one_param[,-c(1:2)], 2, bool_CI_contain_SACE))
  bool_mat = data.frame(rbind(bool_mat, apply(bool_mat, 2, mean)))
  rownames(bool_mat) = c(rows, paste0(unique(substr(rows,1,2)), "Coverage"))
  return(bool_mat)
}


list_all_CI_temp_with_coverage = lapply(1:length(list_all_CI_temp), function(l){
  calculate_coverage(list_all_CI_temp[[l]])
})

coverage_mat = NULL
for (i in 1:length(list_all_CI_temp_with_coverage)) {
  coverage_temp = list_all_CI_temp_with_coverage[[i]]
  coverage_temp = coverage_temp[grep("_Coverage", rownames(coverage_temp)),]
  rownames(coverage_temp) = substr(rownames(coverage_temp),1,1)
  coverage_mat = rbind(coverage_mat, coverage_temp)
}
#coverage_mat = data.frame(t(data.frame(Set = rownames(coverage_mat), coverage_mat)))
coverage_mat = data.frame(t(data.frame(coverage_mat)))
coverage_mat$Matching_estimator = rownames(coverage_mat)

# TODO separate to repF and OLS & repT and WLS and "merge" with lst_final_tables_ALL_est_ALL_dataset
#data_set_names_vec = c("_all", "_wout_O_0_0", "_S1")[-2]
#data_set_names_vec = c("all", "wout_O_0_0", "S1")[-2]

coverage_mat$Replacements = ""
coverage_mat$Replacements[grep("SACE|repF|OLS", rownames(coverage_mat))] = "No"
coverage_mat$Replacements[grep("SACE|repT|WLS", rownames(coverage_mat))] = "Yes"
coverage_mat_melt = melt(coverage_mat, measure.vars = c("A","B"), variable.name = "Set",id=c("Matching_estimator","Replacements"))
colnames(coverage_mat_melt)[which(colnames(coverage_mat_melt)=="value")] = "Coverage"

list_coverage_dataset = list()
for(i in 1:length(data_set_names_vec)){
  list_coverage_dataset[[data_set_names_vec[i]]] = 
    coverage_mat_melt[grep(paste0("naive_|Set|",data_set_names_vec[i]), coverage_mat_melt$Matching_estimator),]
}

#df_coverage = list_coverage_dataset[[1]]
change_est_names = function(df_coverage){
  df_coverage$newname = df_coverage$Matching_estimator
  df_coverage$newname[grep("HL",df_coverage$newname)] ="HL"
  df_coverage$newname[grep("rep",df_coverage$newname)] = "Crude"
  df_coverage$newname[grep("WLS_NOint",df_coverage$newname)] ="WLS without interactions"
  df_coverage$newname[grep("WLS_YESint",df_coverage$newname)] = "WLS with interactions"
  df_coverage$newname[grep("OLS_NOint",df_coverage$newname)] ="OLS without interactions"
  df_coverage$newname[grep("OLS_YESint",df_coverage$newname)] ="OLS with interactions"
  df_coverage$newname[grep("^naive_without_matching$",df_coverage$newname)] = "Naive without matching"
  df_coverage$newname[grep("survivors_naive_without_matching",df_coverage$newname)] = "Survivors naive without matching"
  df_coverage = data.frame(Matching_estimator = df_coverage$newname, 
                           subset(df_coverage, select = -c(Matching_estimator, newname)))
  return(df_coverage)
}

list_coverage_dataset = lapply(1:length(list_coverage_dataset), function(l){
  change_est_names(list_coverage_dataset[[l]])
})
names(list_coverage_dataset) = data_set_names_vec

naives_before_matching_coverage = list()
for(i in 1:length(lst_final_tables_ALL_est_ALL_dataset)){
  #d = data.frame(lst_final_tables_ALL_est_ALL_dataset[[i]])
  #d$Matching_estimator = gsub(d["Matching_estimator"], "Naive", "Crude")
  temp = merge(lst_final_tables_ALL_est_ALL_dataset[[i]], list_coverage_dataset[[i]], 
        by = c("Set", "Replacements", "Matching_estimator")
        #, all.y = TRUE
        )
  lst_final_tables_ALL_est_ALL_dataset[[data_set_names_vec[i]]] = temp
  
  naives_before_matching_coverage[[data_set_names_vec[i]]] = filter(list_coverage_dataset[[i]],
                      Matching_estimator %in% c("Naive without matching", "Survivors naive without matching") )
    
}
naives_before_matching_coverage = naives_before_matching_coverage[[1]]

#lst_final_tables_ALL_est_ALL_dataset
#save(lst_final_tables_ALL_est_ALL_dataset, file = "lst_final_tables_ALL_est_ALL_dataset.RData")
#############################################################################


#############################################################################
# TODO COVERAGE RATE separate lists after main script
# lst_of_all_est_lst = list(list_all_mat_SACE_estimators, 
# list_all_WLS_NOint_regression_estimators, list_all_WLS_YESint_regression_estimators, 
# list_all_OLS_NOint_regression_estimators, list_all_OLS_YESint_regression_estimators)
# 
# 
# list_all_estimators = lst_of_all_est_lst_new_names[[1]]
# change_rownames_to_LETTERS_by_param_set_all_rows = function(list_all_estimators, name_sets_of_param=nrow(mat_gamma)){
#   for(i in 1:length(sets_of_param)){
#     rownames(list_all_estimators[[i]]) = paste0(name_sets_of_param[i], "_", rownames(list_all_estimators[[i]]))
#   }
#   return(list_all_estimators)
# }
# 
# sets_of_param = LETTERS[1:length(list_all_WLS_NOint_regression_estimators)] # nrow(mat_gamma)
# lst_of_all_est_lst_new_names = lst_of_all_est_lst # list_all_CI_temp lst_of_all_est_lst
# for(i in 1:length(lst_of_all_est_lst_new_names)){
#   lst_of_all_est_lst_new_names[[i]] =
#     change_rownames_to_LETTERS_by_param_set_all_rows(lst_of_all_est_lst_new_names[[i]], sets_of_param)
# }
#############################################################################

# TODO add DING and naives
naive_and_ding_est = subset(mat_all_estimators, select = grep("SACE|naive|DING", colnames(mat_all_estimators)))
naive_and_ding_est = subset(naive_and_ding_est, select = -grep("conditional|less0", colnames(naive_and_ding_est)))
naive_and_ding_est = naive_and_ding_est[-grep("med", rownames(naive_and_ding_est)),]
naive_and_ding_est_separate_sets = separate_to_sets(naive_and_ding_est)
naive_and_ding_est_separate_sets_adj = list()
for(i in 1:length(naive_and_ding_est_separate_sets)){
  temp = (round(t(naive_and_ding_est_separate_sets[[i]]),4)) %>% data.frame()
  Set_temp = unique(substring(colnames(temp),1,1))
  colnames(temp) = substring(colnames(temp),3)
  SACE_temp = filter(lst_final_tables_ALL_est_ALL_dataset[[i]], Set==Set_temp)$SACE %>% unique
  temp = data.frame(Set = Set_temp,  Replacements = "", Matching_estimator = rownames(temp), SACE = SACE_temp, temp)
  colnames(temp) = mgsub::mgsub(colnames(temp), c("mean", "sd"), c("Mean", "Emp_SD"))
  #cols = setdiff(colnames(lst_final_tables_ALL_est_ALL_dataset[[1]]), colnames(temp))
  
  # adding estimated se of naive estimators
  ind_rows_se = which(temp$Matching_estimator %in% c("most_naive_est_se", "sur_naive_est_se")==T)
  temp$Est_SE = NA; temp$Est_SE[ind_rows_se - 1] = temp$Mean[ind_rows_se]
  temp = temp[-ind_rows_se,]
  naive_and_ding_est_separate_sets_adj[[i]] = temp
}
naive_and_ding_est_all_sets = list.rbind(naive_and_ding_est_separate_sets_adj)
naive_and_ding_est_all_sets = calculate_MSE_per_row(naive_and_ding_est_all_sets)
naive_and_ding_est_all_sets$Matching_estimator = mgsub::mgsub(as.character(naive_and_ding_est_all_sets$Matching_estimator), 
  c("most_naive_est", "sur_naive_est"), c("Naive without matching", "Survivors naive without matching"))
naive_and_ding_est_all_sets = 
  merge(naive_and_ding_est_all_sets, naives_before_matching_coverage, by = c("Set", "Matching_estimator", "Replacements"), all.x = TRUE)

final_tables = list()
for(i in 1:length(data_set_names_vec)){
  temp = data.frame(bind_rows(lst_final_tables_ALL_est_ALL_dataset[[i]], naive_and_ding_est_all_sets)) %>% arrange(Set)
  temp = filter(temp, Matching_estimator!="SACE")
  temp = temp %>% mutate_if(is.numeric, round, digits=4)
  temp$Matching_estimator = mgsub::mgsub(temp$Matching_estimator, c("DING_est", "DING_model_assisted_est_ps", "most_naive_est", "sur_naive_est"),
                c("DING & Lu", "DING & Lu model assisted", "Naive without matching", "Survivors naive without matching"))
  final_tables[[data_set_names_vec[i]]] = temp
}


save(final_tables, file = "final_tables.RData")
print(final_tables[[2]]  %>% xtable(digits=c(3)
  ,caption = "Matching on mahalanobis distance with caliper on PS, using Abadie Imbens SE estimator from Match function from Matching packgae. It is important to mention that usually AI SE estimator is biased in propensity score macthing (try 2016 paper?), and unbiased when matching on the mahalanobis"),
  size="\\fontsize{10pt}{10pt}\\selectfont")




###############################################
# COVERAGE SIMULATION
rv_binom = c()
rv_binom_lower_0.92 = c()
for(i in 1:10000){
  x = rbinom(1, 200, 0.95)
  rv_binom[i] = x
  # mean of bern 200 Bernoullis
  rv_binom_lower_0.92[i] = ifelse((x / 200) <= 0.92, 1, 0)
}
hist(rv_binom)
# binom is sum of Bernoullis
length(rv_binom[rv_binom <= 0.92*200]) / length(rv_binom)
mean(rv_binom_lower_0.92)

x2 = rbinom(10000, 200, 0.95)
length(x2[x2 <= 0.92*200]) / length(x2)

theoretical_p_lower_than_92 = pbinom(0.92*200, 200, 0.95, lower.tail = TRUE, log.p = FALSE)
pbinom(qbinom(0.5, 200, 0.95, lower.tail = TRUE, log.p = FALSE), 200, 0.95, lower.tail = TRUE, log.p = FALSE)
###############################################
