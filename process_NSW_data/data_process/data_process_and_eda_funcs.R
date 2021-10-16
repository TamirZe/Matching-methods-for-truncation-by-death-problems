# adjust data
######################################################################## 
adjust_data = function(data, divide_salary = 1000, log_salary = F, data_bool="DW"){
  #data$nodegree = 1 - data$nodegree
  data$emp74 = ifelse(data$re74 > 0, 1, 0)
  data$emp75 = ifelse(data$re75 > 0, 1, 0)
  data$re74 = data$re74 / divide_salary
  data$re75 = data$re75 / divide_salary
  data$id = c(1:nrow(data))
  data$intercept = 1
  colnames(data) = gsub("treat[^. ]*", "A", colnames(data))
  colnames(data)[which(colnames(data) == "re78")] = "Y"
  data$S = ifelse(data$Y == 0, 0, 1)
  data$emp_74_75 = paste0("emp(", data$emp74, ",", data$emp75 ,")")
  data = data.frame(subset(data, select = c(id, intercept)),
                    subset(data, select = -c(u74, A, S, Y, id, intercept)), 
                    subset(data, select = c(A, S, Y)))
  if(data_bool=="DW"){
    data = subset(data, select = -data_id)
  }else if(data_bool=="LL"){
    print(identical(data$u75, 1-data$emp75))
    data = subset(data, select = -u75)
    data = subset(data, select = c("id","intercept","age","education","black","hispanic",
                                   "married","nodegree", "re74","re75","emp74","emp75","emp_74_75","A","S","Y"))   
    # remove earnings and employment in 74
    #data = subset(data, select = -c(re74, emp74))
  }
  data$OBS = paste0("O(", data$A, ",", data$S, ")")
  data = data.table(data)
  return(data)
}
######################################################################## 

# covarites descriptive ####
######################################################################## 
covarites_descriptive_table = function(dat, cov_descr){ # cov_descr = variables
  N1 = nrow(filter(dat, A==1)); N0 = nrow(filter(dat, A==0))
  treatment = data.frame(sapply(subset(filter(dat, A==1), select = c("EMest_p_as", "EMest_p_pro", "e_1_as", cov_descr[-1], "S")),
                                function(x) c(mean = round(mean(x),2), sd = sd(x))) %>% t) 
  control = data.frame(sapply(subset(filter(dat, A==0), select = c("EMest_p_as", "EMest_p_pro", "e_1_as", cov_descr[-1], "S")), 
                              function(x) c(mean = round(mean(x),2), sd = sd(x))) %>% t) 
  # control$sd / sqrt(N0) # treatment$sd / sqrt(N1)
  treatment$sd = round(treatment$sd, 2); control$sd = round(control$sd, 2)
  #  vertically
  variables_table_ver = data.frame(Variable = rownames(treatment),
                                   untrt = paste0(control$mean, ' (', control$sd, ')'), 
                                   trt = paste0(treatment$mean, ' (', treatment$sd, ')')) %>% filter(Variable != "age2")
  variables_table_ver = rbind(c("Metric", rep("Full dataset", 3)), 
                              c("N", N0, N1), c("N_match", rep("", 3)), c("N_unq", rep("", 3)), variables_table_ver)
  balance_before_match_disc$SMD = round( (control$prop - treatment$prop) / control$sd, 3 )
  #print(variables_table_ver %>% xtable(caption = "Sample means and standard errors for male."), size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F)
  
  # horizontally
  variables_table_hor = subset(dat, select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", cov_descr[-1], "S"))[,
                                                                                                                lapply(.SD, mean), by="A"] %>% round(2)
  variables_table_hor = arrange(variables_table_hor, A) %>% data.frame()
  est_var_x0 = apply(subset(filter(dat, A==0), select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", cov_descr[-1], "S")), 2, sd)
  sd_smd = sqrt(est_var_x0)
  SMD = (variables_table_hor[1,-1] - variables_table_hor[2,-1]) / sd_smd[-1]
  variables_table_hor = rbind(variables_table_hor, c(rep(" ",1), unlist(mutate_if(SMD, is.numeric, round, 3))))
  variables_table_hor[3,-1] = mutate_if(SMD, is.numeric, round, 3)
  variables_table_hor = subset(variables_table_hor, select = -c(age2))
  #print(variables_table_hor %>% xtable(size="\\fontsize{12pt}{12pt}\\selectfont", include.rownames=F))
  ############################################################  
  return(list(variables_table_ver=variables_table_ver, variables_table_hor=variables_table_hor))
}

# cont and discrete
covarites_descriptive_table_cont_disc = function(dat, cov_descr, metric = "Full dataset", rnd = 1){
  disc_var = c(names(which(apply(subset(dat, select = cov_descr[-1]), 2, function(x) { all(x %in% 0:1)}))), "S")
  cont_var = setdiff(cov_descr[-1], disc_var)
  # cont
  m_data_cont = subset(dat, select = c("A", "EMest_p_as", "EMest_p_pro", "e_1_as", cont_var))
  trt = rbind(c( filter(m_data_cont, A==1)[, .N, by="A"]$N, 0), 
              data.frame(sapply(filter(m_data_cont, A==1), function(x) c(mean = mean(x), sd = sd(x))) %>% t))
  untrt = rbind(c( filter(m_data_cont, A==0)[, .N, by="A"]$N, 0), 
                data.frame(sapply(filter(m_data_cont, A==0), function(x) c(mean = mean(x), sd = sd(x))) %>% t))
  variables_verical_cont = data.frame(Variable = c("N", rownames(trt)[-1]),
                                      untrt = paste0(round(untrt$mean, rnd), ' (', round(untrt$sd, rnd), ')'),
                                      trt = paste0(round(trt$mean, rnd), ' (', round(trt$sd, rnd), ')'))
  variables_verical_cont$SMD = round( (untrt$mean - trt$mean) / untrt$sd, (rnd+1) )
  variables_verical_cont = filter(variables_verical_cont, !Variable %in% c("A", "age2")) 
  # discrete
  m_data_disc = subset(dat, select = c("A", disc_var))
  trt = data.frame(sapply(filter(m_data_disc, A==1), function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
  untrt = data.frame(sapply(filter(m_data_disc, A==0), function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
  variables_verical_disc = data.frame(Variable = rownames(trt),
                                      untrt = paste0(round(untrt$count, (rnd+1)), ' (', 100 * round(untrt$prop, (rnd+1)), '%)'),
                                      trt = paste0(round(trt$count, (rnd+1)), ' (', 100 * round(trt$prop, (rnd+1)), '%)'))
  variables_verical_disc$SMD = round( (untrt$prop - trt$prop) / untrt$sd, (rnd+1) )
  variables_verical = rbind(c("Metric", rep(metric, 3)), variables_verical_cont[1,],
                            c("N_match", rep("", 3)), c("N_unq", rep("", 3)), variables_verical_cont[-1,], variables_verical_disc)
  variables_verical = filter(variables_verical, !Variable=="A")
  return(variables_verical)
}
######################################################################## 

#DESCRIPTION OF THE DATA
######################################################################## 
descriptive = function(dd, by_obs="OBS"){
  obs_table = table(data.frame(dd)[,by_obs])
  names = names(obs_table)
  if(by_obs=="OBS"){
    obs_table = data.frame(rbind(obs_table[c("O(0,0)", "O(0,1)")], obs_table[c("O(1,0)", "O(1,1)")]))
    rownames(obs_table) = c("A=0","A=1"); colnames(obs_table) = c("S=0","S=1")
  }else{ #by_obs="em_74_75"
    obs_table = data.frame(rbind(obs_table[c("emp(0,0)", "emp(0,1)")], obs_table[c("emp(1,0)", "emp(1,1)")]))
    rownames(obs_table) = c("emp74=0","emp74=1"); colnames(obs_table) = c("emp75=0","emp75=1")
  }
  OBS_table = cbind(obs_table, Total = apply(obs_table, 1, sum))
  OBS_table = rbind(OBS_table, Total = apply(OBS_table, 2, sum))
  OBS_table_prop = round(OBS_table / sum(obs_table), 2)
  
  # SUMMARY STATISTICS
  #my.summary = function(x) list(mean = mean(x), N = length(x))
  agg_by_cell = data.table(dd[,lapply(.SD, mean), by=by_obs, .SDcols=c("re74", "re75", "Y")], N=dd[, .N, by = by_obs]$N)
  agg_by_cell = data.table(agg_by_cell, prop = round(agg_by_cell$N / sum(agg_by_cell$N), 2))
  agg_by_cell = agg_by_cell[c(3,4,2,1), ]
  agg_by_cell = mutate_if(agg_by_cell, is.numeric, round, 2)
  agg_by_A = data.table(dd[,lapply(.SD, mean), by=A, .SDcols=c("re74", "re75", "S", "Y")], N=dd[, .N, by = A]$N)
  agg_by_A = mutate_if(agg_by_A, is.numeric, round, 2)
  agg_by_S = data.table(dd[,lapply(.SD, mean), by=S, .SDcols=c("re74", "re75", "Y")], N=dd[, .N, by = S]$N)
  agg_by_S = mutate_if(agg_by_S, is.numeric, round, 2)
  return(list(OBS_table=OBS_table, OBS_table_prop=OBS_table_prop,
              agg_by_cell=agg_by_cell, agg_by_A=agg_by_A, agg_by_S=agg_by_S))
}