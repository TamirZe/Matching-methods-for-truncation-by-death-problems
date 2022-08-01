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

######################################################################## 
# covarites descriptive ####
# cont and discrete
covarites_descriptive_table_cont_disc = function(dat, cov_descr, metric = "Full dataset", rnd = 1){
  disc_var = c(names(which(apply(subset(dat, select = cov_descr[-1]), 2, function(x) { all(x %in% 0:1)}))), "S")
  cont_var = setdiff(cov_descr[-1], disc_var)
  # cont
  m_data_cont = subset(dat, select = c("A", "EMest_p_as", "EMest_p_pro", "pi_tilde_as1", cont_var))
  trt = rbind(c( filter(m_data_cont, A==1)[, .N, by="A"]$N, 0), 
              data.frame(sapply(filter(m_data_cont, A==1), function(x) c(mean = mean(x), sd = sd(x))) %>% t))
  untrt = rbind(c( filter(m_data_cont, A==0)[, .N, by="A"]$N, 0), 
                data.frame(sapply(filter(m_data_cont, A==0), function(x) c(mean = mean(x), sd = sd(x))) %>% t))
  variables_verical_cont = data.frame(Variable = c("N", rownames(trt)[-1]),
                                      untrt = paste0(round(untrt$mean, rnd), ' (', round(untrt$sd, rnd), ')'),
                                      trt = paste0(round(trt$mean, rnd), ' (', round(trt$sd, rnd), ')'))
  variables_verical_cont$SMD = round( (trt$mean - untrt$mean) / untrt$sd, (rnd+1) )
  variables_verical_cont$SDR = round( (trt$sd / untrt$sd), (rnd+1) )
  variables_verical_cont = filter(variables_verical_cont, !Variable %in% c("A", "age2")) 
  # discrete
  m_data_disc = subset(dat, select = c("A", disc_var))
  trt = data.frame(sapply(filter(m_data_disc, A==1), function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
  untrt = data.frame(sapply(filter(m_data_disc, A==0), function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
  variables_verical_disc = data.frame(Variable = rownames(trt),
              untrt = paste0(round(untrt$count, (rnd+1)), ' (', 100 * round(untrt$prop, (rnd+1)), '%)'),
              trt = paste0(round(trt$count, (rnd+1)), ' (', 100 * round(trt$prop, (rnd+1)), '%)'))
  variables_verical_disc$SMD = round( (trt$prop - untrt$prop) / untrt$sd, (rnd+1) )
  variables_verical_disc$SDR = round( trt$sd / untrt$sd, (rnd+1) )
  #variables_verical_disc$SDR2 = round( (sqrt((trt$prop)*(1-trt$prop))) / (sqrt((untrt$prop)*(1-untrt$prop))), (rnd+1) )
  variables_verical = rbind(c("Metric", rep(metric, 3)), variables_verical_cont[1,],
              c("N_match", rep("", 3)), c("N_unq", rep("", 3)), variables_verical_cont[-1,], variables_verical_disc)
  variables_verical = rbind(c("Metric", rep(metric, (ncol(variables_verical_cont)-1))),
                            variables_verical_cont[1,],
                            c("N_match", rep("", ncol(variables_verical_cont)-1)),
                            c("N_unq", rep("", ncol(variables_verical_cont)-1)),
                            variables_verical_cont[-1,],
                            variables_verical_disc)
  
  variables_verical = filter(variables_verical, !Variable=="A")
  return(variables_verical)
}
######################################################################## 

######################################################################## 
balance_after_matching_newWF = function(m_data, match_obj, dt_match, X_sub_cols, metric, rnd=1){
  disc_var = names(which(apply(subset(m_data, select = X_sub_cols[-1]), 2, function(x) { all(x %in% 0:1)})))
  cont_var = setdiff(X_sub_cols[-1], disc_var)
  m_data_trt_cont = dt_match %>% select_if(names(.) %in% c("A", "EMest_p_as", "EMest_p_pro", "pi_tilde_as1", cont_var))
  m_data_untrt_cont = dt_match %>% select_if(names(.) %in% paste0("A0_", c("A", "EMest_p_as", "EMest_p_pro", "pi_tilde_as1", cont_var)))
  
  trt = rbind(c( m_data[A==1, .N, by="A"]$N, 0), 
              data.frame(sapply(m_data_trt_cont, function(x) c(mean = mean(x), sd = sd(x))) %>% t))
  untrt = rbind(c( m_data[A==0, .N, by="A"]$N, 0), 
                data.frame(sapply(m_data_untrt_cont, function(x) c(mean = mean(x), sd = sd(x))) %>% t))
  balance_match_cont = data.frame(Variable = c("N", rownames(trt)[-1]),
                                  untrt = paste0(round(untrt$mean, rnd), ' (', round(untrt$sd, rnd), ')'),
                                  trt = paste0(round(trt$mean, rnd), ' (', round(trt$sd, rnd), ')'))
  balance_match_cont$SMD = round( (trt$mean - untrt$mean) / untrt$sd, (rnd+1) )
  balance_match_cont$SDR = round( trt$sd / untrt$sd, (rnd+1) )
  balance_match_cont = filter(balance_match_cont, !Variable %in% c("A", "age2")) 
  balance_match_cont = rbind(balance_match_cont[1,], c("N_match", nrow(dt_match), nrow(dt_match), ""),
                             c("N_unq", length(unique(match_obj$index.control)), length(unique(match_obj$index.treated)), ""), balance_match_cont[-1,])
  # discrete
  m_data_trt_disc = subset(dt_match, select = c("A", disc_var))
  m_data_untrt_disc = subset(dt_match, select = paste0("A0_", c("A", disc_var)))
  trt = data.frame(sapply(m_data_trt_disc, function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
  untrt = data.frame(sapply(m_data_untrt_disc, function(x) c(count = sum(x), prop = mean(x), sd = sd(x))) %>% t)
  balance_match_disc = data.frame(Variable = rownames(trt),
                                  untrt = paste0(round(untrt$count, (rnd+1)), ' (', 100 * round(untrt$prop, (rnd+1)), '%)'),
                                  trt = paste0(round(trt$count, (rnd+1)), ' (', 100 * round(trt$prop, (rnd+1)), '%)'))
  balance_match_disc$SMD = round( (trt$prop - untrt$prop) / untrt$sd, (rnd+1) )
  balance_match_disc$SDR = round( trt$sd / untrt$sd, (rnd+1) )
  #balance_match_disc$SDR2 = round( (sqrt((trt$prop)*(1-trt$prop))) / (sqrt((untrt$prop)*(1-untrt$prop))), (rnd+1) )
  balance_match = rbind(c("Metric", rep(metric, (ncol(balance_match_cont)-1))),
                        balance_match_cont, balance_match_disc)
  balance_match = filter(balance_match, !Variable=="A")
  
  return(balance_match=balance_match)
}
######################################################################## 
