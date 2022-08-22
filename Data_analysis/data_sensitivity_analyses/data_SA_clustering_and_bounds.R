library(factoextra); library(NbClust); library(pracma)

clustering_func = function(dataset, arm, variables_clstr, num_clstrs=NULL){
  data_arm = dataset %>% filter(A == arm) %>% subset(select = variables_clstr)
  data_arm = scale(data_arm)
  
  library(factoextra)
  # Elbow method
  fviz_nbclust(data_arm, kmeans, method = "wss") +
    geom_vline(xintercept = c(4,7), linetype = 2)+
    labs(subtitle = "Elbow method")
  # Silhouette method
  fviz_nbclust(data_arm, kmeans, method = "silhouette")+
    labs(subtitle = "Silhouette method")
  # Gap statistic
  # nboot = 50 to keep the function speedy. 
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  set.seed(123)
  fviz_nbclust(data_arm, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
    labs(subtitle = "Gap statistic method")
  
  if(is.null(num_clstrs)){
    opt_num_clstrs <- NbClust(data_arm, distance = "euclidean", min.nc = 2,
                  max.nc = 10, method = "kmeans")
    opt_num_clstrs_tab = fviz_nbclust(opt_num_clstrs)
    num_clstrs = 
      as.numeric(as.character(opt_num_clstrs_tab$data[which.max(opt_num_clstrs_tab$data$freq), "Number_clusters"]))
  }
  
  k_means_obj = kmeans(x = data_arm, centers = num_clstrs, iter.max = 20, nstart = 25)
  clusters = k_means_obj$cluster
  return(clusters)
}

calculate_upper_alpha_per_clstr = function(data_per_clstr, Y_sort, as_amount){
  #round(as_amount)
  #second_g_amount = length(Y_sort) - as_amount
  
  # true bounds, with the fraction located in the appropriate location
  up_bound_second_g = mean(Y_sort[ceil(as_amount):length(Y_sort)])
  low_bound_as = mean(Y_sort[1:floor(as_amount)])
  
  # the least restrictive - less wide bounds
  #up_bound_second_g = mean(Y_sort[floor(as_amount):length(Y_sort)]) # floor(as_amount)
  #low_bound_as = mean(Y_sort[1:floor(as_amount)-1])
  
  # round
  #up_bound_second_g = mean(Y_sort[round(as_amount):length(Y_sort)])
  #low_bound_as = mean(Y_sort[1:(round(as_amount) - 1)])
  
  possible_alpha_up = up_bound_second_g / low_bound_as
  return(possible_alpha_up = possible_alpha_up)
}

calculate_lower_alpha_per_clstr = function(data_per_clstr, Y_sort, as_amount){
  
  # true bounds, with the fraction located in the appropriate location
  low_bound_second_g = mean(Y_sort[1:(floor(length(Y_sort) - as_amount))])
  up_bound_as = mean(Y_sort[(ceil(length(Y_sort) - as_amount)):length(Y_sort)])
  
  # the least restrictive - less wide bounds
  #low_bound_second_g = mean(Y_sort[1:(ceil(length(Y_sort) - as_amount))])
  #up_bound_as = mean(Y_sort[(ceil(length(Y_sort) - as_amount) + 1):length(Y_sort)])
  
  # round
  #low_bound_second_g = mean(Y_sort[1:(length(Y_sort) - round(as_amount))])
  #up_bound_as = mean(Y_sort[(length(Y_sort) - round(as_amount) + 1):length(Y_sort)])
  
  possible_alpha_low = low_bound_second_g / up_bound_as
  return(possible_alpha_low = possible_alpha_low)
}



dataset = data_with_PS[data_with_PS$S == 1,]
arm = 0
xi = 0
variables_clstr = c("age", "education", "re75", "black", "hispanic", "married", "nodegree", "emp75")
model_per_cluster = TRUE

dataset_arm = dataset %>% filter(A==arm)
Y_model = lm(as.formula(paste0("Y ~ ", paste(variables_clstr, collapse = " + "))), data = dataset_arm)
Y_hat = predict(Y_model)
print(min(Y_hat))
Y_hat_sort = sort(Y_hat) 
max_alpha = max(Y_hat_sort)
min_alpha = min(Y_hat_sort)
alpha_upper_bound = max_alpha / min_alpha
alpha_lower_bound = min_alpha / max_alpha
hist(Y_hat, breaks = 30); min(Y_hat)

clusters = clustering_func(dataset=dataset, arm=arm, variables_clstr=variables_clstr, num_clstrs=5)
data_with_clstrs = data.frame(cluster = clusters, dataset %>% filter(A==arm))
clstrs_vec = sort(unique(clusters))

max(clusters) == length(clstrs_vec)
max_alpha_per_clstr = min_alpha_per_clstr = mean_pi_tilde = 
  as_amount_vec = second_g_amount_vec = vector(length = max(clusters))
all_Y_hat_vec = c()
for (i in 1:length(clstrs_vec)){
  data_per_clstr = filter(data_with_clstrs, cluster == i)
  hist(data_per_clstr$pi_tilde_as1, breaks = 30)
  mean_pi_tilde[i] = mean(data_per_clstr$pi_tilde_as1)
  
  #mean(data_per_clstr$EMest_p_as) / ( mean(data_per_clstr$EMest_p_as) + mean(data_per_clstr$EMest_p_pro) )
  as_prop_est = ifelse(arm == 1, mean(data_per_clstr$pi_tilde_as1), 1 / (1 + xi))
  Y_sort = sort(data_per_clstr$Y)
    
  if(model_per_cluster == TRUE){
    # regression model per each cluster
    Y_model_clstr = lm(as.formula(paste0("Y~", paste(variables_clstr, collapse="+"))), data=data_per_clstr)
    Y_hat = predict(Y_model_clstr)
    print(min(Y_hat))
    Y_hat = Y_hat[Y_hat>0]
    Y_hat_sort = sort(Y_hat)  
    ##Y_hat_sort = Y_hat[order(data_per_clstr$Y, decreasing = F)]
  }else{ # regression model per all the dataset of treated/untreated survivors
  Y_hat = predict(Y_model, newdata = data_per_clstr)
  print(min(Y_hat))
  Y_hat = Y_hat[Y_hat>0]
  Y_hat_sort = sort(Y_hat)
  }
  
  all_Y_hat_vec = c(all_Y_hat_vec, Y_hat)
  
  as_amount = length(Y_hat_sort) * as_prop_est
  as_amount_vec[i] = as_amount
  second_g_amount_vec[i] = length(Y_hat_sort) - as_amount
  
  max_alpha_per_clstr[i] = calculate_upper_alpha_per_clstr(data_per_clstr=data_per_clstr,
                                                           Y_sort=Y_hat_sort, as_amount=as_amount)
  min_alpha_per_clstr[i] = calculate_lower_alpha_per_clstr(data_per_clstr=data_per_clstr,
                                                           Y_sort=Y_hat_sort, as_amount=as_amount)
  
}

hist(all_Y_hat_vec, breaks = 30)
alpha_upper_bound_clstr = max(max_alpha_per_clstr)
alpha_lower_bound_clstr = min(min_alpha_per_clstr)
1 / min_alpha_per_clstr



