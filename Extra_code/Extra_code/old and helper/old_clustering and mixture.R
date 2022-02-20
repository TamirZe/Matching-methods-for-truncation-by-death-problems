# 445 obs
library(Matching)
data(lalonde)

dim(lalonde)
data = lalonde
data$u74 = 1 - data$u74
data$u75 = 1 - data$u75
colnames(data)[which(colnames(data) %in% c("u74", "u75"))] = c("employed74", "employed75")
data$employed78 = ifelse(data$re78 == 0, 0, 1)

treat_employed = filter(data, treat == 1, employed78 == 1) # GG GD
treat_unemployed = filter(data, treat == 1, employed78 == 0 ) # DD DG
control_employed = filter(data, treat == 0, employed78 == 1) # GG DG
control_unemployed = filter(data, treat == 0, employed78 == 0) # DD GD
dim(treat_employed)

##########################################################################################################3
# k mean clustering
# TODO a better idea is Gaussian Mixture (EM algo)
df = control_employed
df_scaled <- subset( scale(df), select = -c(treat, employed78, re78, employed74, employed75, re74, re75) )
# df_scaled <- subset( scale(df), select = -c(treat, employed78, re78) )
theoretical_size_monotonicity = (p_CG / p_TG) * nrow(df_scaled)

set.seed(101)
fit <- kmeans(df_scaled, 2)
size = fit$size # 99 is less than theoretical_size_monotonicity
#NaRV.omit(df_scaled)
# get cluster means 
aggregate(df_scaled, by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(df, fit$cluster)
a = filter(mydata, fit.cluster == 1)
b = filter(mydata, fit.cluster == 2)
options(scipen = 2)
#format(x, scientific=TRUE)
apply(a, 2, mean); apply(b, 2, mean)

# Ward Hierarchical Clustering
d <- dist(df, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=5, border="red")
##########################################################################################################



##########################################################################################################
install.packages("Matching")
library(Matching)
# model for PS
glm1  <- glm(treat ~  age + educ +  black + 
               hisp + married + re74 + re75 + (1-employed74) +
               I(age^2) + I(educ^2) + I(re74^2) + I(re75^2) + 
               black*(1-employed74), family=binomial, data = data)

############# Match data  #############################################
# X is the estimated (by model) to PS
X  <- glm1$fitted
Y  <- data$re78
Tr  <- data$treat
########################################################################

########################################################################
# TOY example
ATE_MATCH_PS  <- Match(Y=Y, Tr=Tr, X=X, M=1, replace=TRUE, estimand = "ATE")
ATE_MATCH_PS$est
ATE_MATCH_PS$se
ATE_MATCH_PS$nobs
cbind(ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control)
ATE_MATCH_PS$ndrops
ps_comparison = rbind(X[ATE_MATCH_PS$index.treated],X[ATE_MATCH_PS$index.control])
##########################################################################


##########################################################################################################
ATE_MATCH_PS_exact  <- Match(Y=Y, Tr=Tr, X=X, estimand = "ATT", replace = TRUE, exact = T)
matched_pairs_exact = cbind(ATE_MATCH_PS_exact$index.treated, ATE_MATCH_PS_exact$index.control)
matched_pairs_employment_status = 
  data.frame(cbind(matched_pairs_exact, data$employed78[ATE_MATCH_PS_exact$index.treated],
        data$employed78[ATE_MATCH_PS_exact$index.control]))
colnames(matched_pairs_employment_status) = c("trt", "ctr", "S_trt", "S_ctr")
first_way_matched_and_employed78 = filter(matched_pairs_employment_status, S_trt == S_ctr & S_trt == 1)
diff_scale_estimator = mean(data$re78[first_way_matched_and_employed78$trt]) -
  mean(data$re78[first_way_matched_and_employed78$ctr])

'''ATE_MATCH_PS_exact$ndrops
length(unique(ATE_MATCH_PS_exact$index.treated))
ATE_MATCH_PS_exact$index.treated[order(ATE_MATCH_PS_exact$index.treated)]'''

# second way to fimd matche pairs and their S status
ps_comparison_exact = rbind(X[ATE_MATCH_PS_exact$index.treated],X[ATE_MATCH_PS_exact$index.control])

matched_pairs_exact[,1] == as.numeric(as.character(colnames(ps_comparison_exact)))

ps_employed78_comparison_exact = 
  data.frame(rbind(data$employed78[ATE_MATCH_PS_exact$index.treated],
                   data$employed78[ATE_MATCH_PS_exact$index.control]))

ind_S0_S1_equal_1 = which(data$employed78[ATE_MATCH_PS_exact$index.treated] == 
                    data$employed78[ATE_MATCH_PS_exact$index.control] &
                    data$employed78[ATE_MATCH_PS_exact$index.control] == 1)

second_way = matched_and_employed78 = matched_pairs_exact[ind_S0_S1_equal_1,]
View(matched_and_employed78)
##########################################################################################################

##########################################################################################################
# function that match the data by PS model, and estimate the diff Y ONLY for pairs with both S = 1
x = X; y = Y
trt =  Tr; df =  PO_table;  excat_bool = FALSE
replace_bool = FALSE; estimand =  "ATT"
match_by_PS_and_retain_LL = function(x, y, trt, df, excat_bool = FALSE, replace_bool, estimand){
  set.seed(102)
  ATE_MATCH_PS  <- Match(Y=y, Tr=trt, X=x, estimand = estimand, exact = excat_bool, replace = replace_bool)
  matched_pairs = cbind(ATE_MATCH_PS$index.treated, ATE_MATCH_PS$index.control)
  matched_pairs_employment_status = 
    data.frame(cbind(matched_pairs, df$employed78[ATE_MATCH_PS$index.treated],
                                    df$employed78[ATE_MATCH_PS$index.control]))
  colnames(matched_pairs_employment_status) = c("trt", "ctr", "S_trt", "S_ctr")
  first_way_matched_and_employed78 = filter(matched_pairs_employment_status, S_trt == S_ctr & S_trt == 1)
  diff_scale_estimator = mean(df$re78[first_way_matched_and_employed78$trt]) -
    mean(df$re78[first_way_matched_and_employed78$ctr])
  ps_comparison = rbind(X[ATE_MATCH_PS$index.treated],X[ATE_MATCH_PS$index.control])
  return(list(diff_scale_estimator, matched_pairs_employment_status, first_way_matched_and_employed78,
              ATE_MATCH_PS$nobs))
}
##########################################################################################################
x = X; y = Y; trt = Tr; df = data; excat = TRUE; estimand = "ATT"; replace_bool = TRUE
check1 = match_by_PS_and_retain_LL(x, y, trt, df, excat_bool = TRUE, replace_bool, estimand)
check2 = match_by_PS_and_retain_LL(x, y, trt, df, excat_bool = FALSE, replace_bool, estimand)
check3 = match_by_PS_and_retain_LL(x, y, trt, df, excat_bool = TRUE, replace_bool, estimand =  "ATE")
check4 = match_by_PS_and_retain_LL(x, y, trt, df, excat_bool = FALSE, replace_bool, estimand = "ATE")
check5 = match_by_PS_and_retain_LL(x, y, trt, df, excat_bool = TRUE, replace_bool = FALSE, estimand =  "ATT")
check6 = match_by_PS_and_retain_LL(x, y, trt, df, excat_bool = FALSE, replace_bool = FALSE, estimand = "ATT")

check1[[1]]; check2[[1]]
check3[[1]]; check4[[1]]
check5[[1]]; check6[[1]]
