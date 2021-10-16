# m_data from function_my_EM_real_data(data, covariates, iterations)
#data_with_PS = EM_list[[2]]
m_data = data_with_PS
reg_covariates = c("Y", "A", "re74", "age", "married")

# reg_covariates
reg_data =  subset(m_data, select = reg_covariates)
lin_reg = lm(Y~., reg_data_matched)
sum(reg_data$Y==0)
sum = summary(lin_reg)

coeffs = sum$coefficients
X = as.matrix(data.frame(1 ,subset(reg_data, select = reg_covariates[-c(1)])))

# extracting rss
RSS = sum(lin_reg$residuals^2)  # sum(resid(lin_reg)^2)
# more ways for extracting rss
anova(lin_reg); with(summary(lin_reg), df[2] * sigma^2)

S = sqrt(RSS / (nrow(reg_data) - length(reg_covariates)))
S == sigma(lin_reg)
# https://stats.stackexchange.com/questions/57746/what-is-residual-standard-error
# Residual standard error (RSE) from lm is sqrt( RSS / (n-p) ) 
# sigma(lin_reg) is equal to S = RSE.
# rss = S^2 * (n-p) ; p including the intercept
RSS - sigma(lin_reg)^2 * ( nrow(reg_data) - length(reg_covariates) )

# check if in the std err of coeff we need to divide in n or sqrt(n), NOPE!!!!!!!
sqrt(solve(t(X) %*% X)) * S 
sqrt(vcov(lin_reg))
sqrt(diag(vcov(lin_reg))) - sqrt(diag(solve(t(X) %*% X))) * S 
tstats <- coef(lin_reg) / sqrt(diag(vcov(lin_reg)))
tstats - coeffs[,"t value"]
