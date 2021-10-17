# comparing LL and DW

########################################################################
# verify ll contains dw, and missing values in LL
dw_dat = data.frame(subset(nsw, select = -c(data_id, hispanic)), hispanic = nsw$hispanic)
ll_dat = data.frame(subset(LL, select = -c(hispanic, u75)), hispanic = LL$hispanic); colnames(ll_dat)[1] = "treat"
nsw_dup = dw_dat[which(duplicated(dw_dat)),]; ll_dup = ll_dat[which(duplicated(ll_dat)),]
identical(match_df(nsw_dup, ll_dup), nsw_dup)
identical(match_df(ll_dup, nsw_dup), nsw_dup)
match_dw_ll = match_df(dw_dat, ll_dat); identical(dw_dat, match_dw_ll)
duplicated(rbind(dw_dat,ll_dat)) %>% table # extra 15 since there are 15 duplicated rows in the datasets: length(which(duplicated(nsw)))
#duplicated(rbind(dfr_extended_1618,dfr_extended_1617))[(nrow(dfr_extended_1618) + 1): (nrow(dfr_extended_1617) + nrow(dfr_extended_1618))] %>% table
unique(rbind(dw_dat,ll_dat))

# verify the missing values in LL
# in LL and not in DW
library(sqldf)
res <- sqldf('SELECT * FROM ll_dat EXCEPT SELECT * FROM dw_dat')
LL_butnot_DW = anti_join(ll_dat, dw_dat); 722-445 == nrow(LL_butnot_DW); table(LL_butnot_DW$re74, exclude = FALSE)
LL_butnot_DW$emp74 = ifelse(LL_butnot_DW$re74 > 0, 1, 0)
LL_butnot_DW$emp75 = ifelse(LL_butnot_DW$re75 > 0, 1, 0)
LL_butnot_DW$equal_7475 = ifelse(LL_butnot_DW$re74 == LL_butnot_DW$re75, 1, 0)
LL_butnot_DW$equal_7475_earn = ifelse(LL_butnot_DW$equal_7475 == 1 & LL_butnot_DW$re74 !=0, 1, 0)
LL_butnot_DW$equal_7475 == LL_butnot_DW$equal_7475_earn
sum(LL_butnot_DW$equal_7475_earn)
View(subset(LL_butnot_DW, select = -hispanic))

dw_dat$equal_7475 = ifelse(dw_dat$re74 == dw_dat$re75, 1, 0)
dw_dat$equal_7475_earn = ifelse(dw_dat$equal_7475 == 1 & dw_dat$re74 !=0, 1, 0)
sum(dw_dat$equal_7475); sum(dw_dat$equal_7475_earn)
ll_dat$equal_7475 = ifelse(ll_dat$re74 == ll_dat$re75, 1, 0)
ll_dat$equal_7475_earn = ifelse(ll_dat$equal_7475 == 1 & ll_dat$re74 !=0, 1, 0)
sum(ll_dat$equal_7475); sum(ll_dat$equal_7475_earn)

table(nsw$re74); table(LL$re74)
unique(LL$re74)
sum(LL$re74 == LL$re75); sum(LL$re74[LL$re74!=0] == LL$re75[LL$re74!=0])
sum(nsw$re74 == nsw$re75); sum(filter(nsw, re74!=0)$re74 == filter(nsw, re74!=0)$re75)
########################################################################


########################################################################
# READING TXT FILES # https://users.nber.org/~rdehejia/data/.nswdata2.html
# read from LL DW Datasets folder
# These are text files. treatment indicator (1 if treated, 0 if not treated),  age, education, Black (1 if black, 0 otherwise),
# Hispanic (1 if Hispanic, 0 otherwise), married (1 if married, 0 otherwise), nodegree (1 if no degree, 0 otherwise), RE75 (earnings in 1975), and RE78 (earnings in 1978). The last variable is the outcome; other variables are pre-treatment.
LL_cntl = read.table("nsw_control.txt", header = FALSE, sep = "", dec = ".")
LL_trt = read.table("nsw_treated.txt", header = FALSE, sep = "", dec = ".")
DW_cntl = read.table("nswre74_control.txt", header = FALSE, sep = "", dec = ".")
DW_trt = read.table("nswre74_treated.txt", header = FALSE, sep = "", dec = ".")

mean(LL_trt[,ncol(LL_trt)]) - mean(LL_cntl[,ncol(LL_cntl)])
mean(LL_trt[LL_trt$V9>0,ncol(LL_trt)]) - mean(LL_cntl[LL_cntl$V9>0, ncol(LL_cntl)])

LL_data = rbind(LL_cntl, LL_trt)
DW_data = rbind(DW_cntl, DW_trt) %>% subset(select = -V8)
colnames(DW_data) = colnames(LL_data)
LL_butnot_DW = anti_join(LL_data, DW_data)
LL_butnot_DW2 = anti_join(ll_dat, dw_dat); 722-445 == nrow(LL_butnot_DW); table(LL_butnot_DW$re74, exclude = FALSE)
LL_butnot_DW2 = subset(LL_butnot_DW2, select = -c(re74, u74))
LL_butnot_DW2 = data.frame(LL_butnot_DW2[,1:4], subset(LL_butnot_DW2, select = c(hisp, married, nodegree, re75, re78)))
colnames(LL_butnot_DW) = colnames(LL_butnot_DW2)
LL_butnot_DW = data.frame(apply(LL_butnot_DW, 2, round, 2)); LL_butnot_DW2 = data.frame(apply(LL_butnot_DW2, 2, round, 2))
chck = anti_join(LL_butnot_DW, LL_butnot_DW2)
hist(LL_butnot_DW$re75); hist(LL_butnot_DW2$re75)

# regular files
library(readstata13)
nsw <- read.dta13("process NSW data/data/nsw_dw.dta")
library(cem)
data(LL, package = "cem")
# DW check
nsw2 = subset(nsw, select = -c(data_id, u74))
DW_data2 = rbind(DW_cntl, DW_trt)
colnames(DW_data2) = colnames(nsw2)
nsw2 = data.frame(apply(nsw2, 2, round, 2)); DW_data2 = data.frame(apply(DW_data2, 2, round, 2))
hist(nsw2$re75); hist(DW_data2$re75)
DW_chck = anti_join(nsw2, DW_data2)

#LL check
LL2 = subset(LL, select = c("treated","age","education","black","hispanic","married","nodegree","re75","re78")) %>% apply(2, round, 2) %>% data.frame()
LL_data = data.frame(apply(rbind(LL_cntl, LL_trt), 2, round, 2))
colnames(LL_data) = colnames(LL2)
hist(LL2$re78); hist(LL_data$re78)
table(LL2$re78) == table(LL_data$re78)
LL_chck = anti_join(LL2, LL_data)
########################################################################