######################## main script ########################
# load all libraries with apply shortcut...
library(rlist); library(locfit); library(plyr); library(nnet); library(xtable); 
library(matrixStats); library(data.table); library(dplyr); library(plyr);
library(reshape); library(MASS); library(ggplot2); library(rockchalk); library(data.table)
source("simple_estimators.R")

library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
#foreach(j = seq_along(numeric(length(c(1:n_sim)))), .combine=c) %dopar%{
stopCluster(cl)

source("matching_scripts/matching_PS.R")

# TODO I won't use the flu data, although it appears in Ding and lu 2017,
# since this is an encouragement experiment with non-compliance
# flu_data = read.table("fludata.txt")
# str(flu_data)

sg = read.table(file = 'swogdata.txt')

# Ding and Lu code, from swog.r
Z = sg$Z
N = length(Z)
# TODO D is probably the survival status. If score12 is na, the patient didnt survive?
D = ifelse(is.na(sg$score12), 0, 1)
View(data.frame(score12 = sg$score12, D = D))
Y = sg$score12 - sg$score0
# where the subject didnt survive, report Y as 0 for now.
Y[is.na(Y)] = 0
X = as.matrix(cbind(sg$AGE, sg$RACEB, sg$RACEO, sg$score0))
sg_data = data.frame(id = c(1:(nrow(sg))), AGE = sg$AGE, RACEB = sg$RACEB,
    RACEO = sg$RACEO, score0 = sg$score0, A = Z, S = D, Y = Y)

# simple estimators
simple_estimators = some_simple_estimators(sg_data)
