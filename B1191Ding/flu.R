##Preprocessing the data
fl = read.table("~/Desktop/DingLu2015JRSSB_Code/fludata.txt", header=TRUE, quote="\"")
Z = fl$assign 
D = fl$receive
Y = fl$outcome
X = as.matrix(fl[, -c(1, 2, 3)])
N = length(Z)

##Estimation under GPI
require(nnet)
point = PSPS_M_weighting(Z, D, X, Y)
#CACE, AACE and NACE
point$CACE.reg 
point$AACE.reg
point$NACE.reg 

##sensitivity analysis for GPI
#for example let \epsilon_1 = 0.5 and \epsilon_0 = 1.5
point = PSPS_M_weighting(Z, D, X, Y, ep1 = 0.5, ep0 = 1.5)
#CACE, AACE and NACE
point$CACE.reg 
point$AACE.reg
point$NACE.reg

##Notes
#1. For variance estimation, apply function ``PSPS_M_weighting'' on bootstraped samples
#2. For balancing checking, do for example ``PSPS_M_weighting(Z, D, X, X[, 1])

