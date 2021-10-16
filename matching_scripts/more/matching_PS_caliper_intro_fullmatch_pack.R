# matching on Principal score
library(optmatch)
library(DOS)
# examples in DOS cran and in rosenbaum book: design of obs studies
#@@@ this is the DOS cran example @@@
data(costa)
z<-1*(costa$welder=="Y")
aa<-1*(costa$race=="A")
smoker=1*(costa$smoker=="Y")
age<-costa$age
x<-cbind(age,aa,smoker)
controls = length(z[z==0]); treated = length(z[z==1])
dmat<-mahal(z,x)
d = data.frame(z,x)
# Mahalanobis distances
round(dmat[,1:6],2) # Compare with Table 8.5 in Design of Observational Studies (2010)
# Impose propensity score calipers
prop<-glm(z~age+aa+smoker,family=binomial, data = d)$fitted.values # propensity score
# Mahalanobis distanced penalized for violations of a propensity score caliper.
# This version is used for numerical work.
dmat2<-addcaliper(dmat,z,prop,caliper=0.1)
round(dmat[,1:6],2) # Compare with Table 8.5 in Design of Observational Studies (2010)
## Not run:
# Find the minimum distance match within propensity score calipers.
pair_match = optmatch::pairmatch(dmat2,data=costa)
## End(Not run)
# Conceptual versions with infinite distances for violations of propensity caliper.
dmat[dmat>20]<-Inf
round(dmat[,1:6],2) # Compare with Table 8.5 in Design of Observational Studies (2010)

ctrl_omit_frac= ( length(d$z[d$z==0]) - 1 * length(d$z[d$z==1]) ) / length(d$z[d$z==0])
m<-fullmatch(dmat2, min.controls=1, max.controls=1, omit.fraction=ctrl_omit_frac)

#  omit.fraction=ctrl_omit_frac
length(m)
m
matched(m)
sum(matched(m))
#m_data = data_with_PS
m_data = data_with_PS[OBS != "O(0,0)"]
# m_data = data_with_PS[S==1]
# m_data = data_with_PS
#min_PS_weighted_match = 0.5
#min_diff_PS = 0.15
# caliper is in sd

library("Matching")

my_matching_func_rosenbaum = function(X_sub_cols, m_data, weighting = FALSE,
                                      M=1, replace, estimand = "ATC", mahal_match = 2,
                                      min_PS = 0, min_diff_PS, caliper = 0.1, min_PS_weighted_match, OBS_table){
  X_sub_cols = paste0("X", c(1:(dim_x)))
  # mahal_match for Weight = 2 for mahalanobis distance. 1 for inverse of variance
  # TODO find a way to use caliper on the PS in the matching function
  ATE_MATCH_PS  <- Match(Y=m_data[,Y], Tr=m_data[,A], 
                         X=subset(m_data, select = X_sub_cols[-1])
                         #,caliper = 
                         ,M=M, replace = replace, estimand = estimand, Weight = mahal_match)
  
  print(ATE_MATCH_PS$estimand)
  #ATE_MATCH_PS$est; ATE_MATCH_PS$se; ATE_MATCH_PS$nobs; ATE_MATCH_PS$index.dropped
  # TODO matching based on proncipal scores
  
  
  x <- subset(m_data, select = X_sub_cols[-1])
  dmat<-mahal(m_data$A,x)
  # Mahalanobis distanced penalized for violations of a propensity score caliper.
  # This version is used for numerical work.
  dmat2<-addcaliper(dmat, m_data$A, m_data$est_p_as, caliper=0.1
                    #, penalty = 
  )
  controls = ncol(dmat2); treated = nrow(dmat2)
  ctrl_omit_frac = ( controls - 1 * treated ) / controls
  #ctrl_omit_frac = ( length(m_data$A[m_data$A==0]) - 1 * length(m_data$A[m_data$A==1]) ) / 
  #  length(m_data$A[m_data$A==0])
  matching <-fullmatch(dmat2, min.controls=0, max.controls=1
                       # , omit.fraction=ctrl_omit_frac
  )
  
  #  omit.fraction=ctrl_omit_frac
  length(matching)
  table(matching)
  View(matching)
  matched(matching)
  sum(matched(matching)) == controls + treated
  
  
  # TODO seems good @@@
  f = as.formula(paste0("A ~ ", paste(X_sub_cols[-1], collapse = " + ")
                        #, " + scores(est_p_as)"
  ))
  principal_score = m_data$est_p_as
  names(principal_score) <- rownames(m_data)
  mhd <- match_on(f, data = m_data) +
    caliper(match_on(principal_score, z = m_data$A), 1)
  pm2 <- pairmatch(mhd, data = m_data) 
  summary(pm2)
  pm2 = na.omit(pm2)
  matched_data = m_data
  matched_data$id = c(1:nrow(matched_data))
  matched_data = matched_data[as.numeric(names(pm2)) , ]
  matched_data$matched_pair = as.numeric(pm2)
  matched_data = arrange(matched_data, matched_pair)
  
  #@
  mhd1 <- match_on(f, data = m_data, method = "mahalanobis")
  all.equal(mhd1, match_on(f, data = m_data), check.attributes=FALSE)
  mhpc.pm <- pairmatch(mhd1, controls = 1, caliper=1, data = m_data)
  summary(mhpc.pm) 
  
  #@
  x <- subset(m_data, select = X_sub_cols[-1])
  dmat<-mahal(m_data$A,x)
  dmat2<-addcaliper(dmat, m_data$A, m_data$est_p_as, caliper=0.1
                    #, penalty = 
  )
  pair_match <- pairmatch(dmat2, controls = 1, data = m_data)
  summary(pair_match) # oops
  class(pair_match)
  
  
  #@
  mhd <- match_on(dmat2, data = m_data)
  pm2 <- pairmatch(mhd, data = m_data) 
  summary(pm2)
  
  ######################################################################
  ######## calculating the amount of as the matching process excluded
  OBS_table
  as_A0_matched = length(which(dt_match_min_ps$A0_g == "as"))
  as_A0_unmatched = OBS_table[1,2] - as_A0_matched
  as_A1_matched = length(which(dt_match_min_ps$g == "as"))
  pro_A1_matched = length(which(dt_match_min_ps$g == "pro"))
  
  # check the compisition in A1_S1: as and protected
  as_in_A1_S1 = m_data[g =="as" & A == 1,]
  pro_in_A1_S1 = m_data[g =="pro" & A == 1,]
  # senity check
  nrow(as_in_A1_S1) + nrow(pro_in_A1_S1) == OBS_table[2,2]
  as_A1_unmatched = nrow(as_in_A1_S1) - as_A1_matched
  included_excluded_in_matching = 
    data.frame(as_A0_matched, as_A0_unmatched, as_A1_matched, as_A1_unmatched, pro_A1_matched)
  ######################################################################
  
  
  
  
  return(list(SACE_matching_est, SACE_matching_PS_weighted_match_est, included_excluded_in_matching))
}




# Matching in R using the optmatch and RItools packages
# Ben B. Hansen, Mark Fredrickson, Josh Buckner, Josh Errickson, 
# and Peter Solenberger, with embedded Fortran code 
# due to Dimitri P. Bertsekas and Paul Tseng 2019-12-06
library(optmatch)
library(RItools)
data(nuclearplants)
head(nuclearplants)
nuclearplants$pt
table(nuclearplants$pt)
with(nuclearplants, table(pt))
nuke.nopt <- subset(nuclearplants, pt == 0)
head(nuke.nopt)
# (pr == 1) as the treatment group and (pr == 0) as comparisons.
table(nuke.nopt$pr)
pairmatch(pr ~ cap, data = nuke.nopt)
print(pairmatch(pr ~ cap, data = nuke.nopt), grouped = TRUE)
psm <- glm(pr ~ date + t1 + t2 + cap + ne + ct + bw + cum.n + pt,
           family = binomial, data = nuclearplants)
table(nuclearplants$pr)
mhd1 <- match_on(pr ~ date + cap + scores(psm), data=nuclearplants)
mhpc.pm <- pairmatch(mhd1, caliper=1, data=nuclearplants)
mhpc.pm
summary(mhpc.pm) # oops
