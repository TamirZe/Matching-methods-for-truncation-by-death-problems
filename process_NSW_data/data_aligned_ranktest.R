# adjust new groups
adjust_pairs_to_new_grp = function(dataset){
  data_new_grp = arrange(dataset, pair,1-A) %>% data.table()
  # BY PAIR, ADD THE NUMBER OF THE TREATED and REMOVE DUPLICATIONS
  data_new_grp = data_new_grp[, `:=` (trt_grp = id[A==1]), by="pair"] %>% 
    arrange(trt_grp) %>% group_by(id) %>% slice(1) %>% arrange(trt_grp, 1-A)
  data_new_grp = data.frame(subset(data_new_grp, select = trt_grp), subset(data_new_grp, select = -trt_grp))
  data_new_grp$trt_grp = as.numeric(as.factor(data_new_grp$trt_grp))
  return(data_new_grp)
}

# Matching Methods for Observational Microarray Studies, 2008
#Ruth Heller, Elisabetta Manduchi and Dylan Small 
# Function for computing aligned rank test for full matching
# Inputs are outcome vector, a vector which says which matched set 
# each unit belongs to and a vector which says whether the unit is treated (1)
# or control (0)
# Computes two-sided and one sided (larger) p-value

alignedranktest=function(outcome,matchedset,treatment){
  # Compute means in each matched set
  matchedset.mean=tapply(outcome,matchedset,mean);
  # Compute residuals
  matchedset.mean.expand=matchedset.mean[matchedset];
  resids=outcome-matchedset.mean.expand;
  # Rank the residuals
  rankresids=rank(resids);
  
  # Test statistics = Sum of residuals in treatment group
  # sum of the ranks sums after alignment of all matched set within trt 
  teststat=sum(rankresids[treatment==1]);
  
  # Expected value and variance of test statistic
  # Expected value
  mean.matchedset.rankresids=tapply(rankresids,matchedset,mean);
  notreated.matchedset=tapply(treatment,matchedset,sum);
  nocontrol.matchedset=tapply(1-treatment,matchedset,sum);
  no.matchedset=notreated.matchedset+nocontrol.matchedset;
  ev.teststat=sum(mean.matchedset.rankresids*notreated.matchedset);
  
  # var
  mean.matchedset.rankresids.expand=mean.matchedset.rankresids[matchedset];
  rankresids.resid.squared=(rankresids-mean.matchedset.rankresids.expand)^2; 
  squared.resids.sum=tapply(rankresids.resid.squared,matchedset,sum);
  temp = ((notreated.matchedset*nocontrol.matchedset)/(no.matchedset*(no.matchedset-1)))*squared.resids.sum;
  var.teststat=sum(temp,na.rm=T);
  
  # two sided (-abs to limit by 0.5 before multiplying by 2)
  pval_2sided=2*pnorm( -abs((teststat-ev.teststat)/sqrt(var.teststat)) ); # pnorm: distribution function
  pval_2sided;
  
  pval_larger=1-pnorm( (teststat-ev.teststat)/sqrt(var.teststat) ) 
  pval_larger;
  return(list(pval_2sided=pval_2sided, pval_larger=pval_larger))
}

