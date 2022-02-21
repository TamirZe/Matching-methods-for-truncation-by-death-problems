#dfr = list_all_CI[[3]]
calculate_coverage = function(dfr){
  # coverage
  bool_CI_contain_SACE = function(x){
    lower_bound =  as.numeric(sub(",.*", "", x)) #as.numeric(sapply(strsplit(as.character(x), ","), "[", 1)) 
    upper_bound = as.numeric(sub(".*,", "", x))
    bool = dfr$SACE >= lower_bound & dfr$SACE <= upper_bound 
    bool = ifelse(bool==T,1,0)
    return(bool)
  }
  rows = rownames(dfr)
  bool_mat = data.frame(apply(dfr[,-1], 2, bool_CI_contain_SACE))
  colnames(bool_mat) = paste0("Coverage_", colnames(dfr[,-1]))
  bool_mat = data.frame(rbind(bool_mat, apply(bool_mat, 2, mean)))
  bool_mat = data.frame(rbind(dfr,101), bool_mat)
  rownames(bool_mat) = c(paste0("Coverage", rows), paste0("Coverage_Mean"))
  coverage = subset(bool_mat, select = grep("Coverage", colnames(bool_mat)))["Coverage_Mean",]
  return(list(bool_mat=bool_mat, coverage=coverage))
}




