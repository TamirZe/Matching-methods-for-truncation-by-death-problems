########################################################################
# PRINT TO LATEX
print(adjustments_for_tables_before_present(final_tables_general$`_S1`) %>% xtable(digits=c(2)),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
no_mis = adjustments_for_tables_before_present(final_tables_general$`_S1`)
ff_mis = adjustments_for_tables_before_present(final_tables_general$`_S1`) %>% subset(select=-c(1:2))
colnames(ff_mis) = paste0("mis_", colnames(ff_mis))
print(cbind(no_mis, ff_mis) %>% 
        filter(!Estimator %in% c("HL No", "BC Yes", "BC caliper inter Yes", "BC inter Yes", "HL Yes")) %>% 
        xtable(digits=c(2)), size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)

print(adjustments_for_tables_before_present(final_tables_crude$`_S1`) %>% xtable(digits=c(2)),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
print(naive_estimators_sum %>% xtable(digits=c(2)),
      size="\\fontsize{11pt}{11pt}\\selectfont", include.rownames=F)
########################################################################