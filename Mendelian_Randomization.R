#=================================================

# File: Mendelian Randomization
# Author: Joy Fu
# Update Date: 07/24/2019

#=================================================

#=================================================
#   Quick Reference

# Function to imply MR -- get the values for MR

#=================================================

library(MendelianRandomization)

# Get the values
get_MR_value = function(reg_object1, reg_object2, mark) {
  
  summary_reg_1 = summary(reg_object1)
  summary_reg_2 = summary(reg_object2)
  coefficient_1 = summary_reg_1$coefficients
  coefficient_2 = summary_reg_2$coefficients
  
  bx = coefficient_1[2]
  bxse = coefficient_1[4]
  by = coefficient_2[2]
  byse = coefficient_2[4]
  
  MRdata_smk = mr_input(bx = bx, bxse = bxse, by = by, byse = byse)
  
  IVWObject = mr_ivw(MRdata_smk,
                     model = "default",
                     robust = FALSE,
                     penalized = FALSE,
                     correl = FALSE,
                     weights = "simple",
                     psi = 0,
                     distribution = "normal",
                     alpha = 0.05)
  
  p_value = IVWObject$Pvalue
  OR = exp(IVWObject@Estimate)
  lower = exp(IVWObject@CILower)
  upper = exp(IVWObject@CIUpper)
  CI_95 = paste0('(', sprintf('%.2f',lower), ', ', sprintf('%.2f',upper), ')')
  
  # Make to a series
  mr_series = c(mark, sprintf('%.2f',OR), CI_95, sprintf('%.3f',p_value))
  
  return(mr_series)
}