#=================================================

# File: Build OR(95% Ci) Table
# Author: Joy Fu
# Update Date: 07/24/2019

#=================================================


#=================================================
#   Quick Reference

# Function to build OR(95% CI) table:
# extract_OR = function(regression_object, variable_index_lst, mark)
# make_OR_table = function(lst_model, var_position, model_head)


# Function to build coefficient table:
# extract_coeff = function(regression_object, variable_index_lst, mark) 
# make_coeff_table = function(lst_model, model_head, type = c('normal', 'interact'), var_position = NULL)
#=================================================


extract_OR = function(regression_object, variable_index_lst, mark) {
  
  summary_reg = summary(regression_object)
  coefficient = summary_reg$coefficients
  
  or_table = array(NA, dim=c(length(variable_index_lst),8))
  
  for(i in 1:length(variable_index_lst)) {
    
    index = variable_index_lst[i]
    OR = exp(coefficient[index])
    
    confinterval = exp(confint(regression_object, level = 0.95))
    num_variable = length(confinterval)/2
    lower_CI = confinterval[index]
    upper_CI = confinterval[num_variable + index]
    
    CI_95 = paste0('(', sprintf('%.2f',lower_CI), ', ', sprintf('%.2f',upper_CI), ')')
    p = round(coefficient[num_variable*3 + index], 4)
    
    or_table[i,1] = mark
    or_table[i,2] = index
    or_table[i,3] = length(regression_object$y)
    or_table[i,4] = round(OR, 2)
    or_table[i,5] = CI_95
    or_table[i,6] = sprintf('%.3f',p)
    or_table[i,7] = round(lower_CI, 2)
    or_table[i,8] = round(upper_CI, 2)
    
    i = i + 1
  }
  
  or_table_final = as.data.frame(or_table)
  colnames(or_table_final) = c('mark','index', 'N', 'OR', '95% CI', 'p-value', 'lower_CI', 'upper_CI')
  
  return(or_table_final)
}


make_OR_table = function(lst_model, var_position, model_head) {

  # Make to a full table
  or_table = data.frame(
    Mark = character(),
    Index = integer(),
    N = integer(),
    OR = double(),
    CI_95 = character(),
    p_value = character(), 
    lower_CI = double(),
    upper_CI = double(),
    stringsAsFactors = F
  )
  
  for (i in 1:length(lst_model)) {
    
    model_object = get(lst_model[i])
    result = extract_OR(model_object, var_position, paste0(model_head, lst_model[i]))
    
    or_table = rbind(or_table, result)
    
    i = i + 1
  }
  
  return(or_table)
}

#=====================================================
#
# Additional: Extract the coefficients for regression
#
#=====================================================

extract_coeff = function(regression_object, variable_index_lst, mark) {
  
  summary_reg = summary(regression_object)
  coefficient = summary_reg$coefficients
  
  coef_table = array(NA, dim=c(length(variable_index_lst),5))
  
  for(i in 1:length(variable_index_lst)) {
    
    index = variable_index_lst[i]
    coeff = coefficient[index]
    confinterval = exp(confint(regression_object, level = 0.95))
    num_variable = length(confinterval)/2
    p = round(coefficient[num_variable*3 + index], 4)
    
    coef_table[i,1] = mark
    coef_table[i,2] = index
    coef_table[i,3] = length(regression_object$y)
    coef_table[i,4] = sprintf('%.3f',coeff)
    coef_table[i,5] = sprintf('%.3f',p)
    
    i = i + 1
  }
  
  coef_table_final = as.data.frame(coef_table)
  colnames(coef_table_final) = c('mark','index', 'N', 'coeff', 'p-value')
  
  return(coef_table_final)
}


make_coeff_table = function(lst_model, model_head, type = c('normal', 'interact'), var_position = NULL) {
  
  # Make to a full table
  coeff_table = data.frame(
    Mark = character(),
    Index = integer(),
    N = integer(),
    coeff = double(),
    p_value = double(), 
    stringsAsFactors = F
  )
  
  for (i in 1:length(lst_model)) {
    
    model_object = get(lst_model[i])
    try = summary(model_object)
    if (type == "interact") {
      var_position = length(try$coefficients)/4
    } 
    
    result = extract_coeff(model_object, var_position, paste0(model_head, lst_model[i]))
    coeff_table = rbind(coeff_table, result)
    
    i = i + 1
  }
  
  return(coeff_table)
}

