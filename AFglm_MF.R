AFglm_mf = function(object, data, exposure, clusterid, case.control = FALSE){
  call = match.call()
  # Warning if the object is not a glm object
  if(!(as.character(object$call[1]) == "glm"))
    stop("The object is not a glm object", call. = FALSE)
  # Warning if the object is not a logistic regression
  if(!(object$family[1] == "binomial" & object$family[2] == "logit"))
    stop("The object is not a logistic regression", call. = FALSE)
  #### Preparation of dataset ####
  formula = object$formula
  #data <- object$data
  npar = length(object$coef)
  
  ## Delete rows with missing on variables in the model ##
  data = as.data.frame(data)
  rownames(data) = 1:nrow(data)
  m = model.matrix(object = formula, data = data)
  complete = as.numeric(rownames(m))
  data = data[complete, ]
  outcome = as.character(terms(formula)[[2]])
  n = nrow(data)
  n.cases = sum(data[, outcome])
  clusters = data[, clusterid]
  
  if(missing(clusterid)) n.cluster <- 0
  else {
    n.cluster = length(unique(data[, clusterid]))
  }
  
  ## Checks ##
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  
  # Create dataset data0 for counterfactual X = 0
  data0 = data
  data0[, exposure] = 0
  
  ## Design matrices ##
  design = model.matrix(object = delete.response(terms(object)), data = data)
  design0 = model.matrix(object = delete.response(terms(object)), data = data0)
  
  #### Meat: score equations ####
  ## If sampling design is case-control ##
  if (case.control == TRUE){
    ## Create linear predictors to estimate the log odds ratio ##
    diff.design = design0 - design
    linearpredictor = design  %*% coef(object)
    linearpredictor0 = design0 %*% coef(object)
    #log odds ratio#
    log.or = linearpredictor - linearpredictor0
    ## Estimate approximate AF ##
    AF.est   = 1 - sum(data[, outcome] * exp( - log.or)) / sum(data[, outcome])
    #### Meat: score equations ####
    ## Score equation 1 ## individual estimating equations of the estimate of AF
    score.AF = data[, outcome] * (exp( - log.or) - AF.est)
    ## Score equation 2 ## individual estimating equations from conditional logistic reg.
    pred.diff = data[, outcome] - predict(object, newdata = data, type = "response")
    score.beta = design * pred.diff
    score.equations = cbind(score.AF, score.beta)
    if (!missing(clusterid)){
      score.equations = score.equations
      score.equations = aggr(score.equations, clusters = clusters)
    }
    meat <- var(score.equations, na.rm=TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    #### Estimating variance using Sandwich estimator ####
    hessian.AF1 = - data[, outcome]
    hessian.AF2 = (design0 - design) * as.vector(data[, outcome] * exp( - log.or))
    hessian.AF = cbind(mean(hessian.AF1), t(colMeans(hessian.AF2, na.rm = TRUE)))
    hessian.beta = cbind(matrix(rep(0, npar), nrow = npar, ncol = 1), - solve(vcov(object = object)) / n)
    ### Bread ###
    bread = rbind(hessian.AF, hessian.beta)
    #### Sandwich ####
    if (!missing(clusterid))
      sandwich = (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster / n^2 ) [1:2, 1:2]
    else
      sandwich = (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    AF.var = sandwich[1, 1]
    #### Output ####
    out = c(list(AF.est = AF.est, AF.var = AF.var, log.or = log.or,
                  objectcall = object$call, call = call, exposure = exposure, outcome = outcome, object = object,
                  sandwich = sandwich, formula = formula,
                  n = n, n.cases = n.cases, n.cluster = n.cluster))
  }
  ## If sampling design is cross-sectional ##
  else {
    ## Score equation 1 ##
    score.P = data[, outcome]
    pred.Y  = predict(object, newdata = data, type = "response")
    ## Score equation 2 ##
    score.P0 = predict(object, newdata = data0, type = "response")
    ## Score equation 3 ##
    score.beta = design * (score.P - pred.Y)
    ### Meat ###
    score.equations = cbind(score.P, score.P0, score.beta)
    if (!missing(clusterid)){
      score.equations = score.equations
      score.equations = aggr(score.equations, clusters = clusters)
    }
    meat = var(score.equations, na.rm = TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    hessian.P = matrix(c(- 1, 0, rep(0,npar)), nrow = 1, ncol = 2 + npar)
    ## Hessian of score equation 2 ##
    g = family(object)$mu.eta
    dmu.deta = g(predict(object = object, newdata = data0))
    deta.dbeta = design0
    dmu.dbeta = dmu.deta * deta.dbeta
    hessian.P0 = matrix(c(0, - 1, colMeans(dmu.dbeta)), nrow = 1, ncol = 2 + npar)
    ## Hessian of score equation 3 ##
    hessian.beta = cbind(matrix(rep(0, npar * 2), nrow = npar, ncol = 2)
                          , - solve(vcov(object = object)) / n)
    ### Bread ###
    bread = rbind(hessian.P, hessian.P0, hessian.beta)
    #### Sandwich ####
    if (!missing(clusterid))
      sandwich = (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster / n^2 ) [1:2, 1:2]
    else
      sandwich = (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    #### Point estimate of AF ####
    P.est  = mean(score.P, na.rm = TRUE)
    P0.est = mean(score.P0, na.rm = TRUE)
    AF.est = 1 - P0.est / P.est
    ## Delta method for variance estimate ##
    gradient = as.matrix(c(P0.est / P.est ^ 2, - 1 / P.est), nrow = 2, ncol = 1)
    AF.var = t(gradient) %*% sandwich %*% gradient
    P.var = sandwich[1, 1]
    P0.var = sandwich[2, 2]
    
    objectcall = object$call
    #### Output ####
    out = c(list(AF.est = AF.est, AF.var = AF.var, P.est = P.est, P0.est = P0.est, P.var = P.var,
                  P0.var = P0.var, objectcall = objectcall, call = call, exposure = exposure, outcome = outcome,
                  object = object, sandwich = sandwich, gradient = gradient, formula = formula,
                  n = n, n.cases = n.cases, n.cluster = n.cluster))
  }
  class(out) = "AF"
  return(out)
}