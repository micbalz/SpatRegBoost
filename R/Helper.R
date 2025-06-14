# Generate spatial weight matrices of networks
network = function(N, k) {
  # N: number of units
  # k: number of neighbors ahead (and behind)
  
  W = matrix(0, nrow = N, ncol = N)
  
  for (i in 1:N) {
    for (j in 1:k) {
      ahead = (i + j - 1) %% N + 1
      behind = (i - j - 1) %% N + 1
      
      W[i, ahead] = 1
      W[i, behind] = 1
    }
  }
  
  # Normalize rows to sum to 1
  W = sweep(W, 1, rowSums(W), "/")
  
  return(W)
}

# Extract estimates
params = function(res, method, true_vals) {
  coef = res[[method]]
  wanted = names(true_vals)
  out = coef[wanted]
  out[setdiff(wanted, names(out))] = NA  
  return(out[wanted])
}

# Compute bias, MSE, and SE
bias_mse_empse = function(estimates_df, true_vals) {
  result = sapply(names(true_vals), function(param) {
    estimates = estimates_df[[param]]
    truth = true_vals[param]
    
    bias = mean(estimates, na.rm = TRUE) - truth
    mse  = mean((estimates - truth)^2, na.rm = TRUE)
    se   = sd(estimates, na.rm = TRUE)
    
    c(Bias = unname(bias), MSE = mse, SE = se)
  })
  
  t(result)
}

# Define prediction performance metrics
rmse = function(actual, predicted) sqrt(mean((actual - predicted)^2))

mae = function(actual, predicted) mean(abs(actual - predicted))

nll = function(N, W, Y, X, lambda, beta, sigma) {
  # Compute the determinants of (I - lambda * M)
  R = Diagonal(nrow(W)) - lambda * W
  det_R = determinant(R, logarithm = TRUE)
  log_det_R = as.numeric(det_R$modulus)
  
  # Compute the log-likelihood
  ll = -N / 2 * (log(2 * pi*sigma) + 1) +
    log_det_R  - (1 / (2*sigma)) * t(R %*% (Y - X %*% beta)) %*% R %*% (Y - X %*% beta)
  
  return(as.numeric(-ll))  
}

# Deselect algorithm

DeselectBoost = function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  
  tau = ifelse(is.null(tau), 0.01, tau)
  
  if(any(class(object) %in% 'mboost')){
    DeselectBoost_1(object, data = data, fam = fam, tau = tau, method = method[1])
  }else{
    switch(length(names(object))-1,{
      DeselectBoostLSS_2(object, data = data, fam = fam, tau = tau, method = method[1])},{
        DeselectBoostLSS_3(object, data = data, fam = fam, tau = tau, method = method[1])},{
          DeselectBoostLSS_4(object, data = data, fam = fam, tau = tau, method = method[1])  
        })
  }
}  

DeselectBoost_1 = function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  
  if(is.null(data) && class(object$model.frame()) == 'list'){return(stop("Please enter the data."))
  } else if(!is.null(data)){
    data = data
  }else{
    data = object$model.frame()
  }
  
  nameVar = names(coef(object,which = ''))[-1]
  which.response = which(sapply(1:(dim(data)[2]), function(x){identical(as.numeric(data[,x]), as.numeric(object$response))}))
  name.response = colnames(data)[which.response]
  
  mstop = object$mstop()
  RiskRed = object$risk()
  totalRiskRed = RiskRed[1] - RiskRed[mstop+1] 
  diffRiskRed = sapply(seq(1:(mstop)), function(k){RiskRed[k]-RiskRed[k+1]})
  
  if(any(class(object) %in% "glmboost")){
    select = selected(object) - 1
    diffRiskRed = diffRiskRed[selected(object)-1 != 0]
  }else{
    select = selected(object)
  }
  
  select = select[select != 0]
  Var = plyr::count(select)[[1]]
  Risk.Var = lapply(1:length(Var),function(j){sum(diffRiskRed[which(plyr::count(select)[[1]][j] == select)])})
  
  n.parameter = c(names(object$coef()))
  if('(Intercept)' %in% n.parameter) n.parameter = n.parameter[-which(n.parameter == '(Intercept)')]
  
  
  Risk.order = data.frame(Var,n.parameter, as.numeric(Risk.Var))
  Risk.order = Risk.order[order(Risk.order$as.numeric.Risk.Var.),]
  Risk.order$CumRisk = cumsum(Risk.order$as.numeric.Risk.Var.)
  colnames(Risk.order) = c( 'Var', 'VarName', 'Risk', 'CumRisk')
  
  perc = ifelse(is.null(tau), 0.01, tau) 
  percRiskRed = totalRiskRed * perc
  if(method[1] == 'attributable'){RiskRedOver = Risk.order[which(Risk.order$Risk > percRiskRed),]
  }else{RiskRedOver = Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(plyr::empty(RiskRedOver)){form2 = as.formula(paste(name.response, "~ 1"))
  }else{form2 =as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName, collapse= "+")))
  if(!is.null(environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]])){
    dfbase = environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]]
  }}
  #if(is.null(object$call$family)){ fam =  Gaussian()
  #}else fam = eval(parse(text = object$call$family))
  
  if(any(class(object) %in% "glmboost")){
    model_after = glmboost(form2, data = data, weights = model.weights(object), family = fam, control = boost_control(mstop = mstop, nu = object$control$nu, risk = object$control$risk))
  }else{
    model_after = gamboost(form2, data = data, weights = model.weights(object), family = fam, control = boost_control(mstop = mstop, nu = object$control$nu, risk = object$control$risk))
  }
  
  out = model_after
  out$tau  = tau
  out$deselectmethod = method[1] 
  class(out) = c(class(out))
  
  return(out)
}

gbm = function (Y, X, W, M, start = c("ols", "boost", "des"), trace = FALSE) {
  N = length(Y)
  
  # (1) Estimate beta consistently by ols, boosting or boosting + stability selection
  if(start == "ols") {
    mod_pre = lm(Y ~ ., data = X[,-1])
    beta = coef(mod_pre)
  } else if(start == "boost") {
    mod_pre = glmboost(Y ~ ., data = X, control = boost_control(trace = trace, mstop = M, nu = 0.1))
    cvr = cvrisk(mod_pre, folds = cv(model.weights(mod_pre), type = "subsampling"))
    beta = coef(mod_pre[mstop(cvr)], off2int = TRUE)
  } else if(start == "des") {
    mod_pre = glmboost(Y ~ ., data = X, control = boost_control(trace = trace, mstop = M, nu = 0.1))
    cvr = cvrisk(mod_pre, folds = cv(model.weights(mod_pre), type = "subsampling"))
    beta = coef(mod_pre[mstop(cvr)], off2int = TRUE)
    
    mod_pre = DeselectBoost(mod_pre, fam = Gaussian())  
    beta = coef(mod_pre, off2int = TRUE)
  }
  
  # (2) GMM estimation of lambda and sigma based on residuals of first estimator 
  u_t = Y - as.matrix(X[, names(beta)]) %*% beta
  u_b = W %*% u_t
  u_bb = W %*% W %*% u_t
  
  G = rbind(
    cbind(2/N * t(u_t) %*% u_b, -1/N * t(u_b) %*% u_b, 1),
    cbind(2/N * t(u_bb) %*% u_b, -1/N * t(u_bb) %*% u_bb, 1/N * sum(diag(t(W) %*% W))),
    cbind(1/N * (t(u_t) %*% u_bb + t(u_b) %*% u_b), -1/N * t(u_b) %*% u_bb, 0)
  )
  
  g = rbind(1/N * t(u_t) %*% u_t,
            1/N * t(u_b) %*% u_b,
            1/N * t(u_t) %*% u_b)
  
  obj = function(params) {
    lambda = params[1]
    sigma = params[2]
    theta = c(lambda, lambda^2, sigma^2)
    res = G %*% theta - g 
    t(res) %*% res
  }
  
  opt = nlminb(start = c(0, 0.5), objective = obj)
  lambda = opt$par[1]
  sig2 = opt$par[2]
  
  # (3) Compute R and omega at optimum and utilize boosting to estimate final beta
  R = diag(nrow(W)) - lambda * W
  omega = sig2 * t(R) %*% R
  
  mod = glmboost(Y ~ ., data = X, family = SEM(omega), control = boost_control(trace = trace, mstop = M, nu = 0.1))
  cvr = cvrisk(mod, folds = cv(model.weights(mod), type = "subsampling"))
  beta = coef(mod[mstop(cvr)], off2int = TRUE)
  
  res_boost = Y - as.matrix(X[, names(beta)]) %*% beta
  sig2_boost = as.numeric(t(R %*% res_boost) %*% (R %*% res_boost)) / N
  
  res = c(lambda = lambda,
          beta,
          sigma = sqrt(sig2_boost))
  
  return(list(start = start,
              first = mod_pre,
              model = mod,
              coef = res, 
              omega = omega,
              residuals = as.vector(res_boost),
              fitted = as.vector(as.matrix(X[, names(beta)]) %*% beta)))
  
}




