# Generate spatial weight matrices in circular world
network = function(n, k) {
  # n: number of units
  # k: number of neighbors ahead (and behind)
  
  W = matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:k) {
      ahead = (i + j - 1) %% n + 1
      behind = (i - j - 1) %% n + 1
      
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

nll = function(Y, Z, W, lambda, beta, sigma) {
  n = length(Y)
  
  # Compute the determinants of (I - lambda * M)
  R = diag(nrow(W)) - lambda * W
  det_R = determinant(R, logarithm = TRUE)
  log_det_R = as.numeric(det_R$modulus)
  
  # Compute the log-likelihood
  ll = -n / 2 * (log(2 * pi*sigma) + 1) +
    log_det_R  - (1 / (2*sigma)) * t(R %*% (Y - Z %*% beta)) %*% R %*% (Y - Z %*% beta)
  
  return(as.numeric(-ll))  
}

# Spatial error model family for boosting
SEM = function(omega, stabilization = c("none", "MAD", "L2")) {
  
  Family(
    ngradient = function(y, f, w = 1) {
      ngr = omega %*% (y - f)
      ngr = stabilize_ngradient(ngr, w, stabilization)
      return(ngr)
    },
    loss = function(y, f) {
      omega %*% (y - f)^2 
    },
    offset = weighted.mean,
    check_y = function(y) {
      if (!is.numeric(y) || !is.null(dim(y)))
        stop("response is not a numeric vector but ", sQuote("family = SEM"))
      y
    },
    name = "Spatial (Durbin) Error Model",
    fW = function(f) rep(1, length(f)),
    response = function(f) f
  )
}

# Stabilizin ngradients in boosting
stabilize_ngradient = function(ngr, w = 1, stabilization) {
  
  if (stabilization == "none") {
    ngr = ngr
  }
  
  if (stabilization == "MAD") {
    div = weighted.median(abs(ngr - weighted.median(ngr, w = w, na.rm = TRUE)),
                          w = w, na.rm = TRUE)
    div = ifelse(div < 0.0001, 0.0001, div)
    ngr = ngr / div
  }
  
  if (stabilization == "L2") {
    div = sqrt(weighted.mean(ngr^2, w = w,  na.rm = TRUE))
    div = ifelse(div < 1e-04, 1e-04, div)
    div = ifelse(div > 1e+04, 1e+04, div)
    ngr = ngr / div
  }
  
  ngr
}

# Deselect algorithm
DeselectBoost = function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  
  tau = ifelse(is.null(tau), 0.01, tau)
  
  if(is.null(data) && class(object$model.frame()) == 'list'){return(stop("Please enter the data."))
  } else if(!is.null(data)){
    data = data
  } else{
    data = object$model.frame()
  }
  
  nameVar = names(coef(object,which = ''))[-1]
  which.response = which(sapply(1:(dim(data)[2]), function(x){identical(as.numeric(data[,x]), as.numeric(object$response))}))
  name.response = colnames(data)[which.response]
  
  mstop = object$mstop()
  RiskRed = object$risk()
  totalRiskRed = RiskRed[1] - RiskRed[mstop+1] 
  diffRiskRed = sapply(seq(1:(mstop)), function(k){RiskRed[k]-RiskRed[k+1]})
  
  select = selected(object) - 1
  diffRiskRed = diffRiskRed[selected(object)-1 != 0]

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
  
  model_after = glmboost(form2, data = data, weights = model.weights(object), family = fam, control = boost_control(mstop = mstop, nu = object$control$nu, risk = object$control$risk))
 
  out = model_after
  out$tau  = tau
  out$deselectmethod = method[1] 
  class(out) = c(class(out))
  
  return(out)

}  

# Spatial cross-validation
cv_spat = function(map, weights, type = "k-means spatial clustering", B = 10) {
  
  # Calculate centroid coordinates for all geometries
  coords = suppressWarnings(st_coordinates(st_centroid(map)))
  
  # Add centroid coordinates safely
  map = map %>%
    mutate(x = coords[, 1],
           y = coords[, 2]) %>%
    st_set_geometry(NULL) 
  
  # Run spatial partitioning
  resamp = partition_kmeans(
    data = map,
    coords = c("x", "y"),
    nfold = B,
    repetition = 1
  )
  
  # Initialize folds matrix
  folds = matrix(weights, nrow = length(weights), ncol = B)
  
  # Apply test indices
  for (fold in seq_len(B)) {
    test = resamp[[1]][[fold]]$test
    folds[test, fold] = 0
  }
  
  attr(folds, "type") = paste(B, "-fold ", type, sep = "")
  
  return(folds)
}

semboost = function(Y, Z, W, M,
                    start = c("ols", "boost", "des"),
                    type = c("subsampling", "kfold", "bootstrap", "k-means spatial clustering"),
                    stabilization = c("none", "MAD", "L2"),
                    trace = FALSE,
                    map = NULL,
                    B = 10) {
  
  # Match arguments
  start = match.arg(start)
  type = match.arg(type)
  n = length(Y)
  
  # Utility function to get CV folds
  get_folds = function(mod, type, weights, map, B) {
    if (type == "k-means spatial clustering") {
      if (is.null(map)) stop("`map` must be provided for spatial cross-validation.")
      return(cv_spat(map = map, weights = weights, type = type, B = B))
    } else {
      return(cv(weights, type = type))
    }
  }
  
  ### --- Step 1: Initial Estimation of beta ---
  if (start == "ols") {
    mod0 = lm(Y ~ ., data = Z[,-1])
    beta = coef(mod0)
    
  } else {
    mod0 = glmboost(Y ~ ., data = Z, control = boost_control(trace = trace, mstop = M, nu = 0.1))
    folds = get_folds(mod0, type, model.weights(mod0), map, B)
    cvr = cvrisk(mod0, folds = folds)
    mod0 = mod0[mstop(cvr)]
    beta = coef(mod0, off2int = TRUE)
    
    if (start == "des") {
      mod0 = DeselectBoost(mod0, fam = Gaussian())
      beta = coef(mod0, off2int = TRUE)
    }
  }
  
  ### --- Step 2: GMM Estimation of lambda and sigma ---
  u_t = Y - as.matrix(Z[, names(beta)]) %*% beta
  u_b = W %*% u_t
  u_bb = W %*% W %*% u_t
  
  G = rbind(
    cbind(2/n * t(u_t) %*% u_b, -1/n * t(u_b) %*% u_b, 1),
    cbind(2/n * t(u_bb) %*% u_b, -1/n * t(u_bb) %*% u_bb, 1/n * sum(diag(t(W) %*% W))),
    cbind(1/n * (t(u_t) %*% u_bb + t(u_b) %*% u_b), -1/n * t(u_b) %*% u_bb, 0)
  )
  
  g = rbind(
    1/n * t(u_t) %*% u_t,
    1/n * t(u_b) %*% u_b,
    1/n * t(u_t) %*% u_b
  )
  
  obj = function(params) {
    lambda = params[1]
    sigma = params[2]
    theta = c(lambda, lambda^2, sigma^2)
    res = G %*% theta - g
    t(res) %*% res
  }
  
  scorr = c(crossprod(W %*% u_t, u_t) / crossprod(u_t, u_t))
  scorr = scorr / (sum(W) / length(u_t))
  
  pars = c(scorr, var(scale(u_t)))
  
  opt = nlminb(start = pars, objective = obj)
  lambda = opt$par[1]
  sig2 = opt$par[2]
  
  ## --- Step 3: Final SEM-Boost Estimation ---
  R = diag(nrow(W)) - lambda * W
  omega = 1 / sig2^2 * t(R) %*% R
  
  mod = glmboost(Y ~ ., data = Z, family = SEM(omega = omega, stabilization = stabilization),
                  control = boost_control(trace = trace, mstop = M, nu = 0.1))
  
  folds = get_folds(mod, type, model.weights(mod), map, B)
  cvr = cvrisk(mod, folds = folds)
  mod = mod[mstop(cvr)]
  beta = coef(mod, off2int = TRUE)
  
  res_boost = Y - as.matrix(Z[, names(beta)]) %*% beta
  sig2_boost = as.numeric(t(R %*% res_boost) %*% (R %*% res_boost)) / n
  
  res = c(lambda = lambda,
           beta,
           sigma = sqrt(sig2_boost))
  
  return(list(
    start = start,
    first = mod0,
    model = mod,
    coef = res,
    omega = omega,
    residuals = as.vector(res_boost),
    fitted = as.vector(as.matrix(Z[, names(beta)]) %*% beta)
  ))
}





