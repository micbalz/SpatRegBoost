### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, sperrorest, knitr, mboost, gamboostLSS, MASS, spdep, spatialreg, future, future.apply, progressr)

source("R/SEM.R")

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
  
  pars = c(scorr, var(u_t))
  
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
          sigma = sig2)
  
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




plan(multisession, workers = parallel::detectCores() - 1)  

# At the beginning of your script
handlers(global = TRUE)
handlers("cli")

### Simulation Setup
nsim = 100

k = 1
lambda_t = 0
beta_t = c(1, 3.5, -2.5, rep(0,98))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
gamma_t = c(-4, 3, rep(0,98))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sigma_t = 1

# Run the simulation
run = function (n, lambda_t, beta_t, gamma_t, sigma_t, k) {
  p = length(beta_t) + length(gamma_t) - 1
  p_true = sum(beta_t[-1] != 0) + sum(gamma_t != 0)
  
  ### Generate adjacency matrices
  W = network(n, k = k)
  
  ### Generate covariates and error
  X.train = matrix(runif(n * (p / 2), -2, 2),  nrow = n, ncol = p / 2)
  Z.train = cbind(X.train, W %*% X.train)
  Z.train = cbind(rep(1,n), Z.train)
  Z.train = data.frame(Z.train)
  names(Z.train) = c(names(beta_t), names(gamma_t))
  
  X.test = matrix(runif(n * (p / 2), -2, 2),  nrow = n, ncol = p / 2)
  Z.test = cbind(X.test, W %*% X.test)
  Z.test = cbind(rep(1,n), Z.test)
  Z.test = data.frame(Z.test)
  names(Z.test) = c(names(beta_t), names(gamma_t))
  
  eps.train = rnorm(n, mean = 0, sd = sigma_t)
  eps.test = rnorm(n, mean = 0, sd = sigma_t)
  
  u.train = solve(diag(n) - lambda_t * W, eps.train)
  u.test = solve(diag(n) - lambda_t * W, eps.test)
  
  Y.train = as.matrix(Z.train) %*% c(beta_t, gamma_t) + u.train
  Y.test = as.matrix(Z.test) %*% c(beta_t, gamma_t) + u.test
  
  ### Run the models
  mod = semboost(Y.train, Z.train, W, M = 250, start = "boost", type = "kfold", stabilization = "MAD")
  gbgb = mod$coef
  
  est = semboost(Y.train, Z.train, W, M = 250, start = "des", type = "kfold", stabilization = "MAD")
  dsgb = est$coef
  
  ### Performance of variable selection
  nameVar = names(Z.train[-1])[1:p]
  trueVar = nameVar[c(1:(p_true/2), length(beta_t):(length(beta_t) + p_true/2 - 1))]
  falseVar = nameVar[!nameVar %in% trueVar]
  
  selectedVar = names(dsgb)[!names(dsgb) %in% c("lambda", "(Intercept)", "sigma")]
  true.positive = length(which(trueVar %in% selectedVar))
  false.positive = length(which(falseVar %in% selectedVar))
  
  metrics = c(
    TPR = true.positive / length(trueVar),
    TNR = 1 - false.positive / length(falseVar),
    FDR = false.positive / length(selectedVar)
  )
  
  ### Performance of prediction
  models = list(gbgb = gbgb, dsgb = dsgb)
  model_names = names(models)
  
  # Extract components
  lambdas = lapply(models, function(m) m[["lambda"]])
  sigmas2 = lapply(models, function(m) m[["sigma"]]^2)
  
  # Extract delta (excluding "lambda" and "sigma")
  deltas = lapply(models, function(m) m[setdiff(names(m), c("lambda", "sigma"))])
  
  # Align Z.test columns to match deltas
  Z.tests = lapply(deltas, function(delta) Z.test[, names(delta), drop = FALSE])
  
  # Compute predicted values
  Y.preds = Map(function(Z, delta) as.matrix(Z) %*% delta, Z.tests, deltas)
  
  # Define display names (your desired labels for the models)
  model_display_names = c(gbgb = "GBGB", dsgb = "DSGB")
  
  # Compute metrics for each model
  pred = mapply(function(name, y_pred, delta, lambda, sigma2, Z_test) {
    list(
      Model = model_display_names[[name]],
      RMSE = rmse(Y.test, y_pred),
      MAE  = mae(Y.test, y_pred),
      NLL  = nll(Y.test, as.matrix(Z_test), W, lambda = lambda, beta = delta, sigma = sigma2)
    )
  }, 
  names(Y.preds), Y.preds, deltas, lambdas, sigmas2, Z.tests,
  SIMPLIFY = FALSE)
  
  # Convert list of lists to data frame
  pred = do.call(rbind, lapply(pred, as.data.frame))
  rownames(pred) = NULL
  
  ### Add deselection study
  des = DeselectBoost(est$model, fam = SEM(omega = est$omega, stabilization = "none"))
  stabs = coef(des, off2int = TRUE)
  
  selectedVar = names(stabs)[!names(stabs) %in% c("lambda", "(Intercept)", "sigma")]
  true.positive = length(which(trueVar %in% selectedVar))
  false.positive = length(which(falseVar %in% selectedVar))
  
  metrics_stabs = c(
    TPR = true.positive / length(trueVar),
    TNR = 1 - false.positive / length(falseVar),
    FDR = false.positive / length(selectedVar)
  )
  
  
  list("GBGB" = gbgb,
       "DSGB" = dsgb,
       regularization = metrics,
       regularization_stabs = metrics_stabs,
       prediction = pred)
  
}

sims = function(N, nsim) {
  set.seed(123456789)
  
  result = list()
  
  for (n in N) {
    
    cat("Running simulation study for n =", n, "\n")
    
    pb = progressor(along = 1:nsim)
    
    results_for_n = future_lapply(1:nsim, function(i) {
      
      res = run(n, lambda_t, beta_t, gamma_t, sigma_t, k)
      pb(sprintf("n=%.2f, replication %d", n, i))
      res
    }, future.seed = TRUE)
    
    result[[paste0("n=", n)]] = results_for_n
  }
  return(result)
}

# Run for multiple n
N = c(50, 100, 150)
results = sims(N, nsim)

### Evaluate estimation
perf = list()

for (k in names(results)) {
  true_vals = c(lambda = lambda_t, sigma = sigma_t)
  
  # Extract simulations for this n
  res = results[[k]]
  
  # For each method, bind rows and evaluate
  coeffs = lapply(c(GBGB = "GBGB", DSGB = "DSGB"), function(key) {
    as.data.frame(do.call(rbind, lapply(res, function(res) {
      params(res, key, true_vals = true_vals)
    })))
  })
  
  performance = lapply(coeffs, function(df) {
    bias_mse_empse(df, true_vals = true_vals)
  })
  
  perf[[k]] = performance
}

perf = imap_dfr(perf, function(models, n) {
  imap_dfr(models, function(metrics, n) {
    metrics %>%
      as.data.frame() %>%
      rownames_to_column("parameter") %>%
      mutate(model = n)
  }, .id = "model") %>%
    mutate(n = n)
}, .id = "lambda_group") %>%
  dplyr::select(n, model, parameter, Bias, MSE, SE)

perf$n = as.numeric(sub("n=", "", perf$n))

perf = perf %>%
  # keep parameter so pivot_wider has unique id columns
  dplyr::select(n, model, parameter, Bias, MSE, SE) %>%
  pivot_longer(cols = c(Bias, MSE, SE), names_to = "stat", values_to = "value") %>%
  pivot_wider(
    id_cols = c(n, stat, parameter),        # ensure uniqueness
    names_from = model,
    values_from = value
  ) %>%
  arrange(n, match(stat, c("Bias", "MSE", "SE"))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  mutate(
    across(-c(n, stat, parameter), 
           ~ ifelse(stat == "MSE",
                    paste0("(", sprintf("%.4f", .x), ")"),
                    ifelse(stat == "SE",
                           paste0("[", sprintf("%.4f", .x), "]"),
                           sprintf("%.4f", .x))))
  )

# keep the grouping/labeling you had
perf = perf %>%
  group_by(n) %>%
  mutate(n = ifelse(row_number() == 1, paste0("$", n, "$"), "")) %>%
  ungroup()

kable(
  perf,
  align = "lccccc",
  caption = "Estimation performance for the spatial autoregressive parameter lambda in the high-dimensional linear setting"
)
