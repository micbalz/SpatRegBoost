### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, sperrorest, knitr, mboost, gamboostLSS, MASS, spdep, spatialreg, future, future.apply, progressr)

source("R/SEM.R")

plan(multisession, workers = parallel::detectCores() - 1)  

handlers(global = TRUE)
handlers("cli")

### Simulation Setup
nsim = 100

n = 400
lambda_t = 0.4
beta_t = c(1, 3.5, -2.5, rep(0,8))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
gamma_t = c(-4, 3, rep(0,8))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sigma_t = 1

### Compute spatial weight matrix
krs = st_read(dsn = "application/vg5000_ebenen_0101", layer = "VG5000_KRS")
inkar = st_centroid(krs)
knn = knearneigh(inkar, k = 10)
nb = knn2nb(knn, row.names = krs$AGS)
listw = nb2listw(nb, style = "W")
W = listw2mat(listw)  

### Run the simulation
run = function(n, W, lambda_t, beta_t, gamma_t, sigma_t, tau) {
  p = length(beta_t) + length(gamma_t) - 1
  p_true = sum(beta_t[-1] != 0) + sum(gamma_t != 0)
  
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
  mod = semboost(Y.train, Z.train, W, M = 500, start = "des", type = "k-means spatial clustering", stabilization = "none", map = krs)
  dsgb = mod$coef
  
  ### Performance of variable selection
  nameVar = names(Z.train[-1])[1:p]
  trueVar = nameVar[c(1:(p_true/2), length(beta_t):(length(beta_t) + p_true/2 - 1))]
  falseVar = nameVar[!nameVar %in% trueVar]
  
  ### Add deselection study
  des = DeselectBoost(mod$model, tau = tau, fam = SEM(omega = mod$omega, stabilization = "none"))
  stabs = c(mod$coef["lambda"], coef(des, off2int = TRUE),  mod$coef["sigma"])
  

  selectedVar = names(stabs)[!names(stabs) %in% c("lambda", "(Intercept)", "sigma")]
  true.positive = length(which(trueVar %in% selectedVar))
  false.positive = length(which(falseVar %in% selectedVar))
  
  metrics_stabs = c(
    TPR = true.positive / length(trueVar),
    TNR = 1 - false.positive / length(falseVar),
    FDR = false.positive / length(selectedVar)
  )
  
  ### Performance of prediction
  models = list(stabs = stabs)
  model_names = names(models)
  
  lambdas = lapply(models, function(m) m[["lambda"]])
  sigmas2 = lapply(models, function(m) m[["sigma"]]^2)
  
  deltas = lapply(models, function(m) m[setdiff(names(m), c("lambda", "sigma"))])
  
  Z.tests = lapply(deltas, function(delta) Z.test[, names(delta), drop = FALSE])
  
  Y.preds = Map(function(Z, delta) as.matrix(Z) %*% delta, Z.tests, deltas)
  
  model_display_names = c(stabs = "stabs")
  
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
  
  pred = do.call(rbind, lapply(pred, as.data.frame))
  rownames(pred) = NULL
  
  
  list("stabs" = stabs,
       regularization_stabs = metrics_stabs,
       prediction = pred)
  
}


sims = function(tau_values, nsim) {
  set.seed(12345678)
  
  result = list()
  
  for (tau in tau_values) {
    
    cat("Running simulation study for tau =", tau, "\n")
    
    pb = progressor(along = 1:nsim)
    
    results_for_lambda = future_lapply(1:nsim, function(i) {
      
      res = run(n, W, lambda_t, beta_t, gamma_t, sigma_t, tau)
      pb(sprintf("tau=%.2f, replication %d", tau, i))
      res
    }, future.seed = TRUE)
    
    result[[paste0("tau=", tau)]] = results_for_lambda
  }
  return(result)
}

tau_values = c(0.005, 0.0075, 0.01,0.025,0.05,0.075,0.1,0.125)
results = sims(tau_values, nsim)


# ---- Extract into a tidy data.frame ----
rows = list()

tau_names = names(results)
for (tau_name in tau_names) {
  tau_num = as.numeric(str_remove(tau_name, "^tau="))
  runs = results[[tau_name]]
  for (i in seq_along(runs)) {
    run = runs[[i]]
    # safe extraction helpers
    safe_get = function(x, nm) {
      if (is.null(x)) return(NA_real_)
      v = tryCatch(x[nm], error = function(e) NA)
      if (length(v) == 0) return(NA_real_)
      as.numeric(v)
    }
    tpr = safe_get(run$regularization_stabs, "TPR")
    tnr = safe_get(run$regularization_stabs, "TNR")
    rmse = NA_real_
    if (!is.null(run$prediction)) {
      # prediction is a 1-row data.frame in your examples
      pred = run$prediction
      if (is.data.frame(pred) && "RMSE" %in% colnames(pred)) {
        rmse = as.numeric(pred$RMSE[1])
      } else if (!is.null(pred["RMSE"])) {
        rmse = as.numeric(pred["RMSE"])
      }
    }
    rows[[length(rows) + 1]] = data.frame(
      tau = tau_num,
      run = i,
      TPR = tpr,
      TNR = tnr,
      RMSE = rmse,
      stringsAsFactors = FALSE
    )
  }
}

df = bind_rows(rows)

# ---- Tidy for ggplot: pivot longer so stat becomes a facet ----
df_long = df %>%
  pivot_longer(cols = c(TPR, TNR, RMSE), names_to = "stat", values_to = "value") %>%
  mutate(
    tau = factor(tau, levels = sort(unique(tau)))   # preserve ordering on x-axis
  )

# ---- Plot: boxplots per tau, facet by stat ----
ggplot(df_long, aes(x = tau, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ stat, scales = "free", nrow = 1,
             labeller = as_labeller(c(TPR = "TPR", TNR = "TNR", RMSE = "RMSE"))) +
  labs(
    x = expression(tau),   
    y = NULL
  ) +
  theme_bw(base_size = 13)


