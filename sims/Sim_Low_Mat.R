### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, knitr, Matrix, mboost, MASS, spdep, spatialreg, future, future.apply, progressr)

source("R/Helper.R")
source("R/SEM.R")

plan(multisession, workers = parallel::detectCores() - 2)  

# At the beginning of your script
handlers(global = TRUE)
handlers("cli")

### Simulation Setup
nsim = 100

N = 400
lambda_t = 0
beta_t = c(1, 3.5, -2.5, rep(0,8))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
gamma_t = c(-4, 3, rep(0,8))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sigma_t = 1

# Run the simulation
run = function (N, lambda_t, beta_t, gamma_t, sigma_t, k) {
  p = length(beta_t) + length(gamma_t) - 1
  p_true = sum(beta_t[-1] != 0) + sum(gamma_t != 0)
  
  ### Generate adjacency matrices
  W = network(N, k = k)
  
  ### Generate covariates and error
  X.train = matrix(runif(N * (p / 2), -2, 2),  nrow = N, ncol = p / 2)
  Z.train = cbind(X.train, W %*% X.train)
  Z.train = cbind(rep(1,N), Z.train)
  Z.train = data.frame(Z.train)
  names(Z.train) = c(names(beta_t), names(gamma_t))
  
  X.test = matrix(runif(N * (p / 2), -2, 2),  nrow = N, ncol = p / 2)
  Z.test = cbind(X.test, W %*% X.test)
  Z.test = cbind(rep(1,N), Z.test)
  Z.test = data.frame(Z.test)
  names(Z.test) = c(names(beta_t), names(gamma_t))
  
  eps.train = rnorm(N, mean = 0, sd = sigma_t)
  eps.test = rnorm(N, mean = 0, sd = sigma_t)
  
  u.train = solve(Diagonal(N) - lambda_t * W, eps.train)
  u.test = solve(Diagonal(N) - lambda_t * W, eps.test)
  
  Y.train = as.matrix(Z.train) %*% c(beta_t, gamma_t) + u.train
  Y.test = as.matrix(Z.test) %*% c(beta_t, gamma_t) + u.test
  
  ### Run the models
  mod = errorsarlm(Y.train ~ ., data = data.frame(X.train), listw = mat2listw(W, style = "W"), Durbin = TRUE)
  mle = c(coef(mod),sqrt(mod$s2))
  names(mle) = c("lambda", names(Z.train), "sigma")
  
  mod = GMerrorsar(Y.train ~ ., data = Z.train, listw = mat2listw(W, style = "W"))
  gmm = c(coef(mod)[length(coef(mod))], coef(mod)[-length(coef(mod))], sqrt(mod$s2))
  names(gmm) = c("lambda", names(Z.train), "sigma")
  
  mod = gbm(Y.train, Z.train, W, M = 500, start = "ols")
  lsgb = mod$coef
  
  mod = gbm(Y.train, Z.train, W, M = 500, start = "boost")
  gbgb = mod$coef
  
  est = gbm(Y.train, Z.train, W, M = 500, start = "des")
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
  models = list(mle = mle, gmm = gmm, lsgb = lsgb, gbgb = gbgb, dsgb = dsgb)
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
  model_display_names = c(mle = "QML", gmm = "GMM", lsgb = "LSGB", gbgb = "GBGB", dsgb = "DSGB")
  
  # Compute metrics for each model
  pred = mapply(function(name, y_pred, delta, lambda, sigma2, Z_test) {
    list(
      Model = model_display_names[[name]],
      RMSE = rmse(Y.test, y_pred),
      MAE  = mae(Y.test, y_pred),
      NLL  = nll(N, W, Y.test, as.matrix(Z_test), lambda = lambda, beta = delta, sigma = sigma2)
    )
  }, 
  names(Y.preds), Y.preds, deltas, lambdas, sigmas2, Z.tests,
  SIMPLIFY = FALSE)
  
  # Convert list of lists to data frame
  pred = do.call(rbind, lapply(pred, as.data.frame))
  rownames(pred) = NULL
  
  ### Add deselection study
  des = DeselectBoost(est$model, fam = SEM(omega = est$omega))
  stabs = coef(des, off2int = TRUE)
  
  selectedVar = names(stabs)[!names(stabs) %in% c("lambda", "(Intercept)", "sigma")]
  true.positive = length(which(trueVar %in% selectedVar))
  false.positive = length(which(falseVar %in% selectedVar))
  
  metrics_stabs = c(
    TPR = true.positive / length(trueVar),
    TNR = 1 - false.positive / length(falseVar),
    FDR = false.positive / length(selectedVar)
  )
  
  
  list("QML" = mle,
       "GMM" = gmm,
       "LSGB" = lsgb,
       "GBGB" = gbgb,
       "DSGB" = dsgb,
       regularization = metrics,
       regularization_stabs = metrics_stabs,
       prediction = pred)
  
}

sims = function(K, nsim) {
  set.seed(222)
  
  result = list()
  
  for (k in K) {
    
    cat("Running simulation study for k =", k, "\n")
    
    pb = progressor(along = 1:nsim)
    
    results_for_k = future_lapply(1:nsim, function(i) {
      
      res = run(N, lambda_t, beta_t, gamma_t, sigma_t, k)
      pb(sprintf("k=%.2f, replication %d", k, i))
      res
    }, future.seed = TRUE)
    
    result[[paste0("k=", k)]] = results_for_k
  }
  return(result)
}

# Run for multiple k
K = c(1, 2, 3, 5, 10, 20)
results = sims(K, nsim)

### Evaluate variable selection
# Suppose your results list is named results_by_lambda
reg_means = lapply(results, function(sim_list) {
  # Extract regularization metrics for each simulation run within this k
  reg_matrix = do.call(rbind, lapply(sim_list, function(sim) sim$regularization))
  # Compute column means * 100 (if you want percentages), ignoring NA
  colMeans(reg_matrix * 100, na.rm = TRUE)
})

# Convert to data frame for nicer output
reg_means = do.call(rbind, reg_means)
reg_means = data.frame(k = rownames(reg_means), reg_means, row.names = NULL)
reg_means$k = as.numeric(sub("k=", "", reg_means$k))

# Suppose your results list is named results_by_lambda
reg_stabs = lapply(results, function(sim_list) {
  # Extract regularization metrics for each simulation run within this k
  reg_matrix = do.call(rbind, lapply(sim_list, function(sim) sim$regularization_stabs))
  # Compute column means * 100 (if you want percentages), ignoring NA
  colMeans(reg_matrix * 100, na.rm = TRUE)
})

# Convert to data frame for nicer output
reg_stabs = do.call(rbind, reg_stabs)
reg_stabs = data.frame(k = rownames(reg_stabs), reg_stabs, row.names = NULL)
reg_stabs$k = as.numeric(sub("k=", "", reg_stabs$k))


### Evaluate prediction
# Step 1: Loop over each rho-level
pred_means = lapply(names(results), function(lambda_name) {
  # Get all prediction data frames for this rho
  pred_df = do.call(rbind, lapply(results[[lambda_name]], function(res) res$prediction))
  
  # Add the model as a factor and compute mean metrics
  pred_summary = pred_df %>%
    mutate(Model = factor(Model, levels = c("QML", "GMM", "LSGB", "GBGB", "DSGB")),
           k = sub("k=", "", lambda_name)) %>%
    group_by(k, Model) %>%
    summarise(across(c(RMSE, MAE, NLL), mean), .groups = "drop")
  
  pred_summary
})
pred_means = bind_rows(pred_means)
pred_means$k = as.numeric(pred_means$k)

pred_means = pred_means %>%
  rename(RMSEP = RMSE, MAEP = MAE) %>%  
  dplyr::select(k, Model, RMSEP, MAEP, NLL)

pred_means = pred_means %>%
  pivot_longer(cols = c("RMSEP", "MAEP", "NLL"), names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  arrange(k, match(Metric, c("RMSEP", "MAEP", "NLL"))) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.4f", .)))  # Format to 4 decimals

pred = pred_means %>%
  group_by(k) %>%
  mutate(k = ifelse(row_number() == 1, paste0("$", k, "$"), "")) %>%
  ungroup()

### Evaluate estimation
perf = list()

for (k in names(results)) {
  true_vals = c(lambda = lambda_t)
  
  # Extract simulations for this rho
  res = results[[k]]
  
  # For each method, bind rows and evaluate
  coeffs = lapply(c(QML = "QML", GMM = "GMM", LSGB = "LSGB", GBGB = "GBGB", DSGB = "DSGB"), function(key) {
    as.data.frame(do.call(rbind, lapply(res, function(res) {
      params(res, key, true_vals = true_vals)
    })))
  })
  
  performance = lapply(coeffs, function(df) {
    bias_mse_empse(df, true_vals = true_vals)
  })
  
  perf[[k]] = performance
}

perf = imap_dfr(perf, function(models, k) {
  imap_dfr(models, function(metrics, k) {
    metrics %>%
      as.data.frame() %>%
      rownames_to_column("parameter") %>%
      mutate(model = k)
  }, .id = "model") %>%
    mutate(k = k)
}, .id = "lambda_group") %>%
  dplyr::select(k, model, parameter, Bias, MSE, SE)

perf$k = as.numeric(sub("k=", "", perf$k))

perf = perf %>%
  dplyr::select(k, model, Bias, MSE, SE) %>%
  pivot_longer(cols = c(Bias, MSE, SE), names_to = "stat", values_to = "value") %>%
  pivot_wider(names_from = model, values_from = value) %>%
  arrange(k, match(stat, c("Bias", "MSE", "SE"))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  mutate(
    across(-c(k, stat), ~ ifelse(stat == "MSE", paste0("(", sprintf("%.4f", .x), ")"),
                                      ifelse(stat == "SE", paste0("[", sprintf("%.4f", .x), "]"),
                                             sprintf("%.4f", .x))))
  )


perf = perf %>%
  group_by(k) %>%
  mutate(k = ifelse(row_number() == 1, paste0("$", k, "$"), "")) %>%
  ungroup()


### Evaluate regularization
# Define the correct influential variables in order
inf_var = c(beta_t[2:3], gamma_t[which(gamma_t != 0)], beta_t[4])

# Get DSGB coefficients for these variables for each rho
regs_list = lapply(names(results), function(k) {
  res_list = results[[k]]
  
  # Extract coefficients for DSGB
  reg = as.data.frame(do.call(rbind, lapply(res_list, function(res) {
    p = params(res, "DSGB", true_vals = inf_var)
    # Ensure correct column order
    p = p[names(inf_var)]
    return(p)
  })))
  
  # Replace NA with 0
  reg[is.na(reg)] = 0
  colnames(reg)[5] = "non-inf"
  reg
})

names(regs_list) = names(results)

# Step 1: Long-format version of your regression data
reg_combined = bind_rows(
  lapply(names(regs_list), function(name) {
    df = regs_list[[name]]
    df$k = name
    df
  }),
  .id = NULL
)

# Convert to expression for LaTeX-style labels
reg_combined$k = gsub("=", "==", reg_combined$k)

reg_long = pivot_longer(
  reg_combined,
  cols = -k,
  names_to = "Variable",
  values_to = "Coefficient"
)

# Step 2: Make sure 'non-inf' is last
reg_long$Variable = factor(reg_long$Variable, levels = c(setdiff(unique(reg_long$Variable), "non-inf"), "non-inf"))

# Step 1: Extract numeric k value
reg_long$lambda_num = as.numeric(sub("k==", "", as.character(reg_long$k)))

# Step 2: Define desired order
K_labels = paste0("k==", K)

# Step 3: Re-factor with correct order
reg_long$k = factor(
  reg_long$k,
  levels = K_labels,
  labels = K_labels
)

# Step 3: Convert inf_var to data frame in same long format
names(inf_var)[5] = "non-inf"
true.df = data.frame(
  Variable = names(inf_var),
  value = as.numeric(inf_var)
)

# Final plot with LaTeX facet titles
ggplot(reg_long, aes(x = Variable, y = Coefficient)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_boxplot(data = true.df, aes(x = Variable, y = value), color = "red") +
  ylab("Coefficients") +
  xlab("") +
  facet_wrap(~ k, labeller = label_parsed, ncol = 3, nrow = 4) +
  theme_bw(base_size = 15)


### Print the results
kable(reg_means,
      align = "lccc", 
      caption = "Average selection rates in the low-dimensional linear setting"
)

kable(reg_stabs, 
      align = "lccc",
      caption = "Average selection rates for deselection in the low-dimensional linear setting"
)

kable(
  perf,
  align = "lccccc",
  caption = "Estimation performance for the spatial autoregressive parameter lambda in the low-dimensional linear setting"
)


kable(
  pred,
  align = "llccccc",
  caption = "Prediction performance on independent test data for the low-dimensional linear setting"
)
