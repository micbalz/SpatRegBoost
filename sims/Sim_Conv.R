### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, Matrix, mboost, MASS, spdep, spatialreg)

source("R/Helper.R")
source("R/SEM.R")

gbm = function (Y, X, W, M, start = c("ols", "boost", "des")) {
  N = length(Y)
  
  # (1) Estimate beta consistently by ols, boosting or boosting + stability selection
  if(start == "ols") {
    mod = lm(Y ~ ., data = X[,-1])
    beta = coef(mod)
  } else if(start == "boost") {
    mod = glmboost(Y ~ ., data = X, control = boost_control(trace = FALSE, mstop = M, nu = 0.1))
    cvr = cvrisk(mod, folds = cv(model.weights(mod), type = "subsampling"))
    beta = coef(mod[mstop(cvr)], off2int = TRUE)
  } else if(start == "des") {
    mod = glmboost(Y ~ ., data = X, control = boost_control(trace = FALSE, mstop = M, nu = 0.1))
    cvr = cvrisk(mod, folds = cv(model.weights(mod), type = "subsampling"))
    beta = coef(mod[mstop(cvr)], off2int = TRUE)
    
    mod = DeselectBoost(mod, fam = Gaussian())  
    beta = coef(mod, off2int = TRUE)
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
  
  mod = glmboost(Y ~ ., data = X, family = SEM(omega), control = boost_control(trace = FALSE, mstop = M, nu = 0.1))
  beta = coef(mod, off2int = TRUE)
  
  res_boost = Y - as.matrix(X[, names(beta)]) %*% beta
  res_boost = res_boost - lambda * (W %*% res_boost)
  sig2_boost = sqrt(c(crossprod(res_boost)) / N)
  
  res = c(lambda = lambda,
          beta,
          sigma = sig2_boost)
  
  return(res)
  
}

set.seed(2222)

### Simulation Setup
N = 400
rho_t = 0
beta_t = c(1, 3.5, -2.5, rep(0,8))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
gamma_t = c(-4, 3, rep(0,8))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sigma_t = 1

p = length(beta_t) + length(gamma_t) - 1
p_true = sum(beta_t[-1] != 0) + sum(gamma_t != 0)

### Generate adjacency matrices
W = network(N, k = 5)

### Generate covariates and error
X = matrix(runif(N * (p / 2), -2, 2),  nrow = N, ncol = p / 2)
Z = cbind(X, W %*% X)
Z = cbind(rep(1,N), Z)
Z = data.frame(Z)
names(Z) = c(names(beta_t), names(gamma_t))

eps = rnorm(N, mean = 0, sd = sigma_t)

u = solve(Diagonal(N) - rho_t * W, eps)

Y = as.matrix(Z) %*% c(beta_t, gamma_t) + u

### Estimate the model
mod = GMerrorsar(Y ~ ., data = Z, listw = mat2listw(W, style = "W"), zero.policy = FALSE, legacy = TRUE)
gmm = c(coef(mod)[length(coef(mod))], coef(mod)[-length(coef(mod))], sqrt(mod$s2))
names(gmm) = c("lambda", names(Z), "sigma")

gbsem = gbm(Y, Z, W, M = 10000, start = "ols")

gmm
gbsem


