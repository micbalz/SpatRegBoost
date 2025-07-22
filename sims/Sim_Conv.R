### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, sperrorest, knitr, mboost, gamboostLSS, MASS, spdep, spatialreg)

source("R/SEM.R")

semboost = function (Y, Z, W, M) {
  n = length(Y)
  
  # (1) Estimate beta consistently by ols
  mod0 = lm(Y ~ ., data = Z[,-1])

  # (2) GMM estimation of lambda and sigma based on residuals of first estimator 
  u_t = mod0$residuals
  u_b = W %*% u_t
  u_bb = W %*% W %*% u_t
  
  G = rbind(
    cbind(2/n * t(u_t) %*% u_b, -1/n * t(u_b) %*% u_b, 1),
    cbind(2/n * t(u_bb) %*% u_b, -1/n * t(u_bb) %*% u_bb, 1/n * sum(diag(t(W) %*% W))),
    cbind(1/n * (t(u_t) %*% u_bb + t(u_b) %*% u_b), -1/n * t(u_b) %*% u_bb, 0)
  )
  
  g = rbind(1/n * t(u_t) %*% u_t,
            1/n * t(u_b) %*% u_b,
            1/n * t(u_t) %*% u_b)
  
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
  omega = 1/sig2^2 * t(R) %*% R
  
  mod = glmboost(Y ~ ., data = Z, family = SEM(omega, stabilization = "MAD"), control = boost_control(trace = TRUE, mstop = M, nu = 0.1))
  beta = coef(mod, off2int = TRUE)
  
  res_boost = Y - as.matrix(Z[, names(beta)]) %*% beta
  res_boost = res_boost - lambda * (W %*% res_boost)
  sig2_boost = sqrt(c(crossprod(res_boost)) / n)
  
  res = c(lambda = lambda,
          beta,
          sigma = sig2_boost)
  
  return(res)
  
}

set.seed(98765432)
### Simulation Setup
n = 400
lambda_t = 0
beta_t = c(1, 3.5, -2.5, rep(0,8))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
gamma_t = c(-4, 3, rep(0,8))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sigma_t = 1

p = length(beta_t) + length(gamma_t) - 1
p_true = sum(beta_t[-1] != 0) + sum(gamma_t != 0)

### Compute spatial weight matrix
krs = st_read(dsn = "application/vg5000_ebenen_0101", layer = "VG5000_KRS")
inkar = st_centroid(st_geometry(krs))
knn = knn2nb(knearneigh(inkar, k = 10), row.names = krs$AGS)
listw = nb2listw(knn, style = "W")
W = listw2mat(listw)  


### Generate covariates and error
X = matrix(runif(n * (p / 2), -2, 2),  nrow = n, ncol = p / 2)
Z = cbind(X, W %*% X)
Z = cbind(rep(1,n), Z)
Z = data.frame(Z)
names(Z) = c(names(beta_t), names(gamma_t))

eps = rnorm(n, mean = 0, sd = sigma_t)

u = solve(diag(n) - lambda_t * W, eps)

Y = as.matrix(Z) %*% c(beta_t, gamma_t) + u


### Estimate the model
mod = GMerrorsar(Y ~ ., data = Z, listw = listw, zero.policy = FALSE, legacy = TRUE)
gmm = c(coef(mod)[length(coef(mod))], coef(mod)[-length(coef(mod))], sqrt(mod$s2))
names(gmm) = c("lambda", names(Z), "sigma")

lsgb = semboost(Y, Z, W, M = 10000)
gmm
lsgb

res = data.frame(
  Variable = names(gmm),
  GMM = as.numeric(gmm),
  LSGB = as.numeric(lsgb)
)

kable(res, digits = 10, caption = "Comparison of GMM and LSGB Coefficients")
