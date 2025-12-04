# Gradient Boosting for Spatial Regression Models with Autoregressive Disturbances

This repository provides access to three-step feasible model-based gradient boosting for spatial regression models with autoregressive disturbances.
It includes:  

- Helper functions and the spatial error family for model-based gradient boosting   
- Simulation studies under a varying spatial autoregressive parameter and varying spatial weight matrices
- Estimates for application settings including modeling the life expectancy in German districts.

The repository serves as a foundation for replication.  

## Technical Details  

For in-depth derivations and explanations of model-based gradient boosting for spatial regression models with autoregressive disturbances, refer to: 

Balzer, M. (2025).
Gradient boosting for spatial regression models with autoregressive disturbances. Networks and Spatial
Economics. https://doi.org/10.1007/s11067-025-09717-8 

## Example 
```
require(Matrix)
require(mboost)

source("R/SEM.R")

set.seed(12345678)

# Simulate artificial data
n = 400
lambda_t = 0
beta_t = c(1, 3.5, -2.5, rep(0,8))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
gamma_t = c(-4, 3, rep(0,8))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sigma_t = 1

p = length(beta_t) + length(gamma_t) - 1
p_true = sum(beta_t[-1] != 0) + sum(gamma_t != 0)

# Generate spatial weight matrix
W = network(n, k = 5)

# Generate variables and error
X = matrix(runif(n * (p / 2), -2, 2),  nrow = n, ncol = p / 2)
Z = cbind(X, W %*% X)
Z = cbind(rep(1,n), Z)
Z = data.frame(Z)
names(Z) = c(names(beta_t), names(gamma_t))

eps = rnorm(n, mean = 0, sd = sigma_t)

u = solve(diag(n) - lambda_t * W, eps)

Y = as.matrix(Z) %*% c(beta_t, gamma_t) + u

# Model-based gradient boosting
mod = semboost(Y, Z, W, M = 500, start = "ols", type = "kfold", stabilization = "none", trace = TRUE)

coef(mod$model[200], off2int = TRUE)

par(mar = c(5, 4, 4, 6))  
plot(mod$model)
```








