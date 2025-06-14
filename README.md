# Gradient Boosting for Spatial Regression Models with Autoregressive Disturbances

This repository provides access to three-step feasible model-based gradient boosting for spatial regression models with autoregressive disturbances.
It includes:  

- Helper functions and the spatial error family for model-based gradient boosting   
- Simulation studies under a varying spatial autoregressive parameter and varying spatial weight matrices
- Estimates for application settings including modeling the life expectancy in German districts.

The repository serves as a foundation for replication.  

## Technical Details  

For in-depth derivations and explanations of model-based gradient boosting for spatial regression models with autoregressive disturbances, refer to:  

## Example 
```
require(Matrix)
require(mboost)

source("R/Helper.R")
source("R/SEM.R")

set.seed(2222)

# Simulate artificial data
N = 400
rho_t = 0
beta_t = c(1, 3.5, -2.5, rep(0,8))
names(beta_t) = c("(Intercept)", paste0("X", 1:(length(beta_t)-1)))
gamma_t = c(-4, 3, rep(0,8))
names(gamma_t) = paste0("WX", 1:length(gamma_t))
sigma_t = 1

p = length(beta_t) + length(gamma_t) - 1
p_true = sum(beta_t[-1] != 0) + sum(gamma_t != 0)

# Generate adjacency matrices
W = network(N, k = 5)

# Generate variables and error
X = matrix(runif(N * (p / 2), -2, 2),  nrow = N, ncol = p / 2)
Z = cbind(X, W %*% X)
Z = cbind(rep(1,N), Z)
Z = data.frame(Z)
names(Z) = c(names(beta_t), names(gamma_t))

eps = rnorm(N, mean = 0, sd = sigma_t)

u = solve(Diagonal(N) - rho_t * W, eps)

Y = as.matrix(Z) %*% c(beta_t, gamma_t) + u

# Model-based gradient boosting
mod = gbm(Y, Z, W, M = 500, start = "ols", trace = TRUE)

coef(mod$model[200], off2int = TRUE)

plot(mod$model)
```








