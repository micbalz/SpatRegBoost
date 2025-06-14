### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, knitr, mapview, tmap, viridisLite, Matrix, mboost, MASS, spData, sf, spdep, sphet, spatialreg)

source("R/Helper.R")
source("R/SEM.R")

# Load data
data("boston", package = "spData")

# Build the map
boston.tr = sf::st_read(system.file("shapes/boston_tracts.gpkg",
                                    package="spData")[1])
boston_map = left_join(boston.tr, boston.c, by = c("TOWNNO", "TRACT"))

tm_shape(boston_map) +
  tm_polygons("CMEDV.x",
              palette = viridis(n = 100, direction = -1, option = "G"),
              style = "quantile",
              title = "Housing Prices")


# Create spatial weights matrix and design matrix 
Y =  boston.c[, "CMEDV"]
X = cbind(data.frame("(Intercept)" = rep(1,length(Y))), boston.c[, c("CRIM", "ZN", "INDUS", "NOX", "RM", "AGE", "RAD", "DIS", "TAX", "PTRATIO", "B", "LSTAT")])
colnames(X)[1] = "(Intercept)"
nbw = nb2listw(boston.soi, style = "W")  
W = listw2mat(nbw)  
Z = as.data.frame(cbind(X, W %*% data.matrix(X[,-1])))
colnames(Z)[(ncol(X) + 1):ncol(Z)] = paste0("W_", colnames(Z)[(ncol(X)+1):ncol(Z)])

# Estimate the models
qml = errorsarlm(CMEDV ~ ., data = data.frame(CMEDV = Y, X[,-1]), listw = nbw, Durbin = TRUE)
gmm = GMerrorsar(CMEDV ~ ., data = data.frame(CMEDV = Y, Z), listw = nbw)
set.seed(222)
lsgb = gbm(Y, Z, W, M = 500, start = "ols")
gbgb = gbm(Y, Z, W, M = 500, start = "boost")
dsgb = gbm(Y, Z, W, M = 500, start = "des")
des = DeselectBoost(dsgb$model, fam = SEM(omega = dsgb$omega))

# Extract coefficients from all models
coefs_list = list(
  QML   = c(coef(qml), sigma = sqrt(qml$s2)),
  GMM   = c(coef(gmm)[length(coef(gmm))], coef(gmm)[-length(coef(gmm))], sigma = sqrt(gmm$s2)),
  LSGB = lsgb$coef,
  GBGB = gbgb$coef,
  DSGB = dsgb$coef,
  DSDS = c(dsgb$coef[1], coef(des, off2int = TRUE), dsgb$coef[length(dsgb$coef)])
)

coefs_list = lapply(coefs_list, function(coefs) {
  names(coefs) = gsub("^lag\\.", "W_", names(coefs))
  coefs
})


# All variable names used across models
all_vars = unique(unlist(lapply(coefs_list, names)))
all_vars = sort(all_vars)  # For consistent row order

# Create empty table (matrix)
coef_matrix = matrix(NA, nrow = length(all_vars), ncol = length(coefs_list),
                     dimnames = list(all_vars, names(coefs_list)))

# Fill the matrix
for (model in names(coefs_list)) {
  this_coef = coefs_list[[model]]
  coef_matrix[names(this_coef), model] = round(this_coef, 4)
}

# Convert to data.frame for nicer printing
coef_df = as.data.frame(coef_matrix)
coef_df = tibble::rownames_to_column(coef_df, var = "Variable")

coef_df = coef_df %>%
  dplyr::arrange(match(Variable,
                       c("lambda",
                         setdiff(coef_df$Variable, c("lambda", "sigma")),
                         "sigma")))

# Display the table nicely
kable(coef_df, align = "lcccccc", caption = "Coefficient estimates across different estimation strategies for Boston housing prices")

