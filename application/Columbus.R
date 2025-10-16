### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, sperrorest, knitr, mapview, tmap, viridisLite, Matrix, mboost, gamboostLSS, MASS, spData, sf, spdep, sphet, spatialreg, sperrorest)

source("R/SEM.R")

# Load the Columbus shapefile (sf object)
map = st_read(system.file("shapes/columbus.gpkg", package="spData")[1])

# Create spatial weight matrix
nb = poly2nb(map, queen = TRUE)
nbw = spdep::nb2listw(nb, style = "W")

# Estimation
Y = columbus[, "CRIME"]
X = cbind(data.frame("(Intercept)" = rep(1,length(Y))), columbus[, c("INC", "HOVAL", "DISCBD", "PLUMB", "OPEN")])
colnames(X)[1] = "(Intercept)"
nb = poly2nb(map, queen = TRUE)
nbw = nb2listw(nb, style = "W")
W = listw2mat(nbw)  
Z = as.data.frame(cbind(X, W %*% data.matrix(X[,-1])))
colnames(Z)[(ncol(X) + 1):ncol(Z)] = paste0("W_", colnames(Z)[(ncol(X)+1):ncol(Z)])

# Estimate the models
qml = errorsarlm(CRIME ~ ., data = data.frame(CRIME = Y, X[,-1]), listw = nbw, Durbin = TRUE)
gmm = GMerrorsar(CRIME ~ ., data = data.frame(CRIME = Y, Z), listw = nbw)
set.seed(12345678)
lsgb = semboost(Y, Z, W, M = 1000, start = "ols", type = "k-means spatial clustering", stabilization = "MAD", map = map)
gbgb = semboost(Y, Z, W, M = 1000, start = "boost", type = "k-means spatial clustering", stabilization = "MAD", map = map)
dsgb = semboost(Y, Z, W, M = 1000, start = "des", type = "k-means spatial clustering", stabilization = "MAD", map = map)
des = DeselectBoost(dsgb$model, fam = SEM(omega = dsgb$omega, stabilization = "MAD"))


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

# Custom variable ordering: lambda, non-W, W, sigma
non_w_vars = sort(grep("^W_", coef_df$Variable, invert = TRUE, value = TRUE))
non_w_vars = setdiff(non_w_vars, c("lambda", "sigma"))

w_vars = sort(grep("^W_", coef_df$Variable, value = TRUE))

ordered_vars = c("lambda", non_w_vars, w_vars, "sigma")

coef_df = coef_df %>%
  dplyr::arrange(match(Variable, ordered_vars))

# Build the map using CRIME variable
tm_shape(map) +
  tm_polygons("CRIME", 
              palette = viridis(n = 100, direction = -1, option = "G"),
              style = "quantile", 
              title = "Crime rate") 

# Display the table nicely
kable(coef_df, align = "lcccccc", caption = "Coefficient estimates across different estimation strategies for Columbus crime rate")

