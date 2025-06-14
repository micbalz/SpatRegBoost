### Clear workspace and load required libraries
rm(list = ls())
gc()
options(scipen = 900)

pacman::p_load(tidyverse, data.table, mapview, tmap, viridisLite, Matrix, mboost, MASS, spData, sf, spdep, sphet, spatialreg)

source("R/Helper.R")
source("R/SEM.R")


### Load Data
inkar = fread("application/inkar_2021.csv")
krs_spdf = st_read(dsn = "application/vg5000_ebenen_0101", layer = "VG5000_KRS")

### Filter for Kreise and Year of interest for cross sectional data
krs = inkar %>% filter(Raumbezug == "Kreise")
krs = krs %>% filter(Zeitbezug == 2019)

# Create wide format of Kreise for variables 
krs = krs %>%
  distinct(Name, Zeitbezug, Indikator, .keep_all = TRUE) %>%  
  pivot_wider(
    id_cols = c(Name, Kennziffer),                             
    names_from = Indikator,
    values_from = Wert
  )

krs = dplyr::select(krs, where(~ !anyNA(.)))
krs = krs %>% filter(!Name == "Eisenach, Stadt")
krs = krs%>%
  mutate(Kennziffer = sprintf("%05d", Kennziffer))

# Join your data
inkar_spdf = merge(krs_spdf, krs, by.x = "AGS", by.y = "Kennziffer")

# Build the map
tm_shape(inkar_spdf) +
  tm_polygons("Lebenserwartung",
              palette = viridis(n = 100, direction = -1, option = "G"),
              style = "jenks",
              title = "Life Expectancy")

### Generate spatial weight matrix
knn = knearneigh(st_centroid(krs_spdf), k = 10)
nb = knn2nb(knn, row.names = krs_spdf$AGS)
#nb = poly2nb(krs_spdf, queen = TRUE)
listw = nb2listw(nb, style = "W")
W = listw2mat(listw)  

### Filter dependent and independent variables of interest
X = krs[, c("Lebenserwartung",
  "Durchschnittsalter der Bevölkerung",
  "Arbeitslosenquote",
  "Beschäftigtenquote",
  "Erwerbsquote",
  "Selbständigenquote",
  "Schuldnerquote",
  "SGB II - Quote",
  "Beschäftigte mit akademischem Berufsabschluss", 
  "Eheschließungen",
  "Ehescheidungen",
  "Ausländeranteil",
  "Medianeinkommen",
  "Haushaltseinkommen",
  "Verbraucherinsolvenzverfahren",
  "Arbeitsvolumen",
  "Bodenfläche gesamt",
  "Wohnfläche",
  "Siedlungs- und Verkehrsfläche",
  "Waldfläche",
  "Erholungsfläche",
  "Wasserfläche",
  "Mietpreise",
  "Steuereinnahmen",
  "Bruttoinlandsprodukt je Einwohner",
  "Krankenhausversorgung",
  "Ärzte",
  "Pflegebedürftige",
  "Einwohnerdichte",
  "Pkw-Dichte",
  "Pendlersaldo",
  "Straßenverkehrsunfälle",
  "Getötete im Straßenverkehr"
)]

colnames(X) = c(
  "LIFE",
  "AGE",
  "UNEMPLOYMENT",
  "EMPLOYMENT",
  "PART",
  "SELF",
  "DEBT",
  "WELFARE",
  "ACADEMICS",
  "MARRIAGES",
  "DIVORCES",
  "FOREIGN",
  "MEDINC",
  "HHINC",
  "INS",
  "LABOR",
  "LAND",
  "LIVE",
  "URBAN",
  "FOREST",
  "RECR",
  "WATER",
  "RENT",
  "TAX",
  "GDP",
  "HOSP",
  "DR",
  "CARE",
  "POP",
  "CAR",
  "COM",
  "ACC",
  "TRF"
)



for (v in c("AGE", "MARRIAGES", "DIVORCES", "FOREIGN",
  "UNEMPLOYMENT", "EMPLOYMENT",
  "PART", "SELF", "DEBT", "WELFARE",
  "ACADEMICS","MEDINC", "HHINC", "INS", "LABOR",
  "LAND", "LIVE", "URBAN", "FOREST",
  "RECR", "WATER", "RENT", "TAX", "GDP",
  "HOSP", "DR", "CARE",
  "CAR", "POP", "COM", "ACC",
  "TRF"
)) {
  #X[[v]] = asinh(X[[v]])
  #X[[v]] = scale(X[[v]], center = TRUE, scale = FALSE)
  X[[v]] = scale(X[[v]], center = TRUE, scale = TRUE)
}

### Create design matrices
Y = as.matrix(X[,1])
Z = as.data.frame(cbind(X[,-1], W %*% as.matrix(X[,-1])))
colnames(Z)[(ncol(X)):ncol(Z)] = paste0("W_", colnames(Z)[(ncol(X)):ncol(Z)])
Z = data.frame("(Intercept)" = rep(1, length(Y)), Z)
colnames(Z)[1] = "(Intercept)"

gmm = GMerrorsar(Y ~ ., data = Z[,-1], listw = listw)
qml = errorsarlm(Y ~ ., data = X[,-1], listw = listw, Durbin = TRUE)


set.seed(222)
lsgb = gbm(Y, Z, W, M = 12000, start = "ols")
gbgb = gbm(Y, Z, W, M = 5000, start = "boost")
dsgb = gbm(Y, Z, W, M = 1000, start = "des")
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
kable(coef_df, align = "lcccccc", caption = "Coefficient estimates across different estimation strategies for German district life expectancy")

