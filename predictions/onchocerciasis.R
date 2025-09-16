
library (sf) 
library (terra) 
library (raster)
library (ggplot2)
library(viridis)
library (sp)
library (dplyr)
library(readxl)
library(PrevMap)


##### BENIN ####

data_oncho_benin <- read_excel("onchocerciasis.xlsx", sheet = "oncho_benin")

data_oncho_benin$Longitude <- as.numeric(data_oncho_benin$Longitude)
data_oncho_benin$Latitude  <- as.numeric(data_oncho_benin$Latitude)
data_oncho_benin <- subset(data_oncho_benin, !is.na(Longitude) & !is.na(Latitude))

benin_sf_ll <- st_as_sf(data_oncho_benin, coords = c("Longitude", "Latitude"), crs = 4326)

benin_sf_utm <- st_transform(benin_sf_ll, crs = 32631)
utm_mat <- st_coordinates(benin_sf_utm)
data.mod_benin <- cbind(benin_sf_utm, utm_mat)

names(data.mod_benin)[names(data.mod_benin) == "X"] <- "utm_x"
names(data.mod_benin)[names(data.mod_benin) == "Y"] <- "utm_y"

data.mod_benin_df <- cbind(st_drop_geometry(benin_sf_utm), utm_x = utm_mat[,1], utm_y = utm_mat[,2])

cols_to_remove <- c("Country", "ADMIN1_NAME", "ADMIN2_NAME", "IU_NAME", 
                    "LocationName", "LocationType", 
                    "Method_1", "Method_2", 
                    "Age_start", "Age_end", "geometry")

data.mod_benin <- data.mod_benin[, !names(data.mod_benin) %in% cols_to_remove]

data.mod_benin$SurveyYear <- as.numeric(data.mod_benin$SurveyYear)
data.mod_benin_2017 <- subset(data.mod_benin, SurveyYear == 2017)
data.mod_benin_2017 <- st_drop_geometry(data.mod_benin_2017)



## DISTANCE WATER BODIES ## 

wb <- readRDS("benin_water_bundle.rds")
str(wb)  

pts <- st_as_sf(data.mod_benin_2017, coords = c("utm_x", "utm_y"), crs = NA)  
st_crs(pts) <- 32631

if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
  water <- wb$climatology %||% wb$raster %||% wb  
  if (!terra::same.crs(water, st_crs(pts)$wkt)) {
    water <- terra::project(water, st_crs(pts)$wkt)
  }
  
  water_to <- classify(water, rbind(c(-Inf, Inf, 1)))  
  water_to[is.na(water)] <- NA                         
  
  dist_r <- terra::distance(water_to)   
  pts_v <- vect(pts)
  dvals  <- terra::extract(dist_r, pts_v)[,1]  
  
} else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
  rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
  rivers <- st_as_sf(rivers)

  if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set it to the correct UTM.")
  if (st_crs(rivers) != st_crs(pts)) {
    rivers <- st_transform(rivers, st_crs(pts))
  }
  dmat  <- st_distance(pts, rivers)      
  dvals <- apply(dmat, 1, min) |> as.numeric()
} else {
  stop("Don't recognize the contents of benin_water_bundle.rds (not a SpatRaster or sf).")
}

data.mod_benin_2017$dist_to_water_m   <- dvals
data.mod_benin_2017$log_dist_to_water <- scale(log1p(dvals))  


## ANNUAL MEAN TEMPERATURE ##
tmean_2017 <- rast("tmean_2017_benin.tif")

data.mod_benin_2017_sf <- st_as_sf(data.mod_benin_2017,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32631)  
tmean_2017_utm <- project(tmean_2017, "EPSG:32631")
vals_tmean <- terra::extract(tmean_2017_utm, vect(data.mod_benin_2017_sf))
data.mod_benin_2017$mean_temp <- vals_tmean[,2] 
data.mod_benin_2017 <- st_drop_geometry(data.mod_benin_2017)


## PRECIPITATION SEASONALITY ##
prec_season <- rast("wc2.1_2.5m_bio_15.tif") 
data.mod_benin_2017_sf <- st_as_sf(data.mod_benin_2017,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32631)  
prec_season_utm <- project(prec_season, "EPSG:32631")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_benin_2017_sf))
data.mod_benin_2017$precip_seasonality <- vals_prec[,2] 
data.mod_benin_2017 <- st_drop_geometry(data.mod_benin_2017)


## NDVI ##
ndvi_2017_benin <- rast("NDVI_2017_annual_mean_benin.tif")
ndvi_2017_benin_utm <- project(ndvi_2017_benin, "EPSG:32631")
data.mod_benin_2017_sf <- st_as_sf(data.mod_benin_2017,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32631)

vals_ndvi <- terra::extract(ndvi_2017_benin_utm, vect(data.mod_benin_2017_sf))
data.mod_benin_2017$ndvi_mean_benin <- vals_ndvi[,2]
data.mod_benin_2017 <- st_drop_geometry(data.mod_benin_2017)


## SLOPE ##
slope_benin <- rast("slope_benin_1km.tif")
slope_benin_utm <- project(slope_benin, "EPSG:32631")
data.mod_benin_2017_sf <- st_as_sf(data.mod_benin_2017,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32631)

vals_slope <- terra::extract(slope_benin_utm, vect(data.mod_benin_2017_sf))
data.mod_benin_2017$slope_benin <- vals_slope[,2]

data.mod_benin_2017 <- st_drop_geometry(data.mod_benin_2017)


#### PREDICTION BENIN ####
ID.coords.benin <- create.ID.coords(
  data = data.mod_benin_2017,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)
fit.LA.benin <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_benin + slope_benin + log_dist_to_water[,1],
  units.m = ~ Examined,
  coords = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa = 0.5,             
  start.cov.pars = 80,      
  fixed.rel.nugget = 0,     
  data = data.mod_benin_2017,
  ID.coords = ID.coords.benin,
  family = "Binomial"
)
s <- summary(fit.LA.benin)
summary(fit.LA.benin)

cor(data.mod_benin_2017[, c("mean_temp","precip_seasonality","ndvi_mean_benin","slope_benin","log_dist_to_water")])

## MONTE CARLO 
data.mod_benin_2017$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_benin_2017$dist_to_water_m))
)
ID.coords.benin <- create.ID.coords(
  data = data.mod_benin_2017,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

control.mcmc <- control.mcmc.MCML(
  n.sim  = 110000,
  burnin = 10000,
  thin   = 10
)

par0.benin <- coef(fit.LA.benin)
logphi <- s$cov.pars["log(phi)", "Estimate"]
phi_start <- exp(logphi)
phi_start


fit.MCML.benin <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_benin + slope_benin + log_dist_to_water,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.benin,
  kappa          = 0.5,
  start.cov.pars = phi_start,   
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.benin,
  data           = data.mod_benin_2017,
  method         = "nlminb"
)

summary(fit.MCML.benin)


#### PREDICTION PREVALENCE 
benin_0 <- st_read(paste0(benin_dir, "/gadm41_BEN_0.shp"))
benin_0 <- st_transform(benin_0, 32631)


bb   <- st_bbox(benin_0)
nxny <- 250   
xseq <- seq(bb["xmin"], bb["xmax"], length.out = nxny)
yseq <- seq(bb["ymin"], bb["ymax"], length.out = nxny)
grid_df <- expand.grid(x = xseq, y = yseq)


grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 32631)
inside  <- st_within(grid_sf, benin_0, sparse = FALSE)[,1]
grid_sf <- grid_sf[inside, ]
grid_df <- cbind(st_coordinates(grid_sf) |> as.data.frame())
names(grid_df) <- c("x", "y")

extract2 <- function(r, pts, crs = "EPSG:32631") {
  if (inherits(pts, "sf")) {
    v <- terra::vect(pts)  
  } else {
    stopifnot(all(c("x","y") %in% names(pts)))
    v <- terra::vect(pts, geom = c("x","y"), crs = crs)
  }
  terra::extract(r, v)[[2]]  
}

grid_df$mean_temp          <- extract2(tmean_2017_utm,      grid_df)        
grid_df$precip_seasonality <- extract2(prec_season_utm,     grid_df)
grid_df$ndvi_mean_benin    <- extract2(ndvi_2017_benin_utm, grid_df)
grid_df$slope_benin        <- extract2(slope_benin_utm,     grid_df)



`%||%` <- function(x, y) if (is.null(x)) y else x

get_dist_to_water <- function(grid_xy, wb) {
  if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
    water <- wb$climatology %||% wb$raster %||% wb
    if (!terra::same.crs(water, "EPSG:32631")) water <- terra::project(water, "EPSG:32631")
    water_to <- classify(water, rbind(c(-Inf, Inf, 1)))
    water_to[is.na(water)] <- NA
    dist_r <- terra::distance(water_to)  
    terra::extract(dist_r, vect(grid_xy))[[1]]
  } else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
    rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
    if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set to EPSG:32631.")
    if (st_crs(rivers)$epsg != 32631) rivers <- st_transform(rivers, 32631)
    pts_sf <- st_as_sf(grid_xy, coords = c("x", "y"), crs = 32631)
    dmat   <- st_distance(pts_sf, rivers)
    apply(dmat, 1, min) |> as.numeric()
  } else {
    stop("Unrecognized water bundle (not SpatRaster or sf).")
  }
}

grid_df$dist_to_water_m <- get_dist_to_water(grid_df, wb)

train_logd <- log1p(data.mod_benin_2017$dist_to_water_m)
m <- mean(train_logd, na.rm = TRUE)
s <- sd(train_logd,   na.rm = TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s

grid_df <- grid_df |> tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_benin, slope_benin, log_dist_to_water)

grid.pred_benin <- as.matrix(grid_df[, c("x", "y")]) / 1000


pred_benin <- spatial.pred.binomial.MCML(
  object            = fit.MCML.benin,
  grid.pred         = grid.pred_benin,  
  predictors        = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_benin","slope_benin","log_dist_to_water")],
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_benin <- range(pred_benin$prevalence$predictions, na.rm = TRUE)

cat(
  "Predicted prevalence (Benin, 2017) ranges from",
  sprintf("%.3f", rng_benin[1]), "to", sprintf("%.3f", rng_benin[2]), "\n"
)

#Borders 
benin_0 <- st_read(paste0(benin_dir, "/gadm41_BEN_0.shp"))
benin_1 <- st_read(paste0(benin_dir, "/gadm41_BEN_1.shp"))
benin_2 <- st_read(paste0(benin_dir, "/gadm41_BEN_2.shp"))
benin_0 <- st_transform(benin_0, 32631)
benin_1 <- st_transform(benin_1, 32631)
benin_2 <- st_transform(benin_2, 32631)


grid_df$prevalence <- pred_benin$prevalence$predictions

pred_sf_benin <- st_as_sf(
  grid_df,
  coords = c("x", "y"),
  crs = 32631,      
  remove = FALSE
)


## National prevalence map 
ggplot() +
  geom_sf(data = pred_sf_benin, aes(color = prevalence), size = 0.8, alpha = 0.9) +
  geom_sf(data = benin_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent
  ) +
  coord_sf(crs = st_crs(32631)) +
  theme_minimal() +
  labs(
    title = "Predicted Onchocerciasis Prevalence",
    subtitle = "Benin – National Level"
  )

## Average by district 
pred_joined_benin2 <- st_join(pred_sf_benin, benin_2)

oncho_prev_by_district <- pred_joined_benin2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>  
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

benin_avg_oncho_admin2 <- left_join(benin_2, oncho_prev_by_district, by = "NAME_2")

## District level prevalence map 
ggplot() +
  geom_sf(data = benin_avg_oncho_admin2, aes(fill = mean_prev), color = "black", size = 0.3) +
  geom_sf(data = benin_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  labs(
    title = "Onchocerciasis Prevalence by District",
    subtitle = "Benin – District Level"
  )


## Uncertainty map (95% CI width) ##

grid_df$prevalence <- pred_benin$prevalence$predictions
grid_df$lower95    <- pred_benin$prevalence$quantiles[, "2.5%"]
grid_df$upper95    <- pred_benin$prevalence$quantiles[, "97.5%"]
grid_df$ci_width <- grid_df$upper95 - grid_df$lower95

pred_sf_benin <- st_as_sf(
  grid_df,
  coords = c("x", "y"),
  crs = 32631,
  remove = FALSE
)

ggplot() +
  geom_sf(data = pred_sf_benin, aes(color = ci_width), size = 0.8, alpha = 0.9) +
  geom_sf(data = benin_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black",
    high = "white"
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Onchocerciasis Prevalence",
    subtitle = "Benin – National Level"
  )

pred_joined_benin2 <- st_join(pred_sf_benin, benin_2)

oncho_ci_by_district <- pred_joined_benin2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |> 
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

benin_avg_oncho_ci_admin2 <- left_join(benin_2, oncho_ci_by_district, by = "NAME_2")

ggplot() +
  geom_sf(data = benin_avg_oncho_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = benin_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 0.7),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Onchocerciasis Prevalence by District",
    subtitle = "Benin – District Level"
  )


pred_bundle_benin <- list(
  predictions = pred_benin$prevalence$predictions,     
  lower95     = pred_benin$prevalence$quantiles[, "2.5%"],  
  upper95     = pred_benin$prevalence$quantiles[, "97.5%"], 
  ci_width    = pred_benin$prevalence$quantiles[, "97.5%"] -
    pred_benin$prevalence$quantiles[, "2.5%"],  
  grid        = grid_df,        
  grid.pred   = grid.pred_benin,  
  model       = fit.MCML.benin,   
  rng         = rng_benin,      
  date        = Sys.Date(),
  note        = "Onchocerciasis prevalence predictions and uncertainty for Benin, 2017"
)
saveRDS(pred_bundle_benin,
        "benin_oncho_predictions_2017_bundle.rds")

bundle <- readRDS("benin_oncho_predictions_2017_bundle.rds")
str(bundle)

library(writexl)

write_xlsx(data.mod_benin_2017, "benin_oncho_data_2017.xlsx")




#### CIV ####

library (sf) 
library (terra) 
library (raster)
library (ggplot2)
library(viridis)
library (sp)
library (dplyr)
library(readxl)
library(PrevMap)


data_oncho_civ <- read_excel("onchocerciasis.xlsx", sheet = "oncho_civ")
data_oncho_civ$Longitude <- as.numeric(data_oncho_civ$Longitude)
data_oncho_civ$Latitude  <- as.numeric(data_oncho_civ$Latitude)
data_oncho_civ <- subset(data_oncho_civ, !is.na(Longitude) & !is.na(Latitude))

civ_sf_ll <- st_as_sf(data_oncho_civ, coords = c("Longitude", "Latitude"), crs = 4326)
civ_sf_utm <- st_transform(civ_sf_ll, crs = 32630)
utm_mat_civ <- st_coordinates(civ_sf_utm)
data.mod_civ <- cbind(civ_sf_utm, utm_mat_civ)

names(data.mod_civ)[names(data.mod_civ) == "X"] <- "utm_x"
names(data.mod_civ)[names(data.mod_civ) == "Y"] <- "utm_y"

data.mod_civ_df <- cbind(st_drop_geometry(civ_sf_utm),
                         utm_x = utm_mat_civ[, 1],
                         utm_y = utm_mat_civ[, 2])

cols_to_remove <- c("Country", "ADMIN1_NAME", "ADMIN2_NAME", "IU_NAME",
                    "LocationName", "LocationType",
                    "Method_1", "Method_2",
                    "Age_start", "Age_end", "geometry")

data.mod_civ <- data.mod_civ[, !names(data.mod_civ) %in% cols_to_remove]


data.mod_civ$SurveyYear <- as.numeric(data.mod_civ$SurveyYear)
data.mod_civ_2016 <- subset(data.mod_civ, SurveyYear == 2016)
data.mod_civ_2016 <- st_drop_geometry(data.mod_civ_2016)





## DISTANCE WATER BODIES ## 

wb <- readRDS("civ_water_bundle.rds")
str(wb) 

pts <- st_as_sf(data.mod_civ_2016, coords = c("utm_x", "utm_y"), crs = NA)  
st_crs(pts) <- 32630

if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
  water <- wb$climatology %||% wb$raster %||% wb  
  if (!terra::same.crs(water, st_crs(pts)$wkt)) {
    water <- terra::project(water, st_crs(pts)$wkt)
  }
  water_to <- classify(water, rbind(c(-Inf, Inf, 1)))  
  water_to[is.na(water)] <- NA                         
  
  dist_r <- terra::distance(water_to)   
  pts_v <- vect(pts)
  dvals  <- terra::extract(dist_r, pts_v)[,1] 
  
} else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
  rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
  rivers <- st_as_sf(rivers)
  if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set it to the correct UTM.")
  if (st_crs(rivers) != st_crs(pts)) {
    rivers <- st_transform(rivers, st_crs(pts))
  }

    dmat  <- st_distance(pts, rivers)      
  dvals <- apply(dmat, 1, min) |> as.numeric()
} else {
  stop("Don't recognize the contents of civ_water_bundle.rds (not a SpatRaster or sf).")
}

data.mod_civ_2016$dist_to_water_m   <- dvals
data.mod_civ_2016$log_dist_to_water <- scale(log1p(dvals))

## ANNUAL MEAN TEMPERATURE ##
tmean_2016 <- rast("tmean_2016_civ.tif")

data.mod_civ_2016_sf <- st_as_sf(data.mod_civ_2016,
                                 coords = c("utm_x", "utm_y"),
                                 crs = 32630) 
tmean_2016_utm <- project(tmean_2016, "EPSG:32630")
vals_tmean <- terra::extract(tmean_2016_utm, vect(data.mod_civ_2016_sf))
data.mod_civ_2016$mean_temp <- vals_tmean[,2] 
data.mod_civ_2016 <- st_drop_geometry(data.mod_civ_2016)


## PRECIPITATION SEASONALITY ##
prec_season <- rast("wc2.1_2.5m_bio_15.tif") 
data.mod_civ_2016_sf <- st_as_sf(data.mod_civ_2016,
                                 coords = c("utm_x", "utm_y"),
                                 crs = 32630) 
prec_season_utm <- project(prec_season, "EPSG:32630")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_civ_2016_sf))
data.mod_civ_2016$precip_seasonality <- vals_prec[,2] 
data.mod_civ_2016 <- st_drop_geometry(data.mod_civ_2016)


## NDVI ##
ndvi_2016_civ <- rast("NDVI_2016_annual_mean_civ.tif")
ndvi_2016_civ_utm <- project(ndvi_2016_civ, "EPSG:32630")
data.mod_civ_2016_sf <- st_as_sf(data.mod_civ_2016,
                                 coords = c("utm_x", "utm_y"),
                                 crs = 32630)

vals_ndvi <- terra::extract(ndvi_2016_civ_utm, vect(data.mod_civ_2016_sf))
data.mod_civ_2016$ndvi_mean_civ <- vals_ndvi[,2]
data.mod_civ_2016 <- st_drop_geometry(data.mod_civ_2016)



## SLOPE ##
slope_civ <- rast("slope_civ_1km.tif")

slope_civ_utm <- project(slope_civ, "EPSG:32630")
data.mod_civ_2016_sf <- st_as_sf(data.mod_civ_2016,
                                 coords = c("utm_x", "utm_y"),
                                 crs = 32630)

vals_slope <- terra::extract(slope_civ_utm, vect(data.mod_civ_2016_sf))
data.mod_civ_2016$slope_civ <- vals_slope[,2]

data.mod_civ_2016 <- st_drop_geometry(data.mod_civ_2016)


#### PREDICTION CIV ####
ID.coords.civ <- create.ID.coords(
  data = data.mod_civ_2016,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)
fit.LA.civ <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_civ + slope_civ + log_dist_to_water[,1],
  units.m = ~ Examined,
  coords = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa = 0.5,              
  start.cov.pars = 80,      
  fixed.rel.nugget = 0,     
  data = data.mod_civ_2016,
  ID.coords = ID.coords.civ,
  family = "Binomial"
)

s <- summary(fit.LA.civ)
summary(fit.LA.civ)

cor(data.mod_civ_2016[, c("mean_temp","precip_seasonality","ndvi_mean_civ","slope_civ","log_dist_to_water")])

## MONTE CARLO 
data.mod_civ_2016$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_civ_2016$dist_to_water_m))
)
ID.coords.civ <- create.ID.coords(
  data = data.mod_civ_2016,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

control.mcmc <- control.mcmc.MCML(
  n.sim  = 110000,
  burnin = 10000,
  thin   = 10
)

par0.civ <- coef(fit.LA.civ)
logphi <- s$cov.pars["log(phi)", "Estimate"]
phi_start <- exp(logphi)
phi_start


fit.MCML.civ <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_civ + slope_civ + log_dist_to_water,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.civ,
  kappa          = 0.5,
  start.cov.pars = phi_start,   
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.civ,
  data           = data.mod_civ_2016,
  method         = "nlminb"
)

summary(fit.MCML.civ)


#### PREDICTION PREVALENCE 
civ_0 <- st_read(paste0(civ_dir, "/gadm41_CIV_0.shp"))
civ_0 <- st_transform(civ_0, 32630)


bb   <- st_bbox(civ_0)
nxny <- 250   
xseq <- seq(bb["xmin"], bb["xmax"], length.out = nxny)
yseq <- seq(bb["ymin"], bb["ymax"], length.out = nxny)
grid_df <- expand.grid(x = xseq, y = yseq)


grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 32630)
inside  <- st_within(grid_sf, civ_0, sparse = FALSE)[,1]
grid_sf <- grid_sf[inside, ]
grid_df <- cbind(st_coordinates(grid_sf) |> as.data.frame())
names(grid_df) <- c("x", "y")

extract2 <- function(r, pts, crs = "EPSG:32630") {
  if (inherits(pts, "sf")) {
    v <- terra::vect(pts)  
  } else {
    stopifnot(all(c("x","y") %in% names(pts)))
    v <- terra::vect(pts, geom = c("x","y"), crs = crs)
  }
  terra::extract(r, v)[[2]] 
}

grid_df$mean_temp          <- extract2(tmean_2016_utm,      grid_df)        
grid_df$precip_seasonality <- extract2(prec_season_utm,     grid_df)
grid_df$ndvi_mean_civ    <- extract2(ndvi_2016_civ_utm, grid_df)
grid_df$slope_civ        <- extract2(slope_civ_utm,     grid_df)


`%||%` <- function(x, y) if (is.null(x)) y else x

get_dist_to_water <- function(grid_xy, wb) {
  if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
    water <- wb$climatology %||% wb$raster %||% wb
    if (!terra::same.crs(water, "EPSG:32630")) water <- terra::project(water, "EPSG:32630")
    water_to <- classify(water, rbind(c(-Inf, Inf, 1)))
    water_to[is.na(water)] <- NA
    dist_r <- terra::distance(water_to)  
    terra::extract(dist_r, vect(grid_xy))[[1]]
  } else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
    rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
    if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set to EPSG:32630.")
    if (st_crs(rivers)$epsg != 32630) rivers <- st_transform(rivers, 32630)
    pts_sf <- st_as_sf(grid_xy, coords = c("x", "y"), crs = 32630)
    dmat   <- st_distance(pts_sf, rivers)
    apply(dmat, 1, min) |> as.numeric()
  } else {
    stop("Unrecognized water bundle (not SpatRaster or sf).")
  }
}

grid_df$dist_to_water_m <- get_dist_to_water(grid_df, wb)

train_logd <- log1p(data.mod_civ_2016$dist_to_water_m)
m <- mean(train_logd, na.rm = TRUE)
s <- sd(train_logd,   na.rm = TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s
grid_df <- grid_df |> tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_civ, slope_civ, log_dist_to_water)

grid.pred_civ <- as.matrix(grid_df[, c("x", "y")]) / 1000


pred_civ <- spatial.pred.binomial.MCML(
  object            = fit.MCML.civ,
  grid.pred         = grid.pred_civ,  
  predictors        = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_civ","slope_civ","log_dist_to_water")],
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_civ <- range(pred_civ$prevalence$predictions, na.rm = TRUE)

cat(
  "Predicted prevalence (civ, 2016) ranges from",
  sprintf("%.3f", rng_civ[1]), "to", sprintf("%.3f", rng_civ[2]), "\n"
)

## Borders ##
civ_0 <- st_read(paste0(civ_dir, "/gadm41_CIV_0.shp"))
civ_1 <- st_read(paste0(civ_dir, "/gadm41_CIV_1.shp"))
civ_2 <- st_read(paste0(civ_dir, "/gadm41_CIV_2.shp"))
civ_0 <- st_transform(civ_0, 32630)
civ_1 <- st_transform(civ_1, 32630)
civ_2 <- st_transform(civ_2, 32630)


grid_df$prevalence <- pred_civ$prevalence$predictions

pred_sf_civ <- st_as_sf(
  grid_df,
  coords = c("x", "y"),
  crs = 32630,      
  remove = FALSE
)


## National prevalence map ##
ggplot() +
  geom_sf(data = pred_sf_civ, aes(color = prevalence), size = 0.8, alpha = 0.9) +
  geom_sf(data = civ_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent
  ) +
  coord_sf(crs = st_crs(32630)) +
  theme_minimal() +
  labs(
    title = "Predicted Onchocerciasis Prevalence",
    subtitle = "civ – National Level"
  )

## Average by district ##
pred_joined_civ2 <- st_join(pred_sf_civ, civ_2)

oncho_prev_by_district <- pred_joined_civ2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>   
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

civ_avg_oncho_admin2 <- left_join(civ_2, oncho_prev_by_district, by = "NAME_2")

## District level prevalence map ##
ggplot() +
  geom_sf(data = civ_avg_oncho_admin2, aes(fill = mean_prev), color = "black", size = 0.3) +
  geom_sf(data = civ_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  labs(
    title = "Onchocerciasis Prevalence by District",
    subtitle = "civ – District Level"
  )



## Uncertainty map (95% CI width) ##

grid_df$prevalence <- pred_civ$prevalence$predictions
grid_df$lower95    <- pred_civ$prevalence$quantiles[, "2.5%"]
grid_df$upper95    <- pred_civ$prevalence$quantiles[, "97.5%"]

grid_df$ci_width <- grid_df$upper95 - grid_df$lower95

pred_sf_civ <- st_as_sf(
  grid_df_civ,
  coords = c("x", "y"),
  crs = 32630,  
  remove = FALSE
)

# National uncertainty map
ggplot() +
  geom_sf(data = pred_sf_civ, aes(color = ci_width), size = 0.8, alpha = 0.9) +
  geom_sf(data = civ_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black",
    high = "white"
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Onchocerciasis Prevalence",
    subtitle = "Côte d’Ivoire – National Level"
  )

pred_joined_civ2 <- st_join(pred_sf_civ, civ_2)

oncho_ci_by_district_civ <- pred_joined_civ2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |> 
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

civ_avg_oncho_ci_admin2 <- left_join(civ_2, oncho_ci_by_district_civ, by = "NAME_2")

# District level map
ggplot() +
  geom_sf(data = civ_avg_oncho_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = civ_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Onchocerciasis Prevalence by District",
    subtitle = "Côte d’Ivoire – District Level"
  )





# save as bundle
pred_bundle_civ <- list(
  predictions = pred_civ$prevalence$predictions,
  lower95     = pred_civ$prevalence$quantiles[, "2.5%"],
  upper95     = pred_civ$prevalence$quantiles[, "97.5%"],
  ci_width    = pred_civ$prevalence$quantiles[, "97.5%"] -
    pred_civ$prevalence$quantiles[, "2.5%"],
  
  grid        = grid_df,        
  grid.pred   = grid.pred_civ,  
  model       = fit.MCML.civ,    
  rng         = rng_civ,        
  date        = Sys.Date(),
  note        = "Onchocerciasis prevalence predictions + uncertainty for Côte d’Ivoire, 2016"
)

saveRDS(pred_bundle_civ,
        "civ_oncho_predictions_2016_bundle.rds")

bundle <- readRDS("civ_oncho_predictions_2016_bundle.rds")
str(bundle)

library(writexl)

write_xlsx(data.mod_civ_2016, "civ_oncho_data_2016.xlsx")



#### GHANA ####

library (sf) 
library (terra) 
library (raster)
library (ggplot2)
library(viridis)
library (sp)
library (dplyr)
library(readxl)
library(PrevMap)


data_oncho_ghana <- read_excel("onchocerciasis.xlsx", sheet = "oncho_ghana")
data_oncho_ghana$Longitude <- as.numeric(data_oncho_ghana$Longitude)
data_oncho_ghana$Latitude  <- as.numeric(data_oncho_ghana$Latitude)
data_oncho_ghana <- subset(data_oncho_ghana, !is.na(Longitude) & !is.na(Latitude))
ghana_sf_ll <- st_as_sf(data_oncho_ghana, coords = c("Longitude", "Latitude"), crs = 4326)

ghana_sf_utm <- st_transform(ghana_sf_ll, crs = 32630)
utm_mat_ghana <- st_coordinates(ghana_sf_utm)
data.mod_ghana <- cbind(ghana_sf_utm, utm_mat_ghana)

names(data.mod_ghana)[names(data.mod_ghana) == "X"] <- "utm_x"
names(data.mod_ghana)[names(data.mod_ghana) == "Y"] <- "utm_y"

data.mod_ghana_df <- cbind(st_drop_geometry(ghana_sf_utm),
                           utm_x = utm_mat_ghana[, 1],
                           utm_y = utm_mat_ghana[, 2])

cols_to_remove <- c("Country", "ADMIN1_NAME", "ADMIN2_NAME", "IU_NAME",
                    "LocationName", "LocationType",
                    "Method_1", "Method_2",
                    "Age_start", "Age_end", "geometry")

data.mod_ghana <- data.mod_ghana[, !names(data.mod_ghana) %in% cols_to_remove]


data.mod_ghana$SurveyYear <- as.numeric(data.mod_ghana$SurveyYear)
data.mod_ghana_2012 <- subset(data.mod_ghana, SurveyYear == 2012)
data.mod_ghana_2012 <- st_drop_geometry(data.mod_ghana_2012)





## DISTANCE WATER BODIES ## 

wb <- readRDS("ghana_water_bundle.rds")
str(wb)   

pts <- st_as_sf(data.mod_ghana_2012, coords = c("utm_x", "utm_y"), crs = NA)  
st_crs(pts) <- 32630  

if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
  water <- wb$climatology %||% wb$raster %||% wb
  if (!terra::same.crs(water, st_crs(pts)$wkt)) {
    water <- terra::project(water, st_crs(pts)$wkt)
  }
  water_to <- classify(water, rbind(c(-Inf, Inf, 1)))  
  water_to[is.na(water)] <- NA                         
  
  dist_r <- terra::distance(water_to)   
  pts_v <- vect(pts)
  dvals  <- terra::extract(dist_r, pts_v)[,1]  
  
} else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
  rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
  rivers <- st_as_sf(rivers)
  if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set it to the correct UTM.")
  if (st_crs(rivers) != st_crs(pts)) {
    rivers <- st_transform(rivers, st_crs(pts))
  }

    dmat  <- st_distance(pts, rivers)      
  dvals <- apply(dmat, 1, min) |> as.numeric()
} else {
  stop("Don't recognize the contents of ghana_water_bundle.rds (not a SpatRaster or sf).")
}

data.mod_ghana_2012$dist_to_water_m   <- dvals
data.mod_ghana_2012$log_dist_to_water <- scale(log1p(dvals)) 

## ANNUAL MEAN TEMPERATURE ##
tmean_2012 <- rast("tmean_2012_ghana.tif")

data.mod_ghana_2012_sf <- st_as_sf(data.mod_ghana_2012,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32630) 
tmean_2012_utm <- project(tmean_2012, "EPSG:32630")
vals_tmean <- terra::extract(tmean_2012_utm, vect(data.mod_ghana_2012_sf))
data.mod_ghana_2012$mean_temp <- vals_tmean[,2] 
data.mod_ghana_2012 <- st_drop_geometry(data.mod_ghana_2012)


## PRECIPITATION SEASONALITY ##
prec_season <- rast("wc2.1_2.5m_bio_15.tif") 
data.mod_ghana_2012_sf <- st_as_sf(data.mod_ghana_2012,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32630)  
prec_season_utm <- project(prec_season, "EPSG:32630")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_ghana_2012_sf))
data.mod_ghana_2012$precip_seasonality <- vals_prec[,2] 
data.mod_ghana_2012 <- st_drop_geometry(data.mod_ghana_2012)


## NDVI ##
ndvi_2012_ghana <- rast("NDVI_2012_annual_mean_gha.tif")
ndvi_2012_ghana_utm <- project(ndvi_2012_ghana, "EPSG:32630")
data.mod_ghana_2012_sf <- st_as_sf(data.mod_ghana_2012,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32630)

vals_ndvi <- terra::extract(ndvi_2012_ghana_utm, vect(data.mod_ghana_2012_sf))
data.mod_ghana_2012$ndvi_mean_ghana <- vals_ndvi[,2]
data.mod_ghana_2012 <- st_drop_geometry(data.mod_ghana_2012)



## SLOPE ##
slope_ghana <- rast("slope_ghana_1km.tif")

slope_ghana_utm <- project(slope_ghana, "EPSG:32630")
data.mod_ghana_2012_sf <- st_as_sf(data.mod_ghana_2012,
                                   coords = c("utm_x", "utm_y"),
                                   crs = 32630)

vals_slope <- terra::extract(slope_ghana_utm, vect(data.mod_ghana_2012_sf))
data.mod_ghana_2012$slope_ghana <- vals_slope[,2]

data.mod_ghana_2012 <- st_drop_geometry(data.mod_ghana_2012)


#### PREDICTION GHANA ####
ID.coords.ghana <- create.ID.coords(
  data = data.mod_ghana_2012,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)
fit.LA.ghana <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_ghana + slope_ghana + log_dist_to_water[,1],
  units.m = ~ Examined,
  coords = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa = 0.5,              
  start.cov.pars = 20,      
  fixed.rel.nugget = 0,     
  data = data.mod_ghana_2012,
  ID.coords = ID.coords.ghana,
  family = "Binomial"
)

s <- summary(fit.LA.ghana)
summary(fit.LA.ghana)

cor(data.mod_ghana_2012[, c("mean_temp","precip_seasonality","ndvi_mean_ghana","slope_ghana","log_dist_to_water")])

## MONTE CARLO 
data.mod_ghana_2012$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_ghana_2012$dist_to_water_m))
)
ID.coords.ghana <- create.ID.coords(
  data = data.mod_ghana_2012,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

control.mcmc <- control.mcmc.MCML(
  n.sim  = 110000,
  burnin = 10000,
  thin   = 10
)

par0.ghana <- coef(fit.LA.ghana)
logphi <- s$cov.pars["log(phi)", "Estimate"]
phi_start <- exp(logphi)
phi_start


fit.MCML.ghana <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_ghana + slope_ghana + log_dist_to_water,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.ghana,
  kappa          = 0.5,
  start.cov.pars = phi_start,   
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.ghana,
  data           = data.mod_ghana_2012,
  method         = "nlminb"
)


summary(fit.MCML.ghana)


#### PREDICTION PREVALENCE 
ghana_0 <- st_read(paste0(ghana_dir, "/gadm41_gha_0.shp"))
ghana_0 <- st_transform(ghana_0, 32630)


bb   <- st_bbox(ghana_0)
nxny <- 250   
xseq <- seq(bb["xmin"], bb["xmax"], length.out = nxny)
yseq <- seq(bb["ymin"], bb["ymax"], length.out = nxny)
grid_df <- expand.grid(x = xseq, y = yseq)


grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 32630)
inside  <- st_within(grid_sf, ghana_0, sparse = FALSE)[,1]
grid_sf <- grid_sf[inside, ]
grid_df <- cbind(st_coordinates(grid_sf) |> as.data.frame())
names(grid_df) <- c("x", "y")


extract2 <- function(r, pts, crs = "EPSG:32630") {
  if (inherits(pts, "sf")) {
    v <- terra::vect(pts)  
  } else {
    stopifnot(all(c("x","y") %in% names(pts)))
    v <- terra::vect(pts, geom = c("x","y"), crs = crs)
  }
  terra::extract(r, v)[[2]]  
}

grid_df$mean_temp          <- extract2(tmean_2012_utm,      grid_df)        
grid_df$precip_seasonality <- extract2(prec_season_utm,     grid_df)
grid_df$ndvi_mean_ghana    <- extract2(ndvi_2012_ghana_utm, grid_df)
grid_df$slope_ghana        <- extract2(slope_ghana_utm,     grid_df)


`%||%` <- function(x, y) if (is.null(x)) y else x

get_dist_to_water <- function(grid_xy, wb) {
  if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
    water <- wb$climatology %||% wb$raster %||% wb
    if (!terra::same.crs(water, "EPSG:32630")) water <- terra::project(water, "EPSG:32630")
    water_to <- classify(water, rbind(c(-Inf, Inf, 1)))
    water_to[is.na(water)] <- NA
    dist_r <- terra::distance(water_to)  
    terra::extract(dist_r, vect(grid_xy))[[1]]
  } else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
    rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
    if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set to EPSG:32630.")
    if (st_crs(rivers)$epsg != 32630) rivers <- st_transform(rivers, 32630)
    pts_sf <- st_as_sf(grid_xy, coords = c("x", "y"), crs = 32630)
    dmat   <- st_distance(pts_sf, rivers)
    apply(dmat, 1, min) |> as.numeric()
  } else {
    stop("Unrecognized water bundle (not SpatRaster or sf).")
  }
}

grid_df$dist_to_water_m <- get_dist_to_water(grid_df, wb)

train_logd <- log1p(data.mod_ghana_2012$dist_to_water_m)
m <- mean(train_logd, na.rm = TRUE)
s <- sd(train_logd,   na.rm = TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s

grid_df <- grid_df |> tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_ghana, slope_ghana, log_dist_to_water)
grid.pred_ghana <- as.matrix(grid_df[, c("x", "y")]) / 1000


pred_ghana <- spatial.pred.binomial.MCML(
  object            = fit.MCML.ghana,
  grid.pred         = grid.pred_ghana,  
  predictors        = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_ghana","slope_ghana","log_dist_to_water")],
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_ghana <- range(pred_ghana$prevalence$predictions, na.rm = TRUE)

cat(
  "Predicted prevalence (ghana, 2012) ranges from",
  sprintf("%.3f", rng_ghana[1]), "to", sprintf("%.3f", rng_ghana[2]), "\n"
)

ghana_0 <- st_read(paste0(ghana_dir, "/gadm41_gha_0.shp"))
ghana_1 <- st_read(paste0(ghana_dir, "/gadm41_gha_1.shp"))
ghana_2 <- st_read(paste0(ghana_dir, "/gadm41_gha_2.shp"))
ghana_0 <- st_transform(ghana_0, 32630)
ghana_1 <- st_transform(ghana_1, 32630)
ghana_2 <- st_transform(ghana_2, 32630)


grid_df$prevalence <- pred_ghana$prevalence$predictions

pred_sf_ghana <- st_as_sf(
  grid_df,
  coords = c("x", "y"),
  crs = 32630,     
  remove = FALSE
)


## National prevalence map ##
ggplot() +
  geom_sf(data = pred_sf_ghana, aes(color = prevalence), size = 0.8, alpha = 0.9) +
  geom_sf(data = ghana_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent
  ) +
  coord_sf(crs = st_crs(32630)) +
  theme_minimal() +
  labs(
    title = "Predicted Onchocerciasis Prevalence",
    subtitle = "ghana – National Level"
  )

## Average by district ##
pred_joined_ghana2 <- st_join(pred_sf_ghana, ghana_2)

oncho_prev_by_district <- pred_joined_ghana2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>  
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

ghana_avg_oncho_admin2 <- left_join(ghana_2, oncho_prev_by_district, by = "NAME_2")


## District level prevalence map ##
ggplot() +
  geom_sf(data = ghana_avg_oncho_admin2, aes(fill = mean_prev), color = "black", size = 0.3) +
  geom_sf(data = ghana_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  labs(
    title = "Onchocerciasis Prevalence by District",
    subtitle = "ghana – District Level"
  )




## Uncertainty map (95% CI width) – Ghana ##

grid_df$prevalence <- pred_ghana$prevalence$predictions
grid_df$lower95    <- pred_ghana$prevalence$quantiles[, "2.5%"]
grid_df$upper95    <- pred_ghana$prevalence$quantiles[, "97.5%"]

grid_df$ci_width <- grid_df$upper95 - grid_df$lower95
pred_sf_ghana <- st_as_sf(
  grid_df,
  coords = c("x", "y"),
  crs = 32630,  
  remove = FALSE
)


# National uncertainty map

ggplot() +
  geom_sf(data = pred_sf_ghana, aes(color = ci_width), size = 0.8, alpha = 0.9) +
  geom_sf(data = ghana_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black",
    high = "white"
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Onchocerciasis Prevalence",
    subtitle = "Ghana – National Level"
  )


pred_joined_ghana2 <- st_join(pred_sf_ghana, ghana_2)

oncho_ci_by_district_ghana <- pred_joined_ghana2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |> 
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

ghana_avg_oncho_ci_admin2 <- left_join(ghana_2, oncho_ci_by_district_ghana, by = "NAME_2")


# District level map
ggplot() +
  geom_sf(data = ghana_avg_oncho_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = ghana_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 0.7),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Onchocerciasis Prevalence by District",
    subtitle = "Ghana – District Level"
  )



pred_bundle_ghana <- list(
  predictions = pred_ghana$prevalence$predictions,
  lower95     = pred_ghana$prevalence$quantiles[, "2.5%"],
  upper95     = pred_ghana$prevalence$quantiles[, "97.5%"],
  ci_width    = pred_ghana$prevalence$quantiles[, "97.5%"] -
    pred_ghana$prevalence$quantiles[, "2.5%"],
  
  grid        = grid_df,         
  grid.pred   = grid.pred_ghana, 
  model       = fit.MCML.ghana,  
  rng         = rng_ghana,       
  date        = Sys.Date(),
  note        = "Onchocerciasis prevalence predictions + uncertainty for Ghana, 2012"
)

saveRDS(pred_bundle_ghana,
        "ghana_oncho_predictions_2012_bundle.rds")
bundle <- readRDS("ghana_oncho_predictions_2012_bundle.rds")
str(bundle)

library(writexl)
write_xlsx(data.mod_ghana_2012, "ghana_oncho_data_2012.xlsx")


#### TOGO ####

library (sf) 
library (terra) 
library (raster)
library (ggplot2)
library(viridis)
library (sp)
library (dplyr)
library(readxl)
library(PrevMap)


data_oncho_togo <- read_excel("onchocerciasis.xlsx", sheet = "oncho_togo")

data_oncho_togo$Longitude <- as.numeric(data_oncho_togo$Longitude)
data_oncho_togo$Latitude  <- as.numeric(data_oncho_togo$Latitude)
data_oncho_togo <- subset(data_oncho_togo, !is.na(Longitude) & !is.na(Latitude))

togo_sf_ll <- st_as_sf(data_oncho_togo, coords = c("Longitude", "Latitude"), crs = 4326)
togo_sf_utm <- st_transform(togo_sf_ll, crs = 32631)
utm_mat <- st_coordinates(togo_sf_utm)
data.mod_togo <- cbind(togo_sf_utm, utm_mat)

names(data.mod_togo)[names(data.mod_togo) == "X"] <- "utm_x"
names(data.mod_togo)[names(data.mod_togo) == "Y"] <- "utm_y"

data.mod_togo_df <- cbind(st_drop_geometry(togo_sf_utm), utm_x = utm_mat[,1], utm_y = utm_mat[,2])

cols_to_remove <- c("Country", "ADMIN1_NAME", "ADMIN2_NAME", "IU_NAME", 
                    "LocationName", "LocationType", 
                    "Method_1", "Method_2", 
                    "Age_start", "Age_end", "geometry")

data.mod_togo <- data.mod_togo[, !names(data.mod_togo) %in% cols_to_remove]

data.mod_togo$SurveyYear <- as.numeric(data.mod_togo$SurveyYear)
data.mod_togo_2015 <- subset(data.mod_togo, SurveyYear == 2015)
data.mod_togo_2015 <- st_drop_geometry(data.mod_togo_2015)



## DISTANCE WATER BODIES ## 

wb <- readRDS("togo_water_bundle.rds")
str(wb) 

pts <- st_as_sf(data.mod_togo_2015, coords = c("utm_x", "utm_y"), crs = NA)  
st_crs(pts) <- 32631 

if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
  water <- wb$climatology %||% wb$raster %||% wb  

    if (!terra::same.crs(water, st_crs(pts)$wkt)) {
    water <- terra::project(water, st_crs(pts)$wkt)
  }

  water_to <- classify(water, rbind(c(-Inf, Inf, 1))) 
  water_to[is.na(water)] <- NA                         
  
  dist_r <- terra::distance(water_to)   
  pts_v <- vect(pts)
  dvals  <- terra::extract(dist_r, pts_v)[,1]
  
} else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
  rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
  rivers <- st_as_sf(rivers)
  if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set it to the correct UTM.")
  if (st_crs(rivers) != st_crs(pts)) {
    rivers <- st_transform(rivers, st_crs(pts))
  }
  dmat  <- st_distance(pts, rivers)      
  dvals <- apply(dmat, 1, min) |> as.numeric()
} else {
  stop("Don't recognize the contents of togo_water_bundle.rds (not a SpatRaster or sf).")
}

data.mod_togo_2015$dist_to_water_m   <- dvals
data.mod_togo_2015$log_dist_to_water <- scale(log1p(dvals))  


## ANNUAL MEAN TEMPERATURE ##
tmean_2015 <- rast("tmean_2015_togo.tif")

data.mod_togo_2015_sf <- st_as_sf(data.mod_togo_2015,
                                  coords = c("utm_x", "utm_y"),
                                  crs = 32631)  # UTM 31N
tmean_2015_utm <- project(tmean_2015, "EPSG:32631")
vals_tmean <- terra::extract(tmean_2015_utm, vect(data.mod_togo_2015_sf))
data.mod_togo_2015$mean_temp <- vals_tmean[,2] 
data.mod_togo_2015 <- st_drop_geometry(data.mod_togo_2015)


## PRECIPITATION SEASONALITY ##
prec_season <- rast("wc2.1_2.5m_bio_15.tif") 
data.mod_togo_2015_sf <- st_as_sf(data.mod_togo_2015,
                                  coords = c("utm_x", "utm_y"),
                                  crs = 32631) 
prec_season_utm <- project(prec_season, "EPSG:32631")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_togo_2015_sf))
data.mod_togo_2015$precip_seasonality <- vals_prec[,2] 
data.mod_togo_2015 <- st_drop_geometry(data.mod_togo_2015)


## NDVI ##
ndvi_2015_togo <- rast("NDVI_2015_annual_mean_tgo.tif")
ndvi_2015_togo_utm <- project(ndvi_2015_togo, "EPSG:32631")
data.mod_togo_2015_sf <- st_as_sf(data.mod_togo_2015,
                                  coords = c("utm_x", "utm_y"),
                                  crs = 32631)

vals_ndvi <- terra::extract(ndvi_2015_togo_utm, vect(data.mod_togo_2015_sf))
data.mod_togo_2015$ndvi_mean_togo <- vals_ndvi[,2]
data.mod_togo_2015 <- st_drop_geometry(data.mod_togo_2015)


## SLOPE ##
slope_togo <- rast("slope_togo_1km.tif")
slope_togo_utm <- project(slope_togo, "EPSG:32631")
data.mod_togo_2015_sf <- st_as_sf(data.mod_togo_2015,
                                  coords = c("utm_x", "utm_y"),
                                  crs = 32631)

vals_slope <- terra::extract(slope_togo_utm, vect(data.mod_togo_2015_sf))
data.mod_togo_2015$slope_togo <- vals_slope[,2]

data.mod_togo_2015 <- st_drop_geometry(data.mod_togo_2015)
data.mod_togo_2015 <- na.omit(data.mod_togo_2015)

#### PREDICTION togo ####
ID.coords.togo <- create.ID.coords(
  data = data.mod_togo_2015,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)
fit.LA.togo <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_togo + slope_togo + log_dist_to_water[,1],
  units.m = ~ Examined,
  coords = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa = 0.5,              
  start.cov.pars = 50,      
  fixed.rel.nugget = 0,     
  data = data.mod_togo_2015,
  ID.coords = ID.coords.togo,
  family = "Binomial"
)


s <- summary(fit.LA.togo)
summary(fit.LA.togo)
  
  
  cor(data.mod_togo_2015[, c("mean_temp","precip_seasonality","ndvi_mean_togo","slope_togo","log_dist_to_water")])

  

## MONTE CARLO 
data.mod_togo_2015$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_togo_2015$dist_to_water_m))
)
ID.coords.togo <- create.ID.coords(
  data = data.mod_togo_2015,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

control.mcmc <- control.mcmc.MCML(
  n.sim  = 110000,
  burnin = 10000,
  thin   = 10
)

par0.togo <- coef(fit.LA.togo)
logphi <- s$cov.pars["log(phi)", "Estimate"]
phi_start <- exp(logphi)
phi_start


fit.MCML.togo <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_togo + slope_togo + log_dist_to_water,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.togo,
  kappa          = 0.5,
  start.cov.pars = phi_start,   
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.togo,
  data           = data.mod_togo_2015,
  method         = "nlminb"
)

summary(fit.MCML.togo)
  
  
#### PREDICTION PREVALENCE 
togo_0 <- st_read(paste0(togo_dir, "/gadm41_TGO_0.shp"))
togo_0 <- st_transform(togo_0, 32631)


bb   <- st_bbox(togo_0)
nxny <- 250   
xseq <- seq(bb["xmin"], bb["xmax"], length.out = nxny)
yseq <- seq(bb["ymin"], bb["ymax"], length.out = nxny)
grid_df <- expand.grid(x = xseq, y = yseq)

grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 32631)
inside  <- st_within(grid_sf, togo_0, sparse = FALSE)[,1]
grid_sf <- grid_sf[inside, ]
grid_df <- cbind(st_coordinates(grid_sf) |> as.data.frame())
names(grid_df) <- c("x", "y")


extract2 <- function(r, pts, crs = "EPSG:32631") {
  if (inherits(pts, "sf")) {
    v <- terra::vect(pts) 
  } else {
    stopifnot(all(c("x","y") %in% names(pts)))
    v <- terra::vect(pts, geom = c("x","y"), crs = crs)
  }
  terra::extract(r, v)[[2]] 
}

grid_df$mean_temp          <- extract2(tmean_2015_utm,      grid_df)        
grid_df$precip_seasonality <- extract2(prec_season_utm,     grid_df)
grid_df$ndvi_mean_togo    <- extract2(ndvi_2015_togo_utm, grid_df)
grid_df$slope_togo        <- extract2(slope_togo_utm,     grid_df)



`%||%` <- function(x, y) if (is.null(x)) y else x

get_dist_to_water <- function(grid_xy, wb) {
  if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
    water <- wb$climatology %||% wb$raster %||% wb
    if (!terra::same.crs(water, "EPSG:32631")) water <- terra::project(water, "EPSG:32631")
    water_to <- classify(water, rbind(c(-Inf, Inf, 1)))
    water_to[is.na(water)] <- NA
    dist_r <- terra::distance(water_to)
    terra::extract(dist_r, vect(grid_xy))[[1]]
  } else if (inherits(wb, "sf") || inherits(wb$rivers, "sf") || inherits(wb$water, "sf") || inherits(wb$lines, "sf")) {
    rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
    if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set to EPSG:32631.")
    if (st_crs(rivers)$epsg != 32631) rivers <- st_transform(rivers, 32631)
    pts_sf <- st_as_sf(grid_xy, coords = c("x", "y"), crs = 32631)
    dmat   <- st_distance(pts_sf, rivers)
    apply(dmat, 1, min) |> as.numeric()
  } else {
    stop("Unrecognized water bundle (not SpatRaster or sf).")
  }
}

grid_df$dist_to_water_m <- get_dist_to_water(grid_df, wb)
train_logd <- log1p(data.mod_togo_2015$dist_to_water_m)
m <- mean(train_logd, na.rm = TRUE)
s <- sd(train_logd,   na.rm = TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s
grid_df <- grid_df |> tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_togo, slope_togo, log_dist_to_water)
grid.pred_togo <- as.matrix(grid_df[, c("x", "y")]) / 1000


pred_togo <- spatial.pred.binomial.MCML(
  object            = fit.MCML.togo,
  grid.pred         = grid.pred_togo, 
  predictors        = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_togo","slope_togo","log_dist_to_water")],
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_togo <- range(pred_togo$prevalence$predictions, na.rm = TRUE)

cat(
  "Predicted prevalence (togo, 2015) ranges from",
  sprintf("%.3f", rng_togo[1]), "to", sprintf("%.3f", rng_togo[2]), "\n"
)


#Borders 
togo_0 <- st_read(paste0(togo_dir, "/gadm41_TGO_0.shp"))
togo_1 <- st_read(paste0(togo_dir, "/gadm41_TGO_1.shp"))
togo_2 <- st_read(paste0(togo_dir, "/gadm41_TGO_2.shp"))
togo_0 <- st_transform(togo_0, 32631)
togo_1 <- st_transform(togo_1, 32631)
togo_2 <- st_transform(togo_2, 32631)

grid_df$prevalence <- pred_togo$prevalence$predictions

pred_sf_togo <- st_as_sf(
  grid_df,
  coords = c("x", "y"),
  crs = 32631,      
  remove = FALSE
)


##National prevalence map 
ggplot() +
  geom_sf(data = pred_sf_togo, aes(color = prevalence), size = 0.8, alpha = 0.9) +
  geom_sf(data = togo_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent
  ) +
  coord_sf(crs = st_crs(32631)) +
  theme_minimal() +
  labs(
    title = "Predicted Onchocerciasis Prevalence",
    subtitle = "togo – National Level"
  )

## Average by district ##
pred_joined_togo2 <- st_join(pred_sf_togo, togo_2)

oncho_prev_by_district <- pred_joined_togo2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>   
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

togo_avg_oncho_admin2 <- left_join(togo_2, oncho_prev_by_district, by = "NAME_2")


## District level prevalence map ##
ggplot() +
  geom_sf(data = togo_avg_oncho_admin2, aes(fill = mean_prev), color = "black", size = 0.3) +
  geom_sf(data = togo_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  theme_minimal() +
  labs(
    title = "Onchocerciasis Prevalence by District",
    subtitle = "togo – District Level"
  )


## Uncertainty map (95% CI width) – Togo ##

grid_df$prevalence <- pred_togo$prevalence$predictions
grid_df$lower95    <- pred_togo$prevalence$quantiles[, "2.5%"]
grid_df$upper95    <- pred_togo$prevalence$quantiles[, "97.5%"]
grid_df$ci_width <- grid_df$upper95 - grid_df$lower95

pred_sf_togo <- st_as_sf(
  grid_df,
  coords = c("x", "y"),
  crs = 32631,  
  remove = FALSE
)

# --- National uncertainty map ---
ggplot() +
  geom_sf(data = pred_sf_togo, aes(color = ci_width), size = 0.8, alpha = 0.9) +
  geom_sf(data = togo_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black",
    high = "white"
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Onchocerciasis Prevalence",
    subtitle = "Togo – National Level"
  )

ggplot() +
  geom_sf(data = pred_sf_togo, aes(color = ci_width), size = 0.8, alpha = 0.9) +
  geom_sf(data = togo_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(
    name = "Prevalence (%)",
    option = "C",
    direction = 1,
    limits = c(0, 0.7),
    labels = scales::percent
  ) +
  coord_sf(crs = st_crs(32630)) +
  theme_minimal() +
  labs(
    title = "Predicted Onchocerciasis Prevalence",
    subtitle = "Togo – National Level"
  )



togo_2 <- st_read(paste0(togo_dir, "/gadm41_TGO_2.shp")) |> st_transform(32631)

pred_joined_togo2 <- st_join(pred_sf_togo, togo_2)

oncho_ci_by_district_togo <- pred_joined_togo2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |> 
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

togo_avg_oncho_ci_admin2 <- left_join(togo_2, oncho_ci_by_district_togo, by = "NAME_2")

# --- District level uncertainty map ---
ggplot() +
  geom_sf(data = togo_avg_oncho_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = togo_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 0.5),  
    labels = function(x) ifelse(x %in% c(0, 0.5), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Onchocerciasis Prevalence by District",
    subtitle = "Togo – District Level"
  )

pred_bundle_togo <- list(
  predictions = pred_togo$prevalence$predictions,
  lower95     = pred_togo$prevalence$quantiles[, "2.5%"],
  upper95     = pred_togo$prevalence$quantiles[, "97.5%"],
  ci_width    = pred_togo$prevalence$quantiles[, "97.5%"] -
    pred_togo$prevalence$quantiles[, "2.5%"],
  
  grid        = grid_df,         
  grid.pred   = grid.pred_togo, 
  model       = fit.MCML.togo,   
  rng         = rng_togo,        
  date        = Sys.Date(),
  note        = "Onchocerciasis prevalence predictions + uncertainty for Togo, 2015"
)

saveRDS(pred_bundle_togo,
        "togo_oncho_predictions_2015_bundle.rds")
bundle <- readRDS("togo_oncho_predictions_2015_bundle.rds")
str(bundle)

library(writexl)

write_xlsx(data.mod_togo_2015, "togo_oncho_data_2015.xlsx")





