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

data_st_benin <- read_excel("schistosomiasis.xlsx", sheet = "st_benin")
data_st_benin$Longitude <- as.numeric(data_st_benin$Longitude)
data_st_benin$Latitude  <- as.numeric(data_st_benin$Latitude)
data_st_benin <- subset(data_st_benin, !is.na(Longitude) & !is.na(Latitude))

benin_sf_ll <- st_as_sf(data_st_benin, coords = c("Longitude", "Latitude"), crs = 4326)

benin_sf_utm <- st_transform(benin_sf_ll, crs = 32631)
utm_mat <- st_coordinates(benin_sf_utm)
data.mod_benin <- cbind(benin_sf_utm, utm_mat)

names(data.mod_benin)[names(data.mod_benin) == "X"] <- "utm_x"
names(data.mod_benin)[names(data.mod_benin) == "Y"] <- "utm_y"

data.mod_benin_df <- cbind(st_drop_geometry(benin_sf_utm), utm_x = utm_mat[,1], utm_y = utm_mat[,2])

cols_to_remove <- c("Country", "ADMIN1_NAME", "ADMIN2_NAME", "IU_NAME", 
                    "LocationName", "LocationType", "SurveyType",
                    "Method_1", "Method_2", 
                    "Age_start", "Age_end", "geometry")

data.mod_benin <- data.mod_benin[, !names(data.mod_benin) %in% cols_to_remove]

data.mod_benin$SurveyYear <- as.numeric(data.mod_benin$SurveyYear)
data.mod_benin_2015 <- subset(data.mod_benin, SurveyYear == 2015)
data.mod_benin_2015 <- st_drop_geometry(data.mod_benin_2015)


####  BENIN INTESTINAL ####
data.mod_benin_2015_intestinal <- subset(data.mod_benin_2015, SCH_type == "Intestinal")


## DISTANCE WATER BODIES ## 
wb <- readRDS("benin_water_bundle.rds")
str(wb)   

pts <- st_as_sf(data.mod_benin_2015_intestinal, coords = c("utm_x", "utm_y"), crs = NA) 
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

data.mod_benin_2015_intestinal$dist_to_water_m   <- dvals
data.mod_benin_2015_intestinal$log_dist_to_water <- scale(log1p(dvals))


## ANNUAL MEAN TEMPERATURE ##
tmean_2015 <- rast("tmean_2015_benin.tif")
data.mod_benin_2015_intestinal_sf <- st_as_sf(data.mod_benin_2015_intestinal,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631) 
tmean_2015_utm <- project(tmean_2015, "EPSG:32631")
vals_tmean <- terra::extract(tmean_2015_utm, vect(data.mod_benin_2015_intestinal_sf))
data.mod_benin_2015_intestinal$mean_temp <- vals_tmean[,2] 
data.mod_benin_2015_intestinal <- st_drop_geometry(data.mod_benin_2015_intestinal)



## ANNUAL PRECIPITATION ##
ppt_2015 <- rast("ppt_total_2015_benin.tif")
data.mod_benin_2015_intestinal_sf <- st_as_sf(data.mod_benin_2015_intestinal,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)  
ppt_2015_utm <- project(ppt_2015, "EPSG:32631")
vals_ppt <- terra::extract(ppt_2015_utm, vect(data.mod_benin_2015_intestinal_sf))
data.mod_benin_2015_intestinal <- st_drop_geometry(data.mod_benin_2015_intestinal)



## PRECIPITATION SEASONALITY ##
prec_season <- rast("wc2.1_2.5m_bio_15.tif") 
data.mod_benin_2015_intestinal_sf <- st_as_sf(data.mod_benin_2015_intestinal,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631) 
prec_season_utm <- project(prec_season, "EPSG:32631")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_benin_2015_intestinal_sf))
data.mod_benin_2015_intestinal$precip_seasonality <- vals_prec[,2] 
data.mod_benin_2015_intestinal <- st_drop_geometry(data.mod_benin_2015_intestinal)


## ELEVATION ##
elev <- rast("wc2.1_30s_elev.tif") 
data.mod_benin_2015_intestinal_sf <- st_as_sf(data.mod_benin_2015_intestinal,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631) 
elev_utm <- project(elev, "EPSG:32631")
vals_elev <- terra::extract(elev_utm, vect(data.mod_benin_2015_intestinal_sf))
data.mod_benin_2015_intestinal$elevation <- vals_elev[,2] 
data.mod_benin_2015_intestinal <- st_drop_geometry(data.mod_benin_2015_intestinal)



## NDVI ##
ndvi_2015_benin <- rast("NDVI_2015_annual_mean_benin.tif")
ndvi_2015_benin_utm <- project(ndvi_2015_benin, "EPSG:32631")
data.mod_benin_2015_intestinal_sf <- st_as_sf(data.mod_benin_2015_intestinal,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)

vals_ndvi <- terra::extract(ndvi_2015_benin_utm, vect(data.mod_benin_2015_intestinal_sf))
data.mod_benin_2015_intestinal$ndvi_mean_benin <- vals_ndvi[,2]
data.mod_benin_2015_intestinal <- st_drop_geometry(data.mod_benin_2015_intestinal)


## LST ##
lst_2015_benin <- rast("diffLST_2016_annual_mean_benin.tif")
lst_2015_benin_utm <- project(lst_2015_benin, "EPSG:32631")
data.mod_benin_2015_intestinal_sf <- st_as_sf(data.mod_benin_2015_intestinal,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)
vals_ndvi <- terra::extract(lst_2015_benin_utm, vect(data.mod_benin_2015_intestinal_sf))
data.mod_benin_2015_intestinal$lst_mean_benin <- vals_ndvi[,2]
data.mod_benin_2015_intestinal <- st_drop_geometry(data.mod_benin_2015_intestinal)

data.mod_benin_2015_intestinal <- data.mod_benin_2015_intestinal[!is.na(data.mod_benin_2015_intestinal$lst_mean_benin), ]


#### PREDICTION BENIN INTESTINAL ####
ID.coords.benin <- create.ID.coords(
  data = data.mod_benin_2015_intestinal,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)
fit.LA.benin <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_benin +lst_mean_benin + log_dist_to_water + elevation,
  units.m = ~ Examined,
  coords = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa = 0.5,              
  start.cov.pars = 10,      
  fixed.rel.nugget = 0,     
  data = data.mod_benin_2015_intestinal,
  ID.coords = ID.coords.benin,
  family = "Binomial"
)
s <- summary(fit.LA.benin)
summary(fit.LA.benin)



## MONTE CARLO 
data.mod_benin_2015_intestinal$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_benin_2015_intestinal$dist_to_water_m))
)
ID.coords.benin <- create.ID.coords(
  data = data.mod_benin_2015_intestinal,
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
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_benin +lst_mean_benin + log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.benin,
  kappa          = 0.5,
  start.cov.pars = phi_start,  
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.benin,
  data           = data.mod_benin_2015_intestinal,
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

grid_df$mean_temp          <- extract2(tmean_2015_utm,      grid_df)
grid_df$precip_seasonality <- extract2(prec_season_utm,     grid_df)
grid_df$ndvi_mean_benin    <- extract2(ndvi_2015_benin_utm, grid_df)
grid_df$lst_mean_benin     <- extract2(lst_2015_benin_utm,  grid_df)
grid_df$elevation          <- extract2(elev_utm,            grid_df)


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

train_logd <- log1p(data.mod_benin_2015_intestinal$dist_to_water_m)
m <- mean(train_logd, na.rm = TRUE)
s <- sd(train_logd,   na.rm = TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s

grid_df <- grid_df |> tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_benin,lst_mean_benin,log_dist_to_water,elevation)

grid.pred_benin <- as.matrix(grid_df[, c("x", "y")]) / 1000


pred_benin <- spatial.pred.binomial.MCML(
  object            = fit.MCML.benin,
  grid.pred         = grid.pred_benin, 
  predictors = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_benin",
                           "lst_mean_benin","log_dist_to_water","elevation")],  
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_benin <- range(pred_benin$prevalence$predictions, na.rm = TRUE)

cat(
  "Predicted prevalence  St. Intestinal (Benin, 2015) ranges from",
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


##National prevalence map 
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
    title = "Predicted Schistosomiasis (Intestinal) Prevalence",
    subtitle = "Benin – National Level"
  )



## Average by district ##
pred_joined_benin2 <- st_join(pred_sf_benin, benin_2)

oncho_prev_by_district <- pred_joined_benin2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>   
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

benin_avg_oncho_admin2 <- left_join(benin_2, oncho_prev_by_district, by = "NAME_2")


## District level prevalence map ##
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
    title = "Schistosomiasis (Intestinal) Prevalence by District",
    subtitle = "Benin – District Level"
  )




## Uncertainty map (95% CI width) St. Intestinal Benin  ##

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

# National uncertainty map
ggplot() +
  geom_sf(data = pred_sf_benin, aes(color = ci_width), size = 0.8, alpha = 0.9) +
  geom_sf(data = benin_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black",
    high = "white",limits = c(0, 0.9),
    labels = function(x) ifelse(x %in% c(0, 0.9), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Intestinal Schistosomiasis Prevalence",
    subtitle = "Benin – National Level (2015)"
  )


pred_joined_benin2 <- st_join(pred_sf_benin, benin_2)

schisto_ci_by_district <- pred_joined_benin2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |> 
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

benin_avg_schisto_ci_admin2 <- left_join(benin_2, schisto_ci_by_district, by = "NAME_2")

# District level uncertainty map
ggplot() +
  geom_sf(data = benin_avg_schisto_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = benin_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 0.9),
    labels = function(x) ifelse(x %in% c(0, 0.9), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Intestinal Schistosomiasis Prevalence by District",
    subtitle = "Benin – District Level (2015)"
  )


pred_bundle_benin_intestinal <- list(
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
  note        = "Intestinal Schistosomiasis prevalence predictions and uncertainty for Benin, 2015"
)

saveRDS(pred_bundle_benin_intestinal,
        "benin_schisto_intestinal_predictions_2015_bundle.rds")

bundle <- readRDS("benin_st_intestinal_2015_bundle.rds.rds")
str(bundle)

library(writexl)

write_xlsx(data.mod_benin_2015_intestinal, "benin_st_data_2015.xlsx")



##### BENIN – UROGENITAL ####

data.mod_benin_2015_urogenital <- subset(data.mod_benin_2015, SCH_type == "Urogenital")

## DISTANCE WATER BODIES ## 
wb <- readRDS("benin_water_bundle.rds")

pts <- st_as_sf(data.mod_benin_2015_urogenital, coords = c("utm_x", "utm_y"), crs = NA)
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
  stop("Don't recognize the contents of benin_water_bundle.rds")
}

data.mod_benin_2015_urogenital$dist_to_water_m   <- dvals
data.mod_benin_2015_urogenital$log_dist_to_water <- scale(log1p(dvals))


## ANNUAL MEAN TEMPERATURE ##
tmean_2015 <- rast("tmean_2015_benin.tif")
data.mod_benin_2015_urogenital_sf <- st_as_sf(data.mod_benin_2015_urogenital,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)
tmean_2015_utm <- project(tmean_2015, "EPSG:32631")
vals_tmean <- terra::extract(tmean_2015_utm, vect(data.mod_benin_2015_urogenital_sf))
data.mod_benin_2015_urogenital$mean_temp <- vals_tmean[,2]
data.mod_benin_2015_urogenital <- st_drop_geometry(data.mod_benin_2015_urogenital)


## PRECIPITATION SEASONALITY ##
prec_season <- rast("wc2.1_2.5m_bio_15.tif")
data.mod_benin_2015_urogenital_sf <- st_as_sf(data.mod_benin_2015_urogenital,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)
prec_season_utm <- project(prec_season, "EPSG:32631")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_benin_2015_urogenital_sf))
data.mod_benin_2015_urogenital$precip_seasonality <- vals_prec[,2]
data.mod_benin_2015_urogenital <- st_drop_geometry(data.mod_benin_2015_urogenital)


## ELEVATION ##
elev <- rast("wc2.1_30s_elev.tif")
data.mod_benin_2015_urogenital_sf <- st_as_sf(data.mod_benin_2015_urogenital,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)
elev_utm <- project(elev, "EPSG:32631")
vals_elev <- terra::extract(elev_utm, vect(data.mod_benin_2015_urogenital_sf))
data.mod_benin_2015_urogenital$elevation <- vals_elev[,2]
data.mod_benin_2015_urogenital <- st_drop_geometry(data.mod_benin_2015_urogenital)


## NDVI ##
ndvi_2015_benin <- rast("NDVI_2015_annual_mean_benin.tif")
ndvi_2015_benin_utm <- project(ndvi_2015_benin, "EPSG:32631")
data.mod_benin_2015_urogenital_sf <- st_as_sf(data.mod_benin_2015_urogenital,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)
vals_ndvi <- terra::extract(ndvi_2015_benin_utm, vect(data.mod_benin_2015_urogenital_sf))
data.mod_benin_2015_urogenital$ndvi_mean_benin <- vals_ndvi[,2]
data.mod_benin_2015_urogenital <- st_drop_geometry(data.mod_benin_2015_urogenital)


## LST ##
lst_2015_benin <- rast("diffLST_2016_annual_mean_benin.tif")
lst_2015_benin_utm <- project(lst_2015_benin, "EPSG:32631")
data.mod_benin_2015_urogenital_sf <- st_as_sf(data.mod_benin_2015_urogenital,
                                              coords = c("utm_x", "utm_y"),
                                              crs = 32631)
vals_lst <- terra::extract(lst_2015_benin_utm, vect(data.mod_benin_2015_urogenital_sf))
data.mod_benin_2015_urogenital$lst_mean_benin <- vals_lst[,2]
data.mod_benin_2015_urogenital <- st_drop_geometry(data.mod_benin_2015_urogenital)

data.mod_benin_2015_urogenital <- data.mod_benin_2015_urogenital[!is.na(data.mod_benin_2015_urogenital$lst_mean_benin), ]


ID.coords.benin.uro <- create.ID.coords(
  data = data.mod_benin_2015_urogenital,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

fit.LA.benin.uro <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_benin + lst_mean_benin +
    log_dist_to_water + elevation,
  units.m = ~ Examined,
  coords = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa = 0.5,
  start.cov.pars = 10,
  fixed.rel.nugget = 0,
  data = data.mod_benin_2015_urogenital,
  ID.coords = ID.coords.benin.uro,
  family = "Binomial"
)

summary(fit.LA.benin.uro)



## MONTE CARLO refinement ##
data.mod_benin_2015_urogenital$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_benin_2015_urogenital$dist_to_water_m))
)

control.mcmc <- control.mcmc.MCML(
  n.sim  = 110000,
  burnin = 10000,
  thin   = 10
)

par0.benin.uro <- coef(fit.LA.benin.uro)
logphi.uro <- summary(fit.LA.benin.uro)$cov.pars["log(phi)", "Estimate"]
phi_start.uro <- exp(logphi.uro)
phi_start.uro

fit.MCML.benin.uro <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_benin + lst_mean_benin +
    log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.benin.uro,
  kappa          = 0.5,
  start.cov.pars = phi_start.uro,
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.benin.uro,
  data           = data.mod_benin_2015_urogenital,
  method         = "nlminb"
)

summary(fit.MCML.benin.uro)


#### PREDICTION BENIN UROGENITAL ####
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

grid_df$mean_temp          <- extract2(tmean_2015_utm,      grid_df)
grid_df$precip_seasonality <- extract2(prec_season_utm,     grid_df)
grid_df$ndvi_mean_benin    <- extract2(ndvi_2015_benin_utm, grid_df)
grid_df$lst_mean_benin     <- extract2(lst_2015_benin_utm,  grid_df)
grid_df$elevation          <- extract2(elev_utm,            grid_df)


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
train_logd <- log1p(data.mod_benin_2015_urogenital$dist_to_water_m)
m <- mean(train_logd, na.rm = TRUE)
s <- sd(train_logd,   na.rm = TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s

grid_df <- grid_df |> tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_benin,lst_mean_benin,log_dist_to_water,elevation)

grid.pred_benin <- as.matrix(grid_df[, c("x", "y")]) / 1000


pred_benin <- spatial.pred.binomial.MCML(
  object            = fit.MCML.benin.uro,
  grid.pred         = grid.pred_benin,  
  predictors = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_benin",
                           "lst_mean_benin","log_dist_to_water","elevation")],  
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_benin <- range(pred_benin$prevalence$predictions, na.rm = TRUE)

cat(
  "Predicted prevalence  St. urogenital (Benin, 2015) ranges from",
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


##National prevalence map 
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
    title = "Predicted Schistosomiasis (urogenital) Prevalence",
    subtitle = "Benin – National Level"
  )



## Average by district ##
pred_joined_benin2 <- st_join(pred_sf_benin, benin_2)

oncho_prev_by_district <- pred_joined_benin2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>   
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

benin_avg_oncho_admin2 <- left_join(benin_2, oncho_prev_by_district, by = "NAME_2")

## District level prevalence map ##
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
    title = "Schistosomiasis (urogenital) Prevalence by District",
    subtitle = "Benin – District Level"
  )


## Uncertainty map (95% CI width) St Urogenital Benin 2015 ##

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

# National uncertainty map
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
    title = "Uncertainty of Predicted Schistosomiasis (urogenital) Prevalence",
    subtitle = "Benin – National Level (2015)"
  )

pred_joined_benin2 <- st_join(pred_sf_benin, benin_2)

schisto_ci_by_district <- pred_joined_benin2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |> 
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

benin_avg_schisto_ci_admin2 <- left_join(benin_2, schisto_ci_by_district, by = "NAME_2")

# District level uncertainty map
ggplot() +
  geom_sf(data = benin_avg_schisto_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = benin_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Schistosomiasis (urogenital) Prevalence by District",
    subtitle = "Benin – District Level (2015)"
  )


pred_bundle_benin_uro <- list(
  predictions = pred_benin$prevalence$predictions,     
  lower95     = pred_benin$prevalence$quantiles[, "2.5%"],  
  upper95     = pred_benin$prevalence$quantiles[, "97.5%"], 
  ci_width    = pred_benin$prevalence$quantiles[, "97.5%"] -
    pred_benin$prevalence$quantiles[, "2.5%"],  
  grid        = grid_df,          
  grid.pred   = grid.pred_benin,  
  model       = fit.MCML.benin.uro,   
  rng         = rng_benin,        
  date        = Sys.Date(),
  note        = "Schistosomiasis (urogenital) prevalence predictions and uncertainty for Benin, 2015"
)


saveRDS(pred_bundle_benin,
        "benin_st_urogenital_2015_bundle.rds")

bundle <- readRDS("benin_st_urogenital_2015_bundle.rds")
str(bundle)

library(writexl)

write_xlsx(data.mod_benin_2015_urogenital, "benin_st_data_2015_uro.xlsx")




#### CIV ####

library(sf)
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(sp)
library(dplyr)
library(readxl)
library(PrevMap)
library(writexl)


`%||%` <- function(x, y) if (is.null(x)) y else x


data_st_civ <- read_excel("schistosomiasis.xlsx", sheet = "st_civ")

data_st_civ$Longitude <- as.numeric(data_st_civ$Longitude)
data_st_civ$Latitude  <- as.numeric(data_st_civ$Latitude)
data_st_civ <- subset(data_st_civ, !is.na(Longitude) & !is.na(Latitude))

civ_sf_ll  <- st_as_sf(data_st_civ, coords = c("Longitude", "Latitude"), crs = 4326)
civ_sf_utm <- st_transform(civ_sf_ll, 32630)
utm_mat    <- st_coordinates(civ_sf_utm)

data.mod_civ <- cbind(civ_sf_utm, utm_mat)
names(data.mod_civ)[names(data.mod_civ) == "X"] <- "utm_x"
names(data.mod_civ)[names(data.mod_civ) == "Y"] <- "utm_y"

cols_to_remove <- c("Country", "ADMIN1_NAME", "ADMIN2_NAME", "IU_NAME",
                    "LocationName", "LocationType", "SurveyType",
                    "Method_1", "Method_2", "Age_start", "Age_end", "geometry")
data.mod_civ <- data.mod_civ[, !names(data.mod_civ) %in% cols_to_remove]

data.mod_civ$SurveyYear <- as.numeric(data.mod_civ$SurveyYear)
data.mod_civ_2018 <- subset(data.mod_civ, SurveyYear == 2018)
data.mod_civ_2018 <- st_drop_geometry(data.mod_civ_2018)


#### CIV - INTESTINAL ####
data.mod_civ_2018_intestinal <- subset(data.mod_civ_2018, SCH_type == "Intestinal")

wb_civ <- readRDS("civ_water_bundle.rds")
pts <- st_as_sf(data.mod_civ_2018_intestinal, coords = c("utm_x", "utm_y"), crs = NA)
st_crs(pts) <- 32630

if (inherits(wb_civ, "SpatRaster") || inherits(wb_civ$raster, "SpatRaster") || inherits(wb_civ$climatology, "SpatRaster")) {
  water <- wb_civ$climatology %||% wb_civ$raster %||% wb_civ
  if (!terra::same.crs(water, st_crs(pts)$wkt)) {
    water <- terra::project(water, st_crs(pts)$wkt)
  }
  water_to <- classify(water, rbind(c(-Inf, Inf, 1)))
  water_to[is.na(water)] <- NA
  dist_r  <- terra::distance(water_to)
  pts_v   <- vect(pts)
  dvals   <- terra::extract(dist_r, pts_v)[,1]
} else if (inherits(wb_civ, "sf") || inherits(wb_civ$rivers, "sf") || inherits(wb_civ$water, "sf") || inherits(wb_civ$lines, "sf")) {
  rivers <- wb_civ$rivers %||% wb_civ$water %||% wb_civ$lines %||% wb_civ
  rivers <- st_as_sf(rivers)
  if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set it correctly (EPSG:32630).")
  if (st_crs(rivers) != st_crs(pts)) rivers <- st_transform(rivers, st_crs(pts))
  dmat  <- st_distance(pts, rivers)
  dvals <- apply(dmat, 1, min) |> as.numeric()
} else {
  stop("Unrecognized water bundle for CI (not SpatRaster or sf).")
}

data.mod_civ_2018_intestinal$dist_to_water_m   <- dvals
data.mod_civ_2018_intestinal$log_dist_to_water <- scale(log1p(dvals))

# Mean Temperature
tmean_2018_civ <- rast("tmean_2018_civ.tif")  
data.mod_civ_2018_intestinal_sf <- st_as_sf(data.mod_civ_2018_intestinal, coords = c("utm_x","utm_y"), crs = 32630)
tmean_utm <- project(tmean_2018_civ, "EPSG:32630")
vals_tmean <- terra::extract(tmean_utm, vect(data.mod_civ_2018_intestinal_sf))
data.mod_civ_2018_intestinal$mean_temp <- vals_tmean[,2]
data.mod_civ_2018_intestinal <- st_drop_geometry(data.mod_civ_2018_intestinal)



# Precipitation Seasonality 
prec_season <- rast("wc2.1_2.5m_bio_15.tif")
data.mod_civ_2018_intestinal_sf <- st_as_sf(data.mod_civ_2018_intestinal, coords = c("utm_x","utm_y"), crs = 32630)
prec_season_utm <- project(prec_season, "EPSG:32630")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_civ_2018_intestinal_sf))
data.mod_civ_2018_intestinal$precip_seasonality <- vals_prec[,2]
data.mod_civ_2018_intestinal <- st_drop_geometry(data.mod_civ_2018_intestinal)

# Elevation 
elev <- rast("wc2.1_30s_elev.tif")
data.mod_civ_2018_intestinal_sf <- st_as_sf(data.mod_civ_2018_intestinal, coords = c("utm_x","utm_y"), crs = 32630)
elev_utm <- project(elev, "EPSG:32630")
vals_elev <- terra::extract(elev_utm, vect(data.mod_civ_2018_intestinal_sf))
data.mod_civ_2018_intestinal$elevation <- vals_elev[,2]
data.mod_civ_2018_intestinal <- st_drop_geometry(data.mod_civ_2018_intestinal)

# NDVI 
ndvi_2018_civ <- rast("NDVI_2018_annual_mean_civ.tif")
ndvi_utm <- project(ndvi_2018_civ, "EPSG:32630")
data.mod_civ_2018_intestinal_sf <- st_as_sf(data.mod_civ_2018_intestinal, coords = c("utm_x","utm_y"), crs = 32630)
vals_ndvi <- terra::extract(ndvi_utm, vect(data.mod_civ_2018_intestinal_sf))
data.mod_civ_2018_intestinal$ndvi_mean_civ <- vals_ndvi[,2]
data.mod_civ_2018_intestinal <- st_drop_geometry(data.mod_civ_2018_intestinal)

# LST 

lst_2018_civ <- rast("diffLST_2018_annual_mean_civ.tif")
lst_utm <- project(lst_2018_civ, "EPSG:32630")
data.mod_civ_2018_intestinal_sf <- st_as_sf(data.mod_civ_2018_intestinal, coords = c("utm_x","utm_y"), crs = 32630)
vals_lst <- terra::extract(lst_utm, vect(data.mod_civ_2018_intestinal_sf))
data.mod_civ_2018_intestinal$lst_mean_civ <- vals_lst[,2]
data.mod_civ_2018_intestinal <- st_drop_geometry(data.mod_civ_2018_intestinal)
data.mod_civ_2018_intestinal <- data.mod_civ_2018_intestinal[!is.na(data.mod_civ_2018_intestinal$lst_mean_civ), ]

ID.coords.civ <- create.ID.coords(
  data = data.mod_civ_2018_intestinal,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

# Laplace Approximation 
fit.LA.civ <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_civ + lst_mean_civ +
    log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa          = 0.5,
  start.cov.pars = 30,      
  fixed.rel.nugget = 0,
  data           = data.mod_civ_2018_intestinal,
  ID.coords      = ID.coords.civ,
  family         = "Binomial"
)

s.civ <- summary(fit.LA.civ)
print(s.civ)



# MCML
data.mod_civ_2018_intestinal$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_civ_2018_intestinal$dist_to_water_m))
)

ID.coords.civ <- create.ID.coords(
  data = data.mod_civ_2018_intestinal,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

control.mcmc <- control.mcmc.MCML(
  n.sim  = 110000,
  burnin = 10000,
  thin   = 10
)

par0.civ   <- coef(fit.LA.civ)
logphi_civ <- s.civ$cov.pars["log(phi)", "Estimate"]
phi_start_civ <- exp(logphi_civ)
phi_start_civ
fit.MCML.civ <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_civ + lst_mean_civ +
    log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.civ,
  kappa          = 0.5,
  start.cov.pars = phi_start_civ,
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.civ,
  data           = data.mod_civ_2018_intestinal,
  method         = "nlminb"
)

print(summary(fit.MCML.civ))


civ_0 <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))
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
  v <- if (inherits(pts, "sf")) terra::vect(pts) else terra::vect(pts, geom = c("x", "y"), crs = crs)
  terra::extract(r, v)[[2]]
}

grid_df$mean_temp          <- extract2(tmean_utm,          grid_df)
grid_df$precip_seasonality <- extract2(prec_season_utm,    grid_df)
grid_df$ndvi_mean_civ      <- extract2(ndvi_utm,           grid_df)
grid_df$lst_mean_civ       <- extract2(lst_utm,            grid_df)
grid_df$elevation          <- extract2(elev_utm,           grid_df)

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
grid_df$dist_to_water_m <- get_dist_to_water(grid_df, wb_civ)

train_logd <- log1p(data.mod_civ_2018_intestinal$dist_to_water_m)
m <- mean(train_logd, na.rm = TRUE); s <- sd(train_logd, na.rm = TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s

grid_df <- grid_df |>
  tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_civ, lst_mean_civ,
                 log_dist_to_water, elevation)

grid.pred_civ <- as.matrix(grid_df[, c("x","y")]) / 1000


#### PREDICTION CIV INTESTINAL ####
pred_civ <- spatial.pred.binomial.MCML(
  object            = fit.MCML.civ,
  grid.pred         = grid.pred_civ,
  predictors        = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_civ",
                                  "lst_mean_civ","log_dist_to_water","elevation")],
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_civ <- range(pred_civ$prevalence$predictions, na.rm = TRUE)
cat("Predicted prevalence St. Intestinal (Côte d’Ivoire, 2018) ranges from",
    sprintf("%.3f", rng_civ[1]), "to", sprintf("%.3f", rng_civ[2]), "\n")
civ_1 <- st_read(file.path(civ_dir, "gadm41_CIV_1.shp")) |> st_transform(32630)
civ_2 <- st_read(file.path(civ_dir, "gadm41_CIV_2.shp")) |> st_transform(32630)

grid_df$prevalence <- pred_civ$prevalence$predictions
pred_sf_civ <- st_as_sf(grid_df, coords = c("x","y"), crs = 32630, remove = FALSE)

## National prevalence map
ggplot() +
  geom_sf(data = pred_sf_civ, aes(color = prevalence), size = 0.8, alpha = 0.9) +
  geom_sf(data = civ_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(name = "Prevalence (%)", option = "C", direction = 1,
                        limits = c(0, 1), labels = scales::percent) +
  coord_sf(crs = st_crs(32630)) +
  theme_minimal() +
  labs(title = "Predicted Schistosomiasis (Intestinal) Prevalence",
       subtitle = "Côte d’Ivoire – National Level")

## Average by district 
pred_joined_civ2 <- st_join(pred_sf_civ, civ_2)
civ_prev_by_district <- pred_joined_civ2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

civ_avg_prev_admin2 <- left_join(civ_2, civ_prev_by_district, by = "NAME_2")

ggplot() +
  geom_sf(data = civ_avg_prev_admin2, aes(fill = mean_prev), color = "black", size = 0.3) +
  geom_sf(data = civ_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(name = "Prevalence (%)", option = "C", direction = 1,
                       limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  labs(title = "Schistosomiasis (Intestinal) Prevalence by District",
       subtitle = "Côte d’Ivoire – District Level")



## UNCERTAINTY MAP (95% CI width) – CIV INTESTINAL 2018 ##

grid_df$prevalence <- pred_civ$prevalence$predictions
grid_df$lower95    <- pred_civ$prevalence$quantiles[, "2.5%"]
grid_df$upper95    <- pred_civ$prevalence$quantiles[, "97.5%"]
grid_df$ci_width <- grid_df$upper95 - grid_df$lower95

pred_sf_civ <- st_as_sf(
  grid_df,
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
    high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Schistosomiasis (Intestinal) Prevalence",
    subtitle = "Côte d’Ivoire – National Level (2018)"
  )


pred_joined_civ2 <- st_join(pred_sf_civ, civ_2)

schisto_ci_by_district <- pred_joined_civ2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

civ_avg_schisto_ci_admin2 <- left_join(civ_2, schisto_ci_by_district, by = "NAME_2")

# District level uncertainty map
ggplot() +
  geom_sf(data = civ_avg_schisto_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = civ_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Schistosomiasis (Intestinal) Prevalence by District",
    subtitle = "Côte d’Ivoire – District Level (2018)"
  )


pred_bundle_civ_intestinal <- list(
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
  note        = "Schistosomiasis (Intestinal) prevalence predictions and uncertainty for Côte d’Ivoire, 2018"
)

saveRDS(pred_bundle_civ_intestinal,
        "civ_schisto_intestinal_predictions_2018_bundle.rds")



write_xlsx(data.mod_civ_2018_intestinal, "civ_st_data_2018.xlsx")




##### CIV - UROGENITAL ####

library(sf)
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(sp)
library(dplyr)
library(readxl)
library(PrevMap)
library(writexl)


data_st_civ <- read_excel("schistosomiasis.xlsx", sheet = "st_civ")
data_st_civ$Longitude <- as.numeric(data_st_civ$Longitude)
data_st_civ$Latitude  <- as.numeric(data_st_civ$Latitude)
data_st_civ <- subset(data_st_civ, !is.na(Longitude) & !is.na(Latitude))

civ_sf_ll  <- st_as_sf(data_st_civ, coords = c("Longitude", "Latitude"), crs = 4326)
civ_sf_utm <- st_transform(civ_sf_ll, 32630)
utm_mat    <- st_coordinates(civ_sf_utm)

data.mod_civ <- cbind(civ_sf_utm, utm_mat)
names(data.mod_civ)[names(data.mod_civ) == "X"] <- "utm_x"
names(data.mod_civ)[names(data.mod_civ) == "Y"] <- "utm_y"

cols_to_remove <- c("Country", "ADMIN1_NAME", "ADMIN2_NAME", "IU_NAME",
                    "LocationName", "LocationType", "SurveyType",
                    "Method_1", "Method_2", "Age_start", "Age_end", "geometry")
data.mod_civ <- data.mod_civ[, !names(data.mod_civ) %in% cols_to_remove]

data.mod_civ$SurveyYear <- as.numeric(data.mod_civ$SurveyYear)
data.mod_civ_2018 <- subset(data.mod_civ, SurveyYear == 2018)
data.mod_civ_2018 <- st_drop_geometry(data.mod_civ_2018)


# --- UROGENITAL ---
data.mod_civ_2018_uro <- subset(data.mod_civ_2018, SCH_type == "Urogenital")

wb_civ <- readRDS("civ_water_bundle.rds")

pts <- st_as_sf(data.mod_civ_2018_uro, coords = c("utm_x","utm_y"), crs = NA)
st_crs(pts) <- 32630

if (inherits(wb_civ, "SpatRaster") || inherits(wb_civ$raster, "SpatRaster") || inherits(wb_civ$climatology, "SpatRaster")) {
  water <- wb_civ$climatology %||% wb_civ$raster %||% wb_civ
  if (!terra::same.crs(water, st_crs(pts)$wkt)) water <- terra::project(water, st_crs(pts)$wkt)
  water_to <- classify(water, rbind(c(-Inf, Inf, 1))); water_to[is.na(water)] <- NA
  dist_r <- terra::distance(water_to)
  dvals  <- terra::extract(dist_r, vect(pts))[,1]
} else if (inherits(wb_civ, "sf") || inherits(wb_civ$rivers, "sf") || inherits(wb_civ$water, "sf") || inherits(wb_civ$lines, "sf")) {
  rivers <- wb_civ$rivers %||% wb_civ$water %||% wb_civ$lines %||% wb_civ
  if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set EPSG:32630.")
  if (st_crs(rivers) != st_crs(pts)) rivers <- st_transform(rivers, st_crs(pts))
  dvals <- apply(st_distance(pts, rivers), 1, min) |> as.numeric()
} else stop("Unrecognized water bundle (not SpatRaster or sf).")

data.mod_civ_2018_uro$dist_to_water_m   <- dvals
data.mod_civ_2018_uro$log_dist_to_water <- scale(log1p(dvals))


# Mean Temperature
tmean_2018_civ <- rast("tmean_2018_civ.tif")
data.mod_civ_2018_uro_sf <- st_as_sf(data.mod_civ_2018_uro, coords = c("utm_x","utm_y"), crs = 32630)
tmean_utm <- project(tmean_2018_civ, "EPSG:32630")
vals_tmean <- terra::extract(tmean_utm, vect(data.mod_civ_2018_uro_sf))
data.mod_civ_2018_uro$mean_temp <- vals_tmean[,2]
data.mod_civ_2018_uro <- st_drop_geometry(data.mod_civ_2018_uro)


# Precipitation Seasonality
prec_season <- rast("wc2.1_2.5m_bio_15.tif")
data.mod_civ_2018_uro_sf <- st_as_sf(data.mod_civ_2018_uro, coords = c("utm_x","utm_y"), crs = 32630)
prec_season_utm <- project(prec_season, "EPSG:32630")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_civ_2018_uro_sf))
data.mod_civ_2018_uro$precip_seasonality <- vals_prec[,2]
data.mod_civ_2018_uro <- st_drop_geometry(data.mod_civ_2018_uro)

# Elevation
elev <- rast("wc2.1_30s_elev.tif")
data.mod_civ_2018_uro_sf <- st_as_sf(data.mod_civ_2018_uro, coords = c("utm_x","utm_y"), crs = 32630)
elev_utm <- project(elev, "EPSG:32630")
vals_elev <- terra::extract(elev_utm, vect(data.mod_civ_2018_uro_sf))
data.mod_civ_2018_uro$elevation <- vals_elev[,2]
data.mod_civ_2018_uro <- st_drop_geometry(data.mod_civ_2018_uro)

# NDVI 
ndvi_2018_civ <- rast("NDVI_2018_annual_mean_civ.tif")
ndvi_utm <- project(ndvi_2018_civ, "EPSG:32630")
data.mod_civ_2018_uro_sf <- st_as_sf(data.mod_civ_2018_uro, coords = c("utm_x","utm_y"), crs = 32630)
vals_ndvi <- terra::extract(ndvi_utm, vect(data.mod_civ_2018_uro_sf))
data.mod_civ_2018_uro$ndvi_mean_civ <- vals_ndvi[,2]
data.mod_civ_2018_uro <- st_drop_geometry(data.mod_civ_2018_uro)

# LST diff
lst_2018_civ <- rast("diffLST_2018_annual_mean_civ.tif")
lst_utm <- project(lst_2018_civ, "EPSG:32630")
data.mod_civ_2018_uro_sf <- st_as_sf(data.mod_civ_2018_uro, coords = c("utm_x","utm_y"), crs = 32630)
vals_lst <- terra::extract(lst_utm, vect(data.mod_civ_2018_uro_sf))
data.mod_civ_2018_uro$lst_mean_civ <- vals_lst[,2]
data.mod_civ_2018_uro <- st_drop_geometry(data.mod_civ_2018_uro)
data.mod_civ_2018_uro <- data.mod_civ_2018_uro[!is.na(data.mod_civ_2018_uro$lst_mean_civ), ]

ID.coords.civ <- create.ID.coords(
  data = data.mod_civ_2018_uro,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)

fit.LA.civ <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_civ + lst_mean_civ +
    log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  kappa          = 0.5,
  start.cov.pars = 40,
  fixed.rel.nugget = 0,
  data           = data.mod_civ_2018_uro,
  ID.coords      = ID.coords.civ,
  family         = "Binomial"
)
s.civ <- summary(fit.LA.civ)
print(s.civ)



# MCML
data.mod_civ_2018_uro$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_civ_2018_uro$dist_to_water_m))
)
ID.coords.civ <- create.ID.coords(
  data = data.mod_civ_2018_uro,
  ~ I(utm_x / 1000) + I(utm_y / 1000)
)
control.mcmc <- control.mcmc.MCML(n.sim=110000, burnin=10000, thin=10)
par0.civ   <- coef(fit.LA.civ)
logphi_civ <- s.civ$cov.pars["log(phi)", "Estimate"]
phi_start_civ <- exp(logphi_civ)

fit.MCML.civ <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_civ + lst_mean_civ +
    log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x / 1000) + I(utm_y / 1000),
  ID.coords      = ID.coords.civ,
  kappa          = 0.5,
  start.cov.pars = phi_start_civ,
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.civ,
  data           = data.mod_civ_2018_uro,
  method         = "nlminb"
)
print(summary(fit.MCML.civ))
civ_0 <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp")) |> st_transform(32630)

bb   <- st_bbox(civ_0)
nxny <- 250
xseq <- seq(bb["xmin"], bb["xmax"], length.out = nxny)
yseq <- seq(bb["ymin"], bb["ymax"], length.out = nxny)
grid_df <- expand.grid(x = xseq, y = yseq)

grid_sf <- st_as_sf(grid_df, coords = c("x","y"), crs = 32630)
inside  <- st_within(grid_sf, civ_0, sparse = FALSE)[,1]
grid_sf <- grid_sf[inside, ]
grid_df <- cbind(st_coordinates(grid_sf) |> as.data.frame())
names(grid_df) <- c("x","y")

extract2 <- function(r, pts, crs = "EPSG:32630") {
  v <- if (inherits(pts,"sf")) terra::vect(pts) else terra::vect(pts, geom=c("x","y"), crs=crs)
  terra::extract(r, v)[[2]]
}

grid_df$mean_temp          <- extract2(tmean_utm,          grid_df)
grid_df$precip_seasonality <- extract2(prec_season_utm,    grid_df)
grid_df$ndvi_mean_civ      <- extract2(ndvi_utm,           grid_df)
grid_df$lst_mean_civ       <- extract2(lst_utm,            grid_df)
grid_df$elevation          <- extract2(elev_utm,           grid_df)

get_dist_to_water <- function(grid_xy, wb) {
  if (inherits(wb, "SpatRaster") || inherits(wb$raster, "SpatRaster") || inherits(wb$climatology, "SpatRaster")) {
    water <- wb$climatology %||% wb$raster %||% wb
    if (!terra::same.crs(water, "EPSG:32630")) water <- terra::project(water, "EPSG:32630")
    water_to <- classify(water, rbind(c(-Inf, Inf, 1))); water_to[is.na(water)] <- NA
    dist_r <- terra::distance(water_to)
    terra::extract(dist_r, vect(grid_xy))[[1]]
  } else if (inherits(wb,"sf") || inherits(wb$rivers,"sf") || inherits(wb$water,"sf") || inherits(wb$lines,"sf")) {
    rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
    if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set EPSG:32630.")
    if (st_crs(rivers)$epsg != 32630) rivers <- st_transform(rivers, 32630)
    apply(st_distance(st_as_sf(grid_xy, coords=c("x","y"), crs=32630), rivers), 1, min) |> as.numeric()
  } else stop("Unrecognized water bundle.")
}
grid_df$dist_to_water_m <- get_dist_to_water(grid_df, wb_civ)

train_logd <- log1p(data.mod_civ_2018_uro$dist_to_water_m)
m <- mean(train_logd, na.rm=TRUE); s <- sd(train_logd, na.rm=TRUE)
grid_df$log_dist_to_water <- (log1p(grid_df$dist_to_water_m) - m) / s

grid_df <- grid_df |>
  tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_civ, lst_mean_civ,
                 log_dist_to_water, elevation)

grid.pred_civ <- as.matrix(grid_df[, c("x","y")]) / 1000




#### PREDICTION CIV UROGENITAL ####
pred_civ <- spatial.pred.binomial.MCML(
  object            = fit.MCML.civ,
  grid.pred         = grid.pred_civ,
  predictors        = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_civ",
                                  "lst_mean_civ","log_dist_to_water","elevation")],
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_civ <- range(pred_civ$prevalence$predictions, na.rm = TRUE)
cat("Predicted prevalence St. Urogenital (Côte d’Ivoire, 2018) ranges from",
    sprintf("%.3f", rng_civ[1]), "to", sprintf("%.3f", rng_civ[2]), "\n")

civ_1 <- st_read(file.path(civ_dir, "gadm41_CIV_1.shp")) |> st_transform(32630)
civ_2 <- st_read(file.path(civ_dir, "gadm41_CIV_2.shp")) |> st_transform(32630)

grid_df$prevalence <- pred_civ$prevalence$predictions
pred_sf_civ <- st_as_sf(grid_df, coords = c("x","y"), crs = 32630, remove = FALSE)

# National map
ggplot() +
  geom_sf(data = pred_sf_civ, aes(color = prevalence), size = 0.8, alpha = 0.9) +
  geom_sf(data = civ_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(name = "Prevalence (%)", option = "C", direction = 1,
                        limits = c(0, 1), labels = scales::percent) +
  coord_sf(crs = st_crs(32630)) +
  theme_minimal() +
  labs(title = "Predicted Schistosomiasis (Urogenital) Prevalence",
       subtitle = "Côte d’Ivoire – National Level")


# District averages 
pred_joined_civ2 <- st_join(pred_sf_civ, civ_2)
civ_prev_by_district <- pred_joined_civ2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

civ_avg_prev_admin2 <- left_join(civ_2, civ_prev_by_district, by = "NAME_2")

ggplot() +
  geom_sf(data = civ_avg_prev_admin2, aes(fill = mean_prev), color = "black", size = 0.3) +
  geom_sf(data = civ_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(name = "Prevalence (%)", option = "C", direction = 1,
                       limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  labs(title = "Schistosomiasis (Urogenital) Prevalence by District",
       subtitle = "Côte d’Ivoire – District Level")


## UNCERTAINTY MAP (95% CI width) – CIV UROGENITAL ##

grid_df$prevalence <- pred_civ$prevalence$predictions
grid_df$lower95    <- pred_civ$prevalence$quantiles[, "2.5%"]
grid_df$upper95    <- pred_civ$prevalence$quantiles[, "97.5%"]

grid_df$ci_width <- grid_df$upper95 - grid_df$lower95

pred_sf_civ <- st_as_sf(
  grid_df,
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
    high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Schistosomiasis (Urogenital) Prevalence",
    subtitle = "Côte d’Ivoire – National Level (2018)"
  )



pred_joined_civ2 <- st_join(pred_sf_civ, civ_2)

uro_ci_by_district <- pred_joined_civ2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

civ_avg_uro_ci_admin2 <- left_join(civ_2, uro_ci_by_district, by = "NAME_2")


# District level uncertainty map
ggplot() +
  geom_sf(data = civ_avg_uro_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = civ_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0, 1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Schistosomiasis (Urogenital) Prevalence by District",
    subtitle = "Côte d’Ivoire – District Level (2018)"
  )


pred_bundle_civ_uro <- list(
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
  note        = "Schistosomiasis (Urogenital) prevalence predictions and uncertainty for Côte d’Ivoire, 2018"
)

saveRDS(pred_bundle_civ_uro,
        "civ_schisto_urogenital_predictions_2018_bundle.rds")

bundle <- readRDS("civ_schisto_urogenital_predictions_2018_bundle.rds")
str(bundle)

write_xlsx(data.mod_civ_2018_uro, "civ_st_data_2018_urogenital.xlsx")



#### GHANA ####

library(sf)
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(sp)
library(dplyr)
library(readxl)
library(PrevMap)
library(writexl)


`%||%` <- function(x, y) if (is.null(x)) y else x

##### GHANA - UROGENITAL ####

data_st_gha <- read_excel("schistosomiasis.xlsx", sheet = "st_ghana")
data_st_gha$Longitude <- as.numeric(data_st_gha$Longitude)
data_st_gha$Latitude  <- as.numeric(data_st_gha$Latitude)
data_st_gha <- subset(data_st_gha, !is.na(Longitude) & !is.na(Latitude))

gha_sf_ll  <- st_as_sf(data_st_gha, coords = c("Longitude", "Latitude"), crs = 4326)
gha_sf_utm <- st_transform(gha_sf_ll, 32630)
utm_mat    <- st_coordinates(gha_sf_utm)

data.mod_gha <- cbind(gha_sf_utm, utm_mat)
names(data.mod_gha)[names(data.mod_gha) == "X"] <- "utm_x"
names(data.mod_gha)[names(data.mod_gha) == "Y"] <- "utm_y"

cols_to_remove <- c("Country","ADMIN1_NAME","ADMIN2_NAME","IU_NAME",
                    "LocationName","LocationType","SurveyType",
                    "Method_1","Method_2","Age_start","Age_end","geometry")
data.mod_gha <- data.mod_gha[, !names(data.mod_gha) %in% cols_to_remove]

data.mod_gha$SurveyYear <- as.numeric(data.mod_gha$SurveyYear)
data.mod_gha_2021 <- subset(data.mod_gha, SurveyYear == 2021)
data.mod_gha_2021 <- st_drop_geometry(data.mod_gha_2021)

data.mod_gha_2021_urogenital <- subset(data.mod_gha_2021, SCH_type == "Urogenital")

wb_gha <- readRDS("ghana_water_bundle.rds") 

pts <- st_as_sf(data.mod_gha_2021_urogenital, coords = c("utm_x","utm_y"), crs = NA)
st_crs(pts) <- 32630

if (inherits(wb_gha, "SpatRaster") || inherits(wb_gha$raster, "SpatRaster") || inherits(wb_gha$climatology, "SpatRaster")) {
  water <- wb_gha$climatology %||% wb_gha$raster %||% wb_gha
  if (!terra::same.crs(water, st_crs(pts)$wkt)) water <- terra::project(water, st_crs(pts)$wkt)
  water_to <- classify(water, rbind(c(-Inf, Inf, 1))); water_to[is.na(water)] <- NA
  dist_r <- terra::distance(water_to)
  dvals  <- terra::extract(dist_r, vect(pts))[,1]
} else if (inherits(wb_gha, "sf") || inherits(wb_gha$rivers, "sf") || inherits(wb_gha$water, "sf") || inherits(wb_gha$lines, "sf")) {
  rivers <- wb_gha$rivers %||% wb_gha$water %||% wb_gha$lines %||% wb_gha
  rivers <- st_as_sf(rivers)
  if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set EPSG:32630.")
  if (st_crs(rivers) != st_crs(pts)) rivers <- st_transform(rivers, st_crs(pts))
  dvals <- apply(st_distance(pts, rivers), 1, min) |> as.numeric()
} else stop("Unrecognized water bundle for Ghana (not SpatRaster or sf).")

data.mod_gha_2021_urogenital$dist_to_water_m   <- dvals
data.mod_gha_2021_urogenital$log_dist_to_water <- as.numeric(scale(log1p(dvals)))


# Mean Temperature
tmean_2021_gha <- rast("tmean_2021_ghana.tif")
data.mod_gha_2021_urogenital_sf <- st_as_sf(data.mod_gha_2021_urogenital, coords = c("utm_x","utm_y"), crs = 32630)
tmean_utm <- project(tmean_2021_gha, "EPSG:32630")
vals_tmean <- terra::extract(tmean_utm, vect(data.mod_gha_2021_urogenital_sf))
data.mod_gha_2021_urogenital$mean_temp <- vals_tmean[,2]
data.mod_gha_2021_urogenital <- st_drop_geometry(data.mod_gha_2021_urogenital)


# Precipitation Seasonality 
prec_season <- rast("wc2.1_2.5m_bio_15.tif")
data.mod_gha_2021_urogenital_sf <- st_as_sf(data.mod_gha_2021_urogenital, coords = c("utm_x","utm_y"), crs = 32630)
prec_season_utm <- project(prec_season, "EPSG:32630")
vals_prec <- terra::extract(prec_season_utm, vect(data.mod_gha_2021_urogenital_sf))
data.mod_gha_2021_urogenital$precip_seasonality <- vals_prec[,2]
data.mod_gha_2021_urogenital <- st_drop_geometry(data.mod_gha_2021_urogenital)


# Elevation
elev <- rast("wc2.1_30s_elev.tif")
data.mod_gha_2021_urogenital_sf <- st_as_sf(data.mod_gha_2021_urogenital, coords = c("utm_x","utm_y"), crs = 32630)
elev_utm <- project(elev, "EPSG:32630")
vals_elev <- terra::extract(elev_utm, vect(data.mod_gha_2021_urogenital_sf))
data.mod_gha_2021_urogenital$elevation <- vals_elev[,2]
data.mod_gha_2021_urogenital <- st_drop_geometry(data.mod_gha_2021_urogenital)


# NDVI ¡
ndvi_2021_gha <- rast("NDVI_2021_annual_mean_gha.tif")
ndvi_utm <- project(ndvi_2021_gha, "EPSG:32630")
data.mod_gha_2021_urogenital_sf <- st_as_sf(data.mod_gha_2021_urogenital, coords = c("utm_x","utm_y"), crs = 32630)
vals_ndvi <- terra::extract(ndvi_utm, vect(data.mod_gha_2021_urogenital_sf))
data.mod_gha_2021_urogenital$ndvi_mean_gha <- vals_ndvi[,2]
data.mod_gha_2021_urogenital <- st_drop_geometry(data.mod_gha_2021_urogenital)



# LST diff
lst_2021_gha <- rast("diffLST_2021_annual_mean_gha.tif")
lst_utm <- project(lst_2021_gha, "EPSG:32630")
data.mod_gha_2021_urogenital_sf <- st_as_sf(data.mod_gha_2021_urogenital, coords = c("utm_x","utm_y"), crs = 32630)
vals_lst <- terra::extract(lst_utm, vect(data.mod_gha_2021_urogenital_sf))
data.mod_gha_2021_urogenital$lst_mean_gha <- vals_lst[,2]
data.mod_gha_2021_urogenital <- st_drop_geometry(data.mod_gha_2021_urogenital)

data.mod_gha_2021_urogenital <- tidyr::drop_na(
  data.mod_gha_2021_urogenital,
  mean_temp, precip_seasonality, ndvi_mean_gha, lst_mean_gha, elevation
)

ID.coords.gha <- create.ID.coords(
  data = data.mod_gha_2021_urogenital,
  ~ I(utm_x/1000) + I(utm_y/1000)
)

fit.LA.gha <- glgm.LA(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_gha + lst_mean_gha +
    log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x/1000) + I(utm_y/1000),
  kappa          = 0.5,
  start.cov.pars = 10,
  fixed.rel.nugget = 0,
  data           = data.mod_gha_2021_urogenital,
  ID.coords      = ID.coords.gha,
  family         = "Binomial"
)
s.gha <- summary(fit.LA.gha)
print(s.gha)

# MCML
data.mod_gha_2021_urogenital$log_dist_to_water <- as.numeric(
  scale(log1p(data.mod_gha_2021_urogenital$dist_to_water_m))
)
ID.coords.gha <- create.ID.coords(
  data = data.mod_gha_2021_urogenital,
  ~ I(utm_x/1000) + I(utm_y/1000)
)
control.mcmc <- control.mcmc.MCML(n.sim = 110000, burnin = 10000, thin = 10)
par0.gha   <- coef(fit.LA.gha)
logphi_gha <- s.gha$cov.pars["log(phi)", "Estimate"]
phi_start_gha <- exp(logphi_gha)
phi_start_gha

fit.MCML.gha <- binomial.logistic.MCML(
  Positive ~ mean_temp + precip_seasonality + ndvi_mean_gha + lst_mean_gha +
    log_dist_to_water + elevation,
  units.m        = ~ Examined,
  coords         = ~ I(utm_x/1000) + I(utm_y/1000),
  ID.coords      = ID.coords.gha,
  kappa          = 0.5,
  start.cov.pars = phi_start_gha,
  fixed.rel.nugget = 0,
  control.mcmc   = control.mcmc,
  par0           = par0.gha,
  data           = data.mod_gha_2021_urogenital,
  method         = "nlminb"
)
print(summary(fit.MCML.gha))



gha_0 <- st_read(file.path(gha_dir, "gadm41_GHA_0.shp")) |> st_transform(32630)

bb   <- st_bbox(gha_0)
nxny <- 250
xseq <- seq(bb["xmin"], bb["xmax"], length.out = nxny)
yseq <- seq(bb["ymin"], bb["ymax"], length.out = nxny)
grid_df <- expand.grid(x = xseq, y = yseq)

grid_sf <- st_as_sf(grid_df, coords = c("x","y"), crs = 32630)
inside  <- st_within(grid_sf, gha_0, sparse = FALSE)[,1]
grid_sf <- grid_sf[inside, ]
grid_df <- cbind(st_coordinates(grid_sf) |> as.data.frame())
names(grid_df) <- c("x","y")

extract2 <- function(r, pts, crs = "EPSG:32630") {
  v <- if (inherits(pts,"sf")) terra::vect(pts) else terra::vect(pts, geom=c("x","y"), crs=crs)
  terra::extract(r, v)[[2]]
}

grid_df$mean_temp          <- extract2(tmean_utm,       grid_df)
grid_df$precip_seasonality <- extract2(prec_season_utm, grid_df)
grid_df$ndvi_mean_gha      <- extract2(ndvi_utm,       grid_df)
grid_df$lst_mean_gha       <- extract2(lst_utm,        grid_df)
grid_df$elevation          <- extract2(elev_utm,       grid_df)

get_dist_to_water <- function(grid_xy, wb) {
  if (inherits(wb,"SpatRaster") || inherits(wb$raster,"SpatRaster") || inherits(wb$climatology,"SpatRaster")) {
    water <- wb$climatology %||% wb$raster %||% wb
    if (!terra::same.crs(water, "EPSG:32630")) water <- terra::project(water, "EPSG:32630")
    water_to <- classify(water, rbind(c(-Inf, Inf, 1))); water_to[is.na(water)] <- NA
    dist_r <- terra::distance(water_to)
    terra::extract(dist_r, vect(grid_xy))[[1]]
  } else if (inherits(wb,"sf") || inherits(wb$rivers,"sf") || inherits(wb$water,"sf") || inherits(wb$lines,"sf")) {
    rivers <- wb$rivers %||% wb$water %||% wb$lines %||% wb
    if (is.na(st_crs(rivers))) stop("Water sf has no CRS; set EPSG:32630.")
    if (st_crs(rivers)$epsg != 32630) rivers <- st_transform(rivers, 32630)
    apply(st_distance(st_as_sf(grid_xy, coords=c("x","y"), crs=32630), rivers), 1, min) |> as.numeric()
  } else stop("Unrecognized water bundle.")
}
grid_df$dist_to_water_m <- get_dist_to_water(grid_df, wb_gha)

train_logd <- log1p(data.mod_gha_2021_urogenital$dist_to_water_m)
m <- mean(train_logd, na.rm=TRUE); s <- sd(train_logd, na.rm=TRUE)
grid_df$log_dist_to_water <- if (is.finite(s) && s > 0) (log1p(grid_df$dist_to_water_m) - m)/s else 0

grid_df <- grid_df |>
  tidyr::drop_na(mean_temp, precip_seasonality, ndvi_mean_gha, lst_mean_gha,
                 log_dist_to_water, elevation)

grid.pred_gha <- as.matrix(grid_df[, c("x","y")]) / 1000


#### PREDICTION GHANA UROGENITAL ####
pred_gha <- spatial.pred.binomial.MCML(
  object            = fit.MCML.gha,
  grid.pred         = grid.pred_gha,
  predictors        = grid_df[, c("mean_temp","precip_seasonality","ndvi_mean_gha",
                                  "lst_mean_gha","log_dist_to_water","elevation")],
  control.mcmc      = control.mcmc,
  scale.predictions = "prevalence"
)

rng_gha <- range(pred_gha$prevalence$predictions, na.rm = TRUE)
cat("Predicted prevalence St. urogenital (Ghana, 2021) ranges from",
    sprintf("%.3f", rng_gha[1]), "to", sprintf("%.3f", rng_gha[2]), "\n")


gha_1 <- st_read(file.path(gha_dir, "gadm41_GHA_1.shp")) |> st_transform(32630)
gha_2 <- st_read(file.path(gha_dir, "gadm41_GHA_2.shp")) |> st_transform(32630)

grid_df$prevalence <- pred_gha$prevalence$predictions
pred_sf_gha <- st_as_sf(grid_df, coords = c("x","y"), crs = 32630, remove = FALSE)

# National prevalence map
ggplot() +
  geom_sf(data = pred_sf_gha, aes(color = prevalence), size = 0.8, alpha = 0.9) +
  geom_sf(data = gha_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_viridis_c(name = "Prevalence (%)", option = "C", direction = 1,
                        limits = c(0, 1), labels = scales::percent) +
  coord_sf(crs = st_crs(32630)) +
  theme_minimal() +
  labs(title = "Predicted Schistosomiasis (urogenital) Prevalence",
       subtitle = "Ghana – National Level")



pred_joined_gha2 <- st_join(pred_sf_gha, gha_2)
gha_prev_by_district <- pred_joined_gha2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>
  summarize(mean_prev = mean(prevalence, na.rm = TRUE))

gha_avg_prev_admin2 <- left_join(gha_2, gha_prev_by_district, by = "NAME_2")

ggplot() +
  geom_sf(data = gha_avg_prev_admin2, aes(fill = mean_prev), color = "black", size = 0.3) +
  geom_sf(data = gha_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_viridis_c(name = "Prevalence (%)", option = "C", direction = 1,
                       limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  labs(title = "Schistosomiasis (urogenital) Prevalence by District",
       subtitle = "Ghana – District Level")





## UNCERTAINTY MAP (95% CI width) - GHANA UROGENITAL ##
grid_df$prevalence <- pred_gha$prevalence$predictions
grid_df$lower95    <- pred_gha$prevalence$quantiles[, "2.5%"]
grid_df$upper95    <- pred_gha$prevalence$quantiles[, "97.5%"]
grid_df$ci_width <- grid_df$upper95 - grid_df$lower95

pred_sf_gha <- st_as_sf(
  grid_df,
  coords = c("x","y"),
  crs = 32630,
  remove = FALSE
)

# National uncertainty map
ggplot() +
  geom_sf(data = pred_sf_gha, aes(color = ci_width), size = 0.8, alpha = 0.9) +
  geom_sf(data = gha_0, fill = NA, color = "black", linewidth = 0.8) +
  scale_color_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0,1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Predicted Schistosomiasis (Urogenital) Prevalence",
    subtitle = "Ghana – National Level (2021)"
  )


# District level uncertainty
pred_joined_gha2 <- st_join(pred_sf_gha, gha_2)

gha_ci_by_district <- pred_joined_gha2 |>
  st_drop_geometry() |>
  group_by(NAME_2) |>
  summarize(mean_ci_width = mean(ci_width, na.rm = TRUE))

gha_avg_ci_admin2 <- left_join(gha_2, gha_ci_by_district, by = "NAME_2")

ggplot() +
  geom_sf(data = gha_avg_ci_admin2, aes(fill = mean_ci_width), color = "black", size = 0.3) +
  geom_sf(data = gha_0, fill = NA, color = "black", size = 0.6) +
  scale_fill_gradient(
    name = "95% CI width\n(darker = higher certainty)",
    low = "black", high = "white",
    limits = c(0, 1),
    labels = function(x) ifelse(x %in% c(0,1), "", x)
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of Schistosomiasis (Urogenital) Prevalence by District",
    subtitle = "Ghana – District Level (2021)"
  )



pred_bundle_gha_uro <- list(
  predictions = pred_gha$prevalence$predictions,
  lower95     = pred_gha$prevalence$quantiles[, "2.5%"],
  upper95     = pred_gha$prevalence$quantiles[, "97.5%"],
  ci_width    = pred_gha$prevalence$quantiles[, "97.5%"] -
    pred_gha$prevalence$quantiles[, "2.5%"],
  grid        = grid_df,
  grid.pred   = grid.pred_gha,
  model       = fit.MCML.gha,
  rng         = rng_gha,
  date        = Sys.Date(),
  note        = "Schistosomiasis (Urogenital) prevalence predictions and uncertainty for Ghana, 2021"
)

saveRDS(pred_bundle_gha_uro,
        "gha_schisto_urogenital_predictions_2021_bundle.rds")

write_xlsx(data.mod_gha_2021_urogenital, "gha_st_data_2021.xlsx")

