############ SLOPE  #############

library(terra)
library(ggplot2)
library(sf)

#### BENIN ####
slope_benin <- rast("slope_benin.tif")

slope_benin_1km <- aggregate(slope_benin, fact = 10, fun = max, na.rm = TRUE)
writeRaster(slope_benin_1km, "slope_benin_1km.tif", 
            overwrite = TRUE, 
            gdal = c("COMPRESS=LZW"))

slope_benin_1km

#### CÔTE D'IVOIRE ####
slope_civ <- rast("slope_civ.tif")

slope_civ_1km <- aggregate(slope_civ, fact = 10, fun = max, na.rm = TRUE)

writeRaster(slope_civ_1km, "slope_civ_1km.tif", 
            overwrite = TRUE, 
            gdal = c("COMPRESS=LZW"))

slope_civ_1km

#### GHANA ####
slope_ghana <- rast("slope_ghana.tif")

slope_ghana_1km <- aggregate(slope_ghana, fact = 10, fun = max, na.rm = TRUE)

writeRaster(slope_ghana_1km, "slope_ghana_1km.tif", 
            overwrite = TRUE, 
            gdal = c("COMPRESS=LZW"))

slope_ghana_1km


#### TOGO ####
slope_togo <- rast("slope_togo.tif")

slope_togo_1km <- aggregate(slope_togo, fact = 10, fun = max, na.rm = TRUE)

writeRaster(slope_togo_1km, "slope_togo_1km.tif", 
            overwrite = TRUE, 
            gdal = c("COMPRESS=LZW"))

slope_togo_1km

#### BENIN ####
library(sf)
library(terra)
library(ggplot2)


benin_sf <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"), quiet = TRUE) |>
  st_make_valid()
slope_benin <- rast(file.path(data_slope, "slope_benin.tif"))
slope_benin_1km <- aggregate(slope_benin, fact = 10, fun = max, na.rm = TRUE)
writeRaster(slope_benin_1km, file.path(data_slope, "slope_benin_1km.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))
benin_sf_ll <- st_transform(benin_sf, crs(slope_benin_1km))
benin_vect  <- vect(benin_sf_ll)
slope_cm    <- mask(crop(slope_benin_1km, benin_vect), benin_vect)

slope_df <- na.omit(as.data.frame(slope_cm, xy = TRUE))
names(slope_df)[3] <- "slope_deg"

slope_df$slope_class <- cut(
  slope_df$slope_deg,
  breaks  = c(-Inf, 10, 20, 30, Inf),
  labels  = c("0 - 10", "10.1 - 20", "20.1 - 30", "30.1 - 10"),
  right   = TRUE,
  include.lowest = TRUE
)

slope_df$slope_class <- factor(
  slope_df$slope_class,
  levels = c("30.1 - 10", "20.1 - 30", "10.1 - 20", "0 - 10")
)

slope_colors <- c(
  "0 - 10"     = "#006400",  
  "10.1 - 20"  = "#7CFC00",  
  "20.1 - 30"  = "#FFFF00",
  "30.1 - 10"  = "#FFA500"  
)

out_obj <- list(
  data_by_year = list("slope" = slope_df),
  limits  = c(0, 10),
  breaks  = seq(0, 10, 10),
  palette = slope_colors,
  border_sf = benin_sf_ll
)

saveRDS(out_obj,
        file = file.path(data_slope, "benin_slope_bundle.rds"),
        compress = "xz")



#### CÔTE D'IVOIRE ####
library(sf)
library(terra)
library(ggplot2)


civ_sf <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"), quiet = TRUE) |> 
  st_make_valid()

slope_civ <- rast(file.path(data_slope, "slope_civ.tif"))

slope_civ_1km <- aggregate(slope_civ, fact = 10, fun = max, na.rm = TRUE)


civ_sf_ll <- st_transform(civ_sf, crs(slope_civ_1km))
civ_vect  <- vect(civ_sf_ll)
slope_cm  <- mask(crop(slope_civ_1km, civ_vect), civ_vect)

slope_df <- na.omit(as.data.frame(slope_cm, xy = TRUE))
names(slope_df)[3] <- "slope_deg"

slope_df$slope_class <- cut(
  slope_df$slope_deg,
  breaks = c(-Inf, 10, 20, 30, Inf),
  labels = c("0 - 10", "10.1 - 20", "20.1 - 30", "30.1 - 10"),
  right = TRUE,
  include.lowest = TRUE
)

slope_df$slope_class <- factor(
  slope_df$slope_class,
  levels = c("30.1 - 10", "20.1 - 30", "10.1 - 20", "0 - 10")
)

slope_colors <- c(
  "0 - 10"     = "#006400",  
  "10.1 - 20"  = "#7CFC00",  
  "20.1 - 30"  = "#FFFF00", 
  "30.1 - 10"  = "#FFA500"   
)


out_obj <- list(
  data_by_year = list("slope" = slope_df),
  limits  = c(0, 10),
  breaks  = seq(0, 10, 10),
  palette = slope_colors,
  border_sf = civ_sf_ll
)

saveRDS(out_obj,
        file = file.path(data_slope, "civ_slope_bundle.rds"),
        compress = "xz")



#### GHANA ####
library(sf)
library(terra)
library(ggplot2)


ghana_sf <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp"), quiet = TRUE) |> 
  st_make_valid()

slope_ghana <- rast(file.path(data_slope, "slope_ghana.tif"))
slope_ghana_1km <- aggregate(slope_ghana, fact = 10, fun = max, na.rm = TRUE)

ghana_sf_ll <- st_transform(ghana_sf, crs(slope_ghana_1km))
ghana_vect  <- vect(ghana_sf_ll)
slope_cm    <- mask(crop(slope_ghana_1km, ghana_vect), ghana_vect)

slope_df <- na.omit(as.data.frame(slope_cm, xy = TRUE))
names(slope_df)[3] <- "slope_deg"

slope_df$slope_class <- cut(
  slope_df$slope_deg,
  breaks = c(-Inf, 10, 20, 30, Inf),
  labels = c("0 - 10", "10.1 - 20", "20.1 - 30", "30.1 - 10"),
  right = TRUE,
  include.lowest = TRUE
)

slope_df$slope_class <- factor(
  slope_df$slope_class,
  levels = c("30.1 - 10", "20.1 - 30", "10.1 - 20", "0 - 10")
)

slope_colors <- c(
  "0 - 10"     = "#006400",  
  "10.1 - 20"  = "#7CFC00",  
  "20.1 - 30"  = "#FFFF00", 
  "30.1 - 10"  = "#FFA500"   
)


out_obj <- list(
  data_by_year = list("slope" = slope_df),
  limits  = c(0, 10),
  breaks  = seq(0, 10, 10),
  palette = slope_colors,
  border_sf = ghana_sf_ll
)

saveRDS(out_obj,
        file = file.path(data_slope, "ghana_slope_bundle.rds"),
        compress = "xz")



#### TOGO ####
library(sf)
library(terra)
library(ggplot2)


togo_sf <- st_read(file.path(togo_dir, "gadm41_TGO_0.shp"), quiet = TRUE) |>
  st_make_valid()

slope_togo <- rast(file.path(data_slope, "slope_togo.tif"))
slope_togo_1km <- aggregate(slope_togo, fact = 10, fun = max, na.rm = TRUE)
writeRaster(slope_togo_1km, file.path(data_slope, "slope_togo_1km.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=LZW"))

togo_sf_ll <- st_transform(togo_sf, crs(slope_togo_1km))
togo_vect  <- vect(togo_sf_ll)
slope_cm   <- mask(crop(slope_togo_1km, togo_vect), togo_vect)

slope_df <- na.omit(as.data.frame(slope_cm, xy = TRUE))
names(slope_df)[3] <- "slope_deg"
slope_df$slope_class <- cut(
  slope_df$slope_deg,
  breaks  = c(-Inf, 10, 20, 30, Inf),
  labels  = c("0 - 10", "10.1 - 20", "20.1 - 30", "30.1 - 10"),
  right   = TRUE,
  include.lowest = TRUE
)

slope_df$slope_class <- factor(
  slope_df$slope_class,
  levels = c("30.1 - 10", "20.1 - 30", "10.1 - 20", "0 - 10")
)

slope_colors <- c(
  "0 - 10"     = "#006400",
  "10.1 - 20"  = "#7CFC00",
  "20.1 - 30"  = "#FFFF00",
  "30.1 - 10"  = "#FFA500"
)


out_obj <- list(
  data_by_year = list("slope" = slope_df),
  limits  = c(0, 10),
  breaks  = seq(0, 10, 10),
  palette = slope_colors,
  border_sf = togo_sf_ll
)

saveRDS(out_obj,
        file = file.path(data_slope, "togo_slope_bundle.rds"),
        compress = "xz")
