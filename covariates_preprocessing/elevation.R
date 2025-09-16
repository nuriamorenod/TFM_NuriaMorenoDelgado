############ ELEVATION  #############

library(terra)
library(sf)

sf::sf_use_s2(FALSE)
elev <- rast(data_elev)

###### BENIN #####
benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"), quiet = TRUE) |>
  st_make_valid()
benin_elev_ll <- st_transform(benin_sf, crs(elev))
benin_elev_vect <- vect(benin_elev_ll)

elev_benin <- mask(crop(elev, benin_elev_vect), benin_elev_vect)

elev_df_benin <- na.omit(as.data.frame(elev_benin, xy = TRUE))
names(elev_df_benin)[3] <- "elevation"


library(terra)
library(sf)
elev <- rast(elev_path)
benin_sf <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"), quiet = TRUE)
benin_ll <- st_transform(benin_sf, crs(elev))
elev_benin <- mask(crop(elev, vect(benin_ll)), vect(benin_ll))
vals <- values(elev_benin, na.rm = TRUE)
cat("Overall elevation range for Benin: min =", min(vals),
    ", max =", max(vals), "\n")

library(viridisLite) 
library(RColorBrewer)

benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"), quiet = TRUE)

benin_ll <- st_transform(benin_sf, crs(elev))
r_mask   <- mask(crop(elev, vect(benin_ll)), vect(benin_ll))
df       <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "elevation"

elev_palette <- colorRampPalette(
  c("#006400",  
    "#32CD32",  
    "#FFD700",  
    "#CD853F",  
    "#8B4513",  
    "#FFFFFF"   
  )
)(100)

out_obj <- list(
  data_by_year = list("elevation" = df),
  limits  = c(0, 1500),                        
  breaks  = seq(0, 1500, 250),                 
  palette = elev_palette,
  border_sf = benin_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "benin_elevation_bundle.rds"),
        compress = "xz")



###### CÔTE D’IVOIRE #####
library(terra)
library(sf)
elev   <- rast(elev_path)
civ_sf <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"), quiet = TRUE) |>
  st_make_valid()

civ_ll   <- st_transform(civ_sf, crs(elev))
civ_vect <- vect(civ_ll)
elev_civ <- mask(crop(elev, civ_vect), civ_vect)
elev_df_civ <- na.omit(as.data.frame(elev_civ, xy = TRUE))
names(elev_df_civ)[3] <- "elevation"
vals <- values(elev_civ, na.rm = TRUE)
cat("Overall elevation range for Côte d’Ivoire: min =", min(vals),
    ", max =", max(vals), "\n")


civ_sf  <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"), quiet = TRUE)
civ_ll <- st_transform(civ_sf, crs(elev))
r_mask <- mask(crop(elev, vect(civ_ll)), vect(civ_ll))
df     <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "elevation"

elev_palette <- colorRampPalette(
  c("#006400", "#32CD32", "#FFD700", "#CD853F", "#8B4513", "#FFFFFF")
)(100)

out_obj <- list(
  data_by_year = list("elevation" = df),
  limits  = c(0, 1500),
  breaks  = seq(0, 1500, 250),
  palette = elev_palette,
  border_sf = civ_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "civ_elevation_bundle.rds"),
        compress = "xz")



###### GHANA #####
library(terra)
library(sf)
elev   <- rast(elev_path)
gha_sf <- st_read(file.path(gha_dir, "gadm41_GHA_0.shp"), quiet = TRUE) |>
  st_make_valid()

gha_ll   <- st_transform(gha_sf, crs(elev))
gha_vect <- vect(gha_ll)
elev_gha <- mask(crop(elev, gha_vect), gha_vect)
elev_df_gha <- na.omit(as.data.frame(elev_gha, xy = TRUE))
names(elev_df_gha)[3] <- "elevation"

vals <- values(elev_gha, na.rm = TRUE)
cat("Overall elevation range for Ghana: min =", min(vals),
    ", max =", max(vals), "\n")

gha_sf  <- st_read(file.path(gha_dir, "gadm41_GHA_0.shp"), quiet = TRUE)
gha_ll <- st_transform(gha_sf, crs(elev))
r_mask <- mask(crop(elev, vect(gha_ll)), vect(gha_ll))
df     <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "elevation"

elev_palette <- colorRampPalette(
  c("#006400", "#32CD32", "#FFD700", "#CD853F", "#8B4513", "#FFFFFF")
)(100)

out_obj <- list(
  data_by_year = list("elevation" = df),
  limits  = c(0, 1500),
  breaks  = seq(0, 1500, 250),
  palette = elev_palette,
  border_sf = gha_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "ghana_elevation_bundle.rds"),
        compress = "xz")



###### TOGO #####
library(terra)
library(sf)
elev   <- rast(elev_path)
tgo_sf <- st_read(file.path(tgo_dir, "gadm41_TGO_0.shp"), quiet = TRUE) |>
  st_make_valid()

tgo_ll   <- st_transform(tgo_sf, crs(elev))
tgo_vect <- vect(tgo_ll)
elev_tgo <- mask(crop(elev, tgo_vect), tgo_vect)
elev_df_tgo <- na.omit(as.data.frame(elev_tgo, xy = TRUE))
names(elev_df_tgo)[3] <- "elevation"


vals <- values(elev_tgo, na.rm = TRUE)
cat("Overall elevation range for Togo: min =", min(vals),
    ", max =", max(vals), "\n")



tgo_sf  <- st_read(file.path(tgo_dir, "gadm41_TGO_0.shp"), quiet = TRUE)

tgo_ll <- st_transform(tgo_sf, crs(elev))
r_mask <- mask(crop(elev, vect(tgo_ll)), vect(tgo_ll))
df     <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "elevation"

elev_palette <- colorRampPalette(
  c("#006400", "#32CD32", "#FFD700", "#CD853F", "#8B4513", "#FFFFFF")
)(100)

out_obj <- list(
  data_by_year = list("elevation" = df),
  limits  = c(0, 1500),
  breaks  = seq(0, 1500, 250),
  palette = elev_palette,
  border_sf = tgo_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "togo_elevation_bundle.rds"),
        compress = "xz")
