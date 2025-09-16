############ PRECIPITATION SEASONALITY  #############

library(terra)
library(ggplot2)
library(sf)

sf::sf_use_s2(FALSE)
ppt_season <- rast(data_ppt_season)


###### BENIN #####
benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"), quiet = TRUE) |>
  st_make_valid()

benin_ll  <- st_transform(benin_sf, crs(ppt_season))
benin_vect <- vect(benin_ll)
ppt_season_benin <- mask(crop(ppt_season, benin_vect), benin_vect)

ppt_df_benin <- na.omit(as.data.frame(ppt_season_benin, xy = TRUE))
names(ppt_df_benin)[3] <- "ppt_season"

library(terra)
library(sf)

ppt_season <- rast(ppt_path)
benin_sf   <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"), quiet = TRUE)

benin_ll <- st_transform(benin_sf, crs(ppt_season))

ppt_benin <- mask(crop(ppt_season, vect(benin_ll)), vect(benin_ll))

vals <- values(ppt_benin, na.rm = TRUE)

cat("Overall ppt_season range for Benin: min =", min(vals),
    ", max =", max(vals), "\n")


## convert the data 
benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"), quiet = TRUE)

benin_ll <- st_transform(benin_sf, crs(ppt_season))
r_mask <- mask(crop(ppt_season, vect(benin_ll)), vect(benin_ll))
df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "ppt_season"

out_obj <- list(
  data_by_year = list("ppt_season" = df),
  limits  = c(35, 130),   
  breaks  = seq(40, 120, 10),
  palette = viridis(length(seq(40, 120, 10)) - 1, option = "viridis"),
  border_sf = benin_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "benin_ppt_season_bundle.rds"),
        compress = "xz")


#### Cote d'Ivoire ####
library(terra)
library(sf)
library(ggplot2)

sf::sf_use_s2(FALSE)
ppt_season <- rast(ppt_path)
civ_sf <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"), quiet = TRUE) |> st_make_valid()

civ_ll   <- st_transform(civ_sf, crs(ppt_season))
civ_vect <- vect(civ_ll)
ppt_season_civ <- mask(crop(ppt_season, civ_vect), civ_vect)
ppt_df_civ <- na.omit(as.data.frame(ppt_season_civ, xy = TRUE))
names(ppt_df_civ)[3] <- "ppt_season"


library(terra)
library(sf)
ppt_season <- rast(ppt_path)
civ_sf     <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"), quiet = TRUE)
civ_ll <- st_transform(civ_sf, crs(ppt_season))

ppt_civ <- mask(crop(ppt_season, vect(civ_ll)), vect(civ_ll))
vals <- values(ppt_civ, na.rm = TRUE)

cat("Overall ppt_season range for Côte d'Ivoire: min =", min(vals),
    ", max =", max(vals), "\n")


#convert data
civ_sf  <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"), quiet = TRUE)

civ_ll <- st_transform(civ_sf, crs(ppt_season))
r_mask <- mask(crop(ppt_season, vect(civ_ll)), vect(civ_ll))
df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "ppt_season"

out_obj <- list(
  data_by_year = list("ppt_season" = df),
  limits  = c(35, 130),    
  breaks  = seq(40, 120, 10),
  palette = viridis(length(seq(40, 120, 10)) - 1, option = "viridis"),
  border_sf = civ_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "civ_ppt_season_bundle.rds"),
        compress = "xz")


library(terra)
library(sf)
library(viridis)


dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ppt_season <- rast(ppt_path)
civ_sf <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"), quiet = TRUE) |> st_make_valid()

if (is.na(crs(ppt_season))) stop("ppt_season has no CRS; set to EPSG:4326 or correct CRS.")
civ_ll <- st_transform(civ_sf, crs(ppt_season))

r_mask <- mask(crop(ppt_season, vect(civ_ll)), vect(civ_ll))
df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
if (!nrow(df)) stop("Empty dataframe after mask; check overlap/CRS.")

names(df)[3] <- "ppt_season"
lims <- as.numeric(quantile(df$ppt_season, c(0.01, 0.99), na.rm = TRUE))
if (!all(is.finite(lims))) lims <- range(df$ppt_season, na.rm = TRUE, finite = TRUE)

breaks_c <- pretty(lims, n = 7)
breaks_c <- breaks_c[breaks_c >= min(lims) & breaks_c <= max(lims)]
if (length(breaks_c) < 3) breaks_c <- pretty(range(df$ppt_season, na.rm = TRUE), n = 7)

pal_colors <- viridis(length(breaks_c) - 1, option = "viridis")

out_obj <- list(
  data_by_year = list("ppt_season" = df),  # same structure as Benin bundle
  limits  = range(df$ppt_season, na.rm = TRUE, finite = TRUE),
  breaks  = breaks_c,
  palette = pal_colors,
  border_sf = civ_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "civ_ppt_season_bundle.rds"),
        compress = "xz")




#### GHANA ####
library(terra)
library(sf)
library(ggplot2)

sf::sf_use_s2(FALSE)

ppt_season <- rast(ppt_path)
gha_sf <- st_read(file.path(gha_dir, "gadm41_GHA_0.shp"), quiet = TRUE) |> st_make_valid()

gha_ll   <- st_transform(gha_sf, crs(ppt_season))
gha_vect <- vect(gha_ll)

ppt_season_gha <- mask(crop(ppt_season, gha_vect), gha_vect)

ppt_df_gha <- na.omit(as.data.frame(ppt_season_gha, xy = TRUE))
names(ppt_df_gha)[3] <- "ppt_season"


library(terra)
library(sf)

ppt_season <- rast(ppt_path)
gha_sf     <- st_read(file.path(gha_dir, "gadm41_GHA_0.shp"), quiet = TRUE)

gha_ll <- st_transform(gha_sf, crs(ppt_season))
ppt_gha <- mask(crop(ppt_season, vect(gha_ll)), vect(gha_ll))

vals <- values(ppt_gha, na.rm = TRUE)

cat("Overall ppt_season range for Ghana: min =", min(vals),
    ", max =", max(vals), "\n")


# convert data 
gha_sf  <- st_read(file.path(gha_dir, "gadm41_GHA_0.shp"), quiet = TRUE)

gha_ll <- st_transform(gha_sf, crs(ppt_season))
r_mask <- mask(crop(ppt_season, vect(gha_ll)), vect(gha_ll))
df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "ppt_season"

out_obj <- list(
  data_by_year = list("ppt_season" = df),
  limits  = c(35, 130),   
  breaks  = seq(40, 120, 10),
  palette = viridis(length(seq(40, 120, 10)) - 1, option = "viridis"),
  border_sf = gha_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "ghana_ppt_season_bundle.rds"),
        compress = "xz")



#### TOGO ####
library(terra)
library(sf)
library(ggplot2)

sf::sf_use_s2(FALSE)

ppt_season <- rast(ppt_path)
tgo_sf <- st_read(file.path(tgo_dir, "gadm41_TGO_0.shp"), quiet = TRUE) |> st_make_valid()

tgo_ll   <- st_transform(tgo_sf, crs(ppt_season))
tgo_vect <- vect(tgo_ll)

ppt_season_tgo <- mask(crop(ppt_season, tgo_vect), tgo_vect)

ppt_df_tgo <- na.omit(as.data.frame(ppt_season_tgo, xy = TRUE))
names(ppt_df_tgo)[3] <- "ppt_season"

ggplot() +
  geom_tile(data = ppt_df_tgo, aes(x = x, y = y, fill = ppt_season)) +
  geom_sf(data = tgo_ll, fill = NA, color = "black", linewidth = 0.6) +
  coord_sf(xlim = st_bbox(tgo_ll)[c("xmin","xmax")],
           ylim = st_bbox(tgo_ll)[c("ymin","ymax")],
           expand = FALSE) +
  scale_fill_viridis_c(name = "Precip. Seasonality", option = "viridis") +
  labs(title = "Precipitation Seasonality – TOGO",
       x = "Longitude", y = "Latitude") +
  theme_minimal()


ppt_season <- rast(ppt_path)
tgo_sf     <- st_read(file.path(tgo_dir, "gadm41_TGO_0.shp"), quiet = TRUE)

tgo_ll <- st_transform(tgo_sf, crs(ppt_season))
ppt_tgo <- mask(crop(ppt_season, vect(tgo_ll)), vect(tgo_ll))

vals <- values(ppt_tgo, na.rm = TRUE)

cat("Overall ppt_season range for Togo: min =", min(vals),
    ", max =", max(vals), "\n")


tgo_sf  <- st_read(file.path(tgo_dir, "gadm41_TGO_0.shp"), quiet = TRUE)

tgo_ll <- st_transform(tgo_sf, crs(ppt_season))
r_mask <- mask(crop(ppt_season, vect(tgo_ll)), vect(tgo_ll))
df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
names(df)[3] <- "ppt_season"

out_obj <- list(
  data_by_year = list("ppt_season" = df),
  limits  = c(35, 130),    
  breaks  = seq(40, 120, 10),
  palette = viridis(length(seq(40, 120, 10)) - 1, option = "viridis"),
  border_sf = tgo_ll
)

saveRDS(out_obj,
        file = file.path(out_dir, "togo_ppt_season_bundle.rds"),
        compress = "xz")
