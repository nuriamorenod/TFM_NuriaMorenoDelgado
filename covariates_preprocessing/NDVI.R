
############ NDVI #############


#### NDVI ANNUAL MEAN ####
library(terra)
years <- list.dirs(dir_ndvi, recursive = FALSE, full.names = FALSE)

for (year in years) {
  year_path <- file.path(dir_ndvi, year)
  files <- list.files(year_path, pattern = "\\.tif$", full.names = TRUE)
  r_list <- lapply(files, rast)
  
  r_stack <- rast(r_list)
  r_mean <- mean(r_stack, na.rm = TRUE)
  
  out_file <- file.path(dir_ndvi, paste0("NDVI_", year, "_annual_mean.tif"))
  
  writeRaster(r_mean, out_file, overwrite = TRUE)
  
  message("Saved annual mean for ", year, " → ", out_file)
}


library(terra)
library(sf)

############ BENIN #############
benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))
benin_vect <- vect(benin_sf)
ndvi_files <- list.files(ndvi_dir, pattern = "NDVI_\\d{4}_annual_mean\\.tif$", full.names = TRUE)

for (f in ndvi_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  r_crop <- crop(r, benin_vect)
  r_mask <- mask(r_crop, benin_vect)
  
  out_file <- gsub("\\.tif$", "_benin.tif", f)
  writeRaster(r_mask, out_file, overwrite = TRUE)
  
  message("Saved: ", out_file)
}


ndvi_files <- list.files(ndvi_dir, pattern = "_benin\\.tif$", full.names = TRUE)
all_values <- c() 

for (f in ndvi_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall NDVI range for Benin: min =", min(all_values), ", max =", max(all_values), "\n")


library(ggplot2)
library(terra)
library(sf)


df2014_ndvi <- as.data.frame(r2014_ndvi, xy = TRUE, na.rm = TRUE)
names(df2014_ndvi)[3] <- "NDVI"

library(terra)
library(sf)
library(stringr)
library(viridis)

benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))
ndvi_files <- list.files(ndvi_dir, pattern = "_benin\\.tif$", full.names = TRUE)

prep_one_ndvi <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "NDVI"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(ndvi_files, prep_one_ndvi, border_sf = benin_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(-1, 1)
breaks_c <- seq(-1, 1, by = 0.2)
pal_colors <- colorRampPalette(c("red","yellow","darkgreen"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = benin_sf
)

saveRDS(out_obj,
        file = file.path(ndvi_dir, "benin_ndvi_annual_2012_2024_bundle.rds"),
        compress = "xz")





############ CIV #############
Çciv_sf  <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))
civ_vect <- vect(civ_sf)

all_tifs <- list.files(ndvi_dir, pattern = "\\.tif$", full.names = TRUE)
ndvi_files <- all_tifs[
  grepl("NDVI_\\d{4}_annual_mean\\.tif$", all_tifs, ignore.case = TRUE) &
    !grepl("_benin\\.tif$", all_tifs, ignore.case = TRUE) &
    !grepl("_civ\\.tif$",   all_tifs, ignore.case = TRUE)
]

for (f in ndvi_files) {
  r <- rast(f)
  if (is.na(crs(r, describe = TRUE)$code)) crs(r) <- "EPSG:4326"
  
  r_crop <- crop(r, civ_vect)
  r_mask <- mask(r_crop, civ_vect)
  
  out_file <- sub("\\.tif$", "_civ.tif", f)
  writeRaster(r_mask, out_file, overwrite = TRUE)
  message("Saved: ", out_file)
}

civ_ndvi_files <- list.files(ndvi_dir, pattern = "_civ\\.tif$", full.names = TRUE)
all_values <- c()
for (f in civ_ndvi_files) {
  r <- rast(f)
  all_values <- c(all_values, values(r, na.rm = TRUE))
}
cat("Overall NDVI range for Côte d'Ivoire: min =", min(all_values),
    ", max =", max(all_values), "\n")


library(terra)
library(sf)
library(stringr)
library(viridis)

civ_sf   <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))
ndvi_files <- list.files(ndvi_dir, pattern = "_civ\\.tif$", full.names = TRUE)

prep_one_ndvi <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "NDVI"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(ndvi_files, prep_one_ndvi, border_sf = civ_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(-1, 1)
breaks_c <- seq(-1, 1, by = 0.2)
pal_colors <- colorRampPalette(c("red","yellow","darkgreen"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = civ_sf
)

saveRDS(out_obj,
        file = file.path(ndvi_dir, "civ_ndvi_annual_2012_2024_bundle.rds"),
        compress = "xz")



############ GHANA #############
ghana_sf  <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp"))
ghana_vect <- vect(ghana_sf)

all_tifs <- list.files(ndvi_dir, pattern = "\\.tif$", full.names = TRUE)
ndvi_files <- all_tifs[
  grepl("NDVI_\\d{4}_annual_mean\\.tif$", all_tifs, ignore.case = TRUE) &
    !grepl("_(benin|civ|gha)\\.tif$", all_tifs, ignore.case = TRUE)
]

for (f in ndvi_files) {
  r <- rast(f)
  if (is.na(crs(r, describe = TRUE)$code)) crs(r) <- "EPSG:4326"
  
  r_crop <- crop(r, ghana_vect)
  r_mask <- mask(r_crop, ghana_vect)
  
  out_file <- sub("\\.tif$", "_gha.tif", f)
  writeRaster(r_mask, out_file, overwrite = TRUE)
  message("Saved: ", out_file)
}

gha_ndvi_files <- list.files(ndvi_dir, pattern = "_gha\\.tif$", full.names = TRUE)
all_values <- c()
for (f in gha_ndvi_files) {
  r <- rast(f)
  all_values <- c(all_values, values(r, na.rm = TRUE))
}
cat("Overall NDVI range for Ghana: min =", min(all_values),
    ", max =", max(all_values), "\n")



library(terra)
library(sf)
library(stringr)
library(viridis)

ghana_sf  <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp"))

ndvi_files <- list.files(ndvi_dir, pattern = "_gha\\.tif$", full.names = TRUE)

prep_one_ndvi <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "NDVI"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(ndvi_files, prep_one_ndvi, border_sf = ghana_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(-1, 1)
breaks_c <- seq(-1, 1, by = 0.2)
pal_colors <- colorRampPalette(c("red","yellow","darkgreen"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = ghana_sf
)

saveRDS(out_obj,
        file = file.path(ndvi_dir, "ghana_ndvi_annual_2012_2024_bundle.rds"),
        compress = "xz")




############ TOGO #############
library(terra)
library(sf)
library(stringr)
library(viridis)

togo_sf  <- st_read(file.path(togo_dir, "gadm41_TGO_0.shp"))
ndvi_files <- list.files(ndvi_dir, pattern = "_tgo\\.tif$", full.names = TRUE)
prep_one_ndvi <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "NDVI"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(ndvi_files, prep_one_ndvi, border_sf = togo_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(-1, 1)
breaks_c <- seq(-1, 1, by = 0.2)
pal_colors <- colorRampPalette(c("red","yellow","darkgreen"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = togo_sf
)

saveRDS(out_obj,
        file = file.path(ndvi_dir, "togo_ndvi_annual_2012_2024_bundle.rds"),
        compress = "xz")

