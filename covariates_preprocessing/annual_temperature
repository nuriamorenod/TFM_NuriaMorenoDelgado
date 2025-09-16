
############ ANNUAL MEAN TEMPERATURE #############

library(terra)
setwd(data_temp)

years <- 2012:2024
outdir <- file.path(data_temp, "tmean_annual")
dir.create(outdir, showWarnings = FALSE)


for (y in years) {
  fmin <- sprintf("TerraClimate_tmin_%d.nc", y)
  fmax <- sprintf("TerraClimate_tmax_%d.nc", y)
  
  rmin <- rast(file.path(data_temp, fmin)) 
  rmax <- rast(file.path(data_temp, fmax))
  
  
  
  # Monthly mean then annual mean
  tmean_monthly <- (rmin + rmax) / 2
  tmean_annual  <- mean(tmean_monthly)  
  names(tmean_annual) <- sprintf("tmean_%d", y)
  
  outfile <- file.path(outdir, sprintf("tmean_%d.tif", y))
  writeRaster(tmean_annual, outfile, overwrite = TRUE)
  
  message("Saved: ", outfile)
}



############ BENIN #############
benin_0 <- st_read(paste0(benin_dir, "/gadm41_BEN_0.shp"))
benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))
benin_vect <- vect(benin_0)
tmean_files <- list.files(tmean_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)


for (f in tmean_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  r_crop <- crop(r, benin_vect)
  r_mask <- mask(r_crop, benin_vect)
  
  out_file <- gsub("\\.(tif|tiff|img|grd)$", "_benin.tif", f, ignore.case = TRUE)
  writeRaster(r_mask, out_file, overwrite = TRUE)
}

library(terra)

tmean_files <- list.files(tmean_dir, pattern = "_benin\\.tif$", full.names = TRUE)

all_values <- c() 

for (f in tmean_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Benin: min =", min(all_values), ", max =", max(all_values), "\n")


library(terra)
library(sf)
library(stringr)
library(viridis)

benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))
tmean_files <- list.files(tmean_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
tmean_files <- tmean_files[!grepl("_ghana|_togo|_civ", tmean_files, ignore.case = TRUE)]

prep_one_tmean <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "tmean"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(tmean_files, prep_one_tmean, border_sf = benin_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(21, 31)
breaks_c <- seq(21, 31, by = 2) 
pal_colors <- viridis::plasma(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = benin_sf
)

saveRDS(out_obj,
        file = file.path(tmean_dir, "benin_tmean_annual_2012_2024_bundle.rds"),
        compress = "xz")




############ COTE D'IVORE #############
civ_0 <- st_read(paste0(civ_dir, "/gadm41_CIV_0.shp"))
civ_sf  <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))
civ_vect <- vect(civ_0)

tmean_files <- list.files(tmean_dir, pattern = "\\.tif$", full.names = TRUE)
tmean_files <- tmean_files[!grepl("_benin|_ghana|_togo", tmean_files)]

for (f in tmean_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  r_crop <- crop(r, civ_vect)
  r_mask <- mask(r_crop, civ_vect)
  out_file <- gsub("\\.(tif|tiff|img|grd)$", "_civ.tif", f, ignore.case = TRUE)
  writeRaster(r_mask, out_file, overwrite = TRUE)
}


library(terra)

tmean_files <- list.files(tmean_dir, pattern = "_civ\\.tif$", full.names = TRUE)
all_values <- c()  

for (f in tmean_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Civ: min =", min(all_values), ", max =", max(all_values), "\n")

library(terra)
library(sf)
library(stringr)
library(viridis)

civ_sf    <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))

tmean_files <- list.files(tmean_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
tmean_files <- tmean_files[!grepl("_benin|_ghana|_togo", tmean_files, ignore.case = TRUE)]

prep_one_tmean <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "tmean"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(tmean_files, prep_one_tmean, border_sf = civ_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(21, 31)
breaks_c      <- seq(21, 31, by = 2)
pal_colors    <- viridis::plasma(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = civ_sf
)

saveRDS(out_obj,
        file = file.path(tmean_dir, "civ_tmean_annual_2012_2024_bundle.rds"),
        compress = "xz")




############ GHANA #############
ghana_0 <- st_read(paste0(ghana_dir, "/gadm41_GHA_0.shp"))
ghana_sf  <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp"))
ghana_vect <- vect(ghana_0)


tmean_files <- list.files(tmean_dir, pattern = "\\.tif$", full.names = TRUE)
tmean_files <- tmean_files[!grepl("_togo|_benin|_civ", tmean_files)]

for (f in tmean_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  r_crop <- crop(r, ghana_vect)
  r_mask <- mask(r_crop, ghana_vect)
  out_file <- gsub("\\.(tif|tiff|img|grd)$", "_ghana.tif", f, ignore.case = TRUE)
  writeRaster(r_mask, out_file, overwrite = TRUE)
}


library(terra)

tmean_files <- list.files(tmean_dir, pattern = "_ghana\\.tif$", full.names = TRUE)

all_values <- c()  

for (f in tmean_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Ghana: min =", min(all_values), ", max =", max(all_values), "\n")

library(terra)
library(sf)
library(stringr)
library(viridis)

ghana_sf  <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp"))

tmean_files <- list.files(tmean_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
tmean_files <- tmean_files[!grepl("_benin|_togo|_civ", tmean_files, ignore.case = TRUE)]

prep_one_tmean <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "tmean"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(tmean_files, prep_one_tmean, border_sf = ghana_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(21, 31)
breaks_c <- seq(21, 31, by = 2)
pal_colors <- viridis::plasma(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = ghana_sf
)

saveRDS(out_obj,
        file = file.path(tmean_dir, "ghana_tmean_annual_2012_2024_bundle.rds"),
        compress = "xz")


############ TOGO #############
togo_0 <- st_read(paste0(togo_dir, "/gadm41_TGO_0.shp"))
togo_sf  <- st_read(file.path(togo_dir, "gadm41_TGO_0.shp"))
togo_vect <- vect(togo_0)


tmean_files <- list.files(tmean_dir, pattern = "\\.tif$", full.names = TRUE)
tmean_files <- tmean_files[!grepl("_benin|_ghana|_civ", tmean_files)]

for (f in tmean_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  r_crop <- crop(r, togo_vect)
  r_mask <- mask(r_crop, togo_vect)
  out_file <- gsub("\\.(tif|tiff|img|grd)$", "_togo.tif", f, ignore.case = TRUE)
  writeRaster(r_mask, out_file, overwrite = TRUE)
}


library(terra)

tmean_files <- list.files(tmean_dir, pattern = "_togo\\.tif$", full.names = TRUE)

all_values <- c() 

for (f in tmean_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Togo: min =", min(all_values), ", max =", max(all_values), "\n")

library(terra)
library(sf)
library(stringr)
library(viridis)

togo_sf   <- st_read(file.path(togo_dir, "gadm41_TGO_0.shp"))

tmean_files <- list.files(tmean_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
tmean_files <- tmean_files[!grepl("_benin|_ghana|_civ", tmean_files, ignore.case = TRUE)]

prep_one_tmean <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "tmean"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  df
}

dfs <- lapply(tmean_files, prep_one_tmean, border_sf = togo_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(21, 31)
breaks_c <- seq(21, 31, by = 2)
pal_colors <- viridis::plasma(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = togo_sf
)

saveRDS(out_obj,
        file = file.path(tmean_dir, "togo_tmean_annual_2012_2024_bundle.rds"),
        compress = "xz")
