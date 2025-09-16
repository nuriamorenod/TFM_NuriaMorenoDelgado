############ LST #############


#### LAND SURFACE TEMPERATURE ####
## LST_day 
library(terra)
files_all <- list.files(dir_LST_day, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

if (length(files_all) == 0) {
  stop("No .tif files found. Check the folder path or file extensions.")
}

years <- sub(".*_(\\d{4})_.*", "\\1", basename(files_all))
out_dir <- file.path(dir_LST_day, "annual_means")
dir.create(out_dir, showWarnings = FALSE)

for (yr in sort(unique(years))) {
  files_year <- files_all[years == yr]
  message("Year ", yr, ": found ", length(files_year), " files")
  if (length(files_year) == 0) next
  
  r_stack <- rast(files_year)          
  r_mean  <- mean(r_stack, na.rm = TRUE)
  
  out_file <- file.path(out_dir, paste0("LST_day_", yr, "_annual_mean.tif"))
  writeRaster(r_mean, out_file, overwrite = TRUE)
  message("Saved annual mean → ", out_file)
}


## LST_night 
library(terra)
files_all <- list.files(dir_LST_night, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

if (!length(files_all)) stop("No .tif files found in LST_night.")

bn <- basename(files_all)
pat <- "MOD11A2_LST_Night_(\\d{4})"  

extract_year <- function(name) {
  m <- regexec(pat, name, perl = TRUE)
  hits <- regmatches(name, m)[[1]]
  if (length(hits) >= 2) hits[2] else NA_character_
}

years <- vapply(bn, extract_year, character(1))

out_dir <- file.path(dir_LST_night, "annual_means")
dir.create(out_dir, showWarnings = FALSE)

for (yr in sort(unique(years))) {
  files_year <- files_all[years == yr]
  message("Night ", yr, ": found ", length(files_year), " files")
  
  r_stack <- rast(files_year)
  r_mean  <- mean(r_stack, na.rm = TRUE)
  
  out_file <- file.path(out_dir, paste0("LST_night_", yr, "_annual_mean.tif"))
  writeRaster(r_mean, out_file, overwrite = TRUE)
  message("Saved annual mean → ", out_file)
}


### diffLST
library(terra)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f_day   <- list.files(dir_day, pattern = "\\.tif$", full.names = TRUE)
f_night <- list.files(dir_night, pattern = "\\.tif$", full.names = TRUE)

get_year <- function(x) sub(".*_(\\d{4})_.*", "\\1", basename(x))
yr_day   <- get_year(f_day)
yr_night <- get_year(f_night)

years <- sort(intersect(yr_day, yr_night))
if (!length(years)) stop("No overlapping years found between Day and Night annual means.")

for (yr in years) {
  fdy <- f_day[yr_day == yr][1]
  fni <- f_night[yr_night == yr][1]
  r_day   <- rast(fdy)
  r_night <- rast(fni)
  
  if (!compareGeom(r_day, r_night, stopOnError = FALSE)) {
    r_night <- resample(r_night, r_day, method = "bilinear")
    message("Resampled night raster to match day raster for year ", yr)
  }
  
  # Temperature difference
  r_diff_K <- r_day - r_night
  
  out_file <- file.path(out_dir, paste0("diffLST_", yr, "_annual_mean.tif"))
  writeRaster(r_diff_K, out_file, overwrite = TRUE)
  
  message("Saved diurnal LST difference (K) for ", yr, " → ", out_file)
}



library(terra)
diff_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
r_stack <- rast(diff_files)

overall_min <- global(r_stack, "min", na.rm = TRUE)[1,1]
overall_max <- global(r_stack, "max", na.rm = TRUE)[1,1]

cat("Overall min across all years:", overall_min, "\n")
cat("Overall max across all years:", overall_max, "\n")



library(terra)
library(sf)

# ==== BENIN ====
benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))
benin_vect <- vect(benin_sf)

# diffLST 
diff_files <- list.files(diff_dir, pattern = "diffLST_\\d{4}_annual_mean\\.tif$", full.names = TRUE)

for (f in diff_files) {
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



library(terra)
library(sf)
library(stringr)

benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))
diff_files <- list.files(diff_dir, pattern = "_benin\\.tif$", full.names = TRUE)

prep_one_diffLST <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "diffLST"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{4}")
  df
}

dfs <- lapply(diff_files, prep_one_diffLST, border_sf = benin_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))
common_limits <- c(-1.05719, 17.56216)   
breaks_c <- seq(-1, 18, by = 2)          
pal_colors <- colorRampPalette(c("green","yellow","red"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = benin_sf
)

saveRDS(out_obj,
        file = file.path(diff_dir, "benin_diffLST_annual_2012_2024_bundle.rds"),
        compress = "xz")



# ==== Cote D'ivoire ====
library(terra)
library(sf)
civ_sf  <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))
civ_vect <- vect(civ_sf)

# diffLST
diff_files <- list.files(diff_dir, pattern = "diffLST_\\d{4}_annual_mean\\.tif$", full.names = TRUE)

for (f in diff_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  r_crop <- crop(r, civ_vect)
  r_mask <- mask(r_crop, civ_vect)
    out_file <- gsub("\\.tif$", "_civ.tif", f)
  writeRaster(r_mask, out_file, overwrite = TRUE)
  
  message("Saved: ", out_file)
}


library(terra)
library(sf)
library(stringr)

civ_sf  <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))
diff_files <- list.files(diff_dir, pattern = "_civ\\.tif$", full.names = TRUE)

prep_one_diffLST <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "diffLST"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{4}")
  df
}

dfs <- lapply(diff_files, prep_one_diffLST, border_sf = civ_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))
common_limits <- c(-1.05719, 17.56216)   
breaks_c <- seq(-1, 18, by = 2)
pal_colors <- colorRampPalette(c("green","yellow","red"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = civ_sf
)

saveRDS(out_obj,
        file = file.path(diff_dir, "civ_diffLST_annual_2012_2024_bundle.rds"),
        compress = "xz")




library(terra)
library(sf)

# ==== GHANA ====
gha_sf <- st_read(paste0(gha_dir, "/gadm41_GHA_0.shp"))
gha_vect <- vect(gha_sf)

# diffLST
diff_files <- list.files(diff_dir, pattern = "diffLST_\\d{4}_annual_mean\\.tif$", full.names = TRUE)

for (f in diff_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
    r_crop <- crop(r, gha_vect)
  r_mask <- mask(r_crop, gha_vect)
  
  out_file <- gsub("\\.tif$", "_gha.tif", f)
  writeRaster(r_mask, out_file, overwrite = TRUE)
  
  message("Saved: ", out_file)
}

library(terra)
library(sf)
library(stringr)

gha_sf  <- st_read(file.path(gha_dir, "gadm41_GHA_0.shp"))
diff_files <- list.files(diff_dir, pattern = "_gha\\.tif$", full.names = TRUE)

prep_one_diffLST <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "diffLST"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{4}")
  df
}

dfs <- lapply(diff_files, prep_one_diffLST, border_sf = gha_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))
common_limits <- c(-1.05719, 17.56216)  
breaks_c <- seq(-1, 18, by = 2)
pal_colors <- colorRampPalette(c("green","yellow","red"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = gha_sf
)

saveRDS(out_obj,
        file = file.path(diff_dir, "gha_diffLST_annual_2012_2024_bundle.rds"),
        compress = "xz")



library(terra)
library(sf)

# ==== TOGO ====
tgo_sf  <- st_read(file.path(tgo_dir, "gadm41_TGO_0.shp"))
tgo_vect <- vect(tgo_sf)

# diffLST
diff_files <- list.files(diff_dir, pattern = "diffLST_\\d{4}_annual_mean\\.tif$", full.names = TRUE)

for (f in diff_files) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  r_crop <- crop(r, tgo_vect)
  r_mask <- mask(r_crop, tgo_vect)
  
  out_file <- gsub("\\.tif$", "_tgo.tif", f)
  writeRaster(r_mask, out_file, overwrite = TRUE)
  
  message("Saved: ", out_file)
}


library(terra)
library(sf)
library(stringr)


tgo_sf  <- st_read(file.path(tgo_dir, "gadm41_TGO_0.shp"))
diff_files <- list.files(diff_dir, pattern = "_tgo\\.tif$", full.names = TRUE)
prep_one_diffLST <- function(f, border_sf) {
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- crs(vect(border_sf))
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch(r |> crop(vect(border_sf)) |> mask(vect(border_sf)),
                     error = function(e) NULL)
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) return(NULL)
  names(df)[3] <- "diffLST"
  df$year <- str_extract(basename(f), "(?:19|20)\\d{4}")
  df
}

dfs <- lapply(diff_files, prep_one_diffLST, border_sf = tgo_sf)
dfs <- dfs[!sapply(dfs, is.null)]
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(-1.05719, 17.56216)   # global min/max across all years
breaks_c <- seq(-1, 18, by = 2)
pal_colors <- colorRampPalette(c("green","yellow","red"))(length(breaks_c) - 1)

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_c,
  palette      = pal_colors,
  border_sf    = tgo_sf
)

saveRDS(out_obj,
        file = file.path(diff_dir, "tgo_diffLST_annual_2012_2024_bundle.rds"),
        compress = "xz")

