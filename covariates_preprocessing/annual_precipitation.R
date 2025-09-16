############ ANNUAL MEAN PRECIPITATION #############

library(terra)

setwd(data_ppt)

years <- 2012:2024

outdir <- file.path(data_ppt, "ppt_total_annual")
dir.create(outdir, showWarnings = FALSE)

for (y in years) {
  f <- sprintf("TerraClimate_ppt_%d.nc", y)
  r <- rast(file.path(data_ppt, f))    
  stopifnot(nlyr(r) == 12)
    ppt_total_annual <- sum(r)
  names(ppt_total_annual) <- sprintf("ppt_total_%d", y)
  
  outfile <- file.path(outdir, sprintf("ppt_total_%d.tif", y))
  writeRaster(ppt_total_annual, outfile, overwrite = TRUE)
  message("Saved: ", outfile)
}

library(terra)
library(sf)
library(ggplot2)
library(stringr)


############ BENIN #############
benin_0 <- st_read(paste0(benin_dir, "/gadm41_BEN_0.shp"))
benin_sf  <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))
benin_vect <- vect(benin_0)


ppt_files <- list.files(ppt_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)

for (f in ppt_files) {
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

ppt_files <- list.files(ppt_dir, pattern = "_benin\\.tif$", full.names = TRUE)

all_values <- c() 

for (f in ppt_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Benin: min =", min(all_values), ", max =", max(all_values), "\n")




# Prepare data #
library(terra)
library(sf)
library(stringr)
library(dplyr)

benin_sf <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp"))

ppt_files <- list.files(ppt_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
ppt_files <- ppt_files[!grepl("_ghana|_togo|_civ", ppt_files)]  # exclude other countries

prep_one <- function(f, border_sf) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  border_sf <- st_transform(border_sf, crs(r))
  r_mask <- tryCatch({
    r |> crop(vect(border_sf)) |> mask(vect(border_sf))
  }, error = function(e) {
    message("Skipping file due to non-overlapping extents: ", basename(f))
    return(NULL)
  })
  
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (nrow(df) > 0) {
    names(df)[3] <- "precip"
    df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  }
  df
}

dfs <- lapply(ppt_files, prep_one, border_sf = benin_sf)
dfs <- dfs[!sapply(dfs, is.null)]  # remove NULLs
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(0, 2200)             
breaks_mm     <- seq(0, 2200, by = 500)
pal_colors    <- c("#9C3F27", "#E3A322", "#F3E12C", "#68D357", "#35B6B3", "#1E74D1", "#6B30E7") 

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_mm,
  palette      = pal_colors,
  border_sf    = benin_sf
)

saveRDS(out_obj,
        file = file.path(ppt_dir, "benin_ppt_annual_2012_2024_bundle.rds"),
        compress = "xz")




############ GHANA #############
ghana_0 <- st_read(paste0(ghana_dir, "/gadm41_GHA_0.shp"))
ghana_sf  <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp"))
ghana_vect <- vect(ghana_0)

ppt_files <- list.files(ppt_dir, pattern = "\\.tif$", full.names = TRUE)
ppt_files <- ppt_files[!grepl("_benin|_ghana", ppt_files)]

for (f in ppt_files) {
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

ppt_files <- list.files(ppt_dir, pattern = "_ghana\\.tif$", full.names = TRUE)

all_values <- c()  
for (f in ppt_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Ghana: min =", min(all_values), ", max =", max(all_values), "\n")



# Prepare data #
library(terra)
library(sf)
library(stringr)
library(dplyr)


ghana_sf  <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp"))

ppt_files <- list.files(ppt_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
ppt_files <- ppt_files[!grepl("_benin|_togo|_civ", ppt_files)]  # exclude other countries

prep_one <- function(f, border_sf) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch({
    r |> crop(vect(border_sf)) |> mask(vect(border_sf))
  }, error = function(e) {
    message("Skipping file due to non-overlapping extents: ", basename(f))
    return(NULL)
  })
  
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (nrow(df) > 0) {
    names(df)[3] <- "precip"
    df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  }
  df
}

dfs <- lapply(ppt_files, prep_one, border_sf = ghana_sf)
dfs <- dfs[!sapply(dfs, is.null)]  # remove NULLs
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(0, 2500)              
breaks_mm     <- seq(0, 2500, by = 500)
pal_colors    <- c("#9C3F27", "#E3A322", "#F3E12C", "#68D357", "#35B6B3", "#1E74D1", "#6B30E7") 

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_mm,
  palette      = pal_colors,
  border_sf    = ghana_sf
)

saveRDS(out_obj,
        file = file.path(ppt_dir, "ghana_ppt_annual_2012_2024_bundle.rds"),
        compress = "xz")





############ TOGO #############
togo_0 <- st_read(paste0(togo_dir, "/gadm41_TGO_0.shp"))
togo_sf  <- st_read(file.path(togo_dir, "gadm41_TGO_0.shp"))
togo_vect <- vect(togo_0)

ppt_files <- list.files(ppt_dir, pattern = "\\.tif$", full.names = TRUE)
ppt_files <- ppt_files[!grepl("_benin|_ghana", ppt_files)]

for (f in ppt_files) {
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

ppt_files <- list.files(ppt_dir, pattern = "_togo\\.tif$", full.names = TRUE)

all_values <- c() 

for (f in ppt_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Togo: min =", min(all_values), ", max =", max(all_values), "\n")



# Prepare the data  #
library(terra)
library(sf)
library(stringr)
library(dplyr)

togo_sf   <- st_read(file.path(togo_dir, "gadm41_TGO_0.shp"))

ppt_files <- list.files(ppt_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
ppt_files <- ppt_files[!grepl("_benin|_ghana|_civ", ppt_files)]  # exclude other countries

prep_one <- function(f, border_sf) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch({
    r |> crop(vect(border_sf)) |> mask(vect(border_sf))
  }, error = function(e) {
    message("Skipping file due to non-overlapping extents: ", basename(f))
    return(NULL)
  })
  
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (nrow(df) > 0) {
    names(df)[3] <- "precip"
    df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  }
  df
}

dfs <- lapply(ppt_files, prep_one, border_sf = togo_sf)
dfs <- dfs[!sapply(dfs, is.null)]  # remove NULLs
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(0, 2000)              
breaks_mm     <- seq(0, 2000, by = 500)
pal_colors    <- c("#9C3F27", "#E3A322", "#F3E12C", "#68D357", "#35B6B3", "#1E74D1", "#6B30E7") 

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_mm,
  palette      = pal_colors,
  border_sf    = togo_sf
)

saveRDS(out_obj,
        file = file.path(ppt_dir, "togo_ppt_annual_2012_2024_bundle.rds"),
        compress = "xz")




############ COTE D'IVORE #############
civ_0 <- st_read(paste0(civ_dir, "/gadm41_CIV_0.shp"))
civ_sf  <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))
civ_vect <- vect(civ_0)

ppt_files <- list.files(ppt_dir, pattern = "\\.tif$", full.names = TRUE)
ppt_files <- ppt_files[!grepl("_benin|_ghana|_togo", ppt_files)]

for (f in ppt_files) {
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

ppt_files <- list.files(ppt_dir, pattern = "_civ\\.tif$", full.names = TRUE)

all_values <- c()  

for (f in ppt_files) {
  r <- rast(f)
  vals <- values(r, na.rm = TRUE)
  all_values <- c(all_values, vals)
}

cat("Overall range for Civ: min =", min(all_values), ", max =", max(all_values), "\n")



# Prepare the data  #
library(terra)
library(sf)
library(stringr)
library(dplyr)


civ_sf    <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))

ppt_files <- list.files(ppt_dir, pattern = "\\.(tif|tiff|img|grd)$", full.names = TRUE)
ppt_files <- ppt_files[!grepl("_benin|_ghana|_togo", ppt_files)]  # exclude other countries

prep_one <- function(f, border_sf) {
  r <- rast(f)
  
  if (is.na(crs(r, describe=TRUE)$code)) {
    crs(r) <- "EPSG:4326"
  }
  
  border_sf <- st_transform(border_sf, crs(r))
  
  r_mask <- tryCatch({
    r |> crop(vect(border_sf)) |> mask(vect(border_sf))
  }, error = function(e) {
    message("Skipping file due to non-overlapping extents: ", basename(f))
    return(NULL)
  })
  
  if (is.null(r_mask)) return(NULL)
  
  df <- as.data.frame(r_mask, xy = TRUE, na.rm = TRUE)
  if (nrow(df) > 0) {
    names(df)[3] <- "precip"
    df$year <- str_extract(basename(f), "(?:19|20)\\d{2}")
  }
  df
}

dfs <- lapply(ppt_files, prep_one, border_sf = civ_sf)
dfs <- dfs[!sapply(dfs, is.null)]  # remove NULLs
names(dfs) <- sapply(dfs, \(d) unique(d$year))

common_limits <- c(0, 3500)             
breaks_mm     <- seq(0, 3500, by = 500)
pal_colors    <- c("#9C3F27", "#E3A322", "#F3E12C", "#68D357", "#35B6B3", "#1E74D1", "#6B30E7", "#800080") 

out_obj <- list(
  data_by_year = dfs,
  limits       = common_limits,
  breaks       = breaks_mm,
  palette      = pal_colors,
  border_sf    = civ_sf
)

saveRDS(out_obj,
        file = file.path(ppt_dir, "civ_ppt_annual_2012_2024_bundle.rds"),
        compress = "xz")

