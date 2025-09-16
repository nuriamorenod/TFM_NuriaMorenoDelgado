############ POPULATION DENSITY #############
library(terra)
library(sf)
library(stringr)

sf::sf_use_s2(FALSE)


######## BENIN ######################

adm2 <- st_read(file.path(benin_dir, "gadm41_BEN_2.shp"))  
adm1 <- st_read(file.path(benin_dir, "gadm41_BEN_1.shp"))  
adm0 <- st_read(file.path(benin_dir, "gadm41_BEN_0.shp")) 
pd_files <- list.files(
  pd_dir,
  pattern = "^ben_pd_(?:201[2-9]|2020)_1km_UNadj\\.tif$",
  full.names = TRUE
)
brks <- c(0, 50, 100, 250, 500, 1000, Inf)
labs <- c("0–50", "50–100", "100–250", "250–500", "500–1,000", "≥1,000")

prep_one_pd <- function(f, adm2_sf) {
  year <- str_extract(basename(f), "(?:19|20)\\d{4}") |> substr(1, 4)
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  
  if (crs(r) != crs(vect(adm2_sf))) r <- project(r, crs(vect(adm2_sf)))
  
  z <- terra::extract(r, vect(adm2_sf), fun = mean, na.rm = TRUE)
  colname <- names(r)[1]
  
  out <- adm2_sf
  out$density_mean <- z[[colname]]
  out$density_class <- cut(out$density_mean, breaks = brks,
                           include.lowest = TRUE, right = FALSE, labels = labs)
  out$year <- year
  out
}

sf_by_year <- lapply(pd_files, prep_one_pd, adm2_sf = adm2)
names(sf_by_year) <- sapply(sf_by_year, function(x) unique(x$year))

bundle <- list(
  data_by_year = sf_by_year,   
  breaks       = brks,
  labels       = labs,
  palette      = list(type = "brewer", name = "Reds"), 
  overlays     = list(adm1 = adm1, adm0 = adm0)       
)

out_rds <- file.path(pd_dir, "pd_annual", "benin_pd_admin2_2012_2020_bundle.rds")
dir.create(dirname(out_rds), showWarnings = FALSE)
saveRDS(bundle, out_rds, compress = "xz")
message("Saved: ", out_rds)



###### COTE D'IVOIRE ############
library(terra)
library(sf)
library(stringr)
sf::sf_use_s2(FALSE)

civ_adm2 <- st_read(file.path(civ_dir, "gadm41_CIV_2.shp")) 
civ_adm1 <- st_read(file.path(civ_dir, "gadm41_CIV_1.shp"))  
civ_adm0 <- st_read(file.path(civ_dir, "gadm41_CIV_0.shp"))  

civ_files <- list.files(
  pd_dir,
  pattern = "^civ_pd_(?:201[2-9]|2020)_1km_UNadj\\.tif$",
  full.names = TRUE,
  ignore.case = TRUE
)

brks <- c(0, 50, 100, 250, 500, 1000, Inf)
labs <- c("0–50", "50–100", "100–250", "250–500", "500–1,000", "≥1,000")

prep_one_pd <- function(f, adm2_sf) {
  year <- str_extract(basename(f), "(?:19|20)\\d{4}") |> substr(1, 4)
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  if (crs(r) != crs(vect(adm2_sf))) r <- project(r, crs(vect(adm2_sf)))
  
  z <- terra::extract(r, vect(adm2_sf), fun = mean, na.rm = TRUE)
  colname <- names(r)[1]
  
  out <- adm2_sf
  out$density_mean  <- z[[colname]]
  out$density_class <- cut(out$density_mean, breaks = brks,
                           include.lowest = TRUE, right = FALSE, labels = labs)
  out$year <- year
  out
}

sf_by_year <- lapply(civ_files, prep_one_pd, adm2_sf = civ_adm2)
names(sf_by_year) <- sapply(sf_by_year, function(x) unique(x$year))

bundle_civ <- list(
  data_by_year = sf_by_year,                 
  breaks       = brks,
  labels       = labs,
  palette      = list(type = "brewer", name = "Reds"),
  overlays     = list(adm1 = civ_adm1, adm0 = civ_adm0)
)

out_rds <- file.path(pd_dir, "pd_annual", "civ_pd_admin2_2012_2020_bundle.rds")
dir.create(dirname(out_rds), showWarnings = FALSE)
saveRDS(bundle_civ, out_rds, compress = "xz")
message("Saved: ", out_rds)



###### GHANA ############
library(terra)
library(sf)
library(stringr)
sf::sf_use_s2(FALSE)


ghana_adm2 <- st_read(file.path(ghana_dir, "gadm41_GHA_2.shp"))  
ghana_adm1 <- st_read(file.path(ghana_dir, "gadm41_GHA_1.shp"))
ghana_adm0 <- st_read(file.path(ghana_dir, "gadm41_GHA_0.shp")) 

ghana_files <- list.files(
  pd_dir,
  pattern = "^gha_pd_(?:201[2-9]|2020)_1km_UNadj\\.tif$",
  full.names = TRUE,
  ignore.case = TRUE
)

brks <- c(0, 50, 100, 250, 500, 1000, Inf)
labs <- c("0–50", "50–100", "100–250", "250–500", "500–1,000", "≥1,000")

prep_one_pd_gha <- function(f, adm2_sf) {
  year <- str_extract(basename(f), "(?:19|20)\\d{4}") |> substr(1, 4)
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  if (crs(r) != crs(vect(adm2_sf))) r <- project(r, crs(vect(adm2_sf)))
  
  z <- terra::extract(r, vect(adm2_sf), fun = mean, na.rm = TRUE)
  colname <- names(r)[1]
  
  out <- adm2_sf
  out$density_mean  <- z[[colname]]
  out$density_class <- cut(out$density_mean, breaks = brks,
                           include.lowest = TRUE, right = FALSE, labels = labs)
  out$year <- year
  out
}

sf_by_year_gha <- lapply(ghana_files, prep_one_pd_gha, adm2_sf = ghana_adm2)
names(sf_by_year_gha) <- sapply(sf_by_year_gha, function(x) unique(x$year))

bundle_gha <- list(
  data_by_year = sf_by_year_gha,                 
  breaks       = brks,
  labels       = labs,
  palette      = list(type = "brewer", name = "Reds"),
  overlays     = list(adm1 = ghana_adm1, adm0 = ghana_adm0)
)

out_rds_gha <- file.path(pd_dir, "pd_annual", "gha_pd_admin2_2012_2020_bundle.rds")
dir.create(dirname(out_rds_gha), showWarnings = FALSE)
saveRDS(bundle_gha, out_rds_gha, compress = "xz")
message("Saved: ", out_rds_gha)



###### TOGO #####

togo_adm2 <- st_read(file.path(togo_dir, "gadm41_TGO_2.shp")) 
togo_adm1 <- st_read(file.path(togo_dir, "gadm41_TGO_1.shp")) 
togo_adm0 <- st_read(file.path(togo_dir, "gadm41_TGO_0.shp"))  

togo_files <- list.files(
  pd_dir,
  pattern = "^tgo_pd_(?:201[2-9]|2020)_1km_UNadj\\.tif$",
  full.names = TRUE,
  ignore.case = TRUE
)

prep_one_pd_tgo <- function(f, adm2_sf) {
  year <- str_extract(basename(f), "(?:19|20)\\d{4}") |> substr(1, 4)
  r <- rast(f)
  if (is.na(crs(r))) crs(r) <- "EPSG:4326"
  if (crs(r) != crs(vect(adm2_sf))) r <- project(r, crs(vect(adm2_sf)))
  
  z <- terra::extract(r, vect(adm2_sf), fun = mean, na.rm = TRUE)
  colname <- names(r)[1]
  
  out <- adm2_sf
  out$density_mean  <- z[[colname]]
  out$density_class <- cut(out$density_mean, breaks = brks,
                           include.lowest = TRUE, right = FALSE, labels = labs)
  out$year <- year
  out
}

sf_by_year_tgo <- lapply(togo_files, prep_one_pd_tgo, adm2_sf = togo_adm2)
names(sf_by_year_tgo) <- sapply(sf_by_year_tgo, function(x) unique(x$year))

bundle_tgo <- list(
  data_by_year = sf_by_year_tgo,                 
  breaks       = brks,
  labels       = labs,
  palette      = list(type = "brewer", name = "Reds"),
  overlays     = list(adm1 = togo_adm1, adm0 = togo_adm0)
)

out_rds_tgo <- file.path(pd_dir, "pd_annual", "tgo_pd_admin2_2012_2020_bundle.rds")
dir.create(dirname(out_rds_tgo), showWarnings = FALSE)
saveRDS(bundle_tgo, out_rds_tgo, compress = "xz")
message("Saved: ", out_rds_tgo)

