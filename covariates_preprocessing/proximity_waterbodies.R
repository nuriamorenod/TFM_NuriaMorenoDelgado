######## WATER BODIES ################
library(sf)
library(ggplot2)
library(terra)

rivers <- st_read("ne_10m_rivers_lake_centerlines.shp")

#### BENIN ####
library(sf)
library(ggplot2)
library(dplyr)


rivers_benin <- st_read("gis_osm_waterways_free_1.shp")
water_benin <- st_read("gis_osm_water_a_free_1.shp")

benin_sf <- st_make_valid(benin_sf)
rivers_benin <- rivers_benin %>% filter(fclass == "river")

rivers_benin <- st_transform(rivers_benin, st_crs(benin_sf))
water_benin  <- st_transform(water_benin, st_crs(benin_sf))


rivers_benin <- st_intersection(rivers_benin, benin_sf)
water_benin  <- st_intersection(water_benin, benin_sf)
rivers_all_benin <- bind_rows(
  rivers_benin %>% select(geometry) %>% mutate(type = "Rivers"),
  water_benin  %>% select(geometry) %>% mutate(type = "Rivers")
)

ggplot() +
  geom_sf(data = rivers_all_benin, color = "blue", fill = "blue", linewidth = 0.4) +
  geom_sf(data = benin_sf, fill = NA, color = "black", linewidth = 0.6) +
  theme_minimal() +
  labs(title = "Rivers — Benin")


benin_water_bundle <- list(
  rivers = rivers_benin,
  water_bodies = water_benin
)

saveRDS(benin_water_bundle, "benin_water_bundle.rds")



#### CÔTE D’IVOIRE ####
library(sf)
library(ggplot2)
library(dplyr)

rivers_civ <- st_read("gis_osm_waterways_free_1.shp")
water_civ <- st_read("gis_osm_water_a_free_1.shp")

civ_sf <- st_make_valid(civ_sf)
rivers_civ <- rivers_civ %>% filter(fclass == "river")

rivers_civ <- st_transform(rivers_civ, st_crs(civ_sf))
water_civ  <- st_transform(water_civ, st_crs(civ_sf))

rivers_civ <- st_intersection(rivers_civ, civ_sf)
water_civ  <- st_intersection(water_civ, civ_sf)

rivers_all_civ <- bind_rows(
  rivers_civ %>% select(geometry) %>% mutate(type = "Rivers"),
  water_civ  %>% select(geometry) %>% mutate(type = "Rivers")
)

ggplot() +
  geom_sf(data = rivers_all_civ, color = "blue", fill = "blue", linewidth = 0.4) +
  geom_sf(data = civ_sf, fill = NA, color = "black", linewidth = 0.6) +
  theme_minimal() +
  labs(title = "Rivers — Côte d’Ivoire")

civ_water_bundle <- list(
  rivers = rivers_civ,
  water_bodies = water_civ
)

saveRDS(civ_water_bundle, "civ_water_bundle.rds")



#### GHANA ####
library(sf)
library(ggplot2)
library(dplyr)


rivers_ghana <- st_read("gis_osm_waterways_free_1.shp")
water_ghana <- st_read("gis_osm_water_a_free_1.shp")

ghana_sf <- st_read("gadm41_GHA_0.shp")
ghana_sf <- st_make_valid(ghana_sf)

rivers_ghana <- rivers_ghana %>% filter(fclass == "river")
rivers_ghana <- st_transform(rivers_ghana, st_crs(ghana_sf))
water_ghana  <- st_transform(water_ghana, st_crs(ghana_sf))
rivers_ghana <- st_intersection(rivers_ghana, ghana_sf)
water_ghana  <- st_intersection(water_ghana, ghana_sf)

rivers_all_ghana <- bind_rows(
  rivers_ghana %>% select(geometry) %>% mutate(type = "Rivers"),
  water_ghana  %>% select(geometry) %>% mutate(type = "Rivers")
)

ggplot() +
  geom_sf(data = rivers_all_ghana, color = "blue", fill = "blue", linewidth = 0.4) +
  geom_sf(data = ghana_sf, fill = NA, color = "black", linewidth = 0.6) +
  theme_minimal() +
  labs(title = "Rivers — Ghana")


ghana_water_bundle <- list(
  rivers = rivers_ghana,
  water_bodies = water_ghana
)

saveRDS(ghana_water_bundle, "ghana_water_bundle.rds")



#### TOGO ####
library(sf)
library(ggplot2)
library(dplyr)


rivers_togo <- st_read("gis_osm_waterways_free_1.shp")
water_togo <- st_read("gis_osm_water_a_free_1.shp")

togo_sf <- st_read("gadm41_TGO_0.shp")
togo_sf <- st_make_valid(togo_sf)

rivers_togo <- rivers_togo %>% filter(fclass == "river")
rivers_togo <- st_transform(rivers_togo, st_crs(togo_sf))
water_togo  <- st_transform(water_togo, st_crs(togo_sf))

rivers_togo <- st_intersection(rivers_togo, togo_sf)
water_togo  <- st_intersection(water_togo, togo_sf)

rivers_all_togo <- bind_rows(
  rivers_togo %>% select(geometry) %>% mutate(type = "Rivers"),
  water_togo  %>% select(geometry) %>% mutate(type = "Rivers")
)

ggplot() +
  geom_sf(data = rivers_all_togo, color = "blue", fill = "blue", linewidth = 0.4) +
  geom_sf(data = togo_sf, fill = NA, color = "black", linewidth = 0.6) +
  theme_minimal() +
  labs(title = "Rivers — Togo")


togo_water_bundle <- list(
  rivers = rivers_togo,
  water_bodies = water_togo
)

saveRDS(togo_water_bundle, "togo_water_bundle.rds")

