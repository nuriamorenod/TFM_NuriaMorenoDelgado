library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(sf)
library(RColorBrewer)

sf::sf_use_s2(FALSE)

pi_file <- "Poverty Index.xlsx"
pi_raw  <- read_excel(pi_file)

country_col <- "Country"
year_col    <- "Survey Year"
sub_col  <-  "Subnational Region"
pov_col <- "Headcount ratio: Population in multidimensional poverty (% population)"


pi <- pi_raw %>%
  select(
    Country = all_of(country_col),
    Year    = all_of(year_col),
    Region  = all_of(sub_col),
    PovertyIndex = all_of(pov_col)
  ) %>%
  mutate(
    Country      = as.character(Country),
    Region       = as.character(Region),
    Year         = as.integer(Year),
    PovertyIndex = as.numeric(PovertyIndex)
  ) %>%
  filter(!is.na(PovertyIndex), !is.na(Year))



##### BENIN ####
adm1 <- st_read(file.path(benin_dir, "gadm41_BEN_1.shp"), quiet = TRUE) |> st_make_valid()

years <- c(2012, 2014, 2018, 2022)

benin_sf_all <- lapply(years, function(y) {
  benin_year <- pi |>
    filter(Country == "Benin", Year == y) |>
    select(Region, Year, PovertyIndex)
  
  sf_y <- adm1 |>
    left_join(benin_year, by = c("NAME_1" = "Region"))
  
  sf_y$Year <- y  
  sf_y
}) |> dplyr::bind_rows()

benin_sf_all <- benin_sf_all |> filter(!is.na(PovertyIndex))

saveRDS(benin_sf_all, "benin_poverty_adm1_2012_2014_2018_2022.rds")


vmin <- min(benin_sf_all$PovertyIndex, na.rm = TRUE)
vmax <- max(benin_sf_all$PovertyIndex, na.rm = TRUE)

library(ggplot2)
ggplot(benin_sf_all) +
  geom_sf(aes(fill = PovertyIndex), color = "black", linewidth = 0.3) +
  scale_fill_gradientn(colours = brewer.pal(9, "Purples"),
                       limits = c(vmin, vmax),
                       name = "% population") +
  coord_sf(expand = FALSE) +
  facet_wrap(~ Year) +
  labs(
    title = "Benin – Poverty Index (ADM1)",
    subtitle = "Headcount ratio: Population in multidimensional poverty"
  ) +
  theme_minimal()


##### COTE D'IVOIRE ####

civ_adm1 <- st_read(file.path(civ_dir, "gadm41_CIV_1.shp"), quiet = TRUE) |> st_make_valid()
years <- c(2021)

civ_sf_all <- lapply(years, function(y) {
  civ_year <- pi |>
    filter(Country == "Cote d'Ivoire", Year == y) |>
    select(Region, Year, PovertyIndex)
  
  sf_y <- civ_adm1 |>
    left_join(civ_year, by = c("NAME_1" = "Region"))
  
  sf_y$Year <- y
  sf_y
}) |> dplyr::bind_rows()

civ_sf_all <- civ_sf_all |> filter(!is.na(PovertyIndex))

saveRDS(civ_sf_all, "civ_poverty_adm1_2021.rds")


vmin <- min(civ_sf_all$PovertyIndex, na.rm = TRUE)
vmax <- max(civ_sf_all$PovertyIndex, na.rm = TRUE)

ggplot(civ_sf_all) +
  geom_sf(aes(fill = PovertyIndex), color = "black", linewidth = 0.3) +
  scale_fill_gradientn(colours = brewer.pal(9, "Purples"),
                       limits = c(vmin, vmax),
                       name = "% population") +
  coord_sf(expand = FALSE) +
  facet_wrap(~ Year) +
  labs(
    title = "Côte d'Ivoire – Poverty Index (ADM1)",
    subtitle = "Headcount ratio: Population in multidimensional poverty"
  ) +
  theme_minimal()




##### GHANA ####
gha_adm1 <- st_read(file.path(ghana_dir, "gadm41_GHA_1.shp"), quiet = TRUE) |> st_make_valid()

years <- c(2014, 2018, 2022)

gha_sf_all <- lapply(years, function(y) {
  gha_year <- pi |>
    filter(Country == "Ghana", Year == y) |>
    select(Region, Year, PovertyIndex)
  
  sf_y <- gha_adm1 |>
    left_join(gha_year, by = c("NAME_1" = "Region"))
  
  sf_y$Year <- y
  sf_y
}) |> bind_rows()

gha_sf_all <- gha_sf_all |> filter(!is.na(PovertyIndex))
saveRDS(gha_sf_all, "ghana_poverty_adm1_2014_2018_2022.rds")

vmin <- min(gha_sf_all$PovertyIndex, na.rm = TRUE)
vmax <- max(gha_sf_all$PovertyIndex, na.rm = TRUE)

ggplot(gha_sf_all) +
  geom_sf(aes(fill = PovertyIndex), color = "black", linewidth = 0.3) +
  scale_fill_gradientn(colours = brewer.pal(9, "Purples"),
                       limits = c(vmin, vmax),
                       name = "% population") +
  coord_sf(expand = FALSE) +
  facet_wrap(~ Year) +
  labs(
    title = "Ghana – Poverty Index (ADM1)",
    subtitle = "Headcount ratio: Population in multidimensional poverty"
  ) +
  theme_minimal()



#### TOGO ####
tgo_adm1 <- st_read(file.path(togo_dir, "gadm41_TGO_1.shp"), quiet = TRUE) |> st_make_valid()

years <- c(2014, 2017)

tgo_sf_all <- lapply(years, function(y) {
  tgo_year <- pi |>
    filter(Country == "Togo", Year == y) |>
    select(Region, Year, PovertyIndex)
  
  sf_y <- tgo_adm1 |>
    left_join(tgo_year, by = c("NAME_1" = "Region"))
  
  sf_y$Year <- y
  sf_y
}) |> bind_rows()

tgo_sf_all <- tgo_sf_all |> filter(!is.na(PovertyIndex))

saveRDS(tgo_sf_all, "togo_poverty_adm1_2014_2017.rds")

vmin <- min(tgo_sf_all$PovertyIndex, na.rm = TRUE)
vmax <- max(tgo_sf_all$PovertyIndex, na.rm = TRUE)

ggplot(tgo_sf_all) +
  geom_sf(aes(fill = PovertyIndex), color = "black", linewidth = 0.3) +
  scale_fill_gradientn(colours = brewer.pal(9, "Purples"),
                       limits = c(vmin, vmax),
                       name = "% population") +
  coord_sf(expand = FALSE) +
  facet_wrap(~ Year) +
  labs(
    title = "Togo – Poverty Index (ADM1)",
    subtitle = "Headcount ratio: Population in multidimensional poverty"
  ) +
  theme_minimal()


