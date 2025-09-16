
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(sf)
library(RColorBrewer)
library(DT)

sf::sf_use_s2(FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x

# ========= 1) STORE PATHS ONLY (NO big objects in memory) =========
data_bundles <- list(
  "Benin" = list(
    "Annual Precipitation"       = "benin_ppt_annual_2012_2024_bundle.rds",
    "Annual Mean Temperature"    = "benin_tmean_annual_2012_2024_bundle.rds",
    "Population Density"         = "benin_pd_admin2_2012_2020_bundle.rds",
    "NDVI"                       = "benin_ndvi_annual_2012_2024_bundle.rds",
    "Difference LST"             = "benin_diffLST_annual_2012_2024_bundle.rds",
    "Precipitation Seasonality"  = "benin_ppt_season_bundle.rds",
    "Elevation"                  = "benin_elevation_bundle.rds",
    "Slope"                      = "benin_slope_bundle.rds",
    "Water Bodies"               = "benin_water_bundle.rds",
    "Poverty Index"              = "benin_poverty_adm1_2012_2014_2018_2022.rds"
  ),
  "Ghana" = list(
    "Annual Precipitation"       = "ghana_ppt_annual_2012_2024_bundle.rds",
    "Annual Mean Temperature"    = "ghana_tmean_annual_2012_2024_bundle.rds",
    "Population Density"         = "gha_pd_admin2_2012_2020_bundle.rds",
    "NDVI"                       = "ghana_ndvi_annual_2012_2024_bundle.rds",
    "Difference LST"             = "gha_diffLST_annual_2012_2024_bundle.rds",
    "Precipitation Seasonality"  = "ghana_ppt_season_bundle.rds",
    "Elevation"                  = "ghana_elevation_bundle.rds",
    "Slope"                      = "ghana_slope_bundle.rds",
    "Water Bodies"               = "ghana_water_bundle.rds",
    "Poverty Index"              = "ghana_poverty_adm1_2014_2018_2022.rds"
  ),
  "Togo" = list(
    "Annual Precipitation"       = "togo_ppt_annual_2012_2024_bundle.rds",
    "Annual Mean Temperature"    = "togo_tmean_annual_2012_2024_bundle.rds",
    "Population Density"         = "tgo_pd_admin2_2012_2020_bundle.rds",
    "NDVI"                       = "togo_ndvi_annual_2012_2024_bundle.rds",
    "Difference LST"             = "tgo_diffLST_annual_2012_2024_bundle.rds",
    "Precipitation Seasonality"  = "togo_ppt_season_bundle.rds",
    "Elevation"                  = "togo_elevation_bundle.rds",
    "Slope"                      = "togo_slope_bundle.rds",
    "Water Bodies"               = "togo_water_bundle.rds",
    "Poverty Index"              = "togo_poverty_adm1_2014_2017.rds"
  ),
  "Cote d'Ivoire" = list(
    "Annual Precipitation"       = "civ_ppt_annual_2012_2024_bundle.rds",
    "Annual Mean Temperature"    = "civ_tmean_annual_2012_2024_bundle.rds",
    "Population Density"         = "civ_pd_admin2_2012_2020_bundle.rds",
    "NDVI"                       = "civ_ndvi_annual_2012_2024_bundle.rds",
    "Difference LST"             = "civ_diffLST_annual_2012_2024_bundle.rds",
    "Precipitation Seasonality"  = "civ_ppt_season_bundle.rds",
    "Elevation"                  = "civ_elevation_bundle.rds",
    "Slope"                      = "civ_slope_bundle.rds",
    "Water Bodies"               = "civ_water_bundle.rds",
    "Poverty Index"              = "civ_poverty_adm1_2021.rds"
  )
)

# ========= 2) SHAPEFILE PATHS (lazy-read + cache per session) =========
country_outline <- list(
  "Benin"         = "gadm41_BEN_0.shp",
  "Ghana"         = "gadm41_GHA_0.shp",
  "Togo"          = "gadm41_TGO_0.shp",
  "Cote d'Ivoire" = "gadm41_CIV_0.shp"
)
country_borders <- list(
  "Benin"         = "gadm41_BEN_1.shp",
  "Ghana"         = "gadm41_GHA_1.shp",
  "Togo"          = "gadm41_TGO_1.shp",
  "Cote d'Ivoire" = "gadm41_CIV_1.shp"
)
country_districts <- list(
  "Benin"         = "gadm41_BEN_2.shp",
  "Ghana"         = "gadm41_GHA_2.shp",
  "Togo"          = "gadm41_TGO_2.shp",
  "Cote d'Ivoire" = "gadm41_CIV_2.shp"
)
NAME_COL  <- "NAME_1"
NAME2_COL <- "NAME_2"

# ========= 3) OBSERVED, PREDICTIONS, PROJECTS =========
observed_files <- list(
  "Benin" = list(
    "Onchocerciasis" = "benin_oncho_observed_2017_bundle.rds",
    "Schistosomiasis" = list(
      "Intestinal"  = "benin_schisto_intestinal_observed_2015_bundle.rds",
      "Urogenital"  = "benin_schisto_urogenital_observed_2015_bundle.rds"
    )
  ),
  "Ghana" = list(
    "Onchocerciasis" = "ghana_oncho_observed_2012_bundle.rds",
    "Schistosomiasis" = list(
      "Urogenital" = "ghana_schisto_urogenital_observed_2021_bundle.rds"
    )
  ),
  "Togo" = list(
    "Onchocerciasis" = "togo_oncho_observed_2015_bundle.rds"
  ),
  "Cote d'Ivoire" = list(
    "Onchocerciasis" = "civ_oncho_observed_2016_bundle.rds",
    "Schistosomiasis" = list(
      "Intestinal" = "civ_schisto_intestinal_observed_2018_bundle.rds",
      "Urogenital" = "civ_schisto_urogenital_observed_2018_bundle.rds"
    )
  )
)

prediction_files <- list(
  "Benin" = list(
    "Onchocerciasis" = "benin_oncho_predictions_2017_lonlat.rds",
    "Schistosomiasis" = list(
      "Intestinal"  = "benin_schisto_intestinal_predictions_2015_lonlat.rds",
      "Urogenital"  = "benin_schisto_urogenital_predictions_2015_lonlat.rds"
    )
  ),
  "Ghana" = list(
    "Onchocerciasis" = "ghana_oncho_predictions_2012_lonlat.rds",
    "Schistosomiasis" = list(
      "Urogenital" = "gha_schisto_urogenital_predictions_2021_lonlat.rds"
    )
  ),
  "Togo" = list(
    "Onchocerciasis" = "togo_oncho_predictions_2015_lonlat.rds"
  ),
  "Cote d'Ivoire" = list(
    "Onchocerciasis" = "civ_oncho_predictions_2016_lonlat.rds",
    "Schistosomiasis" = list(
      "Intestinal" = "civ_schisto_intestinal_predictions_2018_lonlat.rds",
      "Urogenital" = "civ_schisto_urogenital_predictions_2018_lonlat.rds"
    )
  )
)

anesvad_projects_path <- "anesvad_projects.rds"

# ========= 4) HELPERS =========
.sfcache <- new.env(parent = emptyenv())
get_borders <- function(ctry, level = c("adm1","adm2","outline")) {
  level <- match.arg(level)
  key <- paste(ctry, level, sep = "::")
  if (exists(key, envir = .sfcache)) return(get(key, envir = .sfcache))
  shp <- switch(level,
                "outline" = country_outline[[ctry]],
                "adm1"    = country_borders[[ctry]],
                "adm2"    = country_districts[[ctry]]
  )
  sfobj <- st_read(shp, quiet = TRUE)
  if (any(!st_is_valid(sfobj))) sfobj <- st_make_valid(sfobj)
  assign(key, sfobj, envir = .sfcache)
  sfobj
}

make_pi_bundle <- function(sf_all, borders_adm1) {
  yrs <- sort(unique(sf_all$Year))
  rng <- range(sf_all$PovertyIndex, na.rm = TRUE)
  layers <- lapply(yrs, function(y) {
    ly <- sf_all[sf_all$Year == y, ]
    ly$value <- ly$PovertyIndex
    ly
  })
  names(layers) <- as.character(yrs)
  list(
    data_by_year = layers,
    palette      = brewer.pal(9, "Purples"),
    limits       = rng,
    breaks       = pretty(rng, 5),
    overlays     = list(adm1 = borders_adm1)
  )
}

load_covar_bundle <- function(country, covar_name) {
  f <- data_bundles[[country]][[covar_name]]
  validate(need(!is.null(f) && file.exists(f), paste("Missing file:", f %||% "(NULL)")))
  x <- readRDS(f)
  
  if (is.list(x) && all(c("rivers","water_bodies") %in% names(x)) &&
      inherits(x$rivers, "sf") && inherits(x$water_bodies, "sf")) {
    return(list(
      data_by_year = list(x),
      palette = NULL, limits = NULL, breaks = NULL, overlays = list()
    ))
  }
  
  if (identical(covar_name, "Poverty Index")) {
    adm1 <- get_borders(country, "adm1")
    return(make_pi_bundle(x, adm1))
  }
  
  x
}

get_df_year <- function(obj) {
  if (is.list(obj) && all(c("rivers","water_bodies") %in% names(obj)) &&
      inherits(obj$rivers, "sf") && inherits(obj$water_bodies, "sf")) {
    return(list(kind = "water", data = obj))
  }
  if (inherits(obj, "SpatRaster")) {
    df <- as.data.frame(obj, xy = TRUE, na.rm = TRUE)
    if (!nrow(df)) return(NULL)
    names(df)[1:2] <- c("x","y")
    data_cols <- setdiff(names(df), c("x","y","year"))
    if (length(data_cols) >= 1) names(df)[match(data_cols[1], names(df))] <- "value"
    return(list(kind = "raster", data = df))
  } else if (inherits(obj, "sf")) {
    return(list(kind = "sf", data = obj))
  } else if (is.data.frame(obj)) {
    df <- obj
    if (!nrow(df)) return(NULL)
    names(df)[1:2] <- c("x","y")
    data_cols <- setdiff(names(df), c("x","y","year"))
    if (length(data_cols) >= 1) names(df)[match(data_cols[1], names(df))] <- "value"
    return(list(kind = "raster", data = df))
  } else {
    return(NULL)
  }
}

extract_years <- function(bundle, covar_name = NULL, fallback_start = 2012) {
  n   <- length(bundle$data_by_year)
  nms <- names(bundle$data_by_year)
  
  if (!is.null(covar_name) && grepl("Elevation", covar_name, TRUE)) return("Static Elevation")
  if (!is.null(covar_name) && grepl("Slope", covar_name, TRUE))     return("Static Slope")
  if (!is.null(covar_name) && grepl("Water Bodies|Water", covar_name, TRUE)) return("Static Water")
  
  if (n == 1) return("Climatology over the years")
  if (!is.null(nms)) {
    yr <- stringr::str_extract(nms, "\\b(19|20)\\d{2}\\b")
    if (any(!is.na(yr))) return(yr)
  }
  if (!is.null(bundle$years)) {
    y <- suppressWarnings(as.character(bundle$years))
    if (any(!is.na(y))) return(y)
  }
  seq.int(from = fallback_start, by = 1, length.out = n) |> as.character()
}

get_uncertainty_limits <- function(disease, country) {
  if (disease == "Schistosomiasis") return(c(0,1))
  if (disease == "Onchocerciasis") {
    switch(country,
           "Benin" = c(0,0.8),
           "Ghana" = c(0,0.8),
           "Cote d'Ivoire" = c(0,1),
           "Togo" = c(0,0.5),
           c(0,1)
    )
  } else c(0,1)
}

safe_within <- function(pt_4326, polys) {
  if (is.na(sf::st_crs(polys))) sf::st_crs(polys) <- sf::st_crs(4326)
  if (any(!sf::st_is_valid(polys))) polys <- sf::st_make_valid(polys)
  bb   <- sf::st_bbox(polys)
  lon0 <- (bb["xmin"] + bb["xmax"]) / 2
  zone <- floor((lon0 + 180) / 6) + 1
  epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
  pt_utm    <- sf::st_transform(pt_4326, epsg)
  polys_utm <- sf::st_transform(polys,   epsg)
  sf::st_within(pt_utm, polys_utm, sparse = FALSE)
}

# ========= 5) UI =========
covar_sources <- list(
  "Annual Mean Temperature"    = "TerraClimate",
  "Annual Precipitation"       = "TerraClimate",
  "Difference LST"             = "MODIS",
  "Elevation"                  = "WorldClim",
  "NDVI"                       = "MODIS",
  "Population Density"         = "WorldPop Hub",
  "Poverty Index"              = "OPHI",
  "Precipitation Seasonality"  = "WorldClim",
  "Slope"                      = "WorldPop Hub",
  "Water Bodies"               = "OpenStreetMap"
)

ui <- navbarPage("NTD Dashboard",
                 # ---- Covariates ----
                 tabPanel("Covariates",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("country_pick", "Country", choices = sort(names(data_bundles)),
                                          selected = sort(names(data_bundles))[1]),
                              uiOutput("covar_ui"),
                              uiOutput("year_ui"),
                              tags$hr(),
                              h4("Click information"),
                              verbatimTextOutput("click_info"),
                              uiOutput("covar_source_ui")
                            ),
                            mainPanel(plotOutput("map", height = 600, click = "map_click"))
                          )
                 ),
                 
                 # ---- Observed Data ----
                 tabPanel("Observed Data",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("obs_disease", "Disease",
                                          choices = c("Onchocerciasis","Schistosomiasis"),
                                          selected = "Onchocerciasis"),
                              uiOutput("obs_subtype_ui"),
                              uiOutput("obs_country_ui"),   
                              tags$hr(),
                              h4("Click information"),
                              verbatimTextOutput("obs_click_info"),
                              tags$p("Source: ESPEN (WHO)",
                                     style="color:#555;font-size:90%;font-style:italic;margin-top:8px;")
                            ),
                            mainPanel(
                              plotOutput("obs_map", height = 600, click = "obs_map_click"),
                              tableOutput("obs_table")
                            )
                          )
                 ),
                 
                 # ---- Predictions ----
                 tabPanel("Predictions",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("pred_disease","Disease", choices = c("Onchocerciasis","Schistosomiasis"),
                                          selected = "Onchocerciasis"),
                              uiOutput("pred_subtype_ui"),
                              uiOutput("pred_country_ui"),
                              radioButtons("pred_type", "Map type",
                                           choices = c("Continuous (country level)"="continuous",
                                                       "Average (district level)"="district"),
                                           selected = "continuous"),
                              tags$hr(),
                              h4("Click information"),
                              conditionalPanel("input['tabs'] == 'Prediction'",
                                               verbatimTextOutput("prediction_click_info")),
                              conditionalPanel("input['tabs'] == 'Uncertainty'",
                                               verbatimTextOutput("uncertainty_click_info"))
                            ),
                            mainPanel(
                              tabsetPanel(id = "tabs",
                                          tabPanel("Prediction",  plotOutput("prediction_map",  height=600, click="prediction_map_click")),
                                          tabPanel("Uncertainty", plotOutput("uncertainty_map", height=600, click="uncertainty_map_click"))
                              )
                            )
                          )
                 ),
                 
                 # ---- Projects ----
                 tabPanel("Anesvad Projects",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput("anesvad_country", "Select Country",
                                          choices = c("Benin","Côte d'Ivoire","Ghana","Togo"),
                                          selected = "Benin"),
                              tags$p("Source: Anexo Informes de Impacto 2023, Fundación Anesvad",
                                     style="color:#555;font-size:90%;font-style:italic;margin-top:8px;")
                            ),
                            mainPanel(
                              h3(textOutput("anesvad_title")),
                              DTOutput("anesvad_table")
                            )
                          )
                 )
)

# ========= 6) SERVER =========
server <- function(input, output, session) {
  last_click <- reactiveVal("Click anywhere on the map")
  
  # --- Covariate selector  ---
  output$covar_ui <- renderUI({
    req(input$country_pick)
    covars <- sort(names(data_bundles[[input$country_pick]]))
    sel <- if (!is.null(input$covar_pick) && input$covar_pick %in% covars) input$covar_pick else covars[1]
    selectInput("covar_pick", "Covariate", choices = covars, selected = sel)
  })
  
  cur_bundle <- reactive({
    req(input$country_pick, input$covar_pick)
    load_covar_bundle(input$country_pick, input$covar_pick)
  })
  
  output$year_ui <- renderUI({
    req(cur_bundle())
    years_raw <- extract_years(cur_bundle(), covar_name = input$covar_pick, fallback_start = 2012)
    years_lbl <- suppressWarnings(sort(unique(as.numeric(years_raw[!is.na(years_raw)]))))
    if (!is.numeric(years_lbl) || length(years_lbl) == 0 || is.character(years_raw)) {
      selectInput("year_pick","Year", choices = years_raw, selected = years_raw[1], selectize = FALSE)
    } else {
      full_seq <- seq.int(min(years_lbl), max(years_lbl), by = 1)
      has_gaps <- !identical(unname(years_lbl), unname(full_seq))
      if (has_gaps) {
        selectInput("year_pick", "Year", choices = years_lbl, selected = min(years_lbl), selectize = FALSE)
      } else {
        sliderInput("year_pick", "Year", min=min(years_lbl), max=max(years_lbl),
                    value=min(years_lbl), step=1, sep="")
      }
    }
  })
  
  # Year data extraction
  df_year <- reactive({
    req(input$covar_pick)
    bundle    <- cur_bundle()
    years_raw <- extract_years(bundle, covar_name = input$covar_pick, fallback_start = 2012)
    
    # Static covariates 
    if (grepl("Elevation|Water Bodies|Water", input$covar_pick, ignore.case = TRUE)) {
      return(get_df_year(bundle$data_by_year[[1]]))
    }
    
    req(input$year_pick)
    idx <- match(as.character(input$year_pick), as.character(years_raw))
    validate(need(!is.na(idx) && idx >= 1 && idx <= length(bundle$data_by_year),
                  "Year not found in bundle"))
    get_df_year(bundle$data_by_year[[idx]])
  })
  
  # ---- Map render  ----
  output$map <- renderPlot({
    req(input$country_pick)
    yy     <- df_year(); validate(need(!is.null(yy), "No data for selected year"))
    bundle <- cur_bundle()
    
    b_sf <- get_borders(input$country_pick, "adm1")
    validate(need(!is.null(b_sf), "Missing country borders"))
    bb <- st_bbox(b_sf)
    
    title_txt <- input$covar_pick
    subt_txt <- if (grepl("Elevation", input$covar_pick, TRUE)) {
      paste("Static Elevation | Country:", input$country_pick)
    } else if (grepl("Slope", input$covar_pick, TRUE)) {
      paste("Static Slope | Country:", input$country_pick)
    } else if (grepl("Water Bodies|Water", input$covar_pick, TRUE)) {
      paste("Static Water | Country:", input$country_pick)
    } else if (!is.null(input$year_pick) && grepl("Climatology", input$year_pick, TRUE)) {
      paste("Climatology | Country:", input$country_pick)
    } else {
      paste("Year:", input$year_pick, "| Country:", input$country_pick)
    }
    
    if (yy$kind == "water") {
      wb <- yy$data
      if (is.na(st_crs(wb$water_bodies))) st_crs(wb$water_bodies) <- st_crs(4326)
      if (is.na(st_crs(wb$rivers)))       st_crs(wb$rivers)       <- st_crs(4326)
      ggplot() +
        geom_sf(data = wb$water_bodies, fill = "lightblue", color = NA) +
        geom_sf(data = wb$rivers,       color = "lightblue", linewidth = 0.3) +
        geom_sf(data = b_sf, fill = NA, color = "black", linewidth = 0.6) +
        coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
        theme_minimal() +
        labs(title = title_txt, subtitle = subt_txt, x = NULL, y = NULL)
      
    } else if (yy$kind == "raster") {
      df <- yy$data; validate(need(nrow(df) > 0, "Empty raster dataframe"))
      is_slope_cat <- ("slope_class" %in% names(df)) &&
        is.factor(df$slope_class) &&
        is.character(bundle$palette) && length(bundle$palette) >= 2
      
      pal    <- bundle$palette %||% rev(terrain.colors(20))
      limits <- bundle$limits  %||% range(df$value, na.rm = TRUE)
      breaks <- bundle$breaks  %||% pretty(limits, 5)
      
      if (is_slope_cat) {
        levs <- levels(df$slope_class)
        ggplot() +
          geom_raster(data = df, aes(x = x, y = y, fill = slope_class)) +
          geom_sf(data = b_sf, fill = NA, color = "black", linewidth = 0.6) +
          scale_fill_manual(values = pal, drop = FALSE, limits = levs, name = "Slope (°)") +
          coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
          theme_minimal() +
          labs(title = title_txt, subtitle = subt_txt, x = NULL, y = NULL)
      } else {
        ggplot() +
          geom_raster(data = df, aes(x = x, y = y, fill = value)) +
          geom_sf(data = b_sf, fill = NA, color = "black", linewidth = 0.6) +
          scale_fill_gradientn(colours = pal, limits = limits, breaks = breaks, name = input$covar_pick) +
          coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
          theme_minimal() +
          labs(title = title_txt, subtitle = subt_txt, x = NULL, y = NULL)
      }
      
    } else if (yy$kind == "sf") {
      s <- yy$data
      if (is.na(st_crs(s))) st_crs(s) <- st_crs(4326)
      if (any(!st_is_valid(s))) s <- st_make_valid(s)
      
      has_classes <- "density_class" %in% names(s)
      p <- ggplot() +
        geom_sf(data = s, aes(fill = if (has_classes) density_class else value),
                color = "black", linewidth = 0.2)
      
      if (!is.null(bundle$overlays$adm1)) {
        p <- p + geom_sf(data = bundle$overlays$adm1, fill = NA, color = "black", linewidth = 0.6)
      } else {
        p <- p + geom_sf(data = b_sf, fill = NA, color = "black", linewidth = 0.6)
      }
      
      if (has_classes) {
        pal_name <- (bundle$palette$name %||% "Reds")
        p <- p + scale_fill_brewer(palette = pal_name, drop = FALSE, name = "People / km²",
                                   guide = guide_legend(reverse = TRUE))
      } else {
        pal    <- bundle$palette %||% rev(terrain.colors(20))
        limits <- bundle$limits  %||% range(s$value, na.rm = TRUE)
        breaks <- bundle$breaks  %||% pretty(limits, 5)
        p <- p + scale_fill_gradientn(colours = pal, limits = limits, breaks = breaks,
                                      name = input$covar_pick)
      }
      
      p +
        coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
        theme_minimal() +
        labs(title = title_txt, subtitle = subt_txt, x = NULL, y = NULL)
    }
  })
  
  # ---- Click info  ----
  observeEvent(input$map_click, {
    req(input$country_pick)
    yy <- df_year(); req(!is.null(yy))
    
    if (is.null(input$map_click$x) || is.null(input$map_click$y) ||
        !is.numeric(input$map_click$x) || !is.numeric(input$map_click$y) ||
        !is.finite(input$map_click$x) || !is.finite(input$map_click$y)) {
      last_click("Click the map to see details.")
      return()
    }
    
    p_ll <- st_sfc(st_point(c(input$map_click$x, input$map_click$y)), crs = 4326)
    adm1_sf <- get_borders(input$country_pick, "adm1")
    adm2_sf <- get_borders(input$country_pick, "adm2")
    if (is.na(st_crs(adm1_sf))) st_crs(adm1_sf) <- st_crs(4326)
    if (is.na(st_crs(adm2_sf))) st_crs(adm2_sf) <- st_crs(4326)
    
    inside1 <- safe_within(p_ll, adm1_sf)
    region  <- if (any(inside1)) adm1_sf[[NAME_COL]][which(inside1)[1]] else "(outside regions)"
    
    inside2 <- safe_within(p_ll, adm2_sf)
    district <- if (any(inside2)) adm2_sf[[NAME2_COL]][which(inside2)[1]] else "(outside districts)"
    
    if (yy$kind == "water") {
      wb <- yy$data
      if (is.na(st_crs(wb$water_bodies))) st_crs(wb$water_bodies) <- st_crs(4326)
      inside_water <- safe_within(p_ll, wb$water_bodies)
      last_click(paste0(
        "Region: ", region, "\n",
        "District: ", district, "\n",
        "Water polygon: ", if (any(inside_water)) "Yes" else "No", "\n",
        "Lon: ", round(input$map_click$x, 4), ", Lat: ", round(input$map_click$y, 4)
      ))
      return(invisible())
    }
    
    if (yy$kind == "raster") {
      df <- yy$data; req(nrow(df) > 0)
      i <- which.min((df$x - input$map_click$x)^2 + (df$y - input$map_click$y)^2)
      if (grepl("Slope", input$covar_pick, TRUE)) {
        fmt_val <- paste0(round(as.numeric(df$value[i]), 1), " °")
      } else if (grepl("Seasonality", input$covar_pick, TRUE)) {
        fmt_val <- as.character(round(df$value[i], 1))
      } else if (grepl("Temperature|LST", input$covar_pick, TRUE)) {
        fmt_val <- paste0(round(df$value[i], 1), " °C")
      } else if (grepl("Precipitation", input$covar_pick, TRUE)) {
        fmt_val <- paste0(round(df$value[i], 0), " mm")
      } else if (grepl("Poverty Index", input$covar_pick, TRUE)) {
        fmt_val <- paste0(round(df$value[i], 1), "%")
      } else if (grepl("Population Density", input$covar_pick, TRUE)) {
        fmt_val <- paste0(round(df$value[i], 1), " people/km²")
      } else {
        fmt_val <- as.character(round(df$value[i], 2))
      }
      year_line <- if (!is.null(input$year_pick) && grepl("Static", input$year_pick, TRUE)) "" else paste0("Year: ", input$year_pick, "\n")
      last_click(paste0(
        "Region: ", region, "\n",
        "District: ", district, "\n",
        input$covar_pick, ": ", fmt_val, "\n",
        year_line,
        "Lon: ", round(input$map_click$x, 4), ", Lat: ", round(input$map_click$y, 4)
      ))
    } else if (yy$kind == "sf") {
      s <- yy$data
      if (is.na(st_crs(s))) st_crs(s) <- st_crs(4326)
      inside_poly <- safe_within(p_ll, s)
      if (!any(inside_poly)) {
        last_click("Clicked outside polygons.")
        return()
      }
      row <- which(inside_poly)[1]
      has_classes <- "density_class" %in% names(s)
      val_col <- if ("value" %in% names(s)) {
        "value"
      } else if ("density_mean" %in% names(s)) {
        "density_mean"
      } else if ("density" %in% names(s)) {
        "density"
      } else if ("pop_density" %in% names(s)) {
        "pop_density"
      } else if ("PovertyIndex" %in% names(s)) {
        "PovertyIndex"
      } else NA_character_
      
      fmt_val <- if (!is.na(val_col)) {
        v <- s[[val_col]][row]
        if (grepl("Poverty Index", input$covar_pick, TRUE)) {
          if (is.na(v)) "NA" else paste0(round(as.numeric(v), 1), "%")
        } else if (grepl("Population Density", input$covar_pick, TRUE)) {
          if (is.na(v)) "NA" else paste0(round(as.numeric(v), 1), " people/km²")
        } else if (grepl("Temperature|LST", input$covar_pick, TRUE)) {
          if (is.na(v)) "NA" else paste0(round(as.numeric(v), 1), " °C")
        } else if (grepl("Precipitation", input$covar_pick, TRUE)) {
          if (is.na(v)) "NA" else paste0(round(as.numeric(v), 0), " mm")
        } else {
          if (is.na(v)) "NA" else as.character(round(as.numeric(v), 2))
        }
      } else if (has_classes) {
        paste0(s$density_class[row], " (class)")
      } else "(no value column found)"
      
      year_line <- if (!is.null(input$year_pick) && grepl("Static", input$year_pick, TRUE)) "" else paste0("Year: ", input$year_pick, "\n")
      last_click(paste0(
        "Region: ", region, "\n",
        "District: ", district, "\n",
        input$covar_pick, ": ", fmt_val, "\n",
        year_line,
        "Lon: ", round(input$map_click$x, 4), ", Lat: ", round(input$map_click$y, 4)
      ))
    }
  })
  output$click_info <- renderText({ last_click() })
  
  output$covar_source_ui <- renderUI({
    req(input$covar_pick)
    src <- covar_sources[[input$covar_pick]]
    if (!is.null(src)) {
      tags$p(paste("Source:", src),
             style = "color:#555; font-size:90%; font-style:italic; margin-top:8px;")
    } else NULL
  })
  
  # ===================== Observed Data =====================
  output$obs_country_ui <- renderUI({
    req(input$obs_disease)
    
    if (input$obs_disease == "Schistosomiasis") {
      req(input$obs_subtype)
      avail <- names(observed_files)[sapply(observed_files, function(ctry_list) {
        "Schistosomiasis" %in% names(ctry_list) &&
          input$obs_subtype %in% names(ctry_list[["Schistosomiasis"]])
      })]
    } else { # Onchocerciasis
      avail <- names(observed_files)[sapply(observed_files, function(ctry_list) {
        "Onchocerciasis" %in% names(ctry_list)
      })]
    }
    
    sel <- isolate(input$obs_country)
    if (is.null(sel) || !(sel %in% avail)) sel <- sort(avail)[1]
    selectInput("obs_country", "Country", choices = sort(avail), selected = sel)
  })
  
  observe({
    if (input$obs_disease == "Schistosomiasis") {
      subtypes <- names(observed_files[[input$obs_country]][["Schistosomiasis"]])
      output$obs_subtype_ui <- renderUI({
        selectInput("obs_subtype", "Subtype", choices = subtypes)
      })
    } else {
      output$obs_subtype_ui <- renderUI({ NULL })
    }
  })
  
  obs_data <- reactive({
    req(input$obs_disease, input$obs_country)
    if (input$obs_disease == "Schistosomiasis") {
      req(input$obs_subtype)
      fname <- observed_files[[input$obs_country]][["Schistosomiasis"]][[input$obs_subtype]]
    } else {
      fname <- observed_files[[input$obs_country]][["Onchocerciasis"]]
    }
    validate(need(file.exists(fname), paste("File not found:", fname)))
    readRDS(fname)
  })
  
  output$obs_map <- renderPlot({
    bundle <- obs_data()
    sites  <- bundle$observed_sites
    b0     <- bundle$borders_lvl0
    b1     <- bundle$borders_lvl1
    ggplot() +
      geom_sf(data = b0, fill = "grey95", color = "black", linewidth = 0.7) +
      geom_sf(data = b1, fill = NA, color = "grey40", linewidth = 0.6) +
      geom_sf(data = sites, aes(size = Prevalence), color = "red", alpha = 0.6) +
      theme_minimal() +
      labs(title = paste("Observed", input$obs_disease, "in", input$obs_country),
           subtitle = paste("Year:", bundle$year),
           size = "Prevalence (%)")
  })
  
  # ===== CLICK INFO for OBSERVED DATA =====
  obs_last_click <- reactiveVal("Click anywhere on the map")
  observeEvent(input$obs_map_click, {
    bundle <- obs_data()
    req(!is.null(bundle))
    
    sites <- bundle$observed_sites
    req(!is.null(sites), nrow(sites) > 0)
    p_ll <- st_sfc(st_point(c(input$obs_map_click$x, input$obs_map_click$y)), crs = 4326)
    
    adm1_sf <- get_borders(input$obs_country, "adm1")
    adm2_sf <- get_borders(input$obs_country, "adm2")
    if (is.na(st_crs(adm1_sf))) st_crs(adm1_sf) <- st_crs(4326)
    if (is.na(st_crs(adm2_sf))) st_crs(adm2_sf) <- st_crs(4326)
    
    inside1  <- safe_within(p_ll, adm1_sf)
    region   <- if (any(inside1)) adm1_sf[[NAME_COL]][which(inside1)[1]] else "(outside regions)"
    inside2  <- safe_within(p_ll, adm2_sf)
    district <- if (any(inside2)) adm2_sf[[NAME2_COL]][which(inside2)[1]] else "(outside districts)"
    
    bb   <- st_bbox(sites)
    lon0 <- (bb["xmin"] + bb["xmax"]) / 2
    zone <- floor((lon0 + 180) / 6) + 1
    epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
    
    p_utm     <- st_transform(p_ll, epsg)
    sites_utm <- st_transform(sites, epsg)
    
    dists <- tryCatch(as.numeric(st_distance(p_utm, sites_utm)), error = function(e) NA)
    if (all(is.na(dists))) {
      obs_last_click("No nearby sites"); return()
    }
    i <- which.min(dists)
    dist_km <- round(dists[i] / 1000, 1)
    
    prev <- sites$Prevalence[i] %||% NA
    exam <- sites$Examined[i]   %||% NA
    pos  <- sites$Positive[i]   %||% NA
    
    if (dist_km > 10) {
      obs_last_click(paste0(
        "Region: ", region, "\n",
        "District: ", district, "\n",
        "No observed site within 10 km.\n",
        "Lon: ", round(input$obs_map_click$x, 4),
        ", Lat: ", round(input$obs_map_click$y, 4)
      ))
    } else {
      obs_last_click(paste0(
        "Region: ", region, "\n",
        "District: ", district, "\n",
        "Prevalence: ", ifelse(is.na(prev), "NA", paste0(round(prev, 1), "%")), "\n",
        "Examined: ", exam, "\n",
        "Positive: ", pos, "\n",
        "Distance to click: ", dist_km, " km\n",
        "Lon: ", round(input$obs_map_click$x, 4),
        ", Lat: ", round(input$obs_map_click$y, 4)
      ))
    }
  })
  
  output$obs_click_info <- renderText({ obs_last_click() })
  
  
  
  
  # ===================== Predictions =====================
  output$pred_country_ui <- renderUI({
    req(input$pred_disease)
    if (input$pred_disease == "Schistosomiasis") {
      req(input$pred_subtype)
      avail <- names(prediction_files)[sapply(prediction_files, function(ctry_list) {
        "Schistosomiasis" %in% names(ctry_list) &&
          input$pred_subtype %in% names(ctry_list[["Schistosomiasis"]])
      })]
    } else {
      avail <- names(prediction_files)[sapply(prediction_files, function(ctry_list) {
        "Onchocerciasis" %in% names(ctry_list)
      })]
    }
    sel <- isolate(input$pred_country)
    if (is.null(sel) || !(sel %in% avail)) sel <- sort(avail)[1]
    selectInput("pred_country", "Country", choices = sort(avail), selected = sel)
  })
  
  observe({
    if (input$pred_disease == "Schistosomiasis") {
      req(input$pred_country)
      subtypes <- names(prediction_files[[input$pred_country]][["Schistosomiasis"]])
      output$pred_subtype_ui <- renderUI({
        selectInput("pred_subtype", "Subtype",
                    choices = sort(subtypes),
                    selected = isolate({
                      if (!is.null(input$pred_subtype) && input$pred_subtype %in% subtypes) input$pred_subtype else sort(subtypes)[1]
                    })
        )
      })
    } else {
      output$pred_subtype_ui <- renderUI({ NULL })
    }
  })
  
  pred_data <- reactive({
    req(input$pred_disease, input$pred_country)
    if (input$pred_disease == "Schistosomiasis") {
      req(input$pred_subtype)
      fname <- prediction_files[[input$pred_country]][["Schistosomiasis"]][[input$pred_subtype]]
    } else {
      fname <- prediction_files[[input$pred_country]][["Onchocerciasis"]]
    }
    validate(need(file.exists(fname), paste("Prediction file not found:", fname)))
    readRDS(fname)
  })
  
  # Prediction map
  output$prediction_map <- renderPlot({
    bundle <- pred_data()
    df <- as.data.frame(bundle$grid_ll)
    b_sf <- get_borders(input$pred_country, "adm1")
    d_sf <- get_borders(input$pred_country, "adm2")
    
    if (input$pred_type == "continuous") {
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      bb   <- st_bbox(b_sf); lon0 <- (bb["xmin"] + bb["xmax"]) / 2
      zone <- floor((lon0 + 180)/6) + 1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      b0_sf <- b_sf |> st_transform(epsg) |> st_union() |> st_transform(4326)
      
      ggplot() +
        geom_sf(data = pts, aes(color = predictions), size = 0.8, alpha = 0.9) +
        geom_sf(data = b0_sf, fill = NA, color = "black", linewidth = 0.8) +
        scale_color_viridis_c(option = "C", limits = c(0,1), labels = scales::percent, name = "Prevalence") +
        theme_minimal() +
        labs(title = paste("Predicted", input$pred_disease, "Prevalence"),
             subtitle = paste(input$pred_country,
                              if (input$pred_disease == "Schistosomiasis") input$pred_subtype else ""))
    } else {
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      bb   <- st_bbox(d_sf); lon0 <- (bb["xmin"] + bb["xmax"])/2
      zone <- floor((lon0 + 180)/6)+1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      pts_utm <- st_transform(pts, epsg); d_sf_utm <- st_transform(d_sf, epsg)
      joined <- st_join(pts_utm, d_sf_utm)
      agg <- joined |> st_drop_geometry() |>
        group_by(!!sym(NAME2_COL)) |> summarise(pred_mean = mean(predictions, na.rm = TRUE))
      d_sf2 <- d_sf |> left_join(agg, by = setNames(NAME2_COL, NAME2_COL))
      
      ggplot() +
        geom_sf(data = d_sf2, aes(fill = pred_mean), color = "black", size = 0.3) +
        geom_sf(data = b_sf, fill = NA, color = "black", size = 0.6) +
        scale_fill_viridis_c(option = "C", limits = c(0,1), labels = scales::percent,
                             na.value = "white", name = "Prevalence") +
        theme_minimal() +
        labs(title = paste("Predicted", input$pred_disease, "Prevalence by District"),
             subtitle = input$pred_country, fill = "Mean prevalence")
    }
  })
  
  # Uncertainty map
  output$uncertainty_map <- renderPlot({
    bundle <- pred_data()
    df <- as.data.frame(bundle$grid_ll)
    b_sf <- get_borders(input$pred_country, "adm1")
    d_sf <- get_borders(input$pred_country, "adm2")
    
    lims <- get_uncertainty_limits(input$pred_disease, input$pred_country)
    
    if (input$pred_type == "continuous") {
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      bb   <- st_bbox(b_sf); lon0 <- (bb["xmin"] + bb["xmax"])/2
      zone <- floor((lon0 + 180)/6)+1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      b0_sf <- b_sf |> st_transform(epsg) |> st_union() |> st_transform(4326)
      
      ggplot() +
        geom_sf(data = pts, aes(color = ci_width), size = 0.8, alpha = 0.9) +
        geom_sf(data = b0_sf, fill = NA, color = "black", linewidth = 0.8) +
        scale_color_gradient(name = "95% credible interval width\n(darker = higher certainty)",
                             low = "black", high = "white",
                             limits = lims,
                             breaks = seq(lims[1], lims[2], length.out = 5),
                             labels = function(x) ifelse(x %in% lims, "", x)) +
        theme_minimal() +
        labs(
          title = paste("Uncertainty of Predicted", input$pred_disease, "Prevalence"),
          subtitle = paste(
            input$pred_country,
            if (input$pred_disease == "Schistosomiasis") input$pred_subtype else ""
          )
        )
    } else {
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      bb   <- st_bbox(d_sf); lon0 <- (bb["xmin"] + bb["xmax"])/2
      zone <- floor((lon0 + 180)/6)+1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      pts_utm <- st_transform(pts, epsg); d_sf_utm <- st_transform(d_sf, epsg)
      joined <- st_join(pts_utm, d_sf_utm)
      agg <- joined |> st_drop_geometry() |>
        group_by(!!sym(NAME2_COL)) |> summarise(ci_mean = mean(ci_width, na.rm = TRUE))
      d_sf2 <- d_sf |> left_join(agg, by = setNames(NAME2_COL, NAME2_COL))
      
      ggplot() +
        geom_sf(data = d_sf2, aes(fill = ci_mean), color = "black", size = 0.3) +
        geom_sf(data = b_sf, fill = NA, color = "black", size = 0.6) +
        scale_fill_gradient(name = "95% credible interval width\n(darker = higher certainty)",
                            low = "black", high = "white",
                            limits = lims, na.value = "grey90",
                            breaks = seq(lims[1], lims[2], length.out = 5),
                            labels = function(x) ifelse(x %in% lims, "", x)) +
        theme_minimal() +
        labs(
          title = paste("District-level Uncertainty of Predicted", input$pred_disease, "Prevalence"),
          subtitle = paste(
            input$pred_country,
            if (input$pred_disease == "Schistosomiasis") input$pred_subtype else ""
          )
        )
    }
  })
  
  # Click info for predictions
  prediction_last_click <- reactiveVal("Click anywhere on the map")
  observeEvent(input$prediction_map_click, {
    bundle <- pred_data()
    df <- as.data.frame(bundle$grid_ll)
    p_ll <- st_sfc(st_point(c(input$prediction_map_click$x, input$prediction_map_click$y)), crs = 4326)
    adm1_sf <- get_borders(input$pred_country, "adm1")
    adm2_sf <- get_borders(input$pred_country, "adm2")
    
    inside1 <- safe_within(p_ll, adm1_sf)
    region  <- if (any(inside1)) adm1_sf[[NAME_COL]][which(inside1)[1]] else "(outside regions)"
    inside2 <- safe_within(p_ll, adm2_sf)
    district <- if (any(inside2)) adm2_sf[[NAME2_COL]][which(inside2)[1]] else "(outside districts)"
    
    if (input$pred_type == "continuous") {
      bb   <- st_bbox(adm2_sf); lon0 <- (bb["xmin"] + bb["xmax"])/2
      zone <- floor((lon0 + 180)/6)+1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      p_utm <- st_transform(p_ll, epsg); pts_utm <- st_transform(pts, epsg)
      dists <- as.numeric(st_distance(p_utm, pts_utm))
      i <- which.min(dists)
      pred_val <- df$predictions[i]
      fmt_val <- ifelse(is.na(pred_val), "NA", paste0(round(100*pred_val, 1), "%"))
    } else {
      bb   <- st_bbox(adm2_sf); lon0 <- (bb["xmin"] + bb["xmax"])/2
      zone <- floor((lon0 + 180)/6)+1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      pts_utm <- st_transform(pts, epsg)
      adm2_utm <- st_transform(adm2_sf, epsg)
      joined <- st_join(pts_utm, adm2_utm)
      agg <- joined |> st_drop_geometry() |>
        group_by(!!sym(NAME2_COL)) |>
        summarise(pred_mean = mean(predictions, na.rm = TRUE))
      row <- agg[agg[[NAME2_COL]] == district, ]
      fmt_val <- ifelse(nrow(row) == 0, "NA", paste0(round(100*row$pred_mean, 1), "%"))
    }
    
    prediction_last_click(paste0(
      "Region: ", region, "\n",
      "District: ", district, "\n",
      "Predicted prevalence: ", fmt_val, "\n",
      "Lon: ", round(input$prediction_map_click$x, 4), ", Lat: ", round(input$prediction_map_click$y, 4)
    ))
  })
  output$prediction_click_info <- renderText({ prediction_last_click() })
  
  uncertainty_last_click <- reactiveVal("Click anywhere on the map")
  observeEvent(input$uncertainty_map_click, {
    bundle <- pred_data()
    df <- as.data.frame(bundle$grid_ll)
    p_ll <- st_sfc(st_point(c(input$uncertainty_map_click$x, input$uncertainty_map_click$y)), crs = 4326)
    adm1_sf <- get_borders(input$pred_country, "adm1")
    adm2_sf <- get_borders(input$pred_country, "adm2")
    
    inside1 <- safe_within(p_ll, adm1_sf)
    region  <- if (any(inside1)) adm1_sf[[NAME_COL]][which(inside1)[1]] else "(outside regions)"
    inside2 <- safe_within(p_ll, adm2_sf)
    district <- if (any(inside2)) adm2_sf[[NAME2_COL]][which(inside2)[1]] else "(outside districts)"
    
    if (input$pred_type == "continuous") {
      bb   <- st_bbox(adm2_sf); lon0 <- (bb["xmin"] + bb["xmax"])/2
      zone <- floor((lon0 + 180)/6)+1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      p_utm <- st_transform(p_ll, epsg); pts_utm <- st_transform(pts, epsg)
      dists <- as.numeric(st_distance(p_utm, pts_utm))
      i <- which.min(dists)
      ci_val <- df$ci_width[i]
      fmt_val <- ifelse(is.na(ci_val), "NA", round(ci_val, 3))
    } else {
      bb   <- st_bbox(adm2_sf); lon0 <- (bb["xmin"] + bb["xmax"])/2
      zone <- floor((lon0 + 180)/6)+1; epsg <- if (lon0 >= 0) 32600 + zone else 32700 + zone
      pts <- st_as_sf(df, coords = c("Longitude","Latitude"), crs = 4326)
      pts_utm <- st_transform(pts, epsg); adm2_utm <- st_transform(adm2_sf, epsg)
      joined <- st_join(pts_utm, adm2_utm)
      agg <- joined |> st_drop_geometry() |>
        group_by(!!sym(NAME2_COL)) |>
        summarise(ci_mean = mean(ci_width, na.rm = TRUE))
      row <- agg[agg[[NAME2_COL]] == district, ]
      fmt_val <- ifelse(nrow(row) == 0, "NA", round(row$ci_mean, 3))
    }
    
    uncertainty_last_click(paste0(
      "Region: ", region, "\n",
      "District: ", district, "\n",
      "95% CI width: ", fmt_val, "\n",
      "Lon: ", round(input$uncertainty_map_click$x, 4), ", Lat: ", round(input$uncertainty_map_click$y, 4)
    ))
  })
  output$uncertainty_click_info <- renderText({ uncertainty_last_click() })
  
  # ===================== Projects =====================
  output$anesvad_title <- renderText({ paste("Projects in", input$anesvad_country) })
  anesvad_filtered <- reactive({
    req(input$anesvad_country)
    df_all <- readRDS(anesvad_projects_path)
    df <- df_all[[input$anesvad_country]]
    df |>
      arrange(`Initiative Name`) |>
      select(`Program Name`,`Initiative Name`,`Initiative Code`,Status,`Main Partner`)
  })
  output$anesvad_table <- renderDT({
    datatable(anesvad_filtered(), options = list(pageLength = -1, dom = 'ft'))
  })
}

# ========= 7) RUN =========
shinyApp(ui, server)

