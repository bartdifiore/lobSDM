#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(snakecase)
library(shiny)
library(leaflet)
library(lubridate)
library(bslib)
library(scales)

glorys_df <- readRDS(here::here("Data/Derived/glorys_grid.rds"))

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = "United States of America", returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-75.8, -56.2)

#-------------------------------------------------------------
## Merge biological and covariate data for lobster
#-------------------------------------------------------------
catch_data<- readRDS(here::here("Data/Derived/combined_lobsters_by_life_class.rds"))
str(catch_data)

# Need total biomass by tow/life class
lob_df_bio <- catch_data |>
  mutate(total_weight_at_length = number_at_length * weight_at_length) |>
  group_by(scientific_name, trawl_id, longitude, latitude, season, year, survey, date, life_class) |>
  summarize("total_biomass" = sum(total_weight_at_length)) |>
  ungroup()
summary(lob_df_bio)

env_data <- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
str(env_data)
env_tows <- unique(env_data$ID)
summary(env_data)

all(lob_df_bio$trawl_id %in% env_tows)

all_mod_data_juvenile <- lob_df_bio |>
  filter(life_class == "juvenile") |>
  # Adjust the "true" NA's before adding implicit NA values
  replace_na(list(total_biomass = 99999)) |>
  full_join(env_data, by = c("longitude" = "DECDEG_BEGLON", "latitude" = "DECDEG_BEGLAT", "trawl_id" = "ID", "year" = "EST_YEAR", "season" = "season", "date" = "DATE", "survey" = "survey")) |>
  # Make NAs 0, and then change 99999 fill to NAs
  replace_na(list(total_biomass = 0)) |>
  mutate(total_biomass = na_if(total_biomass, 99999))
summary(all_mod_data_juvenile)

write_rds(all_mod_data_juvenile, "Data/Derived/all_model_data_juvenile.rds", compress = "gz")

all_mod_data_adult <- lob_df_bio |>
  filter(life_class == "adult") |>
  # Adjust the "true" NA's before adding implicit NA values
  replace_na(list(total_biomass = 99999)) |>
  full_join(env_data, by = c("longitude" = "DECDEG_BEGLON", "latitude" = "DECDEG_BEGLAT", "trawl_id" = "ID", "year" = "EST_YEAR", "season" = "season", "date" = "DATE", "survey" = "survey")) |>
  # Make NAs 0, and then change 99999 fill to NAs
  replace_na(list(total_biomass = 0)) |>
  mutate(total_biomass = na_if(total_biomass, 99999))
summary(all_mod_data_adult)

write_rds(all_mod_data_adult, "Data/Derived/all_model_data_adult.rds", compress = "gz")

#---------------------------------------------------
## Merge biological and covariate data for predators
#---------------------------------------------------
catch_data<- readRDS(here::here("Data/Derived/combined_and_filtered_predators.rds"))
str(catch_data)

# Need total biomass by tow/life class
pred_df_bio <- catch_data |>
  mutate(total_weight_at_length = number_at_length * weight_at_length) |>
  group_by(trawl_id, longitude, latitude, season, year, survey, date) |>
  summarize("total_biomass" = sum(total_weight_at_length, na.rm = T))
summary(pred_df_bio)

env_data <- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
str(env_data)

all_mod_data <- pred_df_bio |>
  # Adjust the "true" NA's before adding implicit NA values
  replace_na(list(total_biomass = 99999)) |>
  full_join(env_data, by = c("longitude" = "DECDEG_BEGLON", "latitude" = "DECDEG_BEGLAT", "trawl_id" = "ID", "year" = "EST_YEAR", "season" = "season", "date" = "DATE", "survey" = "survey")) |>
  # Make NAs 0, and then change 99999 fill to NAs
  replace_na(list(total_biomass = 0)) |>
  mutate(total_biomass = na_if(total_biomass, 99999))
summary(all_mod_data)

write_rds(all_mod_data, "Data/Derived/all_model_data_predators.rds", compress = "gz")

#----------------------------------
## Pre-model data exploration
#----------------------------------
# For this analysis, first big cut is going to be based on availability of GLORYs data (1993-2023)
year_min <- 1993
year_max <- 2023

red_mod_data <- all_mod_data |>
    filter(between(year, year_min, year_max))

summary(red_mod_data)

# Question for Bart -- What is going on with the NAs for weight at length? This seems to then naturally chuck the nA for total biomass?
t<- is.na(catch_data$weight_at_length)
t2<- catch_data[t,]
View(t2)
na_wt_len <- catch_data |>
    filter(is.na(weight_at_length))
table(na_wt_len$year, na_wt_len$survey, na_wt_len$scientific_name)

na_bio <- red_mod_data |>
    filter(is.na(total_biomass))
head(na_bio)

#----------------------------------
## Interactive map of the lobster biomass data
#----------------------------------
# Convert data to sf object

# Create SF object
points_sf<- st_as_sf(red_mod_data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

ui <- page_sidebar(
  title = "Spatio-temporal Point Data Explorer",
  theme = bs_theme(bootswatch = "flatly"),
  
  sidebar = sidebar(
    title = "Controls",
    
    # Add dataset selector
    selectInput("selected_survey", "Select Survey",
      choices = unique(points_sf$survey),
      multiple = TRUE,
      selected = unique(points_sf$survey)[1]
    ),
    
    # selectInput("selected_year", "Select Year",
    #   choices = unique(points_sf$year),
    #   multiple = TRUE,
    #   selected = unique(points_sf$year)[1]
    # ),

    # Replace dateInput with dateRangeInput
    sliderInput("year_range", "Select Year Range",
      min = min(points_sf$year, na.rm = TRUE),
      max = max(points_sf$year, na.rm = TRUE),
      step = 1,
      value = c(min(points_sf$year, na.rm = TRUE), max(points_sf$year, na.rm = TRUE)),
      sep = ""
    ),

    selectInput("selected_season", "Select Season",
      choices = unique(points_sf$season),
      multiple = TRUE,
      selected = unique(points_sf$season)[1]
    ),

    selectInput("selected_life_class", "Select Life Class",
      choices = unique(points_sf$life_class),
      multiple = TRUE,
      selected = unique(points_sf$life_class)[1]
    ),

    # Add basemap selector
    selectInput("basemap", "Select Basemap",
      choices = c(
        "Ocean" = "ocean",
        "Satellite" = "satellite",
        "Terrain" = "terrain",
        "Light" = "light"
      ),
      selected = "ocean"
    ),

    sliderInput("opacity", "Point Opacity",
      min = 0, max = 1, value = 0.8
    )
  ),
  
  layout_columns(
    fill = FALSE,
    value_box(
      title = "Number of Points",
      value = textOutput("n_points"),
      showcase = bsicons::bs_icon("geo-alt")
    ),
    value_box(
      title = "Total biomass",
      value = textOutput("total_biomass"),
      showcase = bsicons::bs_icon("calculator")
    ),
    value_box(
      title = "Survey Shown",
      value = textOutput("n_surveys"),
      showcase = bsicons::bs_icon("database")
    ),
    value_box(
      title = "Year Range",
      value = textOutput("year_range"),
      showcase = bsicons::bs_icon("calendar-range")
    )
  ),
  
  card(
    full_screen = TRUE,
    card_header("Interactive Map"),
    leafletOutput("map", height = "600px")
  )
)

server <- function(input, output, session) {

  #  # Create a reactive expression for the scaled sizes
  # point_sizes <- reactive({
  #   # Rescale values to a reasonable range for circle markers
  #   # rescale() maps the values to a range between 5 and 30,
  #   # then multiply by the user's size multiplier
  #   rescale(locations$value, to = c(5, 30)) * input$size_multiplier
  # })

  # Filter data for selected date and dataset
  filtered_data <- reactive({
    points_sf %>%
      filter(
        between(year, input$year_range[1], input$year_range[2]),
        survey %in% input$selected_survey,
        # year %in% input$selected_year,
        season %in% input$selected_season,
        life_class %in% input$selected_life_class
      )
  })
  
  # Create map
  output$map <- renderLeaflet({
    data <- filtered_data()
    
    pal <- colorFactor(
        palette = c('#1b9e77','#d95f02','#7570b3','#e7298a'),
        domain = NULL
    )
    
    # Base map with different tile options
    map <- leaflet() %>%
      setView(lng = mean(st_coordinates(points_sf)[,1]), 
              lat = mean(st_coordinates(points_sf)[,2]), 
              zoom = 5)
    
    # Add appropriate basemap based on selection
    map <- switch(input$basemap,
      "ocean" = map %>% addProviderTiles(providers$CartoDB.DarkMatter),
      "satellite" = map %>% addProviderTiles(providers$Esri.WorldImagery),
      "terrain" = map %>% addProviderTiles(providers$Stamen.Terrain),
      "light" = map %>% addProviderTiles(providers$CartoDB.Positron)
    )
    
    map %>%
      addTiles() %>%
      addCircleMarkers(
        data = data,
        radius = rescale(filtered_data()$total_biomass, to = c(5, 30)),
        color = ~pal(survey),
        fillOpacity = input$opacity,
        popup = ~paste(
          "Survey:", survey, "<br>",
          "Biomass:", round(total_biomass, 2), "<br>",
          "Year:", year, "<br>",
          "Season", season, "<br>",
          "Life Class", life_class
        )
      ) %>%
      addLegend(
        position = "bottomright",
        pal = pal,
        values = unique(points_sf$survey),
        title = "Surveys"
      )
  })
  
  # Summary statistics
  output$n_points <- renderText({
    nrow(filtered_data())
  })
  
  output$total_biomass <- renderText({
    total_biomass <- sum(filtered_data()$total_biomass, na.rm = TRUE)
    round(total_biomass, 2)
  })
  
  output$n_surveys <- renderText({
    length(input$selected_survey)
  })

  output$year_range <- renderText({
    paste(input$year_range[1], "to", input$year_range[2])
  })
}

# When launching the app
app_object <- shinyApp(ui, server)
runApp(app_object)

# To stop it programmatically
stopApp()
