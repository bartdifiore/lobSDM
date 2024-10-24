#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

glorys_df <- readRDS(here::here("Data/Derived/glorys_grid.rds"))

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = "United States of America", returnclass = "sf")


lat_lims <- c(35.2, 48)
lon_lims<- c(-75.8, -56.2)

#----------------------------------
## Merge biological and covariate data
#----------------------------------
catch_data<- readRDS(here::here("Data/Derived/combined_lobsters_by_life_class.rds"))
str(catch_data)

env_data <- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
str(env_data)

all_mod_data <- catch_data |>
    left_join(env_data, by = c("longitude" = "DECDEG_BEGLON", "latitude" = "DECDEG_BEGLAT", "trawl_id" = "ID", "year" = "EST_YEAR", "season" = "season", "date" = "DATE", "survey" = "survey"))
summary(all_mod_data)

write_rds(all_mod_data, "Data/Derived/all_model_data.rds", compress = "gz")

#----------------------------------
## Pre-model data exploration
#----------------------------------
# For this analysis, first big cut is going to be based on availability of GLORYs data (1993-2023)
year_min <- 1993
year_max <- 2023

red_mod_data <- all_mod_data |>
    filter(between(year, year_min, year_max))

summary(red_mod_data)

# What is going on with the NAs??
nas <- red_mod_data[is.na(red_mod_data$BT_seasonal),]
summary(nas)

# Specific surveys??
table(nas$year, nas$survey)

# Weird...plot the points and the GLORYs grid?
nas_sf <- st_as_sf(nas, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

ggplot() +
    geom_sf(data = region, fill = "white") +
    geom_sf(data = states, color = "gray", fill = NA) +
    geom_raster(data = glorys_df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = nas_sf) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
    scale_fill_viridis_c(name = "SST example", na.value = NA) +
    theme_minimal()

ggplot() +
    geom_sf(data = region, fill = "white") +
    geom_sf(data = states, color = "gray", fill = NA) +
    geom_raster(data = glorys_df, aes(x = x, y = y, fill = value)) +
    geom_sf(data = nas_sf) +
    coord_sf(xlim = c(-71, -69), ylim = c(41, 45), expand = FALSE) +
    scale_fill_viridis_c(name = "SST example", na.value = NA) +
    theme_minimal()

# As expected, some issues with the GLORYs grid and the location of some of the very nearshore samples. It'd be great to not have to drop these observations. 
