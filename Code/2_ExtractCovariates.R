library(devtools)
library(tidyverse)
library(stringr)
library(sf)
library(raster)
library(lubridate)
library(gmRi)
source("Code/enhance_r_funcs.R")


# Paths to Box folder
proj_box_path <- cs_path(box_group = "Mills Lab", subfolder = "Projects/sdm_workflow/data 2")
glorys_dir_path <- cs_path(box_group = "RES_Data", subfolder = "GLORYs/NW_Atl_MonthlyTemps")

# Extra helper functions
file_path_sans_ext <- function(x, compression = FALSE) {
  if (compression)
    x <- sub("[.](gz|bz2|xz)$", "", x)
  sub("([^.]+.+)\\.[[:alnum:]]+$", "\\1", x)
}

points_to_sf <- function(points) {
  sf_out <- st_as_sf(points, coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
  return(sf_out)
}

# Helper function to convert season names to middle-of-season dates
season_to_date <- function(year_season) {
  # Split all year_season values at once
  parts <- do.call(rbind, strsplit(year_season, "-"))
  years <- as.numeric(parts[, 1])
  seasons <- tolower(parts[, 2])
  
  # Define middle dates for each season
  season_dates <- list(
    winter = "02-15",
    spring = "05-15",
    summer = "08-15",
    fall = "11-15"
  )
  
  # Create the date strings
  date_str <- paste0(years, "-", season_dates[seasons])
  return(date_str)
}

#---------------------------------------------
## Generate dataframe with all unique tow locations across surveys
#---------------------------------------------
ma_meta <- read_rds("Data/TrawlSurvey_data/mass_weight_at_length.rds") %>%
  mutate(date = lubridate::date(date))

dfo_meta <- read_rds("Data/TrawlSurvey_data/dfo_weight_at_length.rds")

nefsc_meta <- read_rds("Data/TrawlSurvey_data/nefsc_both_weight_at_length.rds") %>% 
  mutate(scientific_name = str_to_sentence(scientific_name))

me_meta <- read_rds("Data/TrawlSurvey_data/me_both_weight_at_length.rds") %>%
  mutate(date = lubridate::date(date))

all_tows <- rbind(ma_meta, dfo_meta, nefsc_meta, me_meta) %>%
  mutate(season = str_to_sentence(season)) %>%
  ungroup() %>% 
  distinct(longitude, latitude, trawl_id, survey, season, year, date) %>%
  rename(ID = trawl_id, DATE = date, EST_YEAR = year, DECDEG_BEGLAT = latitude, DECDEG_BEGLON = longitude)

# Save it
write_rds(all_tows, "Data/Derived/all_tows.rds", compress = "gz")

# Check observations
plot(all_tows$DECDEG_BEGLON, all_tows$DECDEG_BEGLAT)

#---------------------------------------------
## Extract static covariates
#---------------------------------------------
tow_name_check<- c("ID", "DATE", "EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON")
all(tow_name_check %in% names(all_tows))

# Read in static raster layer files from a directory and store them in a list. 
static_dir<- paste0(proj_box_path, "/covariates/static")
static_files<- list.files(static_dir, pattern = ".grd$", full.names = TRUE)
static_layers<- vector("list", length(static_files))

for(i in seq_along(static_files)){
  static_layers[[i]]<- raster(static_files[[i]])
}

# Names?
names(static_layers)[i]<- "Depth"

# Get all_tows into sf object
all_tows_sf<- points_to_sf(all_tows)

# Run static_extract_wrapper
all_tows_with_static_covs<- static_extract_wrapper(static_covariates_list = static_layers, sf_points = all_tows_sf, interp_nas = TRUE, date_col_name = "DATE", df_sf = "sf", out_dir = "Data/Derived/")

# Check it out
summary(all_tows_with_static_covs)

#---------------------------------------------
## Extract dynamic
#---------------------------------------------
# Read in dynamic raster stack files from a directory and store them in a list.
dynamic_files <- list.files(glorys_dir_path, pattern = ".nc$", full.names = TRUE)[c(1, 5)]
dynamic_stacks<- vector("list", length(dynamic_files))

for(i in seq_along(dynamic_files)){
  dynamic_stacks[[i]]<- raster::stack(dynamic_files[[i]])
}

# Names?
names(dynamic_stacks)<- c("BT", "SST")

# Run dynamic_2d_extract_wrapper function
all_tows_with_all_covs<- dynamic_2d_extract_wrapper(dynamic_covariates_list = dynamic_stacks, t_summ = "seasonal", t_position = NULL, interp_nas = TRUE, sf_points = all_tows_with_static_covs, date_col_name = "DATE", df_sf = "df", out_dir = "Data/Derived/")

# Check it out
summary(all_tows_with_all_covs)

#---------------------------------------------
## Get a GLORYs grid for later use
#---------------------------------------------
glorys1 <- raster::stack(dynamic_files[[1]])[[1]] |>
  as.data.frame(xy = TRUE)
colnames(glorys1)[3] <- "value"
write_rds(glorys1, here::here("Data/Derived/glorys_grid.rds"), compress = "gz")

#---------------------------------------------
## Get a GLORYs prediction dataframe
#---------------------------------------------
glorys1<- readRDS(here::here("Data/Derived/glorys_grid.rds"))
glorys_sf <- st_as_sf(glorys1, coords = c("x", "y"), crs = 4326, remove = FALSE)

# We really only need points within the study domain. A lot of different ways of doing this...I don't think we have a shapefile for all of the areas, particularly coastal state surveys. So, trying convex hull instead?
t<- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
t_sf<- st_as_sf(t, coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326)

create_hull_buffer <- function(data_sf, buffer_dist = 0.1, concavity = 2) {
  # Extract coordinates
  coords <- st_coordinates(data_sf)
  # Create concave hull
  hull <- concaveman(coords, concavity = concavity)
  hull_sf <- st_polygon(list(hull)) %>% st_sfc(crs = 4326)
  if (buffer_dist > 0) {
    hull_sf <- st_buffer(hull_sf, dist = buffer_dist)
  }
  return(hull_sf)
}

boundary_sf<- create_hull_buffer(t_sf, buffer_dist = 0, concavity = 2)
plot(boundary_sf)

# Keep prediction points within the mesh
# Method 1: Using st_filter (faster for large datasets)
glorys_df <- st_filter(glorys_sf, boundary_sf) |>
  st_drop_geometry() |>
  dplyr::select(-value)

points(glorys_df$x, glorys_df$y)

# We are going to want to have this for every time step...
env_data <- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
str(env_data)

# Make prediction dataframe
seasons <- c("Spring", "Fall", "Summer")
years <- sort(unique(env_data$EST_YEAR))

glorys_df<- replicate_df(glorys_df, time_name = "season", time_values = seasons) |>
  replicate_df(time_name = "year", time_values = years) |>
  mutate(
    "Year_Season" = paste(year, season, sep = "-"),
    "Date" = as.Date(season_to_date(Year_Season))
  )

# Read in static raster layer files from a directory and store them in a list. 
static_dir<- paste0(proj_box_path, "/covariates/static")
static_files<- list.files(static_dir, pattern = ".grd$", full.names = TRUE)
static_layers<- vector("list", length(static_files))

for(i in seq_along(static_files)){
  static_layers[[i]]<- raster(static_files[[i]])
}

# Names?
names(static_layers)[i]<- "Depth"

# Get all_tows into sf object
all_tows_sf<- st_as_sf(glorys_df, coords = c("x", "y"), crs = 4326, remove = FALSE)

# Run static_extract_wrapper
all_tows_with_static_covs <- static_extract_wrapper(static_covariates_list = static_layers, sf_points = all_tows_sf, interp_nas = FALSE, date_col_name = "Date", df_sf = "sf", out_dir = NULL)

# Now dynamic...
# Not sure the best process, but maybe calculating the "seasonal" covariates first will be easiest
# To do the matching for seasons, need some type of look up table.
month_season_table <- data.frame("Month" = str_pad(seq(from = 1, to = 12, by = 1), 2, "left", 0), "Season" = c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer", "Summer", "Summer", "Fall", "Fall", "Fall", "Winter"))

dynamic_files <- list.files(glorys_dir_path, pattern = ".nc$", full.names = TRUE)[c(1, 5)]
dynamic_stacks<- vector("list", length(dynamic_files))

for (i in seq_along(dynamic_files)) {
  # Read it in
  rast_stack_temp <- raster::stack(dynamic_files[[i]])

  # Names to dates, and season match to get indices
  stack_dates <- as.Date(gsub("X", "", gsub("[.]", "-", paste(names(rast_stack_temp), ".16", sep = ""))))
  stack_season <- month_season_table$Season[match(format(stack_dates, "%m"), month_season_table$Month)]
  stack_season_years <- paste(format(stack_dates, format = "%Y"), stack_season, sep = "-")
  indices <- as.numeric(factor(stack_season_years, levels = unique(stack_season_years)))

  # Calculate mean across indices
  rast_agg_temp <- stackApply(rast_stack_temp, indices = indices, fun = mean, na.rm = TRUE)

  # Going to need the names of the stack to be actual dates
  dates <- as.Date(as.character(sapply(unique(stack_season_years), season_to_date)))
  names(rast_agg_temp) <- dates

  # Save it
  dynamic_stacks[[i]] <- rast_agg_temp
}

names(dynamic_stacks)<- c("BT_seasonal", "SST_seasonal")

# Run dynamic_2d_extract_wrapper function
all_tows_with_all_covs<- dynamic_2d_extract_wrapper(dynamic_covariates_list = dynamic_stacks, t_summ = "match", t_position = NULL, interp_nas = FALSE, sf_points = all_tows_with_static_covs, date_col_name = "Date", df_sf = "df", out_dir = NULL)

all_tows_with_all_covs <- all_tows_with_all_covs |>
  rename(c(BT_seasonal = BT_seasonal_match, SST_seasonal = SST_seasonal_match))

write_rds(all_tows_with_all_covs, here::here("Data/Derived/pred_glorys_with_covs.rds"), compress = "gz")



#---------------------------------------------
## Get covariates for a given mesh to help with predictions
#---------------------------------------------

# Mesh coordinates first
mesh_coords <- data.frame(sdmTMB_mesh$mesh$loc[, 1:2])
names(mesh_coords)<- c("DECDEG_BEGLON", "DECDEG_BEGLAT")

# We are going to want to have this for every time step...
env_data <- readRDS(here::here("Data/Derived/all_tows_all_covs.rds"))
str(env_data)

seasons <- c("Spring", "Fall", "Summer")
years <- sort(unique(env_data$EST_YEAR))

mesh_pred <- replicate_df(mesh_coords, time_name = "season", time_values = seasons) |>
  replicate_df(time_name = "year", time_values = years) |>
  mutate(
    "Year_Season" = paste(year, season, sep = "-"),
    "Date" = as.Date(season_to_date(Year_Season))
  )

# Read in static raster layer files from a directory and store them in a list. 
static_dir<- paste0(proj_box_path, "/covariates/static")
static_files<- list.files(static_dir, pattern = ".grd$", full.names = TRUE)
static_layers<- vector("list", length(static_files))

for(i in seq_along(static_files)){
  static_layers[[i]]<- raster(static_files[[i]])
}

# Names?
names(static_layers)[i]<- "Depth"

# Get all_tows into sf object
all_tows_sf<- points_to_sf(mesh_pred)

# Run static_extract_wrapper
all_tows_with_static_covs <- static_extract_wrapper(static_covariates_list = static_layers, sf_points = all_tows_sf, interp_nas = FALSE, date_col_name = "Date", df_sf = "sf", out_dir = NULL)

# Now dynamic...
# Not sure the best process, but maybe calculating the "seasonal" covariates first will be easiest
# To do the matching for seasons, need some type of look up table.
month_season_table <- data.frame("Month" = str_pad(seq(from = 1, to = 12, by = 1), 2, "left", 0), "Season" = c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer", "Summer", "Summer", "Fall", "Fall", "Fall", "Winter"))

dynamic_files <- list.files(glorys_dir_path, pattern = ".nc$", full.names = TRUE)[c(1, 5)]
dynamic_stacks<- vector("list", length(dynamic_files))

for (i in seq_along(dynamic_files)) {
  # Read it in
  rast_stack_temp <- raster::stack(dynamic_files[[i]])

  # Names to dates, and season match to get indices
  stack_dates <- as.Date(gsub("X", "", gsub("[.]", "-", paste(names(rast_stack_temp), ".16", sep = ""))))
  stack_season <- month_season_table$Season[match(format(stack_dates, "%m"), month_season_table$Month)]
  stack_season_years <- paste(format(stack_dates, format = "%Y"), stack_season, sep = "-")
  indices <- as.numeric(factor(stack_season_years, levels = unique(stack_season_years)))

  # Calculate mean across indices
  rast_agg_temp <- stackApply(rast_stack_temp, indices = indices, fun = mean, na.rm = TRUE)

  # Going to need the names of the stack to be actual dates
  dates <- as.Date(as.character(sapply(unique(stack_season_years), season_to_date)))
  names(rast_agg_temp) <- dates

  # Save it
  dynamic_stacks[[i]] <- rast_agg_temp
}

names(dynamic_stacks)<- c("BT_seasonal", "SST_seasonal")

# Run dynamic_2d_extract_wrapper function
all_tows_with_all_covs<- dynamic_2d_extract_wrapper(dynamic_covariates_list = dynamic_stacks, t_summ = "match", t_position = NULL, interp_nas = FALSE, sf_points = all_tows_with_static_covs, date_col_name = "Date", df_sf = "df", out_dir = NULL)

all_tows_with_all_covs <- all_tows_with_all_covs |>
  rename(c(BT_seasonal = BT_seasonal_match, SST_seasonal = SST_seasonal_match))

write_rds(all_tows_with_all_covs, here::here("Data/Derived/pred_mesh_with_covs.rds"), compress = "gz")
