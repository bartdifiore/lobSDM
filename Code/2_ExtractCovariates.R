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
