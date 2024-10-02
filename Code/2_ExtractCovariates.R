library(devtools)
library(tidyverse)
library(stringr)
library(sf)
library(raster)
library(lubridate)
library(gmRi)

points_to_sf<- function(points){
  sf_out<- st_as_sf(points, coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
  return(sf_out)
}

all_tows <- readRDS("Data/Derived/for_covariate_extraction.rds")

tow_name_check<- c("ID", "DATE", "EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON")
all(tow_name_check %in% names(all_tows))

# Paths to Box folder
proj_box_path<- cs_path(box_group = "Mills Lab", subfolder = "Projects/sdm_workflow/data 2")

#####
## Static covariates
#####
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
all_tows_with_static_covs<- static_extract_wrapper(static_covariates_list = static_layers, sf_points = all_tows_sf, date_col_name = "DATE", df_sf = "sf", out_dir = "Data/Derived/")

# Check it out
summary(all_tows_with_static_covs)

#####
## Dynamic covariates
#####
# Read in dynamic raster stack files from a directory and store them in a list. 
dynamic_dir<- paste0(proj_box_path, "covariates/dynamic")
dynamic_files<- list.files(dynamic_dir, pattern = "BT.grd$", full.names = TRUE)
dynamic_stacks<- vector("list", length(dynamic_files))

for(i in seq_along(dynamic_files)){
  dynamic_stacks[[i]]<- raster::stack(dynamic_files[[i]])
}

# Names?
names(dynamic_stacks)<- c("BT")

# Run dynamic_2d_extract_wrapper function
all_tows_with_all_covs<- dynamic_2d_extract_wrapper(dynamic_covariates_list = dynamic_stacks, t_summ = "seasonal", t_position = NULL, sf_points = all_tows_with_static_covs, date_col_name = "DATE", df_sf = "df", out_dir = "Data/Derived/")

# Check it out
summary(all_tows_with_all_covs)

all_tows_with_all_covs %>% 
  filter(is.na(BT_seasonal)) %>% 
  group_by(EST_YEAR) %>%
  summarize(n = n()) %>% View() # Going to need to update the BT product. 
  

unique(all_tows_with_all_covs$EST_YEAR[is.na(all_tows_with_all_covs$BT_seasonal) == T])
