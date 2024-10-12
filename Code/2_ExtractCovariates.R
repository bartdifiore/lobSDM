library(devtools)
library(tidyverse)
library(stringr)
library(sf)
library(raster)
library(lubridate)
library(gmRi)
source("Code/enhance_r_funcs.R")

file_path_sans_ext <- function(x, compression = FALSE) {
  if (compression)
    x <- sub("[.](gz|bz2|xz)$", "", x)
  sub("([^.]+.+)\\.[[:alnum:]]+$", "\\1", x)
}

points_to_sf<- function(points){
  sf_out<- st_as_sf(points, coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
  return(sf_out)
}

all_tows <- readRDS("Data/Derived/for_covariate_extraction.rds")

tow_name_check<- c("ID", "DATE", "EST_YEAR", "DECDEG_BEGLAT", "DECDEG_BEGLON")
all(tow_name_check %in% names(all_tows))

# Paths to Box folder
proj_box_path <- cs_path(box_group = "Mills Lab", subfolder = "Projects/sdm_workflow/data 2")
glorys_dir_path<- cs_path(box_group = "RES_Data", subfolder = "GLORYs/NW_Atl_MonthlyTemps")

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
dynamic_files <- list.files(glorys_dir_path, pattern = ".nc$", full.names = TRUE)[c(1, 5)]
dynamic_stacks<- vector("list", length(dynamic_files))

for(i in seq_along(dynamic_files)){
  dynamic_stacks[[i]]<- raster::stack(dynamic_files[[i]])
}

# Names?
names(dynamic_stacks)<- c("BT", "SST")

# Run dynamic_2d_extract_wrapper function
all_tows_with_all_covs<- dynamic_2d_extract_wrapper(dynamic_covariates_list = dynamic_stacks, t_summ = "seasonal", t_position = NULL, sf_points = all_tows_with_static_covs, date_col_name = "DATE", df_sf = "df", out_dir = "Data/Derived/")

# Check it out
summary(all_tows_with_all_covs)

# Filter to GLORYs years (1993-2024)
tow_data_out <- all_tows_with_all_covs |>
  filter(EST_YEAR > 1993 & EST_YEAR < 2024)
summary(tow_data_out)

# Not entirely clear to me why we are missing some during the "good" years...potentially insure points?
t <- tow_data_out[which(is.na(tow_data_out$SST_seasonal) == T), ] |>
  distinct()

plot(t$DECDEG_BEGLON, t$DECDEG_BEGLAT) #

saveRDS(t, here::here("Data/Derived/all_tows_all_covs.rds"))
