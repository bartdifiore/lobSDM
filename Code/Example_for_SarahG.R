#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)

# Scaling/unscaling function to facilitate model convergence
set.seed(13)
x <- rnorm(100)

scaled_x <- scale(x)
unscaled_x <- as.numeric((scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center"))
all.equal(x, unscaled_x)

unscale <- function(scaled_x, center, scale){
  if(is.null(attr(scaled_x, "scaled:scale")) == F){
    #(scaled_x * sd) + m
    (scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center")
  }
  if(is.null(attr(scaled_x, "scaled:scale")) == T){
    (scaled_x * scale) + center
  }
}

#----------------------------------
## Get data
#----------------------------------

#during GLORYs time series
year_min <- 1993
year_max <- 2023

lobster<- readRDS(here::here("Data/Derived/all_model_data_juvenile.rds"))|>
  filter(survey == "NEFSC")|>
  filter(between(year, year_min, year_max)) |>
  drop_na(total_biomass)
summary(lobster)

predators <- readRDS(here::here("Data/Derived/all_model_data_predators.rds"))|>
  filter(survey == "NEFSC")|>
  filter(between(year, year_min, year_max)) |>
  drop_na(total_biomass)
summary(predators)



#----------------------------------------------------------------
## Model Data Prep for fitting seasonal model to lobster
#----------------------------------------------------------------

# Going to want to have a continuous time column
all_years<- seq(from = min(lobster$year), to = max(lobster$year))
seasons<- c("Spring", "Fall")
time_fac_levels<- paste(rep(all_years, each = length(unique(seasons))), seasons, sep = "_")
time_ints<- as.numeric(factor(time_fac_levels, levels = time_fac_levels))

lobster <- lobster |>
  mutate(season = factor(season, levels = seasons),
         year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
         year_season_int = as.numeric(year_season_fac)) %>%
  arrange(year_season_int)

# Now for the model matrix trickery
mm_season <- model.matrix(~ 0 + factor(season), data = lobster)
# mm_year <- model.matrix(~ 0 + factor(est_year), data = dat)

lobster <- lobster |>
  select(!contains("factor")) |>
  cbind(mm_season) |>
  #   cbind(mm_year) |>
  as_tibble()

fa <- names(lobster)[grepl("factor", names(lobster))]
fo <- paste0("`", paste(fa, collapse = "` + `"), "`")
svc <- as.formula(paste("~", fo))

# Check
svc

#----------------------------------
## Make mesh
#----------------------------------
# This is definitely something we will want to come back to after we have gotten things up and running through to model inferences. 
lobster <- lobster %>%
  sdmTMB::add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") 

# Create sdmTMB mesh -- check out Owen Liu's code here (https://github.com/owenrliu/eDNA_eulachon/blob/63a1b4d21fa4ffbc629cbb0657bc032998565f17/scripts/eulachon_sdms.Rmd#L217) for more ideas? This was taking forever, and no idea why...
# sdmTMB_mesh<- sdmTMB::make_mesh(lobster, xy_cols = c("longitude", "latitude"), type = "cutoff", cutoff = 100, fmesher_func = fmesher::fm_mesh_2d_inla)
sdmTMB_mesh <- sdmTMB::make_mesh(lobster, xy_cols = c("longitude", "latitude"), n_knots = 200, type = "kmeans")
# sdmTMB_mesh<- readRDS("~/Desktop/mesh_20241026_170037.rds")
plot(sdmTMB_mesh)


tictoc::tic()
fit5 <- sdmTMB(
  total_biomass ~ factor(season) + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  data = lobster,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = TRUE,
  share_range = FALSE,
  spatiotemporal = "ar1",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1",
  extra_time = c(55, 56),
  mesh = sdmTMB_mesh,
  family = tweedie(link = "log"),
  silent = FALSE,
  do_fit = TRUE
)
tictoc::toc()

fit5 <- readRDS("Data/Derived/model_fit5.rds")

tidy(fit5, effects = "fixed")
tidy(fit5, effects = "ran_pars")
fit5

write_rds(fit5, "Data/Derived/model_fit5_lobster.rds")


tictoc::tic()
fit_bd_lob <- sdmTMB(
  total_biomass ~ season + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  data = lobster,
  spatial = "on",
  # offset = dat_mod$swept,
  #anisotropy = TRUE,
  #share_range = FALSE,
  spatiotemporal = "iid",
  # spatial_varying = svc,
  time = "year_season_int",
  # time_varying = ~ 0 + year_season_int,
  # time_varying_type = "ar1",
  # extra_time = c(55, 56),
  mesh = sdmTMB_mesh,
  family = tweedie(link = "log"),
  silent = FALSE,
  do_fit = TRUE
)
tictoc::toc()

# fit_bd_lob <- readRDS("Data/Derived/model_fit_bd_lob.rds")

sanity(fit_bd_lob)
tidy(fit_bd_lob, effects = "fixed")
tidy(fit_bd_lob, effects = "ran_pars")
fit_bd_lob

write_rds(fit_bd_lob, "Data/Derived/model_fit_bd_lob.rds")



#----------------------------------------------------------------
## Model Data Prep for fitting seasonal model to predators
#----------------------------------------------------------------

# Going to want to have a continuous time column
all_years<- seq(from = min(predators$year), to = max(predators$year))
seasons<- c("Spring", "Fall")
time_fac_levels<- paste(rep(all_years, each = length(unique(seasons))), seasons, sep = "_")
time_ints<- as.numeric(factor(time_fac_levels, levels = time_fac_levels))

predators <- predators |>
  mutate(season = factor(season, levels = seasons),
         year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
         year_season_int = as.numeric(year_season_fac)) %>%
  arrange(year_season_int)

# Now for the model matrix trickery
mm_season <- model.matrix(~ 0 + factor(season), data = predators)
# mm_year <- model.matrix(~ 0 + factor(est_year), data = dat)

predators <- predators |>
  select(!contains("factor")) |>
  cbind(mm_season) |>
  #   cbind(mm_year) |>
  as_tibble()

fa <- names(predators)[grepl("factor", names(predators))]
fo <- paste0("`", paste(fa, collapse = "` + `"), "`")
svc <- as.formula(paste("~", fo))

# Check
svc

#---------------------------------------
## Make mesh and fit model to predators
#---------------------------------------
# This is definitely something we will want to come back to after we have gotten things up and running through to model inferences. 
predators <- predators %>%
  sdmTMB::add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") 

sdmTMB_mesh <- sdmTMB::make_mesh(predators, xy_cols = c("longitude", "latitude"), n_knots = 200, type = "kmeans")
plot(sdmTMB_mesh)


tictoc::tic()
fit_pred <- sdmTMB(
  total_biomass ~ season + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  data = predators,
  spatial = "on",
  # offset = dat_mod$swept,
  #anisotropy = TRUE,
  #share_range = FALSE,
  spatiotemporal = "iid",
  # spatial_varying = svc,
  time = "year_season_int",
  # time_varying = ~ 0 + year_season_int,
  # time_varying_type = "ar1",
  # extra_time = c(55, 56),
  mesh = sdmTMB_mesh,
  family = tweedie(link = "log"),
  silent = FALSE,
  do_fit = TRUE
)
tictoc::toc()

sanity(fit_pred)
tidy(fit_pred, effects = "fixed")
tidy(fit_pred, effects = "ran_pars")
fit_pred

write_rds(fit_pred, "Data/Derived/model_fit_predators.rds")


#---------------------------------------
## Estimate overlap metrics and plot
#---------------------------------------

glorys <- readRDS("Data/Derived/glorys_grid.rds")
source("Code/Carroll_2019_functions.R")


