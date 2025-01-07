#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)
library(sf)

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
sdmTMB_mesh <- sdmTMB::make_mesh(lobster, xy_cols = c("X", "Y"), n_knots = 200, type = "kmeans")
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

# fit5 <- readRDS("Data/Derived/model_fit5.rds")

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

sdmTMB_mesh <- sdmTMB::make_mesh(predators, xy_cols = c("X", "Y"), n_knots = 200, type = "kmeans")
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

lob_mod <- readRDS("Data/Derived/model_fit_bd_lob.rds")
pred_mod <- readRDS("Data/Derived/model_fit_predators.rds")

source("Code/Carroll_2019_functions.R")

temp <- lobster %>%
  ungroup() %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_geometry()

plot(st_buffer(st_concave_hull(st_union(temp), ratio = 0.01),dist = 1/12), col = "red") 
plot(temp, cex = 0.1, add = T)

clipper <- st_buffer(st_concave_hull(st_union(temp), ratio = 0.1),dist = 1/12) %>%
  st_as_sf() %>%
  mutate(formerge = "in")

df_predictions <- readRDS("Data/Derived/pred_glorys_with_covs.rds") %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_set_crs(4326) %>%
  st_transform(crs = 32619) %>%
  filter(year >= 1993 & year < 2020, season %in% c("Spring", "Fall")) %>%
  mutate(year_season_int = as.numeric(as.factor(Year_Season)), 
         X = as.numeric(st_coordinates(.)[,1])/1000, 
         Y = as.numeric(st_coordinates(.)[,2])/1000) %>%
  drop_na(Depth) %>%
  drop_na(BT_seasonal) %>%
  st_drop_geometry()
  
grid <- readRDS("Data/Derived/pred_glorys_with_covs.rds") %>%
  select(x,y) %>%
  distinct() %>%
  sf::st_as_sf(coords = c("x", "y")) %>%
  st_set_crs(4326) %>%
  st_make_grid(cellsize = 1/12) %>%
  st_as_sf() %>%
  st_join(clipper, join = st_within) %>%
  filter(formerge == "in") %>%
  mutate(id = 1:n())

plot(grid)

lob_predictions <- predict(lob_mod, newdata = df_predictions, type = "response") %>%
  mutate(X = X*1000,
         Y = Y*1000) %>%
  sf::st_as_sf(coords = c("X", "Y")) %>%
  st_set_crs(32619) %>%
  st_transform(crs = 4326)

grid_lob_predictions <- st_join(grid, lob_predictions)

grid_lob_predictions %>%
  filter(year %in% c(1995, 2000, 2005, 2010, 2015, 2019)) %>%
  drop_na(year, est) %>%
  mutate(log_est = log(est)) %>%
  ggplot()+
  geom_sf(aes(fill = est), color = "transparent")+
  scale_fill_gradient(high = "yellow", low = "blue", trans = "log")+
  facet_grid(season~year)

ggsave("Figures/predicted_lobster.png")

predator_predictions <- predict(pred_mod, newdata = df_predictions, type = "response") %>%
  mutate(X = X*1000,
         Y = Y*1000) %>%
  sf::st_as_sf(coords = c("X", "Y")) %>%
  st_set_crs(32619) %>%
  st_transform(crs = 4326)

grid_predator_predictions <- st_join(grid, predator_predictions)

grid_predator_predictions %>%
  filter(year %in% c(1995, 2000, 2005, 2010, 2015, 2019)) %>%
  drop_na(year, est) %>%
  mutate(log_est = log(est)) %>%
  ggplot()+
  geom_sf(aes(fill = est), color = "transparent")+
  scale_fill_gradient(high = "yellow", low = "blue", trans = "log")+
  facet_grid(season~year)

ggsave("Figures/predicted_predators.png")


grid_predator_predictions %>%
  left_join(grid_lob_predictions %>% select(id, season, year, Year_Season, Date, Depth, BT_seasonal, SST_seasonal, year_season_int, est) %>% rename(lob_est = est) %>% st_drop_geometry()) %>%
  mutate(lob_to_pred_ratio = lob_est/ est) %>%
  filter(year %in% c(1995, 2000, 2005, 2010, 2015, 2019)) %>%
  drop_na(year, lob_to_pred_ratio) %>%
  ggplot()+
  geom_sf(aes(fill = lob_to_pred_ratio), color = "transparent")+
  scale_fill_gradient2(high = "red", low = "blue", mid = "white", trans = "log")+
  facet_grid(season~year)

ggsave("Figures/prey_to_pred_ratio.png")

df_overlap <- grid_predator_predictions %>%
  select(id, season, year, Depth, BT_seasonal, SST_seasonal, est) %>%
  rename(predator_est = est) %>%
  st_drop_geometry() %>%
  left_join(grid_lob_predictions %>% 
              select(id, season, year, Depth, BT_seasonal, SST_seasonal, est) %>%
              rename(lobster_est = est) %>%
              st_drop_geometry)

write_rds(df_overlap, "Data/Derived/overlap_df.rds", compress = "gz")

#-------------------------------------------
## Calculate overlap metrics
#-------------------------------------------

df <- readRDS("Data/Derived/overlap_df.rds") %>%
  left_join(grid) %>%
  st_as_sf() %>%
  st_join(ecodata::epu_sf %>% select(EPU) %>% st_transform(crs = 4326) %>% st_make_valid())


# temp <- df_overlap %>% 
#   filter(year == 2012, season == "Spring")
# 
# b_overlap <- biomass_overlapfn(temp$lobster_est, temp$predator_est)
# hist(b_overlap)
# 
# loc_col <- loc_collocfn(temp$lobster_est, temp$predator_est)


df %>%
  group_by(season, year, EPU) %>%
  st_drop_geometry() %>%
  drop_na(season) %>%
  drop_na(EPU) %>%
  # filter(EPU == "GOM") %>%
  summarize(b_overlap = biomass_overlapfn(lobster_est, predator_est), 
            loc_col = loc_collocfn(lobster_est, predator_est), 
            a_alpha = asymmalpha_overlapfn(lobster_est, predator_est), 
            ab = AB_overlapfn(lobster_est, predator_est), 
            bhatta = bhatta_coeffn(lobster_est, predator_est), 
            schoeners = schoeners_overlapfn(lobster_est, predator_est)) %>%
  pivot_longer(cols = b_overlap:schoeners) %>%
  ggplot(aes(x = year, y = value))+
  geom_line(aes(color = season))+
  facet_grid(name~EPU, scales = "free")
ggsave("Figures/overlap_metrics.png")


























