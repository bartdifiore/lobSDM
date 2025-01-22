#----------------------------------
## Libraries
#----------------------------------

library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)


#----------------------------------
## Get data
#----------------------------------

# during GLORYs time series
year_min <- 1993
year_max <- 2023

predators <- readRDS(here::here("Data/Derived/all_model_data_predators.rds")) %>%
  dplyr::filter(between(year, year_min, year_max)) %>% 
  rename(predator_biomass = total_biomass) %>% 
  select(trawl_id, longitude, latitude, season, year, survey, date, predator_biomass)

lobster <- readRDS(here::here("Data/Derived/all_model_data_juvenile.rds")) %>%
  dplyr::filter(between(year, year_min, year_max))%>% 
  rename(lobster_biomass = total_biomass)

grid.2 <- predators %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>% 
  st_bbox() %>%
  st_as_sfc() %>%
  st_make_grid(cellsize = 2) %>%
  st_sf() %>%
  mutate(cell_id = 1:n()) %>%
  st_make_valid()

mod_data <- lobster %>% 
  left_join(predators) %>% 
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>%
  st_join(grid.2) %>% 
  st_drop_geometry() %>% 
  left_join(predators %>% select(trawl_id, longitude, latitude, season, year, survey, date))
  
  

# Scale/center covariates
# Get means and sds
column_means <- colMeans(mod_data[, c("Depth", "BT_seasonal", "predator_biomass")], na.rm = TRUE)
column_sds <- apply(mod_data[, c("Depth", "BT_seasonal", "predator_biomass")], 2, sd, na.rm = TRUE)

# Scale the data
mod_data <- mod_data |>
  ungroup() %>% # So in this version of the code i had to ungroup otherwise I got NA's. I didn't do this when prepping the juvenile lobster data. We need to make sure that the juvenile lobster code isn't only scaling within groups....
  mutate(across(
    c(Depth, BT_seasonal, predator_biomass),
    ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),
    .names = "{.col}_scaled"
  ))
summary(mod_data)

#----------------------------------
## Model Data Prep for fitting seasonal model
#----------------------------------

# Going to want to have a continuous time column
all_years<- seq(from = min(mod_data$year), to = max(mod_data$year))
seasons<- c("Spring", "Summer", "Fall")
time_fac_levels<- paste(rep(all_years, each = length(unique(seasons))), seasons, sep = "_")
time_ints<- as.numeric(factor(time_fac_levels, levels = time_fac_levels))

mod_data <- mod_data |>
  mutate(season = factor(season, levels = seasons),
         year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
         year_season_int = as.numeric(year_season_fac)) %>%
  arrange(year_season_int)

# Now for the model matrix trickery
mm_season <- model.matrix(~ 0 + factor(season), data = mod_data)
# mm_year <- model.matrix(~ 0 + factor(est_year), data = dat) 

mod_data <- mod_data |>
  dplyr::select(!contains("factor")) |>
  cbind(mm_season) |>
  #   cbind(mm_year) |>
  as_tibble()

fa <- names(mod_data)[grepl("factor", names(mod_data))]
fo <- paste0("`", paste(fa, collapse = "` + `"), "`")
svc <- as.formula(paste("~", fo))

# Check
svc

#----------------------------------
## Make mesh
#----------------------------------
# This is definitely something we will want to come back to after we have gotten things up and running through to model inferences. 
mod_data <- mod_data %>%
  sdmTMB::add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") 

sdmTMB_mesh <- sdmTMB::make_mesh(mod_data, xy_cols = c("X", "Y"), n_knots = 200, type = "kmeans")

sdmTMB_mesh <- sdmTMB::make_mesh(mod_data, xy_cols = c("X", "Y"),cutoff = 100)
# sdmTMB_mesh<- readRDS("~/Desktop/mesh_20241026_170037.rds")
plot(sdmTMB_mesh)

#----------------------------------
## Fit model 6
#----------------------------------
#| label: Exclude space-invariant temporal autocorrelation.
#| include: false
#| echo: false
#| warning: false



n_seasons <- length(unique(mod_data$season)) 
# tau_Z_map <- factor(cbind(rep(1, n_seasons), rep(2, n_seasons))) # Delta model, two LPs
tau_Z_map <- factor(cbind(rep(1, n_seasons)))

tictoc::tic()
fit_lp <- sdmTMB(
  lobster_biomass ~ factor(season) + factor(survey) + predator_biomass_scaled + s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
  # control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = TRUE,
  share_range = FALSE,
  spatiotemporal = "ar1",
  spatial_varying = ~predator_biomass_scaled,
  time = "year_season_int",
  # time_varying = ~ 0 + year_season_int,
  # time_varying_type = "ar1",
  extra_time = time_ints,
  mesh = sdmTMB_mesh,
  family = tweedie(),
  silent = FALSE,
  do_fit = TRUE
)
tictoc::toc()

tidy(fit_lp, effects = "fixed")
tidy(fit_lp, effects = "ran_pars")
sanity(fit_lp)

write_rds(fit_lp, here::here("Pred_Fit_lp.rds"), compress = "gz")


mod_data2 <- mod_data %>% 
  mutate(predator_biomass_log_scaled = scale(log(predator_biomass + 0.0001)))

tictoc::tic()
fit_lp2 <- sdmTMB(
  lobster_biomass ~ factor(season) + factor(survey) + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data2,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = TRUE,
  share_range = FALSE,
  spatiotemporal = "ar1",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ cell_id*predator_biomass_log_scaled,
  time_varying_type = "ar1",
  extra_time = time_ints,
  mesh = sdmTMB_mesh,
  family = tweedie(),
  silent = FALSE,
  do_fit = TRUE
)
tictoc::toc()

tidy(fit_lp2, effects = "fixed")
tidy(fit_lp2, effects = "ran_pars")
sanity(fit_lp2)

write_rds(fit_lp2, here::here("Pred_Fit_lp2.rds"), compress = "gz")


nd <- expand.grid(season = "Fall", 
                  survey = "NEFSC", 
                  predator_biomass_log_scaled = seq(min(mod_data2$predator_biomass_log_scaled, na.rm = T), max(mod_data2$predator_biomass_log_scaled, na.rm = T), length.out = 100),
                  # predator_biomass_log_scaled = quantile(mod_data2$predator_biomass_log_scaled, c(0.025, 0.5, 0.975)), 
                  Depth = mean(mod_data2$Depth), 
                  BT_seasonal = mean(mod_data2$BT_seasonal), 
                  year_season_int = unique(mod_data2$year_season_int), 
                  cell_id = unique(mod_data2$cell_id)) %>%
  mutate(`factor(season)Spring` = 0, 
         `factor(season)Summer` = 0, 
         `factor(season)Fall` = 1)

p <- predict(fit_lp2, newdata = nd, re_form = NA)

p %>% 
  filter(year_season_int %in% c(1, 10, 20, 30, 40, 50, 70)) %>%
  # filter(year_season_int %in% c(50)) %>%
  ggplot(aes(x = predator_biomass_log_scaled, y = exp(est)))+
  geom_line(aes(color = as.factor(year_season_int)), show.legend = T)+
  scale_y_log10()+
  facet_wrap(~cell_id, scales = "free")+
  labs(y = "Predicted lobster biomass", x = "predator_biomass_log_scaled")




n_seasons <- length(unique(mod_data$season)) 
# tau_Z_map <- factor(cbind(rep(1, n_seasons), rep(2, n_seasons))) # Delta model, two LPs
tau_Z_map <- factor(c(rep(1, n_seasons), 2))

mod_data3 <- mod_data2 %>% 
  ungroup() %>%
  mutate(time = as.numeric(case_when(season == "Spring" ~ paste(year, "25", sep = "."), 
                                 season == "Summer" ~ paste(year, "50", sep = "."), 
                                 season == "Fall" ~ paste(year, "75", sep = "."))), 
         scaled_time = (time - mean(time))/10) %>% 
  as.data.frame()

tictoc::tic()
fit_lp3 <- sdmTMB(
  lobster_biomass ~ factor(season) + factor(survey) + predator_biomass_log_scaled + s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
  # control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data3,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = TRUE,
  share_range = FALSE,
  spatiotemporal = "ar1",
  spatial_varying = ~ predator_biomass_log_scaled:scaled_time,
  time = "year_season_int",
  # time_varying = ~ 0 + year_season_int,
  # time_varying_type = "ar1",
  extra_time = time_ints,
  mesh = sdmTMB_mesh,
  family = tweedie(),
  silent = FALSE,
  do_fit = TRUE
)
tictoc::toc()


tidy(fit_lp3, effects = "fixed")
tidy(fit_lp3, effects = "ran_pars")
sanity(fit_lp3)


#---------------------------------
## Predict
#---------------------------------

# In order to have a continuous prediction surface, I need to use the predictions from the predator model to build the prediction grid for this model.

pred_mod <- readRDS("Data/Derived/Pred_Fit6.rds")

prediction_grid <- readRDS(here::here("Data/Derived/pred_glorys_with_covs.rds")) %>%
  filter(between(year, year_min, year_max) & season %in% c("Spring", "Summer", "Fall")) %>%
  rename(longitude = x, latitude = y)

# Going to want to have a continuous time column
all_years<- seq(from = min(prediction_grid$year), to = max(prediction_grid$year))
seasons<- c("Spring", "Summer", "Fall")
time_fac_levels<- paste(rep(all_years, each = length(unique(seasons))), seasons, sep = "_")
time_ints<- as.numeric(factor(time_fac_levels, levels = time_fac_levels))

# Now for the model matrix trickery
mm_season <- model.matrix(~ 0 + factor(season), data = prediction_grid)
# mm_year <- model.matrix(~ 0 + factor(est_year), data = dat) 

column_means <- colMeans(prediction_grid[, c("Depth", "BT_seasonal")], na.rm = TRUE)
column_sds <- apply(prediction_grid[, c("Depth", "BT_seasonal")], 2, sd, na.rm = TRUE)

prediction_grid <- prediction_grid |>
  mutate(
    Depth_scaled = (Depth - column_means["Depth"]) / column_sds["Depth"],
    BT_seasonal_scaled = (BT_seasonal - column_means["BT_seasonal"]) / column_sds["BT_seasonal"]
  ) |>
  mutate(season = factor(season, levels = seasons),
         year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
         year_season_int = as.numeric(year_season_fac)) %>%
  dplyr::select(!contains("factor")) %>%
  cbind(mm_season) %>%
  as_tibble() %>%
  add_utm_columns(ll_names = c("longitude", "latitude")) %>%
  mutate(survey = "NEFSC") %>%
  drop_na()

pred_yhat <- predict(pred_mod, newdata = prediction_grid)

prediction_grid$pred_yhat <- pred_yhat$est 
prediction_grid$predator_biomass_log_scaled <- scale(log(prediction_grid$pred_yhat + 0.0001))

prediction_grid <- prediction_grid %>%
  ungroup() %>%
  mutate(time = as.numeric(case_when(season == "Spring" ~ paste(year, "25", sep = "."), 
                                     season == "Summer" ~ paste(year, "50", sep = "."), 
                                     season == "Fall" ~ paste(year, "75", sep = "."))), 
         scaled_time = (time - mean(time))/10) %>% 
  as.data.frame() %>% 
  drop_na()

fitlp_preds <- predict(fit_lp3, newdata = prediction_grid)

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)



fitlp_preds %>%
  as_tibble() %>% 
  rename(interaction = `zeta_s_predator_biomass_log_scaled:scaled_time`) %>%
  filter(year == 2023, season == "Fall") %>% 
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = interaction)) +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  scale_fill_viridis_c() +
  theme_minimal()+
  labs(
    fill = "Interaction", 
    title = "Temporal effect of predators on juvenile lobster", x = "", y = ""
  )
ggsave("Figures/its_something.png")


nd <- expand.grid(season = "Fall", 
                  survey = "NEFSC", 
                  predator_biomass_log_scaled = seq(min(mod_data3$predator_biomass_log_scaled, na.rm = T), max(mod_data3$predator_biomass_log_scaled, na.rm = T), length.out = 100),
                  # predator_biomass_log_scaled = quantile(mod_data2$predator_biomass_log_scaled, c(0.025, 0.5, 0.975)), 
                  Depth_scaled = mean(mod_data3$Depth_scaled), 
                  BT_seasonal_scaled = mean(mod_data3$BT_seasonal_scaled), 
                  year_season_int = unique(mod_data3$year_season_int),
                  scaled_time = unique(mod_data3$scaled_time)
                  ) %>%
  mutate(`factor(season)Spring` = 0, 
         `factor(season)Summer` = 0, 
         `factor(season)Fall` = 1)

p <- predict(fit_lp3, newdata = nd, re_form = NA)

p %>% 
  filter(year_season_int == c(1, 50, 70)) %>%
ggplot(aes( x = predator_biomass_log_scaled, y = exp(est)))+
  geom_line(aes(color = as.factor(year_season_int)))







