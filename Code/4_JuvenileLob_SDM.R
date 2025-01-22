#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sdmTMB)
library(sdmTMBextra)
library(ggeffects)

# Base map land info
region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

#----------------------------------
## Get juvenile data
#----------------------------------

# Focus on just lobster, and during GLORYs time series
year_min <- 1993
year_max <- 2023

df<- readRDS(here::here("Data/Derived/all_model_data_juvenile.rds")) %>%
  dplyr::filter(between(year, year_min, year_max))

# all(all_mod_data$trawl_id %in% env_tows)
# all(env_tows %in% all_mod_data$trawl_id)
# Every tow should have two observations and the unique tows should be the same as the ones in the environmental data
# t<- table(all_mod_data$trawl_id)
# length(unique(all_mod_data$trawl_id)) == length(unique(env_data$ID))

df %>% 
  group_by(survey) %>% 
  count()


# Still some weird NA biomass values to figure out, dropping those for now
mod_data<- df |>
  drop_na(total_biomass)
summary(mod_data)

# What the heck is going on with that 2014 value?
t<- mod_data[which.max(mod_data$total_biomass),]
plot(mod_data$total_biomass)

# Doesn't seem like there is anyway that point is real.
mod_data <- mod_data[-which.max(mod_data$total_biomass), ] |>
  ungroup()
plot(mod_data$total_biomass)

# Scale/center covariates
# Get means and sds
column_means <- colMeans(mod_data[, c("Depth", "BT_seasonal")], na.rm = TRUE)
column_sds <- apply(mod_data[, c("Depth", "BT_seasonal")], 2, sd, na.rm = TRUE)

# Scale the data
mod_data <- mod_data |>
  mutate(across(
    c(Depth, BT_seasonal),
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

# Create sdmTMB mesh -- check out Owen Liu's code here (https://github.com/owenrliu/eDNA_eulachon/blob/63a1b4d21fa4ffbc629cbb0657bc032998565f17/scripts/eulachon_sdms.Rmd#L217) for more ideas? This was taking forever, and no idea why...
# sdmTMB_mesh<- sdmTMB::make_mesh(mod_data, xy_cols = c("longitude", "latitude"), type = "cutoff", cutoff = 100, fmesher_func = fmesher::fm_mesh_2d_inla)
sdmTMB_mesh <- sdmTMB::make_mesh(mod_data, xy_cols = c("X", "Y"), n_knots = 200, type = "kmeans")

sdmTMB_mesh <- sdmTMB::make_mesh(mod_data, xy_cols = c("X", "Y"),cutoff = 100)
# sdmTMB_mesh<- readRDS("~/Desktop/mesh_20241026_170037.rds")
plot(sdmTMB_mesh)

#----------------------------------
## Fit model 0
#----------------------------------
#| label: Basic environment-only SDM
#| include: false
#| echo: false
#| warning: false

# Starting with the most simple model that just includes the environmental covariates.
fit0<- sdmTMB(
  total_biomass ~ s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
  data = mod_data,
  spatial = "off",
  # offset = dat_mod$swept, 
  # anisotropy = FALSE, 
  # share_range = TRUE,
  spatiotemporal = "off",
  spatial_varying = NULL,
  time = NULL,
  time_varying = NULL,
  # time_varying_type = "ar1",
  extra_time = NULL,
  mesh = sdmTMB_mesh,
  family = tweedie(), 
  # family = gengamma(link = "log"),
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit0, model = 1, effects = "fixed") # Intercept

###
# Diagnostics
###
sanity(fit0)

# Residuals (https://pbs-assess.github.io/sdmTMB/articles/residual-checking.html)
fit0_sims <- simulate(fit0, nsim = 500, type = "mle-mvn")
fit0_resids <- dharma_residuals(fit0_sims, fit0, return_DHARMa = TRUE)
plot(fit0_resids)

DHARMa::testResiduals(fit0_resids)
DHARMa::testZeroInflation(fit0_resids)

###
# Inferences
###
# Where are the smooths?
names(fit0)
names(fit0$smoothers) # All in here...
str(fit0$smoothers$basis_out) # List, one element for each smooth term, within each list element all the information of the fitted smooth is included

# Plot them
plot <- plot_smooths(mod_fit = fit0, n_preds = 200, y_lab = "Predicted biomass", rescale_means = column_means, rescale_sds = column_sds)
plot

###
# Predictions
###
# I have already done this in the `2_extract` code, though if time, probably a good idea to figure out how to make some helper functions to do this quickly. For now, just reading in that prediction dataframe that I already created. Again, this will have depth/BT_seasonal for every year-season time step.
pred_data <- readRDS(here::here("Data/Derived/pred_glorys_with_covs.rds"))

# Scale appropriately...want to use the mean/sd from the original scaleing
pred_data <- pred_data |>
  filter(between(year, year_min, year_max) & season %in% c("Spring", "Summer", "Fall")) |>
  mutate(
    Depth_scaled = (Depth - column_means["Depth"]) / column_sds["Depth"],
    BT_seasonal_scaled = (BT_seasonal - column_means["BT_seasonal"]) / column_sds["BT_seasonal"]
  ) |>
  drop_na()

# Predict
fit0_preds <- predict(fit0, newdata = pred_data, type = "response", se = FALSE, return_tmb_object = TRUE)$data
str(fit0_preds)

# Nest and map
fit0_preds <- fit0_preds |>
  group_by(season, year, Year_Season, Date) |>
  nest()

map_nested <- function(pred_df, time, region_use = region, states_use = states, xlim_use = lon_lims, ylim_use = lat_lims) {
  ggplot() +
    geom_raster(data = pred_df, aes(x = x, y = y, fill = est)) +
    geom_sf(data = region_use, fill = "#f0f0f0") +
    geom_sf(data = states_use, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_fill_viridis_c(name = "Predicted biomass") +
    theme_minimal() +
    labs(
      fill = "Predicted biomass",
      title = time
    )
}

fit0_preds <- fit0_preds |>
  mutate("Pred_Map" = map2(data, Year_Season, map_nested))

fit0_preds$Pred_Map[[50]]
fit0_preds$Pred_Map[[36]]
fit0_preds$Pred_Map[[93]]

#----------------------------------
## Fit model 1
#----------------------------------

#| label: Env-only with season and catchability factors
#| include: false
#| echo: false
#| warning: false

fit1<- sdmTMB(
  total_biomass ~ factor(season) + factor(survey) + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  data = mod_data,
  spatial = "off",
  # offset = dat_mod$swept, 
  # anisotropy = TRUE, 
  # share_range = TRUE,
  spatiotemporal = "off",
  spatial_varying = NULL,
  time = NULL,
  time_varying = NULL,
  # time_varying_type = "ar1",
  extra_time = NULL,
  mesh = sdmTMB_mesh,
  family = tweedie(), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit1, model = 1, effects = "fixed")


#----------------------------------
## Fit model 2
#----------------------------------

#| label: Base model with timevarying intercept that has an AR1 structure
#| include: false
#| echo: false
#| warning: false

# Importantly, here we are going to switch the "extra_time" argument to be a continuous vector of all year-season combinations. This will help fill in missing gaps.

fit2<- sdmTMB(
  total_biomass ~ + factor(season) + factor(survey) + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  data = mod_data,
  spatial = "off",
  # offset = dat_mod$swept, 
  # anisotropy = TRUE, 
  # share_range = TRUE,
  spatiotemporal = "off",
  spatial_varying = NULL,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1", # RW default
  extra_time = time_ints,
  mesh = sdmTMB_mesh,
  family = tweedie(), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit2, model = 1, effects = "fixed") # 

# Where are the time varying intercepts?
# Make a little helpfer function
get_time_var_est<- function(fit){
  est<- as.list(fit$sd_report, "Estimate")$b_rw_t
  se<- as.list(fit$sd_report, "Std. Error")$b_rw_t
  
  delta<- ifelse(length(dim(est)) == 3 & dim(est)[3] == 2, TRUE, FALSE)
  
  if(delta){
    df1 <- data.frame("Model" = 1, "Time" = fit$time_lu$time_from_data, "Est" = est[,,1], "SE" = se[,,1])
    df2<- data.frame("Model" = 2, "Time" = fit$time_lu$time_from_data, "Est" = est[,,2], "SE" = se[,,2])
    df_out<- df1 |>
      bind_rows(df2) |>
      mutate("pct025" = Est + qnorm(0.025) * SE,
             "pct975" = Est + qnorm(0.975) * SE)
    return(df_out)
  } else {
    df1 <- data.frame("Model" = 1, "Time" = fit$time_lu$time_from_data, "Est" = est[,1,1], "SE" = se[,1,1])
    df_out<- df1 |>
      mutate("pct025" = Est + qnorm(0.025) * SE,
             "pct975" = Est + qnorm(0.975) * SE)
    return(df_out)
  }
}

time_var_est<- get_time_var_est(fit2)
est <- as.list(fit2$sd_report, "Estimate")
se <- as.list(fit2$sd_report, "Std. Error")

ggplot() + 
  geom_errorbar(data = time_var_est, aes(x = Time, y = Est, ymin = Est - SE, ymax = Est + SE)) + 
  facet_wrap( ~ Model) +
  theme_light()


#----------------------------------
## Fit model 3
#----------------------------------
#| label: Adding in seasonal SVC
#| include: false
#| echo: false
#| warning: false

n_seasons <- length(unique(mod_data$season)) 
# tau_Z_map <- factor(cbind(rep(1, n_seasons), rep(2, n_seasons))) # Delta model, two LPs
tau_Z_map <- factor(cbind(rep(1, n_seasons)))

fit3<- sdmTMB(
  total_biomass ~ factor(season) + factor(survey) + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data,
  spatial = "off",
  # offset = dat_mod$swept, 
  # anisotropy = TRUE, 
  # share_range = TRUE,
  spatiotemporal = "off",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1",
  extra_time = time_ints,
  mesh = sdmTMB_mesh,
  family = tweedie(), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit3, effects = "fixed")
tidy(fit3, effects = "ran_pars") # Now have range, and three sigma_Z's, representing seasonal variability.
fit3
sanity(fit3)

write_rds(fit3, here::here("Juve_Fit3.rds"), compress = "gz")

#----------------------------------
## Fit model 4
#----------------------------------
#| label: Adding in persistent spatial variation
#| include: false
#| echo: false
#| warning: false


fit4<- sdmTMB(
  total_biomass ~ factor(season) + factor(survey) + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data,
  spatial = "on",
  # offset = dat_mod$swept, 
  anisotropy = TRUE, 
  # share_range = TRUE,
  # spatiotemporal = "ar1",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1",
  extra_time = time_ints,
  mesh = sdmTMB_mesh,
  family = tweedie(), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit4, effects = "fixed")
tidy(fit4, effects = "ran_pars")
fit4
sanity(fit4) # Smooth sigma issue!

write_rds(fit4, here::here("Juve_Fit4.rds"), compress = "gz")

#----------------------------------
## Fit model 5
#----------------------------------
#| label: Final addition, spatio-temporal variation
#| include: false
#| echo: false
#| warning: false

tictoc::tic()
fit5 <- sdmTMB(
  total_biomass ~ factor(season) + factor(survey) + s(Depth, k = 4) + s(BT_seasonal, k = 4),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = TRUE,
  share_range = FALSE,
  spatiotemporal = "ar1",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1",
  extra_time = time_ints,
  mesh = sdmTMB_mesh,
  family = tweedie(),
  silent = FALSE,
  do_fit = TRUE
)
tictoc::toc()

tidy(fit5, effects = "fixed")
tidy(fit5, effects = "ran_pars")
sanity(fit5)

write_rds(fit5, here::here("Juve_Fit5.rds"), compress = "gz")

#----------------------------------
## Fit model 6
#----------------------------------
#| label: Exclude space-invariant temporal autocorrelation.
#| include: false
#| echo: false
#| warning: false

#|

tictoc::tic()
fit6 <- sdmTMB(
  total_biomass ~ factor(season) + factor(survey) + s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = TRUE,
  share_range = FALSE,
  spatiotemporal = "ar1",
  spatial_varying = svc,
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

tidy(fit6, effects = "fixed")
tidy(fit6, effects = "ran_pars")
sanity(fit6)

write_rds(fit6, here::here("Juve_Fit6.rds"), compress = "gz")

# # Predict
# 
# pred_data <- readRDS(here::here("Data/Derived/pred_glorys_with_covs.rds"))
# 
# # Now for the model matrix trickery
# mm_season <- model.matrix(~ 0 + factor(season), data = pred_data)
# # mm_year <- model.matrix(~ 0 + factor(est_year), data = dat) 
# 
# pred_data <- readRDS(here::here("Data/Derived/pred_glorys_with_covs.rds")) |>
#   filter(between(year, year_min, year_max) & season %in% c("Spring", "Summer", "Fall")) |>
#   mutate(
#     Depth_scaled = (Depth - column_means["Depth"]) / column_sds["Depth"],
#     BT_seasonal_scaled = (BT_seasonal - column_means["BT_seasonal"]) / column_sds["BT_seasonal"]
#   ) |>
#   mutate(season = factor(season, levels = seasons),
#          year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
#          year_season_int = as.numeric(year_season_fac)) %>%
#   dplyr::select(!contains("factor")) %>%
#   cbind(mm_season) %>%
#   as_tibble() %>%
#   add_utm_columns(ll_names = c("x", "y")) %>%
#   mutate(survey = "ME_NH") %>%
#   drop_na()
# 
# fit6_preds <- predict(fit6, newdata = pred_data, type = "response", se = FALSE)
# str(fit6_preds)
# 
# # Nest and map
# fit6_preds <- fit6_preds |>
#   group_by(season, year, Year_Season, Date) |>
#   nest()
# 
# map_nested <- function(pred_df, time, region_use = region, states_use = states, xlim_use = lon_lims, ylim_use = lat_lims) {
#   ggplot() +
#     geom_raster(data = pred_df, aes(x = x, y = y, fill = est)) +
#     geom_sf(data = region_use, fill = "#f0f0f0") +
#     geom_sf(data = states_use, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
#     coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
#     scale_fill_viridis_c(name = "Predicted biomass", trans = "log") +
#     theme_minimal() +
#     labs(
#       fill = "Predicted biomass",
#       title = time
#     )
# }
# 
# fit6_preds <- fit6_preds |>
#   mutate("Pred_Map" = map2(data, Year_Season, map_nested))
# 
# cowplot::plot_grid(fit6_preds$Pred_Map[[1]], 
#                    fit6_preds$Pred_Map[[10]], 
#                    fit6_preds$Pred_Map[[20]], 
#                    fit6_preds$Pred_Map[[30]], 
#                    fit6_preds$Pred_Map[[40]], 
#                    fit6_preds$Pred_Map[[50]], 
#                    fit6_preds$Pred_Map[[60]],
#                    fit6_preds$Pred_Map[[70]], 
#                    fit6_preds$Pred_Map[[80]], 
#                    fit6_preds$Pred_Map[[90]], ncols = 3, nrows = 3)
# 









