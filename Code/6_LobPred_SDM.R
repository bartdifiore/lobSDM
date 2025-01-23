#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sdmTMB)
library(sdmTMBextra)
library(ggeffects)
library(patchwork)
source(here::here("Code/sdmTMB_validation_Functions.R"))

# Scaling/unscaling function to facilitate model convergence
set.seed(13)
x <- rnorm(100)

scaled_x <- scale(x)
unscaled_x <- as.numeric((scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center"))
all.equal(x, unscaled_x)

scaled<- function(x, center, scale){
  (x - center) / scale
}

unscale <- function(scaled_x, center, scale) {
  if (is.null(attr(scaled_x, "scaled:scale")) == F) {
    # (scaled_x * sd) + m
    (scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center")
  }
  if (is.null(attr(scaled_x, "scaled:scale")) == T) {
    (scaled_x * scale) + center
  }
}

unscale_aja <- function(scaled_x, orig_mean, orig_sd) {
  (scaled_x * orig_sd) + orig_mean
}

# Base map land info
region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

#----------------------------------
## Get data
#----------------------------------
all_mod_data<- readRDS(here::here("Data/Derived/all_model_data_predators.rds"))
# all(all_mod_data$trawl_id %in% env_tows)
# all(env_tows %in% all_mod_data$trawl_id)
# Every tow should have two observations and the unique tows should be the same as the ones in the environmental data
t<- table(all_mod_data$trawl_id)
# length(unique(all_mod_data$trawl_id)) == length(unique(env_data$ID))

# Focus on just lobster, and during GLORYs time series
year_min <- 1993
year_max <- 2018

red_mod_data <- all_mod_data |>
  dplyr::filter(between(year, year_min, year_max)) |>
  drop_na(total_biomass)

summary(red_mod_data)

# Cutting depth at 400 m
red_mod_data <- red_mod_data |>
  dplyr::filter(Depth <= 400)
summary(red_mod_data)

mod_data<- red_mod_data
# Scale/center covariates
# Get means and sds
column_means <- colMeans(mod_data[, c("Depth", "BT_seasonal")], na.rm = TRUE)
column_sds <- apply(mod_data[, c("Depth", "BT_seasonal")], 2, sd, na.rm = TRUE)

# Scale the data
mod_data <- mod_data |>
  ungroup() |>
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
  mutate(
    season = factor(season, levels = seasons),
    year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
    year_season_int = as.numeric(year_season_fac)
  ) %>%
  arrange(year_season_int)

# Survey as factor
mod_data <- mod_data |>
  mutate("survey" = factor(survey, levels = c("ME_NH", "MA", "NEFSC", "DFO")))

# Now for the model matrix trickery
mm_season <- model.matrix(~ 0 + season, data = mod_data)
# mm_year <- model.matrix(~ 0 + factor(est_year), data = dat) 

mod_data <- mod_data |>
  # dplyr::select(!contains("factor")) |>
  cbind(mm_season) |>
  #   cbind(mm_year) |>
  as_tibble()

# fa <- names(mod_data)[grepl("factor", names(mod_data))]
fa <- colnames(mm_season)
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
sdmTMB_mesh <- sdmTMB::make_mesh(mod_data, xy_cols = c("longitude", "latitude"), n_knots = 300, type = "kmeans")
# sdmTMB_mesh<- readRDS("~/Desktop/mesh_20241026_170037.rds")
plot(sdmTMB_mesh)

#----------------------------------
## Fit base
#----------------------------------
#| label: Basic environment-only SDM
#| include: false
#| echo: false
#| warning: false
n_seasons <- length(unique(mod_data$season)) 
# tau_Z_map <- factor(cbind(rep(1, n_seasons), rep(2, n_seasons))) # Delta model, two LPs
tau_Z_map <- factor(cbind(rep(1, n_seasons)))

fit_base <- sdmTMB(
  total_biomass ~ season + survey + s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
  # total_biomass ~ season + survey + poly(Depth_scaled, 2) + poly(BT_seasonal_scaled, 2),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data,
  spatial = "off",
  # offset = dat_mod$swept,
  # anisotropy = TRUE,
  # share_range = TRUE,
  spatiotemporal = "off",
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

saveRDS(fit_base, here::here("PredFit_Base.rds"))

sdmTMB::plot_smooth(fit_base, select = 2, ggplot = TRUE, level = 0.1, n = 200)
# sdmTMB::plot_smooth(fit_base, select = 1, ggplot = TRUE, level = 0.1, n = 200)
# plot(ggeffect(fit_base, "BT_seasonal_scaled [-2:4 by=1]"))

plot(mod_data$Depth_scaled, mod_data$total_biomass)
plot(mod_data$BT_seasonal_scaled, mod_data$total_biomass)

#----------------------------------
## Fit model with space
#----------------------------------
#| label: Adding in persistent spatial variation
#| include: false
#| echo: false
#| warning: false
fit_sp <- sdmTMB(
  total_biomass ~ season + survey + s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = mod_data,
  spatial = "on",
  # offset = dat_mod$swept,
  anisotropy = TRUE,
  # share_range = TRUE,
  spatiotemporal = "off",
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
saveRDS(fit_sp, here::here("PredFit_Sp.rds"))

#----------------------------------
## Fit model sp and spatio-temporal
#----------------------------------
#| label: Exclude space-invariant temporal autocorrelation.
#| include: false
#| echo: false
#| warning: false
#|

tictoc::tic()
fit_spst <- sdmTMB(
  total_biomass ~ season + survey + s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
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
saveRDS(fit_spst, here::here("PredFit_SpST.rds"))

sdmTMB::plot_smooth(fit_spst, select = 2, ggplot = TRUE, level = 0.1, n = 200)

#----------------------------------
## Fit model sp and spatio-temporal and AR1 on intercept
#----------------------------------
#| label: Exclude space-invariant temporal autocorrelation.
#| include: false
#| echo: false
#| warning: false
#|

tictoc::tic()
fit_spst_ar1 <- sdmTMB(
  total_biomass ~ season + survey + s(Depth_scaled, k = 4) + s(BT_seasonal_scaled, k = 4),
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

saveRDS(fit_spst_ar1, here::here("PredFit_SpST_AR1T.rds"))
sanity(fit_spst_ar1)
fit_spst_ar1

#----------------------------------
## Model comparison
#----------------------------------
#| label: Evaluate and validate models
#| include: false
#| echo: false
#| warning: false
#|

# Read in fits?
fit_base <- readRDS(here::here("PredFit_Base.rds"))
fit_sp <- readRDS(here::here("PredFit_Sp.rds"))
fit_spst <- readRDS(here::here("PredFit_SpST.rds"))
fit_spst_ar1 <- readRDS(here::here("PredFit_SpST_AR1T.rds"))

# Make a tibble with fitted model as its own column so we can map different functions to the models
all_fits<- list("Base" = fit_base, "Sp" = fit_sp, "SpST" = fit_spst, "SpST_AR1" = fit_spst_ar1)

fits_df <- tibble("Name" = c("fit_base", "fit_sp", "fit_spst", "fit_spst_ar1"), "Mod" = all_fits)

# Evaluation table
eval_df <- data.frame("Model" = c("Base", "Base + Spatial Variation", "Base + Spatial + Spatio-temporal Variation", "Base + Spatial + Spatio-temporal + AR1 Temporal Variation"), "AIC" = c(AIC(fits_df$Mod[[1]]), AIC(fits_df$Mod[[2]]), AIC(fits_df$Mod[[3]]), AIC(fits_df$Mod[[4]])), "Log Likelihood" = c(logLik(fits_df$Mod[[1]]), logLik(fits_df$Mod[[2]]), logLik(fits_df$Mod[[3]]), logLik(fits_df$Mod[[4]])))
eval_df

# Validation metrics and plot
# Get the prediction data frame
hold_out_ints<- as.numeric(factor(time_fac_levels, levels = time_fac_levels))[(length(time_ints)-15): length(time_ints)]
hold_out_yr_seas <- time_fac_levels[hold_out_ints]

pred_data_use<- mod_data |>
  filter(year_season_int %in% hold_out_ints)

pred_sdmTMB_nested<- function(fit, pred_data = pred_data_use){
  pred_out<- predict(fit, newdata = pred_data, type = "response")
  return(pred_out)
}

fits_df<- fits_df |>
  mutate(Preds = map(Mod, pred_sdmTMB_nested))

pred_table <- fits_df |>
  mutate(
    Corr_Coeff = pmap_dbl(list(df = Preds, obs_col = "total_biomass", mod_col = "est"), pearson_corr_coeff_func),
    RMSE = pmap_dbl(list(df = Preds, obs_col = "total_biomass", mod_col = "est"), rmse_func)
  ) |>
  dplyr::select(Name, Corr_Coeff, RMSE)
pred_table
# Data prep, need one big data frame with a "group" column for each mode
td_dat<- fits_df |>
  dplyr::select(Name, Preds) |>
  unnest(cols = c(Preds)) |>
  mutate(Model = factor(Name, levels = c("fit_base", "fit_sp", "fit_spst", "fit_spst_ar1"), labels = c("Base", "Base + Spatial Variation", "Base + Spatial + Spatio-temporal Variation", "Base + Spatial + Spatio-temporal + AR1 Temporal Variation")))

# Make the plot
td_plot<- taylor_diagram_func(dat = td_dat, obs = "total_biomass", mod = "est", group = "Model", out.file = "~/GitHub/lobSDM/Figures/Pred_Model_TD.png", grad.corr.lines = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), pcex = 1, cex.axis = 1, normalize = TRUE, mar = c(5, 4, 6, 6), sd.r = 1, fill.cols = c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'), color.cols = rep("black", 5), shapes = c(21, 22, 23, 24, 25), alpha = 0.5, example = FALSE)
td_plot


#----------------------------------
## Model predictions - DOUBLE CHECK MODEL
#----------------------------------
#| label: Make predictions from fitted models
#| include: false
#| echo: false
#| warning: false

# Grab best model
mod<- fit_spst_ar1 #!!!UPDATE THIS IF NEEDED!!#
pred_data <- readRDS(here::here("Data/Derived/pred_glorys_with_covs.rds")) |>
  mutate("season" = factor(season, levels = c("Spring", "Summer", "Fall"))) |>
  filter(Depth <= 400)

pred_data <- pred_data |>
  dplyr::filter(between(year, year_min, year_max) & season %in% c("Spring", "Summer", "Fall")) |>
  mutate(
    Depth_scaled = (Depth - column_means["Depth"]) / column_sds["Depth"],
    BT_seasonal_scaled = (BT_seasonal - column_means["BT_seasonal"]) / column_sds["BT_seasonal"]
  ) |>
  mutate(year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
         year_season_int = as.numeric(year_season_fac))
         
# Bind season stuff
mm_season <- model.matrix(~ 0 + season, data = pred_data)
pred_data<- pred_data |>
  cbind(mm_season) %>%
  as_tibble() %>%
  rename("longitude" = "x", "latitude" = "y") |>
  add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") %>%
  mutate(survey = factor("ME_NH", levels = levels(mod$data$survey)))%>%
  drop_na() 

fit_preds <- predict(mod, newdata = pred_data, type = "response", se = FALSE)
str(fit_preds)
quants<- quantile(fit_preds$est, probs = seq(from = 0, to = 1, by = 0.005), na.rm = TRUE)

max_pred<- ceiling(as.numeric(quants)[length(quants)-1])
fit_preds<- fit_preds |>
  mutate(est = ifelse(est >= max_pred, max_pred, est))
pred_lims<- c(min(fit_preds$est), max_pred)

# Nest and map
fit_preds <- fit_preds |>
  group_by(season, year, Year_Season, Date) |>
  nest()

map_nested <- function(pred_df, time, region_use = region, states_use = states, xlim_use = lon_lims, ylim_use = lat_lims, fill_lims = pred_lims, pred_max) {
  ggplot() +
    geom_raster(data = pred_df, aes(x = longitude, y = latitude, fill = est)) +
    geom_sf(data = region_use, fill = "#f0f0f0") +
    geom_sf(data = states_use, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_fill_viridis_c(name = "Predicted biomass", limits = fill_lims) +
    theme_minimal() +
    labs(
      fill = "Predicted biomass",
      title = time
    )
}

fit_preds <- fit_preds |>
  mutate("Pred_Map" = map2(data, Year_Season, map_nested)) |>
  arrange(Date)

start<- 1 # Spring, 2 = Summer, 3 = Fall
pred_out<- fit_preds$Pred_Map[[start]] + 
  fit_preds$Pred_Map[[start+15]] +
  fit_preds$Pred_Map[[start+30]] +
  fit_preds$Pred_Map[[start+45]] +
  fit_preds$Pred_Map[[start+60]] +
  fit_preds$Pred_Map[[start+75]] +
  patchwork::plot_layout(guides = "collect")
ggsave(here::here("Figures/PredExamplePredDistAbundance.jpg"), adult_out, dpi = 300)


#----------------------------------
## Model projections
#----------------------------------
#| label: Evaluate and validate models
#| include: false
#| echo: false
#| warning: false

# First, make sure we have the right sdmTMB version
# pak::pkg_install("pbs-assess/sdmTMB@project")
# Read in projection dataframes and rescale values
# proj_dat <- readRDS(here::here("Data/Derived/proj_glorys_with_covs.rds"))
proj_dat<- readRDS(here::here("Data/Derived/proj_025grid_with_covs.rds"))

covariate_rescale_func<- function(x, type, center = NULL, scale = NULL){
  if(type == "JT"){
    x.out<- x/100
    return(x.out)
  } else {
    x.out<- as.numeric(scale(abs(x), center = center, scale = scale))
    return(x.out)
  }
}

cols_keep<- names(proj_dat)[-c(9, 11)]

proj_dat_sub<- proj_dat |>
  filter(Depth <= 400) |>
  dplyr::select(one_of(cols_keep)) |>
  mutate(BT_seasonal_SSP1_26_scaled = as.numeric(scale(abs(SSP1_26_bot_temp_mean_match), center = column_means["BT_seasonal"], scale = column_sds["BT_seasonal"])),
  BT_seasonal_SSP5_85_scaled = as.numeric(scale(abs(SSP5_85_bot_temp_mean_match), center = column_means["BT_seasonal"], scale = column_sds["BT_seasonal"])),
  Depth_scaled = as.numeric(scale(abs(Depth), center = column_means["Depth"], scale = column_sds["Depth"])))
summary(proj_dat_sub)

# Add in factor columns needed to leverage SVC
# Bind season stuff
mm_season <- model.matrix(~ 0 + season, data = proj_dat_sub)
proj_dat_use<- proj_dat_sub |>
  cbind(mm_season) %>%
  as_tibble() %>%
  rename("longitude" = "x", "latitude" = "y") |>
  add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") %>%
  mutate(survey = factor("ME_NH", levels = levels(mod$data$survey)))%>%
  drop_na()

# Also need the right year_season stuff...
proj_year_max<- 2100
all_years <- seq(from = min(mod$data$year), proj_year_max, by = 1)
all_seasons <- levels(mod$data$season)
year_season_set <- expand.grid("season" = all_seasons, "year" = all_years)
all_year_season_levels <- apply(year_season_set[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")
proj_dat_use<- proj_dat_use |>
  mutate(year_season_fac = factor(paste(year, season, sep = "_"), levels = all_year_season_levels),
         year_season_int = as.numeric(year_season_fac)) |>
         drop_na(year_season_int)

# Make projections -- 
bt_vars<- c("BT_seasonal_SSP1_26_scaled", "BT_seasonal_SSP5_85_scaled")
bt_vars<- c("BT_seasonal_SSP5_85_scaled")
mod # Sp + SpST
n_sims<- 1

proj_dat_temp<- proj_dat_use |>
    rename(BT_seasonal_scaled = bt_vars)
proj_out <- sdmTMB::project(fit_spst, newdata = proj_dat_temp, nsim = 50, uncertainty = "random")

# Summarize, save, visualize
est_mean <- apply(proj_out$est, 1, mean) # summarize however you'd like
est_se <- apply(proj_out$est, 1, sd)
proj_dat_temp$proj_biomass_mean <- est_mean
proj_dat_temp$proj_biomass_se<- est_se

write_rds(proj_dat_temp, here::here("Data/Derived/pred_projected_biomass.rds"), compress = "gz")

# Visualize
sdmTMB_ind_list <- proj_out[["est"]]
sdmTMB_ind_mat <- as.matrix(data.frame(sdmTMB_ind_list))
dimnames(sdmTMB_ind_mat)[[1]] <- proj_dat_temp$year_season_int
dimnames(sdmTMB_ind_mat)[[2]] <- NULL
attr(sdmTMB_ind_mat, "time") <- "year"
attr(sdmTMB_ind_mat, "link") <- "log"

# Summarize
sdmTMB_ind_mvn <- get_index_sims(sdmTMB_ind_mat,
  level = 0.9,
  # return_sims = FALSE,
  area = rep(1, nrow(sdmTMB_ind_mat))
  # est_function = stats::median,
  # area_function = function(x, Area) x + log(Area),
  # agg_function = function(x) sum(exp(x))
)
# A lot going on there...simplify?
str(sdmTMB_ind_mvn)
season_info<- proj_dat_temp |>
  dplyr::select(year, season, year_season_fac, year_season_int) |>
  distinct()  |>
  rename(year_plot = year)
proj_dat_all<- sdmTMB_ind_mvn|>
  left_join(season_info, by = c("year" = "year_season_int")) 
summary(proj_dat_all)

p <- qnorm(0.9)
proj_dat_all$log_lwr0.1 <- proj_dat_all$log_est - p * proj_dat_all$se
proj_dat_all$log_upr0.9 <- proj_dat_all$log_est + p * proj_dat_all$se

# Reference...
bio_clim<- proj_dat_all |>
  group_by(season) |>
  filter(between(year_plot, 2010, max(mod$data$year))) |>
  summarize("mean" = mean(log_est))

proj_dat_all$season<- factor(proj_dat_all$season, levels = c("Spring", "Summer", "Fall"))

proj_ind_plot <- ggplot() +
  geom_ribbon(data = proj_dat_all, aes(x = year_plot, ymin = log_lwr0.1, ymax = log_upr0.9, fill = season), alpha = 0.3, linewidth = 0) +
  geom_line(data = proj_dat_all, aes(x = year_plot, y = log_est, color = season), lwd = 2) +
  # geom_hline(data = bio_clim, aes(yintercept = mean)) +
  scale_color_manual(name = "season", values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  scale_fill_manual(name = "season", values = c("#66c2a5", "#fc8d62", "#8da0cb")) + 
  xlab("Year") +
  ylab("Juvenile Relative Biomass Index") +
  facet_wrap(~season, nrow = 3) +
  theme_bw()
proj_ind_plot

