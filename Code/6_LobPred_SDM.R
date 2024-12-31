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

df<- readRDS(here::here("Data/Derived/all_model_data_predators.rds")) %>%
  dplyr::filter(between(year, year_min, year_max))

df %>% 
  group_by(survey) %>% 
  count()

mod_data <- df

# Scale/center covariates
# Get means and sds
column_means <- colMeans(mod_data[, c("Depth", "BT_seasonal")], na.rm = TRUE)
column_sds <- apply(mod_data[, c("Depth", "BT_seasonal")], 2, sd, na.rm = TRUE)

# Scale the data
mod_data <- mod_data |>
  ungroup() %>% # So in this version of the code i had to ungroup otherwise I got NA's. I didn't do this when prepping the juvenile lobster data. We need to make sure that the juvenile lobster code isn't only scaling within groups....
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

sdmTMB_mesh <- sdmTMB::make_mesh(mod_data, xy_cols = c("X", "Y"), n_knots = 200, type = "kmeans")
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
fit6 <- sdmTMB(
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

write_rds(fit6, here::here("Pred_Fit6.rds"), compress = "gz")

