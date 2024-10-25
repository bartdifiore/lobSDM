#----------------------------------
## Libraries
#----------------------------------

library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)


x <- rnorm(1000)

scaled_x <- scale(x)
unscaled_x <- (scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center")
  
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

covariates <- read_rds("Data/Derived/all_tows_all_covs.rds") %>% 
  rename(trawl_id = ID, year = EST_YEAR, latitude = DECDEG_BEGLAT, longitude = DECDEG_BEGLON) %>%
  select(trawl_id, latitude, longitude, season, year, survey, Depth, BT_seasonal)

df <- read_rds("Data/Derived/combined_and_filtered_lobsters.rds") %>%
  mutate(total_weight_at_length = number_at_length*weight_at_length) %>%
  group_by(trawl_id, longitude, latitude, season, year, survey, date) %>% 
  summarize(biomass = sum(total_weight_at_length)) %>%
  full_join(covariates) %>%
  janitor::clean_names() %>%
  drop_na(bt_seasonal, depth, year) %>%
  replace_na(list(biomass=0)) # Zero fill tows in which juvenile lobster where NOT found

names(df)

#----------------------------------
## AA data prep trickery
#----------------------------------

# Going to want to have a continuous time column
all_years<- seq(from = min(df$year), to = max(df$year))
seasons<- c("Spring", "Summer", "Fall")
time_fac_levels<- paste(rep(all_years, each = length(unique(seasons))), seasons, sep = "_")
time_ints<- as.numeric(factor(time_fac_levels, levels = time_fac_levels))

df <- df %>%
  mutate(season = factor(season, levels = seasons),
         year_season_fac = factor(paste(year, season, sep = "_"), levels = time_fac_levels),
         year_season_int = as.numeric(year_season_fac)) %>%
  arrange(year_season_int)

# Now for the model matrix trickery
mm_season <- model.matrix(~ 0 + factor(season), data = df)
# mm_year <- model.matrix(~ 0 + factor(est_year), data = dat)

df_mod <- df |>
  select(!contains("factor")) |>
  cbind(mm_season) |>
  #   cbind(mm_year) |>
  as_tibble()

fa <- names(df_mod)[grepl("factor", names(df_mod))]
fo <- paste0("`", paste(fa, collapse = "` + `"), "`")
svc <- as.formula(paste("~", fo))

# Check
svc



#----------------------------------
## Make mesh
#----------------------------------

df_mod <- df_mod %>%
  sdmTMB::add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") 

# Create sdmTMB mesh -- check out Owen Liu's code here (https://github.com/owenrliu/eDNA_eulachon/blob/63a1b4d21fa4ffbc629cbb0657bc032998565f17/scripts/eulachon_sdms.Rmd#L217) for more ideas?
sdmTMB_mesh<- sdmTMB::make_mesh(df_mod, xy_cols = c("X", "Y"), type = "cutoff", cutoff = 30, fmesher_func = fmesher::fm_mesh_2d_inla)
plot(sdmTMB_mesh)


df_mod %>% 
  group_by(survey) %>% 
  summarize(total = sum(biomass))


#----------------------------------
## Fit model 0
#----------------------------------

#| label: Basic environment-only SDM
#| include: false
#| echo: false
#| warning: false

# Starting with the most simple model that just includes the environmental covariates.
fit0<- sdmTMB(
  biomass ~ s(depth, k = 4) + s(bt_seasonal, k = 4),
  data = df_mod,
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
  family = delta_lognormal(type = "poisson-link"), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit0, model = 1, effects = "fixed") # Intercept

# Where are the smooths?
names(fit0)
names(fit0$smoothers) # All in here...
str(fit0$smoothers$basis_out) # List, one element for each smooth term, within each list element all the information of the fitted smooth is included



#----------------------------------
## Fit model 1
#----------------------------------

#| label: Env-only with season and catchability factors
#| include: false
#| echo: false
#| warning: false

fit1<- sdmTMB(
  biomass ~ factor(season) + factor(survey) + s(depth, k = 4) + s(bt_seasonal, k = 4),
  data = df_mod,
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
  family = delta_lognormal(type = "poisson-link"), 
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

fit2<- sdmTMB(
  biomass ~ + factor(season) + factor(survey) + s(depth, k = 4) + s(bt_seasonal, k = 4),
  data = df_mod,
  spatial = "off",
  # offset = dat_mod$swept, 
  # anisotropy = TRUE, 
  # share_range = TRUE,
  spatiotemporal = "off",
  spatial_varying = NULL,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1", # RW default
  extra_time = NULL,
  mesh = sdmTMB_mesh,
  family = delta_lognormal(type = "poisson-link"), 
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
    df1 <- data.frame("Model" = 1, "Time" = fit$time_lu$time_from_data, "Est" = est[,], "SE" = se[,])
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

n_seasons <- length(unique(df_mod$season)) 
tau_Z_map <- factor(cbind(rep(1, n_seasons), rep(2, n_seasons)))

fit3<- sdmTMB(
  biomass ~ factor(season) + factor(survey) + s(depth, k = 4) + s(bt_seasonal, k = 4),
  control = sdmTMBcontrol(map = list(ln_tau_Z = tau_Z_map)),
  data = df_mod,
  spatial = "off",
  # offset = dat_mod$swept, 
  # anisotropy = TRUE, 
  # share_range = TRUE,
  spatiotemporal = "off",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1",
  extra_time = NULL,
  mesh = sdmTMB_mesh,
  family = delta_lognormal(type = "poisson-link"), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit3, effects = "fixed")
tidy(fit3, effects = "ran_pars") # Now have range, and three sigma_Z's, representing seasonal variability.
fit3
sanity(fit3)



#----------------------------------
## Fit model 4
#----------------------------------

fit4<- sdmTMB(
  biomass ~ factor(season) + factor(survey) + s(depth, k = 4) + s(bt_seasonal, k = 4),
  data = df_mod,
  spatial = "on",
  # offset = dat_mod$swept, 
  anisotropy = TRUE, 
  # share_range = TRUE,
  # spatiotemporal = "ar1",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1",
  extra_time = NULL,
  mesh = sdmTMB_mesh,
  family = delta_lognormal(type = "poisson-link"), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit4, effects = "fixed")
tidy(fit4, effects = "ran_pars")
fit4
sanity(fit4)


#----------------------------------
## Fit model 5
#----------------------------------

fit5<- sdmTMB(
  biomass ~ factor(season) + factor(survey) + s(depth, k = 4) + s(bt_seasonal, k = 4),
  data = df_mod,
  spatial = "on",
  # offset = dat_mod$swept, 
  anisotropy = TRUE, 
  share_range = FALSE,
  spatiotemporal = "ar1",
  spatial_varying = svc,
  time = "year_season_int",
  time_varying = ~ 0 + year_season_int,
  time_varying_type = "ar1",
  extra_time = NULL,
  mesh = sdmTMB_mesh,
  family = delta_lognormal(type = "poisson-link"), 
  silent = FALSE,
  do_fit = TRUE
)

tidy(fit5, effects = "fixed")
tidy(fit5, effects = "ran_pars")
fit5

write_rds(fit5, "Data/Derived/model_fitbd.rds")


