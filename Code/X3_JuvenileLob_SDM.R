#----------------------------------
## Libraries
#----------------------------------

library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)

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