library(tidyverse)
library(sdmTMB)
library(sdmTMBextra)
library(sf)

#-----------------------------------------
## Get data
#-----------------------------------------

pred_mod <- readRDS("Data/Derived/Pred_Fit6.rds")
lob_mod <- readRDS("Data/Derived/Juve_Fit6.rds")


#-----------------------------------------
## Get prediction grid
#-----------------------------------------

year_min <- 1993
year_max <- 2023

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

#-------------------------------------------------------
## Generate spatiotemporal predictions for each model
#-------------------------------------------------------

lob_yhat <- predict(lob_mod, newdata = prediction_grid)

pred_yhat <- predict(pred_mod, newdata = prediction_grid)

source("Code/Carroll_2019_functions.R")

# Combine into one df

df_overlap <- pred_yhat %>% 
  arrange(season, year, Year_Season, Date, longitude, latitude) %>%
  ungroup() %>%
  mutate(id = 1:n()) %>%
  select(id, season, year, Depth, BT_seasonal, Depth_scaled, BT_seasonal_scaled, est, longitude, latitude) %>%
  rename(predator_est = est) %>%
  st_drop_geometry() %>%
  left_join(lob_yhat %>% 
              arrange(season, year, Year_Season, Date, longitude, latitude) %>%
              ungroup() %>%
              mutate(id = 1:n()) %>%
              select(id,season, year, Depth, BT_seasonal, Depth_scaled, BT_seasonal_scaled, est, longitude, latitude) %>%
              rename(lobster_est = est) %>%
              st_drop_geometry)

#------------------------------------------------------------------------------
## Estimate overlap metrics for each cell (e.g. do not integrate across space)
#------------------------------------------------------------------------------

grid <- df_overlap %>%
  distinct(latitude, longitude) %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>%
  st_buffer(dist = 1/12) %>%
  mutate(cell_id = 1:n(), 
         area = st_area(geometry))

plot(grid)

df_overlap2 <- df_overlap %>% 
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>% 
  st_join(grid) %>%
  mutate(area_km2 = as.numeric(area)/1000000)

overlap <- df_overlap2 %>%
  group_by(id) %>%
  mutate(lobster_est_response = exp(lobster_est), 
         predator_est_response = exp(predator_est))

overlap$lobster_bin <- ifelse(overlap$lobster_est_response <= quantile(overlap$lobster_est_response, 0.05), 0, 1)
overlap$predator_bin <- ifelse(overlap$predator_est_response <= quantile(overlap$predator_est_response, 0.05), 0, 1)

overlap$co_occurance <- ifelse(overlap$lobster_bin == 1 & overlap$predator_bin == 1, 1, 0)

overlaps <- overlap %>% 
  group_by(year, season) %>%
  mutate(p_lobster = lobster_est_response/sum(lobster_est_response), 
         p_predator = predator_est_response/sum(predator_est_response), 
         asymmalpha = p_lobster*p_predator/sum(p_lobster^2), 
         loc_colloc = p_lobster*p_predator/(sqrt(sum(p_predator^2)*sum(p_lobster^2))), 
         biomass_overlap = (lobster_est_response/max(lobster_est_response))*(predator_est_response/max(predator_est_response))/sum(lobster_est_response/max(lobster_est_response)), 
         schoeners = 1 - 0.5 * (abs(p_lobster - p_predator)), 
         bhatta = sqrt(p_lobster*p_predator), 
         AB = (predator_est_response - mean(predator_est_response))*(lobster_est_response - mean(lobster_est_response))/(mean(predator_est_response)*mean(lobster_est_response)), 
         pred_prey_ratio = predator_est_response/lobster_est_response) %>% 
  pivot_longer(cols = c(co_occurance, asymmalpha:pred_prey_ratio), names_to = "overlap_metric", values_to = "value" )



write_rds(overlaps, "Data/Derived/overlap_metrics.rds", compress = "gz")


#------------------------------------------------------------------------------
## Estimate overlap metrics for larger spatial areas
#------------------------------------------------------------------------------

library(rnaturalearth)
library(rnaturalearthdata)

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

# Build a 1 degree grid, overlay the observations, estimate metrics 

grid.1 <- st_bbox(overlaps) %>%
  st_as_sfc() %>%
  st_make_grid(cellsize = 1) %>%
  st_sf() %>%
  mutate(cell_id = 1:n(), 
         area = st_area(geometry), 
         area_km2 = as.numeric(area/1000000)) %>%
  select(-area) %>%
  st_make_valid()

ggplot() +
  geom_sf(data = grid.1) +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE)

# Merge the data to the 1 degree grid

df_overlap3 <- df_overlap %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_join(grid.1)

# Estimate the metrics

overlaps_1deg <- df_overlap3 %>% 
  st_drop_geometry() %>%
  group_by(season, year, cell_id, area_km2) %>%
  mutate(p_res = exp(predator_est), # Here I just use the predicted estimates on the response scale. However, should we use get_index() to estimate relative biomas in each spatial unit???
         l_res = exp(lobster_est)
  )

overlaps_1deg$lobster_bin <- ifelse(overlaps_1deg$l_res <= quantile(overlaps_1deg$l_res, 0.05), 0, 1)
overlaps_1deg$predator_bin <- ifelse(overlaps_1deg$p_res <= quantile(overlaps_1deg$p_res, 0.05), 0, 1)

out_1_deg <- overlaps_1deg %>% 
  group_by(season, year, cell_id, area_km2) %>%
  summarize(area_overlap = area_overlapfn(l_res, p_res, area_km2),
         range_overlap = range_overlapfn(l_res, p_res, area_km2),
         asymmalpha = asymmalpha_overlapfn(l_res, p_res), 
         loc_colloc = loc_collocfn(l_res, p_res), 
         biomass_overlap = biomass_overlapfn(l_res, p_res),
         hurlbert = hurlbert_overlapfn(l_res, p_res, area_km2),
         schoeners = schoeners_overlapfn(l_res, p_res), 
         bhatta = bhatta_coeffn(l_res, p_res), 
         AB = AB_overlapfn(l_res, p_res)) %>%
  pivot_longer(cols = c(area_overlap:AB), names_to = "overlap_metric", values_to = "value" ) %>%
  left_join(grid.1)

write_rds(out_1_deg, "Data/Derived/overlap_metrics_1deg.rds")  
  




















