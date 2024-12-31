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

#------------------------------------------------------------
## Estimate overlap metrics for each cell
#------------------------------------------------------------

grid <- df_overlap %>%
  distinct(latitude, longitude) %>%
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>%
  st_buffer(dist = 1/12) %>%
  mutate(cell_id = 1:n(), 
         area = st_area(geometry))

plot(grid)

df_overlap <- df_overlap %>% 
  sf::st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4326) %>% 
  st_join(grid) %>%
  mutate(area_km2 = as.numeric(area)/1000000)

overlap <- df_overlap %>%
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


# overlaps <- df_overlap %>%
#   group_by(id) %>%
#   mutate(lobster_est_response = exp(lobster_est), 
#          predator_est_response = exp(predator_est), 
#          lobster_bin = ifelse(lobster_est_response <= quantile(lobster_est_response, 0.05), 0, 1), 
#          predator_bin = ifelse(predator_est_response <= quantile(predator_est_response, 0.05), 0, 1)) %>%
#   mutate(area_overlap = area_overlapfn(prey = lobster_bin, pred = predator_bin, area = area_km2), 
#          loc_colloc = loc_collocfn(prey = lobster_est_response, pred = predator_est_response), 
#          asymmalpha_overlap = asymmalpha_overlapfn(prey = lobster_est_response, pred = predator_est_response), 
#          biomass_overlap = biomass_overlapfn(prey = lobster_est_response, pred = predator_est_response), 
#          hurlbert_overlap = hurlbert_overlapfn(prey = lobster_est_response, pred = predator_est_response, area = area_km2), 
#          schoeners_overlap = schoeners_overlapfn(prey = lobster_est_response, pred = predator_est_response), 
#          bhatta_coef = bhatta_coeffn(prey = lobster_est_response, pred = predator_est_response), 
#          AB_overlap = AB_overlapfn(prey = lobster_est_response, pred = predator_est_response)) #%>%
#   # pivot_longer(cols = area_overlap:AB_overlap, names_to = "overlap_metric", values_to = "value")

write_rds(overlaps, "Data/Derived/overlap_metrics.rds", compress = "gz")

# summary(area_overlapfn_local(overlaps$lobster_bin, overlaps$predator_bin, area = overlaps$area_km2))
# 
# overlaps$loc_colloc <- loc_collocfn(prey = overlaps$lobster_est_response, pred = overlaps$predator_est_response)
# 
# summary(overlaps$loc_colloc)
# 
# AB_overlapfn2 <- function(prey, pred) { 
#   ((pred - mean(pred, na.rm = T)) * (prey - mean(prey, na.rm = T)))/(mean(pred, na.rm = T) * mean(prey, na.rm = T)) 
# }
# 
# temp <- AB_overlapfn2(overlaps$lobster_est_response, overlaps$predator_est_response)
# 
# 
# p_prey = overlaps$lobster_est_response/sum(overlaps$lobster_est_response)
# p_pred <- overlaps$predator_est_response/sum(overlaps$predator_est_response)
# sum(p_prey*p_pred, na.rm = T)/(sqrt(sum(p_prey^2, na.rm = T)*sum(p_pred^2, na.rm = T)))
# 
# temp <- df_overlap %>%
#   filter(id == 1) %>%
#   mutate(lobster_est_response = exp(lobster_est), 
#          predator_est_response = exp(predator_est))
# 
# area_overlapfn(prey = temp$lobster_est_response, pred = temp$predator_est_response, area = temp$area_km2)
# loc_collocfn(prey = temp$lobster_est_response, pred = temp$predator_est_response)


































