library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#------------------------------------------------------------
## Get data
#------------------------------------------------------------

df <- readRDS("Data/Derived/overlap_metrics.rds") %>%
  mutate(x = st_coordinates(geometry)[,1], 
         y = st_coordinates(geometry)[,2]) %>%
  st_drop_geometry()

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

#------------------------------------------------------------
## Visualize lobster and predators
#------------------------------------------------------------

map_biomass <- function(data, time_year, time_season, region_use = region, states_use = states, xlim_use = lon_lims, ylim_use = lat_lims, response = "lobster_est_response"){
  
  temp_df <- data %>%
    filter(year == time_year, 
           season == time_season) %>%
    rename(predicted = response)

    ggplot() +
      geom_raster(data = temp_df, aes(x = x, y = y, fill = predicted)) +
      geom_sf(data = region_use, fill = "#f0f0f0") +
      geom_sf(data = states_use, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
      # scale_fill_viridis_c() +
      scale_fill_viridis_c(trans = "log") +
      theme_minimal() +
      labs(
        fill = "Predicted biomass",
        title = paste(time_year, time_season, sep = " ")
      )
}


map_biomass(data = df, 
            time_year = 2013, 
            time_season = "Fall", 
            response = "predator_est_response")

df %>%
  filter(year %in% c(1993, 1995, 2000, 2005, 2010, 2015, 2020, 2023), 
         season == "Fall") %>% 
  select(year, lobster_est_response,x, y) %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = lobster_est_response)) +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  # scale_fill_viridis_c() +
  # scale_fill_viridis_c(transform = "log") +
  scale_fill_viridis_c(transform = "log") +
  theme_minimal() +
  facet_wrap(~year)+
  labs(
    fill = "Predicted biomass", 
    title = "Juvenile lobster biomass", x = "", y = ""
  )
ggsave("Figures/predicted_juvenile_lobster.png")

hist(df$lobster_est_response)

df %>%
  filter(year %in% c(1993, 1995, 2000, 2005, 2010, 2015, 2020, 2023), 
         season == "Fall") %>% 
  select(year, predator_est_response,x, y) %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = predator_est_response)) +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  # scale_fill_viridis_c() +
  scale_fill_viridis_c(trans = "log") +
  theme_minimal() +
  facet_wrap(~year)+
  labs(
    fill = "Predicted biomass", 
    title = "Aggregated predator biomass", x = "", y = ""
  )
ggsave("Figures/predicted_predator.png")
  


#------------------------------------------------------------
## Visualize overlap metrics
#------------------------------------------------------------

map_overlap <- function(data, time_year, time_season, region_use = region, states_use = states, xlim_use = lon_lims, ylim_use = lat_lims, metric, log_scale = T){
  temp_df <- data %>%
    filter(year == time_year, 
           season == time_season, 
           overlap_metric == metric)
  # return(temp_df)
  if(log_scale == T){
    ggplot() +
      geom_raster(data = temp_df, aes(x = x, y = y, fill = value)) +
      geom_sf(data = region_use, fill = "#f0f0f0") +
      geom_sf(data = states_use, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
      # scale_fill_viridis_c() +
      scale_fill_viridis_c(trans = "log") +
      theme_minimal() +
      labs(
        fill = metric,
        title = paste(time_year, time_season, sep = " ")
      )
    
  }else{
    ggplot() +
      geom_raster(data = temp_df, aes(x = x, y = y, fill = value)) +
      geom_sf(data = region_use, fill = "#f0f0f0") +
      geom_sf(data = states_use, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
      scale_fill_viridis_c() +
      # scale_fill_viridis_c(trans = "log") +
      theme_minimal() +
      labs(
        fill = metric,
        title = paste(time_year, time_season, sep = " ")
      )
  }


}

map_overlap(data = df, 
            time_year = 2023, 
            time_season = "Fall", 
            metric = "pred_prey_ratio", 
            log_scale = T)

map_overlap(data = df, 
            time_year = 2023, 
            time_season = "Fall", 
            metric = "AB", 
            log_scale = F)

#------------------------------------------------------------
## Visualize log pred/prey ratio
#------------------------------------------------------------

forplot <- df %>%
  filter(year %in% c(1993, 1995, 2000, 2005, 2010, 2015, 2020, 2023), 
         season == "Fall", 
         overlap_metric == "pred_prey_ratio")

quantile(log(forplot$value), c(0.025, 0.5, 0.975))

breaks <- as.vector(quantile(log(forplot$value), seq(0.025, 0.975, by = 0.025*2)))
# labels <- seq(0, 1, by= 0.25)

cols <- colorRampPalette(c("white", "darkred"))

cols_use <- cols(length(breaks))


ggplot(forplot) +
    geom_raster(aes(x = x, y = y, fill = log(value))) +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
    # scale_fill_viridis_c() +
    scale_fill_stepsn(colors = cols_use, breaks = breaks, labels = NULL) +
    # scale_fill_steps(n.breaks = 10, nice.breaks = T, low = "white", high = "darkred", transform = "log2")+
    facet_wrap(~year)+
    labs(fill = "log(Predator-prey ratio)", x = "", y = "", title = "log Predator-prey ratio")+
    theme_minimal()
ggsave("Figures/pred_prey_ratio_byyear.png")


