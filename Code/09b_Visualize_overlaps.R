library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#------------------------------------------------------------
## Get data
#------------------------------------------------------------

df <- readRDS("Data/Derived/overlap_metrics_1deg.rds") %>%
  st_as_sf()

region <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")

lat_lims <- c(35.2, 48)
lon_lims<- c(-76, -56.2)

df %>% 
  st_drop_geometry() %>%
  group_by(overlap_metric) %>% 
  summarize(low = quantile(value, 0.025, na.rm = T), 
            mean = mean(value, na.rm = T), 
            high = quantile(value, 0.975, na.rm = T) 
            )


#------------------------------------------------------------
## Visualize overlap metrics
#------------------------------------------------------------


df %>%
  st_make_valid() %>%
  filter(year %in% c(1993, 1995, 2000, 2005, 2010, 2015, 2020, 2023), 
         season == "Fall", 
         overlap_metric == "biomass_overlap") %>%
  ggplot() +
  geom_sf(aes(fill = value), color = "transparent") +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  # scale_fill_viridis_c() +
  scale_fill_viridis_c(trans = "log") +
  facet_wrap(~year)+
  theme_minimal()


df %>%
  st_make_valid() %>%
  filter(year %in% c(1993, 1995, 2000, 2005, 2010, 2015, 2020, 2023), 
         season == "Fall", 
         overlap_metric == "asymmalpha") %>%
  ggplot() +
  geom_sf(aes(fill = value), color = "transparent") +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  # scale_fill_viridis_c() +
  scale_fill_gradient(low = "white", high = scales::muted("red")) +
  facet_wrap(~year)+
  labs(title = "Asymetrical alpha")+
  theme_minimal()
ggsave("Figures/asymetricala_paneled_by_year.png")

df %>%
  filter(cell_id == 149, 
         season == "Fall", 
         overlap_metric == "asymmalpha") %>%
  ggplot() +
  geom_line(aes(x = year, y = value))

#------------------------------------------------------------
## Focus on asymmetrical alpha
#------------------------------------------------------------

# According to Carroll et al. 2019, asymetric alpha is a usefule metric if we are interested in: 
  # 1. Trophic transfer from prey to predator
  # 2. Relative abundance more than absolute abundance. 
# Therefore, I'm going to focus first on asymmetrical alpha (for little reason more than it seems interesting and well justified). 

df %>%
  st_make_valid() %>%
  filter(year %in% c(1993, 1995, 2000, 2005, 2010, 2015, 2020, 2023), 
         season == "Fall", 
         overlap_metric == "asymmalpha") %>%
  ggplot() +
  geom_sf(aes(fill = value), color = "transparent") +
  geom_sf(data = region, fill = "#f0f0f0") +
  geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
  coord_sf(xlim = lon_lims, ylim = lat_lims, expand = FALSE) +
  # scale_fill_viridis_c() +
  scale_fill_viridis_c(trans = "log") +
  facet_wrap(~year)+
  theme_minimal()

# What I'm interested in visualizing is the relative change through time. So for each grid cell has the value of asymmetrical alpha increased or decreased since 1993? To do this we will need to fit a lm() to each time series (or simply estimate a correlation coefficient) and then visualize the slope of this relationship. 

grid.1 <- df %>% 
  select(geometry, cell_id) %>%
  distinct()


forplot <- df %>% 
  drop_na(cell_id) %>%
  filter(!overlap_metric %in% c("area_overlap", "range_overlap")) %>%
  st_drop_geometry() %>%
  group_by(cell_id, season, overlap_metric) %>% 
  nest() %>% 
  mutate(lm_obj = map(data, ~lm(value ~ year, data = .))) %>%
  mutate(lm_tidy = map(lm_obj, broom::tidy)) %>%
  ungroup() %>%
  transmute(cell_id, season, overlap_metric, lm_tidy) %>%
  unnest(cols = c(lm_tidy)) %>%
  filter(term == "year") %>% 
  select(cell_id, season, overlap_metric, estimate) %>%
  left_join(grid.1)


metrics <- unique(forplot$overlap_metric)
labels <- c("Asymetrical alpha", "Local colocation", "Biomass overlap", "Hurlbert's overlap", "Schoener's D", "Bhattacharyya's coefficient", "AB ratio")
out <- list()
for(i in 1:length(metrics)){
  out[[i]] <- forplot %>%
    st_as_sf() %>%
    st_make_valid() %>%
    filter(season == "Fall", 
           overlap_metric == metrics[i]) %>%
    ggplot() +
    geom_sf(aes(fill = estimate), color = "transparent") +
    geom_sf(data = region, fill = "#f0f0f0") +
    geom_sf(data = states, color = "dark gray", lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = lon_lims, ylim = lat_lims, expand = T) +
    # scale_fill_viridis_c() +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue")) +
    labs(title = labels[i])+
    theme_minimal()
}
cowplot::plot_grid(out[[1]], out[[2]], out[[3]], out[[4]],
                   out[[5]], out[[6]], out[[7]])
ggsave("Figures/Temporaltrends_in_overlap.png")  






















