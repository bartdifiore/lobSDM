library(sf)

sites <- read_rds("Data/Derived/for_covariate_extraction.rds") %>%
  st_as_sf(coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT")) %>% 
  st_set_crs(4326)

library(giscoR)
coast <- gisco_get_coastallines(resolution = 3)

# boundary <-  df_raw %>% 
#   st_union() %>%
#   st_concave_hull(ratio = 0.01, allow_holes = F) %>% 
#   st_buffer(dist = 1000*5) # 5 km buffer


# plot(boundary, col = "red")
plot(sites, add = T, cex = 0.01)

g <- sites %>%
  st_make_grid(cellsize =0.25, square = T) %>%
  st_intersection(coast)

ggplot(g) + 
  geom_sf()+
  geom_sf(data = coast, fill = "grey")

plot(boundary, col = "red")
plot(points, add = T, cex = 0.01)
plot(pred_grid, add = T)

# big_map <- ggplot()+
#   geom_sf(data = df_raw, size = 0.01, color = "gray")+
#   geom_sf(data = df_spatial, fill = "transparent", color = "blue", lwd = 1)+
#   geom_sf(data = epu, fill = "transparent", color = "darkred", lwd = 1)+
#   geom_sf(data = boundary, fill = "transparent")+
#   geom_sf(data = pred_grid, fill = "transparent", color = "orange", lwd = 1)+
#   coord_sf(xlim = c(-100000, 800000), ylim = c(1000*4000, 1000*5000))
# ggsave("Figures/statarea_epu_grid.png", plot = big_map)








sites %>% 
  filter(survey == "DFO") %>%
  ggplot()+
  geom_sf(size = 0.1)+
  geom_sf(data = coast, fill = "grey") +
  coord_sf(
    xlim = c(-55, -70),
    ylim = c(40, 50)
  )



sites %>% 
  filter(survey == "DFO", EST_YEAR > 1999) %>%
  ggplot()+
  geom_sf(size = 0.1)+
  geom_sf(data = coast, fill = "grey") +
  coord_sf(
    xlim = c(-55, -70),
    ylim = c(40, 50)
  )



dfo_tows <- sites %>% filter(survey == "DFO") %>% 
  select(ID)

coca_tows <- read_csv("~/Desktop/COCA19_DFO_Tows.csv")

coca_tows <- coca_tows$x
dfo_tows <- dfo_tows$ID

setdiff(coca_tows, dfo_tows)
setdiff(dfo_tows, coca_tows)


x <- c(1,2,3)
y <- c(3,4,5)
setdiff(x,y)
