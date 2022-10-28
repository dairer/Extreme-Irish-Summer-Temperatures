rm(list = ls())
library(tidyverse)
setwd("~/Extreme-Irish-Summer-Temperatures/")


obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Long, Lat, Station)

extreme_date = read_csv("data/processed/clim_data_full.csv") %>%
  filter(date == "1991-07-15") 



my_pal = c(
  '#062c30', # extra dark 
  '#003f5c',
  '#2f4b7c',
  '#665191',
  '#a05195',
  '#d45087',
  '#f95d6a',
  '#ff7c43',
  '#ffa600')

# map for plotting
ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object, change this to "sp" if needed, this will break the plot though
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]






my_plt_clim = extreme_date %>%
  ggplot()+
  geom_sf(data = ireland_sf, alpha = 0, col = "black")+
  geom_point(aes(Long, Lat, col = maxtp), size = 2)+
  ggplot2::scale_color_gradientn(colours = my_pal)+
  theme_minimal(12)+
  labs(col = "\n°C",
       title = "",
       x = "Longitude",
       y = "Latitude")


my_plt_obs = obs_sites %>%
  group_by(Station, Long, Lat) %>%
  summarise(num_obs = round(n()/92)) %>%
  ggplot()+
  geom_sf(data = ireland_sf, alpha = 0, col = "black")+
  geom_point(aes(Long, Lat, col = num_obs, size = num_obs))+
  scale_size_continuous(limits=c(1, 90), breaks=seq(0, 90, by = 20))+
  scale_color_continuous(limits=c(1, 90), breaks=seq(0, 90, by = 20)) +
  guides(color= guide_legend(), size=guide_legend())+
  ggplot2::scale_color_gradientn(limits=c(1, 90), breaks=seq(0, 90, by = 20),colours = my_pal)+
  geom_point(data = (obs_sites %>%
                       group_by(Station, Long, Lat) %>%
                       summarise(num_obs = round(n()/92)) %>%
                       filter(Station %in% c("malin_head", "roches_point", "phoenixpark", "claremorris", "mullingar")) %>%
                       dplyr::select(Station, Long, Lat, num_obs) %>%
                       unique()), aes(Long, Lat), shape = 4, col = 'black', size = 5)+
  theme_minimal(12)+
  labs(col = "Years\nof data",
       size = "Years\nof data",
       title = "",
       x = "Longitude",
       y = "Latitude")




my_plt = gridExtra::grid.arrange(my_plt_obs, my_plt_clim, nrow = 1)

ggsave(plot = my_plt, filename = "output/figs/clim_and_obs_sites.pdf",  height = 4.5, width = 9)



