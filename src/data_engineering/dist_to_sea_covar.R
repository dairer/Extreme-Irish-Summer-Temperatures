# add distance to see as a covariate

rm(list = ls())
library(tidyverse)
library(rnaturalearth)
library(sf)
library(raster)

#obs_data = read_csv("data/processed_data/obs_data_w_clim.csv")
obs_data = read_csv("~/Extreme-Irish-Summer-Temperatures/data/processed/obs_data.csv")

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


irel <-rnaturalearth::ne_countries(type = "map_units",
                                   scale= 10, # resolution of map
                                   returnclass = 'sf', # spatial object, change this to "sp" if needed, this will break the plot though
                                   continent = "europe") %>%
  .[.$name %in% c('Ireland'),]

nirel <-rnaturalearth::ne_countries(type = "map_units",
                                    scale= 10, # resolution of map
                                    returnclass = 'sf', # spatial object, change this to "sp" if needed, this will break the plot though
                                    continent = "europe") %>%
  .[.$name %in% c( 'N. Ireland'),]


ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale= 10, # resolution of map
                                         returnclass = 'sf', # spatial object, change this to "sp" if needed, this will break the plot though
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland','N. Ireland'),]


irel =sf::st_union(irel, nirel) 

irel <- st_transform(irel, 29902)
ireland_sf <- st_transform(ireland_sf, 29902)

grid <- st_make_grid(irel, cellsize = 100, what = "centers")
grid <- st_intersection(grid, irel) 
plot(grid)

dist_sea <- st_distance((irel %>% st_cast( "MULTILINESTRING")), grid)

df <- data.frame(
  dist_sea = as.vector(dist_sea)/1000,
  st_coordinates(grid))

sites = obs_data %>%
  dplyr::select(Long, Lat, Station, Long.projected, Lat.projected) %>%
  unique()

dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

ids_closest = c()
for(s in seq(nrow(sites))){
  print(s)
  smallest_dist = 99999999
  smallest_dist_id = 99999999
  
  for(i in seq(nrow(df))){
    this_dist = dist_h(sites[s,]$Long, sites[s,]$Lat, df[i,]$X, df[i,]$Y)
    
    if(this_dist < smallest_dist) {
      smallest_dist = this_dist
      smallest_dist_id = i
    }
  }
  ids_closest = c(ids_closest, smallest_dist_id)
}

sites$dist_sea = df[ids_closest,]$dist_sea

sites %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = dist_sea))

sites %>%
  write_csv("~/Extreme-Irish-Summer-Temperatures/data/processed/obs_data_dist_to_sea.csv")



# -------- Add these covariates to tibble of clim data
clim_data = read_csv("~/Extreme-Irish-Summer-Temperatures/data/processed_data/ireland_clim_extreme.csv")

irel_grid  = st_transform(grid,4326)
df <- data.frame(dist_sea = as.vector(dist_sea)/1000,
                 st_coordinates(irel_grid))

sites = clim_data %>%
  dplyr::select(Long, Lat, id, Long.projected, Lat.projected) %>%
  unique()

dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

ids_closest = c()
for(s in seq(nrow(sites))){
  print(s)
  smallest_dist = 99999999
  smallest_dist_id = 99999999
  
  for(i in seq(nrow(df))){
    this_dist = dist_h(sites[s,]$Long, sites[s,]$Lat, df[i,]$X, df[i,]$Y)
    
    if(this_dist < smallest_dist) {
      smallest_dist = this_dist
      smallest_dist_id = i
    }
  }
  print(smallest_dist_id)
  ids_closest = c(ids_closest, smallest_dist_id)
}

sites$dist_sea = df[ids_closest,]$dist_sea


sites %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = dist_sea))

sites %>%
  write_csv("~/Extreme-Irish-Summer-Temperatures/data/processed/clim_data_dist_to_sea.csv")
