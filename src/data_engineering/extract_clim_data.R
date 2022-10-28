rm(list =ls())
library(raster) # package for netcdf manipulation
library(tidyverse)
library(rnaturalearth)
library(magrittr)
library(ncdf4)
library(raster)



trans_to_longlat = function(rotated_long = 0, rotated_lat = 0, axis_long = -162, axis_lat = 39.25){
  alpha = (axis_lat + 90)
  beta = axis_long
  
  rotated_long = (rotated_long*pi)/180
  rotated_lat = (rotated_lat*pi)/180
  
  alpha = (alpha*pi)/180
  beta = (beta*pi)/180
  
  #  Convert from spherical to cartesian coordinates
  x = cos(rotated_long)*cos(rotated_lat)
  y = sin(rotated_long)*cos(rotated_lat)
  z = sin(rotated_lat)
  
  x_new = cos(alpha)*cos(beta)*x + sin(beta)*y + sin(alpha)*cos(beta)*z
  y_new = -1*cos(alpha)*sin(beta)*x + cos(beta)*y - sin(alpha)*sin(beta)*z
  z_new = -1*sin(alpha)*x + cos(alpha)*z
  
  
  # Convert cartesian back to spherical coordinates
  lon_new = atan2(y_new,x_new) 
  lat_new = asin(z_new)
  
  #Convert radians back to degrees
  lon_new = (lon_new*180)/pi
  lat_new = (lat_new*180)/pi;
  
  -c(lon_new,lat_new)
}


file.loc = "~/CORDEX-DATA-HISTORIC/r12i1p1-ICHEC-EC-EARTH-(Ireland)CLMcom-CLM-CCLM4-8-17 (EU)/"
file.save = "~/CORDEX-DATA-HISTORIC/processed/r12i1p1-ICHEC-EC-EARTH-(Ireland)CLMcom-CLM-CCLM4-8-17 (EU).csv"

my_files = list.files(file.loc)
for(f in paste0(file.loc, my_files)){
  
  print(f)
  my_dat = raster::brick(f, varname = 'tasmmax') # read in netcdf file
  
  # add id to each loc, pick which ids are in irealnd
  # we only want irish data
  ROI <- extent(-18, -12, 3, 7.5)
  
  # find coords in rotated?
  my_dat <- crop(my_dat,ROI) # takes a minute
  
  my_data = raster::coordinates(my_dat) %>% as_tibble()
  names(my_data) = c('Long', 'Lat')
  
  my_data = my_data %>%
    mutate(Long = signif(Long, 4),
           Lat = signif(Lat, 4)) %>%
    dplyr::mutate(id = group_indices(., Long, Lat))
  
  
  all_data = my_dat %>%
    values() %>%
    as_tibble()
  all_data$id = my_data$id
  names(all_data) = names(all_data) %>% str_remove_all("X")
  
  # ----- down size data
  transformed_coords = apply(my_data, 1, function(x) trans_to_longlat(x[1], x[2])) %>% t()
  
  transformed_coords = tibble(Long = transformed_coords[,1],
                              Lat = transformed_coords[,2],
                              orig_lon = my_data$Long,
                              orig_lat = my_data$Lat, 
                              id = my_data$id)
  
  ireland_sf <- ne_countries(type = "map_units",
                             scale='large', # resolution of map
                             returnclass = 'sf', # spatial object, change this to "sp" if needed, this will break the plot though
                             continent = "europe") %>%
    .[.$name %in% c('Ireland', 'N. Ireland'),]
  
  # create a spatial object
  sp::coordinates(transformed_coords)<- c("Long", "Lat")
  
  # add Coordinate Reference System, which is WGS84
  sp::proj4string(transformed_coords) <- sp::CRS("+proj=longlat +datum=WGS84")
  intersected = sf::st_intersection(sf::st_as_sf(transformed_coords), ireland_sf)
  
  locs = tibble(id = intersected$id,
                Long = sf::st_coordinates(intersected)[,1],
                Lat = sf::st_coordinates(intersected)[,2])
  
  all_data %>%
    right_join(locs, by = 'id') %>%
    dplyr::select(-id) %>%
    tidyr::pivot_longer(c(-Long, -Lat), names_to = "date", values_to = "maxtp") %>% # Making each row an observation
    mutate(date = date %>% lubridate::ymd()) %>%
    mutate(maxtp = maxtp - 273.15) %>%
    dplyr::mutate(id = group_indices(., Long, Lat)) %>%
    write.table(file.save, 
                sep = ",", 
                col.names = !file.exists(file.save), 
                append = T,
                row.names = F)
  
}


