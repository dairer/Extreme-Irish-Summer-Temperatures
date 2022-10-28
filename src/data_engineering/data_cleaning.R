# Author: Dáire Healy
# Date: 17 Feb 2022

# Description: This script prepares observational temperature data for marginal 
# modelling. This involves incorporating all covariates as descriped in the 
# paper associated to this work.
# Note: bulk modelling in preformed in this script. 
gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(mgcv)
library(raster) # package for netcdf manipulation

setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")

# i. ========  ========  Global parameters ========  ========  ========
# 
# # marginal threshold 
marg_thresh_quantile = 0.8 # equiv to 95th quantile of highest 1/4 of data

# months to include in the analysis
months_to_study = c(6,7,8)

potential_shape_values_climate = seq(-0.2, -0.15, length.out = 250)
#quantiles_to_estimate_bulk = seq(0.01, 0.99, length.out = 300)

refit_clim_shape_param = F


# # 1. ========  ========  Read in all data  ========  ========  ======== 
# 1a. ---- for plotting 
# colour palatte 
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

# map of ireland sf
ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]

# 1b. ---- observational data
obs_data = read_csv("data/processed/obs_all_data.csv") %>%
  mutate(year = lubridate::year(date)) %>%
  filter(lubridate::month(date) %in% months_to_study) # just take summer months

# project coordinates
remove.stations = obs_data %>% group_by(Station) %>% summarise(n = n()) %>% arrange(n) %>% filter(n < 100) %>% pull(Station)
obs_data = obs_data %>% filter(!Station %in% remove.stations)


# # 1b. ---- climate data
historic = vroom("~/JRSS_organised_code/data/raw_data/clim/r12i1p1-ICHEC-EC-EARTH-(Ireland)CLMcom-CLM-CCLM4-8-17 (EU).csv")
clim_data = historic %>%
  filter(lubridate::month(date) %in% months_to_study) # only take summer months
rm(list = c('historic')) # clear memory

# project coordinates
my_coords = clim_data %>% dplyr::select(Long, Lat)
sp::coordinates(my_coords)<- c("Long", "Lat")
sp::proj4string(my_coords) <- sp::CRS("+proj=longlat +datum=WGS84")
proj_cords <- sp::spTransform(my_coords,
                              sp::CRS(paste0("+proj=utm +zone=29 ellps=WGS84")))
proj_cords = proj_cords %>% sp::coordinates()/100000
clim_data$Long.projected = proj_cords[,1]
clim_data$Lat.projected = proj_cords[,2]
clim_data$id = clim_data %>% group_indices(Long, Lat)
clim_data %>% write_csv("data/processed/full_clim_data.csv")
# --- delete --- clim_data= read_csv("~/JRSS_organised_code/corrections/data/processed_data/full_clim_data.csv")

# 1c. Add temporal covariates (HADCRUT5)
hadcrut_data = raster::brick("data/raw/miscl/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc")

hadcrut_loc_dat = raster::coordinates(hadcrut_data) %>% as_tibble()
names(hadcrut_loc_dat) = c('Long', 'Lat')
hadcrut_data = hadcrut_data %>% values() %>% as_tibble()
hadcrut_data$Long = hadcrut_loc_dat$Long
hadcrut_data$Lat = hadcrut_loc_dat$Lat
names(hadcrut_data) = names(hadcrut_data) %>% str_remove_all("X")

hadcrut_data = hadcrut_data %>%
  tidyr::pivot_longer(c(-Long, -Lat), names_to = "date", values_to = "temp_anom") %>% # Making each row an observation
  mutate(date = date %>% lubridate::ymd())


irealnd_hadcrut = hadcrut_data %>%
  filter(Long > -11, Long < -6,   Lat > 50, Lat < 56) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% c(6,7,8)) 

global_hadcrut = hadcrut_data %>% 
  drop_na() %>% 
  mutate(year = lubridate::year(date))%>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% c(6,7,8)) %>% 
  group_by(year, month) %>%
  summarise(temp_anom = mean(temp_anom)) %>% ungroup() %>% unique()


loess_fit_irel = predict(loess(temp_anom ~ year, data = irealnd_hadcrut), se=T)
loess_fit_glob = loess(temp_anom ~ year, data = global_hadcrut)

irealnd_hadcrut$loess_temp_anom = loess_fit_irel$fit
irealnd_hadcrut$loess_temp_anom_u = loess_fit_irel$fit+1.96*(loess_fit_irel$se.fit)
irealnd_hadcrut$loess_temp_anom_l = loess_fit_irel$fit-1.96*(loess_fit_irel$se.fit)

loess_glob_preds = predict(loess_fit_glob, irealnd_hadcrut, se = T)
irealnd_hadcrut$loess_glob_temp_anom = loess_glob_preds$fit
irealnd_hadcrut$loess_glob_temp_anom_u = loess_glob_preds$fit+1.96*(loess_glob_preds$se.fit)
irealnd_hadcrut$loess_glob_temp_anom_l = loess_glob_preds$fit-1.96*(loess_glob_preds$se.fit)




# ----- PLOTS 
had_crut_yearly_dat_plot = irealnd_hadcrut %>%
  ggplot()+
  geom_ribbon(aes(x = year,ymin = loess_temp_anom_l, ymax = loess_temp_anom_u, fill = 'Ireland'), alpha = 0.25)+
  geom_ribbon(aes(x = year,ymin = loess_glob_temp_anom_l, ymax = loess_glob_temp_anom_u, fill = 'Global'), alpha = 0.25)+
  geom_line(aes(year, loess_temp_anom, col = 'Ireland'), size = 1)+
  geom_line(aes(year, loess_glob_temp_anom, col = 'Global'), linetype = 'dashed', size = 1)+
  theme_minimal()+
  theme_minimal(12)+
  theme(legend.position = 'none')+
  labs(x = 'Year',
       y = "Temperature\nanomaly (°C)")+
  scale_x_continuous(breaks = c(1931, 1960, 1990, 2022), limits = c(1931,2022))+
  ylim(c(-0.2,1.15))+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))



europe_sf <-rnaturalearth::ne_countries(type = "map_units",
                                        scale='large', # resolution of map
                                        returnclass = 'sf')

hadcrut_map = hadcrut_data %>%
  drop_na() %>%
  filter(Long > -15,Long< 40,
         Lat >35, Lat <70) %>%
  filter(date == "1991-07-16") %>%
  ggplot() +
  geom_tile(aes(Long, Lat, fill = temp_anom)) +
  #geom_tile(aes(Long, Lat), col = 'black') +
  viridis::scale_fill_viridis() +
  geom_sf(data = europe_sf, alpha = 0, col = 'black', size = 0.5)+
  geom_point(data =  hadcrut_data %>%
               filter(Long > -11, Long < -6,   Lat > 50, Lat < 56) %>% unique(),
             aes(Long, Lat), shape = 4, size = 7, col = 'red')+
  xlim(-15, 40) +
  ylim(35,70)+
  theme_minimal(12)+
  labs(fill = "°C",
       x = "Longitude", y = "Latitude")
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

hadcrut_plts = gridExtra::grid.arrange(had_crut_yearly_dat_plot, hadcrut_map, nrow = 1)
ggsave(hadcrut_plts, filename = "~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/output/figs/hadcrut_plots.pdf",
       height = 3.5, width = 9)

obs_data = obs_data %>%
  left_join((irealnd_hadcrut %>%
               dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
               unique()), by = 'year')

# 1d. get climate subset that matches up with observational stations 
#  ------- calculate the distance between two points
dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

obs_sites = obs_data %>%
  dplyr::select(Station, Long, Lat, Long.projected, Lat.projected) %>%
  unique()

clim_sites = clim_data %>%
  dplyr::select(id, Long, Lat, Long.projected, Lat.projected) %>%
  unique()

# itterate through all stations and find climate grid point
closest_irel = c()
for(i in seq(nrow(obs_sites))){
  smallest_dist = 9999999999
  id_of_smallest_dist_irel = 9999999999
  
  for(j in seq(nrow(clim_sites))){
    this_dist = dist_h(obs_sites[i,]$Long.projected,
                       obs_sites[i,]$Lat.projected,
                       clim_sites[j,]$Long.projected,
                       clim_sites[j,]$Lat.projected)
    
    if(this_dist < smallest_dist){
      smallest_dist = this_dist
      id_of_smallest_dist_irel = clim_sites$id[j]
    }
  }
  closest_irel = c(closest_irel, id_of_smallest_dist_irel)
}

# add column "id" to obs data which is the id of the closest clim grid point
obs_sites$id = closest_irel
obs_data = obs_data %>%
  left_join(obs_sites)


# # ---- Observational data threshold model 
# --- get thresold at each site
clim_thresh_values = clim_data %>%
  filter(id %in% obs_data$id) %>%
  group_by(id) %>%
  summarise(clim_thresh_value_9 = quantile(maxtp, 0.9))

obs_data = obs_data %>%
  left_join(clim_thresh_values)




qunatile_model <-  maxtp ~ clim_thresh_value_9
qunatile_model_fit <- evgam(qunatile_model, obs_data, family = "ald", ald.args = list(tau = 0.9))
obs_data$threshold_9 = qunatile_model_fit$location$fitted
qunatile_model_fit %>% saveRDS("output/threshold_model_9")





# 3. ========  ========  Tail Model        ========  ========  ======== 
print("model climate tail")


# 3a. ---- Climate model
clim_data_extreme_9 = clim_data %>%
  group_by(id) %>%
  mutate(threshold = quantile(maxtp, 0.9),
         excess = maxtp - threshold) %>%
  filter(excess > 0) %>%
  ungroup()



ngll = function(par){
  if(par <= 0) return(2^30)
  if(par > -1/shape_param) return(2^30)
  if(any((1+shape_param*this.dat/par)< 0)) return(2^30)
  
  -sum(evd::dgpd(x = this.dat, loc=0, scale = par, shape=shape_param, log=T))
}

estiamte_scale_fixed_shape = function(x,shape_c){
  this.dat <<- x
  shape_param <<- shape_c
  optim(par = c(0), fn = ngll, method = 'Brent', lower=0, upper = 5)
}

# fit scale parameter to each location with constant shap
num_sites = clim_data_extreme_9$id %>% unique()


if(refit_clim_shape_param){
  scales = c()
  loglik_sum = c()
  
  for(potential_shape in potential_shape_values_climate){
    print(potential_shape)
    
    loglik = c()

    for(i in (clim_data_extreme_9$id %>% unique())){
      this_clim_extm_irel = clim_data_extreme_9 %>% filter(id == i) %>% pull(excess)

      model_fit = estiamte_scale_fixed_shape(this_clim_extm_irel, potential_shape)

      scales = c(scales, model_fit$par)
      loglik = c(loglik, model_fit$value)

    }
    loglik_sum = c(loglik_sum, sum(loglik))
  }
  
  optimal_shape = potential_shape_values_climate[which.min(loglik_sum)]

  #optimal_shape =-0.2038076
  
}else{
  # optimal_shape_8 = -0.2038076
  optimal_shape_9 = -0.1889558
  #optimal_shape_925 = -0.1753507
  #optimal_shape =-0.2038076
}



scales_9 = c()
for(i in (clim_data_extreme_9$id %>% unique())){
  this_clim_extm_irel_9 = clim_data_extreme_9 %>% filter(id == i) %>% pull(excess)
  model_fit_9 = estiamte_scale_fixed_shape(this_clim_extm_irel_9, optimal_shape_9)
  scales_9 = c(scales_9, model_fit_9$par)
}


# clim_data_extreme_8 %>% 
#   dplyr::select(Long, Lat, id) %>%
#   unique() %>% 
#   mutate(scale_8 = scales_8,
#          scale_9 = scales_9,
#          scale_925 = scales_925) %>%
#   pivot_longer(-c(Long, Lat, id)) %>%
#   rename(thresh_val = name) %>%
#   ggplot()+
#   geom_point(aes(Long, Lat, col = value), size = 2.5)+
#   geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
#   ggplot2::scale_color_gradientn(colors = my_pal)+
#   facet_wrap(~thresh_val)+
#   theme_minimal(12)+
#   labs(col = expression(sigma))



# 
clim_data_extreme_9 %>%
  dplyr::select(Long, Lat, id, Long.projected, Lat.projected) %>%
  unique() %>%
  mutate(scale_9 = scales_9) %>%
  write_csv("data/processed/clim_scale_grid.csv")
# 
# # 
obs_sites = obs_sites %>%
  left_join(clim_data_extreme_9 %>%
              dplyr::select(Long, Lat, id, Long.projected, Lat.projected) %>%
              unique() %>%
              mutate(scale_9 = scales_9) %>%
              dplyr::select(id, scale_9), by = 'id')


# obs_sites %>%
#   ggplot()+
#   geom_point(aes(Long, Lat, col = scale_8), size = 2.5)+
#   geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
#   ggplot2::scale_color_gradientn(colors = my_pal)+
#   theme_minimal(12)+
#   labs(col = expression(sigma))
# 
# 
obs_data = obs_data %>%
  left_join(obs_sites)


obs_data %>% write_csv("data/processed/obs_data.csv")