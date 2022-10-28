gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(mgcv)
library(raster) # package for netcdf manipulation

setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")


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



# i. ========  ========  Global parameters ========  ========  ========
# 

# months to include in the analysis
months_to_study = c(6,7,8)
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


ngll_free_xi = function(par){
  if(par[1] <= 0) return(2^30)
  if(par[1] > -1/par[2]) return(2^30)
  if(any((1+par[2]*this.dat/par[1])< 0)) return(2^30)
  
  -sum(evd::dgpd(x = this.dat, loc=0, scale = par[1], shape=par[2], log=T))
}

estiamte_scale_free_shape = function(x){
  this.dat <<- x
  optim(par = c(1.5,-0.15), fn = ngll_free_xi)
}






# 3a. ---- Climate model
clim_data_extreme_9 = clim_data %>%
  group_by(id) %>%
  mutate(threshold = quantile(maxtp, 0.9),
         excess = maxtp - threshold) %>%
  filter(excess > 0) %>%
  ungroup()


optimal_shape_9 = -0.1889558


grid = read_csv('data/processed/clim_scale_grid.csv')


scales_fixed_xi = c()
loglik_fixed_xi = c()

scales_free_xi = c()
loglik_free_xi = c()


for(i in grid$id){
  print(i)
  this_clim_extm_irel = clim_data_extreme_9 %>% filter(id == i) %>% pull(excess)
  
  # --- fixed shape
  model_fit = estiamte_scale_fixed_shape(this_clim_extm_irel, optimal_shape_9)
  scales_fixed_xi = c(scales_fixed_xi, model_fit$par)
  loglik_fixed_xi = c(loglik_fixed_xi, model_fit$value)
  
  # --- free shape
  model_fit = estiamte_scale_free_shape(this_clim_extm_irel)
  scales_free_xi = c(scales_free_xi, model_fit$par)
  loglik_free_xi = c(loglik_free_xi, (evd::fpot(this_clim_extm_irel, threshold = 0) %>% logLik %>% as.numeric()))
}





grid$LL_null = loglik_fixed_xi
grid$LL_alt = -loglik_free_xi
grid$LL_ratio = -2*log(grid$LL_alt/grid$LL_null)


plt_clim = grid %>%
  filter(LL_ratio < 0.015) %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = LL_ratio))+
  coord_map()+
  theme_minimal()+
  geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  scale_x_continuous(breaks= -c(10, 8, 6))+
  theme_minimal(12)+
  labs(col = '2ln(LR)', x = '', y = '')








obs_dat = read_csv('data/processed/obs_data.csv') %>%
  mutate(excess = maxtp - threshold_9) %>%
  filter(excess>0)
  





loglik_sum = c()

for(potential_shape in seq(-0.2, -0.15, length.out = 100)){
  print(potential_shape)
  
  loglik = c()
  
  for(i in (obs_dat$Station %>% unique())){
    this_extm = obs_dat %>% filter(Station == i) %>% pull(excess)
    
    model_fit = estiamte_scale_fixed_shape(this_extm, potential_shape)
    
    loglik = c(loglik, model_fit$value)
    
  }
  loglik_sum = c(loglik_sum, sum(loglik))
}

optimal_shape = seq(-0.2, -0.15, length.out = 100)[which.min(loglik_sum)]
optimal_shape = -0.1737374



loglik_fixed_xi = c()
loglik_free_xi = c()

obs_grid = obs_dat %>% dplyr::select(Long, Lat, Station) %>% unique
for(i in obs_grid$Station){
  print(i)
  this_extm_irel = obs_dat %>% filter(Station == i) %>% pull(excess)
  
  # --- fixed shape
  model_fit = estiamte_scale_fixed_shape(this_extm_irel, -0.1737374)
  loglik_fixed_xi = c(loglik_fixed_xi, model_fit$value)
  
  # --- free shape
  loglik_free_xi = c(loglik_free_xi, (evd::fpot(this_extm_irel, threshold = 0) %>% logLik %>% as.numeric()))
}





obs_grid$LL_null = loglik_fixed_xi
obs_grid$LL_alt = -loglik_free_xi
obs_grid$LL_ratio = -2*log(obs_grid$LL_alt/obs_grid$LL_null)


plt_obs = obs_grid %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = LL_ratio))+
  coord_map()+
  theme_minimal()+
  geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  scale_x_continuous(breaks= -c(10, 8, 6))+
  theme_minimal(12)+
  labs(col = '2ln(LR)', x = '', y = '')


plt = gridExtra::grid.arrange(plt_clim, plt_obs, nrow = 1, widths = c(1, 0.95))






ggsave(plot = plt, filename = 'output/figs/LRT.pdf', width = 7, height =3)


