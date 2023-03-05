gc()
rm(list = ls())
library(tidyverse)

marg_mod = "mod_1"
setwd("~/Extreme-Irish-Summer-Temperatures/")
source('src/models/marginal_models/gpd_models.R')

obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
  unique() %>%
  left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))

obs_grid = read_csv("data/processed/thresh_exceedance_lambda_num_quantiles_15.csv") %>%
  left_join(obs_sites) %>%
  left_join(read_csv("data/processed/obs_data.csv") %>%
              dplyr::select(year, loess_temp_anom) %>% unique)

if(marg_mod == "mod_0"){
  # ----- model 0
  obs_grid$scale = my_predict_0(as.numeric(unlist(read_csv("output/gpd_model_fits/model_0_true.csv"))), obs_grid$scale_9)$scale
  obs_grid$shape = my_predict_0(as.numeric(unlist(read_csv("output/gpd_model_fits/model_0_true.csv"))), obs_grid$scale_9)$shape
}else if(marg_mod == "mod_1"){
  # # ----- model 1
  obs_grid$scale = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/model_1_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom)$scale
  obs_grid$shape = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/model_1_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom)$shape
}else if(marg_mod == "mod_2"){
  # # ----- model 2
  obs_grid$scale = my_predict_2(as.numeric(unlist(read_csv("output/gpd_model_fits/model_2_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$scale
  obs_grid$shape = my_predict_2(as.numeric(unlist(read_csv("output/gpd_model_fits/model_2_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$shape
}

extreme_temps = seq(22.9, 36, by = 0.25)

res = c()
for(tmp in extreme_temps){
  res = cbind(res,(-1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - obs_grid$threshold_9),  scale = obs_grid$scale, shape =  obs_grid$shape[1])))))
  }

res = as_tibble(res)
names(res) = as.character(extreme_temps)

res %>%
  mutate(year = obs_grid$year,
         Station = obs_grid$Station) %>%
  pivot_longer(-c(year, Station)) %>%
  rename(temp = name,
         frechet_value = value) %>%
  mutate(temp = as.numeric(temp)) %>%
  write_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_",marg_mod,".csv"))


# # ---- Climate data ----------------------
obs_smoothed_quantiles = readRDS("output/quant_models_clim_num_quantiles_20.csv")

clim_grid = read_csv("data/processed/clim_grid_simulated_on_not_projected.csv") %>%
  left_join(read_csv('~/JRSS_organised_code/data/processed_data/clim_data_dist_to_sea.csv') %>%
              dplyr::select(Long, Lat, dist_sea)) %>%
  left_join(read_csv("data/processed/clim_scale_grid.csv") %>% 
              dplyr::select(Long, Lat, scale_9))

clim_grid = rbind(clim_grid %>% mutate(year = 2020), clim_grid %>% mutate(year = 1942)) %>%
  left_join(read_csv("data/processed/obs_data.csv") %>%
              dplyr::select(year, loess_temp_anom) %>%
              filter(year %in% c(2020, 1942)) %>% 
              unique()) %>%
  left_join(obs_smoothed_quantiles %>%
              dplyr::select(Long, Lat, tau_to_temp, threshold_9, thresh_exceedance_9, year))

marg_mod = 'mod_2'
this_fit_mod = read_csv("output/gpd_model_fits/model_2_true.csv") %>%
  unlist() %>% as.numeric
clim_grid$shape = this_fit_mod[length(this_fit_mod)]
clim_grid$scale = my_predict_2(this_fit_mod, clim_grid$scale_9, clim_grid$loess_temp_anom, clim_grid$dist_sea)$scale
extreme_temps = seq(26, 36, by = 0.25)

res = c()
for(tmp in extreme_temps){
  res = cbind(res,(-1/log(1-clim_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - clim_grid$threshold_9),  scale = clim_grid$scale, shape =  clim_grid$shape[1])))))
}

res = as_tibble(res)
names(res) = as.character(extreme_temps)

res %>%
  mutate(year = clim_grid$year,
         Long = clim_grid$Long,
         Lat = clim_grid$Lat) %>%
  pivot_longer(-c(year, Long, Lat)) %>%
  rename(temp = name,
         frechet_value = value) %>%
  mutate(temp = as.numeric(temp)) %>%
  write_csv(paste0("output/clim_grid_extreme_temps_frechet_scale_",marg_mod,".csv"))

