# find extreme temperatures on frechet scale



gc()
rm(list = ls())
library(tidyverse)

marg_mod = "mod_2"
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland//")
source('src/models/marginal_models/gpd_models.R')


obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
  unique() %>%
  left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))


obs_grid = read_csv("data/processed/thresh_exceedance_lambda_num_quantiles_25.csv") %>%
  left_join(obs_sites) %>%
  left_join(read_csv("data/processed/obs_data.csv") %>%
              dplyr::select(year, loess_temp_anom) %>% unique)

if(marg_mod == "mod_1"){
  # ----- model 1
  obs_grid$scale = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/model_1_true.csv"))), obs_grid$scale_9)$scale
  obs_grid$shape = my_predict_1(as.numeric(unlist(read_csv("output/gpd_model_fits/model_1_true.csv"))), obs_grid$scale_9)$shape
}else if(marg_mod == "mod_2"){
  # # ----- model 2
  obs_grid$scale = my_predict_2(as.numeric(unlist(read_csv("output/gpd_model_fits/model_2_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom)$scale
  obs_grid$shape = my_predict_2(as.numeric(unlist(read_csv("output/gpd_model_fits/model_2_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom)$shape
}else if(marg_mod == "mod_4"){
  # # ----- model 4b
  obs_grid$scale = my_predict_4b(as.numeric(unlist(read_csv("output/gpd_model_fits/model_4_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$scale
  obs_grid$shape = my_predict_4b(as.numeric(unlist(read_csv("output/gpd_model_fits/model_4_true.csv"))), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$shape
}


extreme_temps = seq(27, 36, by = 0.25)

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


# 
# # ---- Climate data
# 
# 
# 
# 
# clim_grid = rbind(read_csv("data/processed_data/clim_data_dist_to_sea.csv") %>% 
#                     dplyr::select(Long, Lat, id, Long.projected, Lat.projected, dist_sea) %>%
#                     mutate(year = 1942),
#                   read_csv("data/processed_data/clim_data_dist_to_sea.csv") %>% 
#                     dplyr::select(Long, Lat, id, Long.projected, Lat.projected, dist_sea) %>%
#                     mutate(year = 2020)) %>%
#   left_join(temp_data %>%
#               dplyr::select(loess_temp_anom,  year) %>%
#               unique()) %>%
#   left_join(readRDS("corrections/marginal_models/bulk/quant_models_clim") %>%
#               dplyr::select(id, Lat, Long, year, Long.projected, Lat.projected, threshold_9, thresh_exceedance_9)) %>%
#   left_join(read_csv("corrections/clim_scale_grid.csv"))
# 
# 
# 
# if(marg_mod == "mod_1"){
#   # ----- model 1
#   gpd_pars = read_csv("model_1_true_pars.csv",
#                       col_names = c('beta.0', 'beta.1', 'beta.2'))
#   clim_grid$scale = my_predict_1(unlist(gpd_pars), clim_grid$scale_9)$scale
#   clim_grid$shape = my_predict_1(unlist(gpd_pars), clim_grid$scale_9)$shape
# }else if(marg_mod == "mod_2"){
#   # # ----- model 2
#   gpd_pars = read_csv("model_2_true_pars.csv",
#                       col_names = c('beta.0', 'beta.1', 'beta.2', 'beta.3'))
#   clim_grid$scale = my_predict_2(unlist(gpd_pars), clim_grid$scale_9, clim_grid$loess_temp_anom)$scale
#   clim_grid$shape = my_predict_2(unlist(gpd_pars), clim_grid$scale_9, clim_grid$loess_temp_anom)$shape
# }else if(marg_mod == "mod_4"){
#   # # ----- model 4b
#   gpd_pars = read_csv("model_4b_true_pars.csv",
#                       col_names = c('beta.0', 'beta.1', 'beta.2', 'beta.3', 'beta.4', 'beta.5'))
#   clim_grid$scale = my_predict_4b(unlist(gpd_pars), clim_grid$scale_9, clim_grid$loess_temp_anom, clim_grid$dist_sea)$scale
#   clim_grid$shape = my_predict_4b(unlist(gpd_pars), clim_grid$scale_9, clim_grid$loess_temp_anom, clim_grid$dist_sea)$shape
# }
# 
# 
# 
# 
# extreme_temps = seq(25, 36, by = 0.25)
# 
# res = c()
# for(tmp in extreme_temps){
#   res = cbind(res,(-1/log(1-clim_grid$thresh_exceedance_9*(1-evd::pgpd(q = (tmp - clim_grid$threshold_9),  scale = clim_grid$scale, shape =  clim_grid$shape[1])))))
# }
# 
# res = as_tibble(res)
# names(res) = as.character(extreme_temps)
# 
# res %>%
#   mutate(year = clim_grid$year,
#          id = clim_grid$id) %>%
#   pivot_longer(-c(year, id)) %>%
#   rename(temp = name,
#          frechet_value = value) %>%
#   mutate(temp = as.numeric(temp)) %>%
#   write_csv(paste0("clim_grid_extreme_temps_frechet_scale_",marg_mod,".csv"))
# 
# 
# 
