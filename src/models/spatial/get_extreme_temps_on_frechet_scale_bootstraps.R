# find extreme temperatures on frechet scale



gc()
rm(list = ls())
library(tidyverse)

marg_mod = "mod_2"
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
source('src/models/marginal_models/gpd_models.R')

obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
  unique() %>%
  left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))

# --- needs to change with the bootstrap
for(this_bts in seq(1, 500)){
  print(this_bts)
  obs_grid = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",15,"_bts_",this_bts,".csv")) %>%
    rename(bts = X1,
           Station = X2,
           year = X3,
           thresh_exceedance_9 = X4) %>% dplyr::select(-bts) %>%
    left_join(obs_sites) %>%
    left_join(read_csv("data/processed/obs_data.csv") %>%
                dplyr::select(year, loess_temp_anom) %>% unique)
  
  if(marg_mod == "mod_1"){
    # ----- model 1
  
    pars = read_csv("output/gpd_model_fits/model_1_pars_corrected.csv",
                    col_names = c("bts", "b0", "b1", "b2")) %>%
      mutate(bts = (as.numeric(str_remove_all(bts, "num_quantiles_50_bts_")))) %>%
      filter(bts == this_bts) %>% unlist %>% .[2:length(.)] %>% as.numeric
    
    obs_grid$scale = my_predict_1(pars, obs_grid$scale_9)$scale
    obs_grid$shape = my_predict_1(pars, obs_grid$scale_9)$shape
    
  }else if(marg_mod == "mod_2"){
    # # ----- model 2
    pars = read_csv("output/gpd_model_fits/model_2_pars_corrected.csv",
                    col_names = c("bts", "b0", "b1", "b2", "b3")) %>%
      mutate(bts = (as.numeric(str_remove_all(bts, "num_quantiles_50_bts_")))) %>%
      filter(bts == this_bts) %>% unlist %>% .[2:length(.)] %>% as.numeric
    
    obs_grid$scale = my_predict_2(pars, obs_grid$scale_9, obs_grid$loess_temp_anom)$scale
    obs_grid$shape = my_predict_2(pars, obs_grid$scale_9, obs_grid$loess_temp_anom)$shape
    
  }else if(marg_mod == "mod_4"){
    # # ----- model 4b
    
    pars = read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
                    col_names = c("bts", "b0", "b1", "b2", "b3", "b4", "b5")) %>% 
      mutate(bts =  gsub("num_quantiles_50_bts_","", bts) %>% as.numeric()) %>%
      filter(bts == this_bts) %>%
      unlist() %>% as.numeric %>% .[2:length(.)]
    
    obs_grid$scale = my_predict_4b(pars, obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$scale
    obs_grid$shape = my_predict_4b(pars, obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$shape
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
           Station = obs_grid$Station, 
           bts = this_bts) %>%
    pivot_longer(-c(year, Station, bts)) %>%
    rename(temp = name,
           frechet_value = value) %>%
    mutate(temp = as.numeric(temp)) %>%
    write_csv(paste0("output/extreme_frechet_values_bts/obs_sites_extreme_temps_bootstraps_frechet_scale_",marg_mod,"_bts_",this_bts,".csv"), append = T)
}





# frechet_vals_old = c()
# frechet_vals_new = c()
# for(bts in seq(300)){
#   dat = vroom::vroom(paste0("output/extreme_frechet_values_bts/obs_sites_extreme_temps_bootstraps_frechet_scale_mod_4_bts_",bts,".csv"),
#                      col_names = c('year', 'Station', 'bts', 'temp', 'frechet_value'))
#   
#   frechet_vals_old = c(frechet_vals_old, dat %>% filter(year == 1940, Station == "casement", temp == 27) %>% pull(frechet_value) )
#   frechet_vals_new = c(frechet_vals_new, dat %>% filter(year == 2020, Station == "casement", temp == 27) %>% pull(frechet_value) )
# }
# 
# 
# 
# 
