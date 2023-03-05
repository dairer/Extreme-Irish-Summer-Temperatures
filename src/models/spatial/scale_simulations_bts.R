# scale up simulations
rm(list=ls())
library(tidyverse)


scale_bts_sim = function(bts_range, marg_mod, yr, num_quantiles, tmp){
  
  setwd("~/Extreme-Irish-Summer-Temperatures/")
  source("src/models/marginal_models/gpd_models.R")

  obs_sites = read_csv("data/processed/obs_data.csv") %>%
    dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
    unique() %>%
    left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))
  
  
  grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
    left_join(obs_sites)
  
  # ---- for each model get scale shape and lambda 

  for(bts in bts_range){

    print(bts)
    extreme_temp_frechet = vroom::vroom(paste0("output/extreme_frechet_values_bts/obs_sites_extreme_temps_bootstraps_frechet_scale_",marg_mod,"_bts_",bts,".csv"),
                                        col_names = c('year', 'Station', 'bts', 'temp', 'frechet_value')) %>%
      dplyr::select(-bts)

    # ---- for model 1 + 2
    this_grid = grid_simulated %>%
      as_tibble() %>%
      left_join(read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",bts,".csv")) %>%
                  rename(bts = X1,
                         Station = X2,
                         year = X3,
                         thresh_exceedance_9 = X4) %>%
                  dplyr::select(-bts) %>% filter(year == yr)) %>%
      left_join(extreme_temp_frechet %>% filter(temp == tmp))
    
    
    
    
    data_sets = list.files('output/simulations/simulations_on_obs_grid/bootstraps/')
    data_sets = data_sets[which(grepl(paste0(marg_mod,"_bootstrap_",bts,"_"), data_sets))]
    
    if(length(data_sets)>0){
      

      my_simulations_extremes = c()
      for(i in seq(1, 100)){
        my_simulations_extremes = c(my_simulations_extremes, 
                                    readRDS(paste0("output/simulations/simulations_on_obs_grid/bootstraps/",marg_mod,"_bootstrap_",bts,"_run_",i)))
      }
      my_simulations_standardised = list()
      for (i in 1:length(my_simulations_extremes)) {
        this_cost = mean(my_simulations_extremes[[i]])
        my_simulations_standardised[[i]] = my_simulations_extremes[[i]]/this_cost
      }
      # --- maximum frechet margin at each of my simulation sites
      max_at_each_site = c()
      for(s in seq(length(my_simulations_standardised[[1]]))){
        max_at_each_site = c(max_at_each_site, lapply(my_simulations_standardised, "[[", s) %>% unlist %>% max)
      }
      
      this_grid$scaler = this_grid$frechet_value/max_at_each_site
      this_grid$scaler[this_grid$scaler == -Inf] = Inf
      br = min(this_grid$scaler)
      my_r = evd::rgpd(n=length(my_simulations_extremes), 1,1,1)
      
      # --- scale simulations back up
      my_simulations_rescaled = my_simulations_standardised
      for (i in 1:length(my_simulations_standardised)) {
        my_simulations_rescaled[[i]] = (my_r[i] * br)*my_simulations_standardised[[i]]
        }
      
      my_simulations_rescaled %>% saveRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bootstraps/bootstrap_",bts,"_model_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp))
    }
  }
}

this_temp_to_scale_upto = 30
# job::job({scale_bts_sim(bts_range = seq(1, 300), "mod_4", 1942,30, this_temp_to_scale_upto)})
# job::job({scale_bts_sim(bts_range = seq(1,300), "mod_4", 2020,30, this_temp_to_scale_upto)})


