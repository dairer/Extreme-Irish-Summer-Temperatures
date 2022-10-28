# scale up simulations
rm(list=ls())
library(tidyverse)
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
source("corrections/marginal_models/tail/gpd_models.R")


read_csv("output/simulations/simulations_on_obs_grid/true/")
scale_true_sim = function(marg_mod, yr, tmp){
  extreme_temp_frechet = read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_",marg_mod,".csv"))
  data_sets = list.files('output/simulations/simulations_on_obs_grid/true/')
  data_sets = data_sets[which(grepl(marg_mod, data_sets))]
  
  obs_sites = read_csv("data/processed/obs_data.csv") %>%
    dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
    unique() %>%
    left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))

  grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
    left_join(obs_sites)
  
  # ---- for each model get scale shape and lambda 
  this_grid = grid_simulated %>%
    left_join(read_csv("data/processed/thresh_exceedance_lambda_num_quantiles_25.csv") %>% filter(year == yr)) %>%
    left_join(extreme_temp_frechet %>% filter(temp == tmp))
  

  my_simulations_extremes = c()
  for(i in seq(1, 100)){ 
    my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true/",marg_mod,"_run_",i)))
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
    #my_simulations_rescaled[[i]] = (my_r[i] * this_grid$scaler)*my_simulations_standardised[[i]]
    my_simulations_rescaled[[i]] = (my_r[i] * br)*my_simulations_standardised[[i]]
  }
  
  # --- model 1 + 2
  my_simulations_rescaled %>% saveRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp))
}

job::job({scale_true_sim("mod_4", 1942, 27)})
job::job({scale_true_sim("mod_4", 2020, 27)})
job::job({scale_true_sim("mod_4", 1942, 28)})
job::job({scale_true_sim("mod_4", 2020, 28)})
job::job({scale_true_sim("mod_4", 1942, 29)})
job::job({scale_true_sim("mod_4", 2020, 29)})
job::job({scale_true_sim("mod_4", 1942, 30)})
job::job({scale_true_sim("mod_4", 2020, 30)})
job::job({scale_true_sim("mod_4", 1942, 31)})
job::job({scale_true_sim("mod_4", 2020, 31)})
job::job({scale_true_sim("mod_4", 1942, 32)})
job::job({scale_true_sim("mod_4", 2020, 32)})
job::job({scale_true_sim("mod_4", 1942, 33)})
job::job({scale_true_sim("mod_4", 2020, 33)})
job::job({scale_true_sim("mod_4", 1942, 34)})
job::job({scale_true_sim("mod_4", 2020, 34)})
job::job({scale_true_sim("mod_4", 1942, 35)})
job::job({scale_true_sim("mod_4", 2020, 35)})

job::job({scale_true_sim("mod_2", 1942, 27)})
job::job({scale_true_sim("mod_2", 2020, 27)})
job::job({scale_true_sim("mod_2", 1942, 28)})
job::job({scale_true_sim("mod_2", 2020, 28)})
job::job({scale_true_sim("mod_2", 1942, 29)})
job::job({scale_true_sim("mod_2", 2020, 29)})
job::job({scale_true_sim("mod_2", 1942, 30)})
job::job({scale_true_sim("mod_2", 2020, 30)})
job::job({scale_true_sim("mod_2", 1942, 31)})
job::job({scale_true_sim("mod_2", 2020, 31)})
job::job({scale_true_sim("mod_2", 1942, 32)})
job::job({scale_true_sim("mod_2", 2020, 32)})
job::job({scale_true_sim("mod_2", 1942, 33)})
job::job({scale_true_sim("mod_2", 2020, 33)})
job::job({scale_true_sim("mod_2", 1942, 34)})
job::job({scale_true_sim("mod_2", 2020, 34)})
job::job({scale_true_sim("mod_2", 1942, 35)})
job::job({scale_true_sim("mod_2", 2020, 35)})
# 
mod_2