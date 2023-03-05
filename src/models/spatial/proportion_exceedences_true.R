rm(list=ls())
library(tidyverse)
marg_mod = 'mod_2'
temp_conditioned_on = 27


setwd("~/Extreme-Irish-Summer-Temperatures/")
obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station, Long.projected, Lat.projected) %>% 
  unique

yr = 2020

if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/",marg_mod,"_yr_",yr,"_min_temp_conditioned_on_", temp_conditioned_on))){
  
  my_simulations  = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/",marg_mod,"_yr_",yr,"_min_temp_conditioned_on_", temp_conditioned_on))
  prop_ex_2020 = c()
    for(tmp in c(temp_conditioned_on)){
      
      frechet_val = read_csv("data/processed/obs_grid_simulated_on.csv") %>%
        left_join(obs_sites) %>%
        left_join(read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_",marg_mod,".csv")) %>% filter(temp == tmp, year == yr)) %>%
        pull(frechet_value)
  
      frechet_val[frechet_val == -Inf] = Inf
      
      num_exceed_tmp = my_simulations %>%
        map(~{ sum(.x > frechet_val)}) %>%
        unlist()
      
      prop_ex_2020 = c(prop_ex_2020, mean(num_exceed_tmp[num_exceed_tmp>0]/182))
    }
  
  yr = 1942
  my_simulations  = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/",marg_mod,"_yr_",yr,"_min_temp_conditioned_on_", temp_conditioned_on))
  prop_ex_1942 = c()
    for(tmp in c(temp_conditioned_on)){
      
    frechet_val = read_csv("data/processed/obs_grid_simulated_on.csv") %>%
      left_join(obs_sites) %>%
      left_join(read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_",marg_mod,".csv")) %>% filter(temp == tmp, year == yr)) %>%
      pull(frechet_value)
    
    frechet_val[frechet_val == -Inf] = Inf
    
    num_exceed_tmp = my_simulations %>% 
      map(~{ sum(.x > frechet_val)}) %>%
      unlist()
    
    prop_ex_1942 = c(prop_ex_1942, mean(num_exceed_tmp[num_exceed_tmp>0]/182))
  }
  
  tibble(temp = temp_conditioned_on, prop_ex_1942, prop_ex_2020) %>%
    write_csv(paste0("output/simulations/simulation_summary/prop_exceedance_model_true_",marg_mod,"_temp_condtioned_on_",temp_conditioned_on,".csv"))
}
