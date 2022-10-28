# simulate under true model
# ---- similate on cliscale grid
rm(list=ls())
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")

library(tidyverse)
library(doParallel)
library(foreach)

#-------- getting grid
# locs_to_pred = as.data.frame(read_csv("data/processed_data/clim_data_dist_to_sea.csv") %>% dplyr::select(Long.projected, Lat.projected))
simulate_this_mod = function(bts_nums, marg_mod){

  locs_to_pred = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv"))
  data_sets = list.files("output/simulations/simulations_on_obs_grid/bootstraps/")
  data_sets = data_sets[which(grepl(marg_mod, data_sets))]
  
  bts_ran = sub(paste0(marg_mod,"_bootstrap_"), "", data_sets) %>%
    map(~{
      pos_start = .x %>% str_locate("_") %>% .[,1]
      bts_num = .x %>% str_sub(end = (pos_start-1)) }) %>% 
    unlist %>% 
    unique %>% 
    as.numeric()
  
  #bts_nums = bts_nums[!(bts_nums%in%bts_ran)]
  variogram_model <- function(h){
    h = sqrt(norm(h,type = "2")^2)
    nu=0.2
    res = rep(NA, length(h))
    res[h == 0] = 0
    res[h != 0] = alpha_var*(1 - ((((h[h != 0] /beta_var)^nu) * besselK(x = (h[h != 0] / beta_var), nu = nu))/(gamma(nu)*2^(nu - 1))))
    res
  }
  
  my_simulate <- function(file_name){
    print(file_name)
    for(i in seq(1,100)){
      nCores <- 4
      cl <- parallel::makeCluster(nCores)
      clusterSetRNGStream(cl)
      
      simulation_score <- mvPot::simulPareto(n = 100,
                                             loc = locs_to_pred,
                                             vario = variogram_model,
                                             nCores = nCores,
                                             cl = cl)
      parallel::stopCluster(cl)
      
      
      file.remove(paste0(file_name,i))
      simulation_score %>%
        saveRDS(paste0(file_name,i))
    }
  }
  
  for(x in bts_nums){
    print(x)
    # stationary_models = read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", marg_mod, ".csv-Mead.csv"),
    #                              col_names = c('bts', "alpha", "beta"))
    
    # --- model 1+2
    stationary_models = read_csv(paste0("output/rpareto_model_fits/bts_rpareto_fits_model_", marg_mod, "_L-BFGS-B.csv"),
                                 col_names = c('bts', "alpha", "beta"))
    
    
    this_rparp_mod = stationary_models %>% filter(bts == x)
    alpha_var <- this_rparp_mod$alpha[nrow(this_rparp_mod)]
    beta_var <- this_rparp_mod$beta[nrow(this_rparp_mod)]
    if(length(alpha_var)>0){
      
      my_simulate(paste0("output/simulations/simulations_on_obs_grid/bootstraps/",marg_mod,"_bootstrap_",this_rparp_mod$bts[nrow(this_rparp_mod)],"_run_"))
    }
  }
}


# job::job({simulate_this_mod(seq(1, 25),  'mod_4')})
# job::job({simulate_this_mod(seq(26, 50),  'mod_4')})
# job::job({simulate_this_mod(seq(51, 75),  'mod_4')})
# job::job({simulate_this_mod(seq(76, 100),  'mod_4')})
# 
# job::job({simulate_this_mod(seq(51, 100),  'mod_4')})
# job::job({simulate_this_mod(seq(101, 150), 'mod_4')})
# job::job({simulate_this_mod(seq(151, 200), 'mod_4')})
# job::job({simulate_this_mod(seq(201, 250), 'mod_4')})
# job::job({simulate_this_mod(seq(251, 300), 'mod_4')})
# 
# job::job({simulate_this_mod(seq(300, 350),  'mod_4')})
# job::job({simulate_this_mod(seq(351, 400), 'mod_4')})
# job::job({simulate_this_mod(seq(401, 450), 'mod_4')})
# job::job({simulate_this_mod(seq(451, 500), 'mod_4')})


job::job({simulate_this_mod(c(376),  'mod_2')})
job::job({simulate_this_mod(c(328, 329, 428, 429),  'mod_2')})

# job::job({simulate_this_mod(seq(301, 350),  'mod_1')})
# job::job({simulate_this_mod(seq(351, 400), 'mod_1')})
# job::job({simulate_this_mod(seq(401, 450), 'mod_1')})
# job::job({simulate_this_mod(seq(451, 500), 'mod_1')})
# 
# job::job({simulate_this_mod(seq(301, 350),  'mod_2')})
# job::job({simulate_this_mod(seq(351, 400), 'mod_2')})
# job::job({simulate_this_mod(seq(401, 450), 'mod_2')})
# job::job({simulate_this_mod(seq(451, 500), 'mod_2')})

