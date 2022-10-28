# simulate under true model
# ---- similate on clim scale grid
rm(list=ls())
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
library(tidyverse)
library(doParallel)
library(foreach)

# --- marginal model names
marg_mod = 'mod_1'
read_csv("data/processed/obs_data.csv") %>% dplyr::select(Long.projected, Lat.projected) %>% unique %>% write_csv("data/processed/obs_grid_simulated_on.csv")
locs_to_pred = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv"))

variogram_model = function(h){
  h = sqrt(norm(h,type = "2")^2)
  nu=0.2
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha_var*(1 - ((((h[h != 0] /beta_var)^nu) * besselK(x = (h[h != 0] / beta_var), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}

simulate = function(file_name){
  for(i in seq(1, 100)){
    print(i)
    nCores <- 10
    cl <- parallel::makeCluster(nCores)
    clusterSetRNGStream(cl)
    
    # this returns samples in the unit frechet margin
    simulation_score <- mvPot::simulPareto(n = 2500,
                                           loc = locs_to_pred,
                                           vario = variogram_model,
                                           nCores = nCores,
                                           cl = cl)
    parallel::stopCluster(cl)
    
    simulation_score %>%
      saveRDS(paste0(file_name,i))
  }
}

fit = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_",marg_mod, ".csv"),
                     col_names = F) %>% 
  .[nrow(.),] %>% # get last row in this file
  unlist %>% as.numeric()


alpha_var <<- fit[1]
beta_var <<- fit[2]

simulate(paste0("output/simulations/simulations_on_obs_grid/true/",marg_mod,"_run_"))