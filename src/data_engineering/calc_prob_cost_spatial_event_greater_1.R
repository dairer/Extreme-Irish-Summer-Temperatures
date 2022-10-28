rm(list = ls())
library(tidyverse)

setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")

marg_mod = 'mod_1'
source('src/models/marginal_models/gpd_models.R')
obs_data = vroom::vroom("data/processed/obs_data.csv") %>%
  left_join(vroom::vroom("data/processed/thresh_exceedance_lambda_num_quantiles_50.csv"))  %>%
  left_join(vroom::vroom("data/processed/obs_data_dist_to_sea.csv") %>% 
              dplyr::select(Station, dist_sea), by = 'Station')


obs_smoothed_quantiles = readRDS("output/quant_models_num_quantiles_50.csv")

# pull date which have 'conditioning' site
dates_to_keep = obs_data %>%
  group_by(date) %>%
  summarise(check_obs_for_s = 'aldergrove' %in% Station) %>%
  filter(check_obs_for_s) %>%
  pull(date)

obs_data = obs_data %>% filter(date %in% dates_to_keep)

if(marg_mod == 'mod_1'){
  this_fit_mod = read_csv("output/gpd_model_fits/model_1_true.csv") %>%
    unlist() %>% as.numeric
  pred <<- my_predict_1(this_fit_mod, obs_data$scale_9)
  
}else if(marg_mod == 'mod_2'){
  
  this_fit_mod = read_csv("output/gpd_model_fits/model_2_true.csv") %>%
    unlist() %>% as.numeric
  pred <<- my_predict_2(this_fit_mod, obs_data$scale_9, obs_data$loess_temp_anom)
  
}else if(marg_mod == 'mod_4'){
  this_fit_mod = read_csv("output/gpd_model_fits/model_4_true.csv") %>%
    unlist() %>% as.numeric
  pred <<- my_predict_4b(this_fit_mod, obs_data$scale_9, obs_data$loess_temp_anom, obs_data$dist_sea)
  
}

obs_data$scale = pred$scale
obs_data$shape = pred$shape

obs_data_standardised = obs_data %>%
  group_by(Station, year) %>%
  group_map(~{
    data = .x$maxtp
    threshold = .x$threshold_9
    res = rep(NA, length(data))
    num_extremes = sum(data>threshold)
    
    if(num_extremes >0){ # if extreme obs in this year at this site
      scle = .x$scale
      shpe = .x$shape
      my_lambda = .x$thresh_exceedance_9
      
      res[data > threshold] = 1 - my_lambda[data > threshold]*(1-evd::pgpd((data[data > threshold] - threshold[data > threshold]), loc = 0,scale = scle[data > threshold], shape = shpe[1]))
    }
    
    this_quant_mod = obs_smoothed_quantiles %>%
      filter(Station == .x$Station[1],
             year == .x$year[1]) %>%
      pull(temp_to_tau) %>%
      .[[1]]
    
    res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])       # alternative ... myecdf(data[data <= threshold])
    res[res<0]=0
    
    .x$unif = res
    .x$frechet_marg = -1/log(.x$unif)
    .x
  },.keep = T) %>%
  plyr::rbind.fill()%>%
  as_tibble() 


# -------- MULTIVARIATE EVENTS
# get cost of each event
extreme_dates = obs_data_standardised %>%
  dplyr::group_by(date) %>%
  dplyr::summarise(cost = mean(frechet_marg)) %>%
  ungroup() %>%
  arrange(desc(cost))


(sum(extreme_dates$cost > 1)/nrow(extreme_dates)) %>%
  saveRDS(paste0("data/processed/prob_event_greater_1_", marg_mod))
