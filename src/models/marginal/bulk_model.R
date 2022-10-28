# This script fits and saves quantile regression model and lambda estimates

gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("~/Extreme-Irish-Summer-Temperatures/")

num_quantiles = 25
obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))

# ---- get covariates for prediction
temporal_covariates = obs_data %>%
  dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
  unique() %>%
  arrange(year)

quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
fit_clim_quants = T # bool, re-estimate clim quantiles?

if(fit_clim_quants){
  file.remove(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"))
  for(q in seq_along(quantiles_to_estimate_bulk)){
    #obs_data_for_quant_reg = obs_data
    obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
    zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
    print(paste0("Fitting  tau = ", zeta))

    obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
    qunatile_model_fit <- evgam(maxtp ~ value + loess_temp_anom, obs_data_for_quant_reg,
                                  family = "ald", ald.args = list(tau = zeta))
    
    # save parameter estimates
    tibble(tau = zeta,
           beta_0 = qunatile_model_fit$location$coefficients[1],
           beta_1 = qunatile_model_fit$location$coefficients[2],
           beta_2 = qunatile_model_fit$location$coefficients[3]) %>%
      write_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"), append = T)
  }
}

# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
         col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))

# # --- creates a tibble with each station and its quantule model
obs_smoothed_quantiles = obs_data %>%
  group_by(Station) %>%
  group_map(~{

    # --- get the climate qunatile estimates closest to current station
    clim_vals <<- obs_data %>%
      filter(Station == .x$Station[1]) %>%
      dplyr::select(quantile, value) %>%
      unique() %>%
      pull(value) %>%
      unlist()

    # predict quantile for each year and site
    quant_reg_pars = quant_reg_model_pars %>%
      arrange(tau)

    res = c()
    for(q in seq_along(quantiles_to_estimate_bulk)){
      qpars = quant_reg_pars[q,]
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates$year,
                         loess_temp_anom = temporal_covariates$loess_temp_anom,
                         quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(temporal_covariates$loess_temp_anom)))
    }

    print(paste0("Interpolating quantile estimates for ", .x$Station[1]))
    # interpolate quantiles over tau for each year
    res %>%
      group_by(year) %>%
      group_map(~{
        tibble(year = .x$year[1],
               tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
               temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      mutate(Station = .x$Station[1])

  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()



# save quantile models
obs_smoothed_quantiles %>% saveRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))
#obs_smoothed_quantiles = readRDS("output/quant_models")



# Calculate lambda
lambda_thresh_ex = obs_data %>%
  group_by(Station) %>%
  group_map(~{

    print(.x$Station[1])

    print(obs_smoothed_quantiles%>%
            filter(Station == .x$Station[1]))

    thresh_exceedance_9 = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_9[1], x))


    tibble(Station = .x$Station[1],
           year = temporal_covariates$year,
           thresh_exceedance_9 = 1-thresh_exceedance_9)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()




lambda_thresh_ex %>%
  write_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))


# ------------ get splines on clim scale

clim_grid = read_csv("data/processed/clim_scale_grid.csv")
clim_dat_full = read_csv("~/JRSS_organised_code/corrections/data/processed_data/full_clim_data.csv")

# estimate empiracle quantules for climate data
clim_quantiles = clim_dat_full %>%
  group_by(id) %>%
  group_map(~{
    res = c()
    for(q in quantiles_to_estimate_bulk){
      res = rbind(res, tibble(id = .x$id[1],
                              Long = .x$Long[1],
                              Lat = .x$Lat[1],
                              Long.projected = .x$Long.projected[1],
                              Lat.projected = .x$Lat.projected[1],
                              quantile = q,
                              value = as.numeric(quantile(.x$maxtp, q))))
    }
    res
  }, .keep = T)%>%
  plyr::rbind.fill() %>%
  as.tibble()

clim_quantiles_subset = clim_quantiles %>%
  # filter(id %in% obs_data$id) %>%
  group_by(id) %>%
  group_map(~{
    
    tibble(id = .x$id[1], quantile = list(.x$quantile), value = list(.x$value))
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


clim_grid = clim_grid %>%
  left_join(clim_quantiles_subset)

clim_date_w_quantile_mod = c()


# here
for(i in seq(nrow(clim_grid))){
  print(clim_grid[i,])

  this_data_for_pred = tibble(year = c(1931, 1942, 2020, 2022)) %>%
    left_join(temporal_covariates) %>%
    mutate(Long = clim_grid[i,]$Long,
           Lat = clim_grid[i,]$Lat,
           Long.projected = clim_grid[i,]$Long.projected,
           Lat.projected = clim_grid[i,]$Lat.projected,
           id = clim_grid[i,]$id,
           # threshold = clim_grid[i,]$threshold,
           # clim_scale = clim_grid[i,]$clim_scale,
           quantile = clim_grid[i,]$quantile,
           value = clim_grid[i,]$value) 
  
  # predict quantile for each year and site
  quant_reg_pars = quant_reg_model_pars %>%
    arrange(tau)

  clim_vals = clim_grid[i,]$value[[1]]
  res = c()
  for(q in seq_along(quantiles_to_estimate_bulk)){
    qpars = quant_reg_pars[q,]
    res = rbind(res,
                tibble(quantile =  qpars$tau,
                       year = this_data_for_pred$year,
                       loess_temp_anom = this_data_for_pred$loess_temp_anom,
                       quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + (qpars$beta_2)*(this_data_for_pred$loess_temp_anom)))
  }
  
  res = res %>%
    group_by(year) %>%
    group_map(~{
      tibble(year = .x$year[1],
             tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
             temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble() 
  
  clim_date_w_quantile_mod = rbind(clim_date_w_quantile_mod, 
                                   this_data_for_pred %>% left_join(res))

}

clim_date_w_quantile_mod %>%
  saveRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))

# thresold at these points 

qunatile_model_fit_9 = readRDS("output/threshold_model_9")
clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))

# get clim threshold values
clim_thresh = clim_dat_full %>%
  group_by(id) %>%
  summarise(clim_thresh_value_9 = quantile(maxtp, 0.9))

clim_date_w_quantile_mod = clim_date_w_quantile_mod %>% left_join(clim_thresh)
clim_date_w_quantile_mod$threshold_9= predict(qunatile_model_fit_9, clim_date_w_quantile_mod)$location

#---- get lambda for climate model
# Calculate lambda
lambda_thresh_ex = clim_date_w_quantile_mod %>%
  group_by(id) %>%
  group_map(~{
    
    thresh_exceedance_9 = clim_date_w_quantile_mod%>%
      filter(id == .x$id[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_9[1], x))
    
    
    tibble(id = .x$id[1],
           year = c(1931, 1942, 2020, 2022),
           thresh_exceedance_9 = 1-thresh_exceedance_9)  
    
    
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()


lambda_thresh_ex %>% 
  write_csv(paste0("data/processed/climate_thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv"))

clim_date_w_quantile_mod %>% 
  left_join(lambda_thresh_ex) %>%
  saveRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv"))
