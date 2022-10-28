gc()
rm(list = ls())
library(tidyverse)
library(spatialsample)
library(evgam)
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")

fit_quant_reg = T
calc_lam = T

fit_quant_regression = function(bts_range, marg_mod, num_quantiles){
  
  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))
  quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
  
  # # ------ Bootstrap quant reg models to get CI for each parameter
  hadcrut_var = obs_data %>%
    dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
    unique()
  
  res = c()
  dat = c()
  for(file_name in bts_range){
    if(marg_mod == "mod_1"){
      dat <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name)) %>%
        rename(maxtp = maxtp_1)
    }
    
    if(marg_mod == "mod_2"){
      dat <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name)) %>%
        rename(maxtp = maxtp_2)
    }
    
    if(marg_mod == "mod_4"){
      dat <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name)) %>%
        rename(maxtp = maxtp_4)
    }
    
    dat = dat %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(hadcrut_var)
    
    print(paste0("Bootstrap number ", file_name))
    for(q in seq_along(quantiles_to_estimate_bulk)){
      
      zeta = obs_data$quantile[[1]][q] # quantile to estimate
      
      print(paste0("quantile ", zeta))
      
      # spatial covariate [q^tau_c(s)]
      obs_data$thisval = obs_data$value %>% lapply(`[[`, q) %>% unlist
      
      dat = dat %>%
        dplyr::select(date,Station,maxtp,year,loess_temp_anom,loess_glob_temp_anom) %>%
        left_join(obs_data %>% dplyr::select(Station, value = thisval) %>% unique) %>% drop_na()
      
      
      qunatile_model_fit <- evgam(maxtp ~ value + loess_temp_anom, dat,
                                  family = "ald", ald.args = list(tau = zeta))
      
      # save parameter estimates w CI
      tibble(bts = file_name,
             tau = zeta,
             beta_0 = qunatile_model_fit$location$coefficients[1],
             beta_1 = qunatile_model_fit$location$coefficients[2],
             beta_2 = qunatile_model_fit$location$coefficients[3])%>%
        write_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"), append = T)
    }
  }
}




# 
  # job::job({fit_quant_regression(seq(1,20),  'mod_1', 15)}) # 22 mins
  # job::job({fit_quant_regression(seq(21,40), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(41,60), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(61,80), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(81,100),'mod_1', 15)})
  # job::job({fit_quant_regression(seq(101,120), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(121,140), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(141,160), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(161,180), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(181,200), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(201,220), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(221,240), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(241,260), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(261,280), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(281,300), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(301,320), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(321,480), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(341,360), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(361,380), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(381,400), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(401,420), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(421,440), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(441,460), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(461,480), 'mod_1', 15)})
  # job::job({fit_quant_regression(seq(481,500), 'mod_1', 15)})

  
  
  
  
  # 
  # job::job({fit_quant_regression(seq(1,20),  'mod_2', 15)}) # 22
  # job::job({fit_quant_regression(seq(21,40), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(41,60), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(61,80), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(81,100),'mod_2', 15)})
  # job::job({fit_quant_regression(seq(101,120), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(121,140), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(141,160), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(161,180), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(181,200), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(201,220), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(221,240), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(241,260), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(261,280), 'mod_2', 15)}) 
  # job::job({fit_quant_regression(seq(281,300), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(301,320), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(321,340), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(341,360), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(361,380), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(381,400), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(401,420), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(421,440), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(441,460), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(461,480), 'mod_2', 15)})
  # job::job({fit_quant_regression(seq(481,500), 'mod_2', 15)})
  
  
  ## ---- jobs commented out are ran 
  # job::job({fit_quant_regression(seq(1,20),    'mod_4', 25)}) 
  # job::job({fit_quant_regression(seq(21,40),   'mod_4', 25)})
  # job::job({fit_quant_regression(seq(41,60),   'mod_4', 25)})
  # job::job({fit_quant_regression(seq(61,80),   'mod_4', 25)})
  # job::job({fit_quant_regression(seq(81,100),  'mod_4', 25)})
  # job::job({fit_quant_regression(seq(101,120), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(121,140), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(141,160), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(161,180), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(181,200), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(201,220), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(221,240), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(241,260), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(261,280), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(281,300), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(301,320), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(321,340), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(341,360), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(361,380), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(381,400), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(401,420), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(421,460), 'mod_4', 25)})
  # job::job({fit_quant_regression(seq(461,500), 'mod_4', 25)})





calc_lambda_bts = function(bts_range, marg_mod, num_quantiles){

  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))
  quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
  
  # # # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
    unique() %>%
    arrange(year)
  
  bootstrapped_models = read_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"),
                                 col_names = c("bts", "tau", "beta_0", "beta_1", "beta_2"))
  
  for(bts_num in bts_range){
    print(paste0("Calculating lambda for bts num ", bts_num))
    this_est_of_qnt_mods = bootstrapped_models %>%
      filter(bts == bts_num)
    
    # # --- creates a tibble with each station and its quantule model
    obs_smoothed_quantiles = obs_data %>%
      group_by(Station) %>%
      group_map(~{
        
        # --- get the climate qunatile estimates closest to current station
        clim_vals = obs_data %>%
          filter(Station == .x$Station[1]) %>%
          dplyr::select(quantile, value) %>%
          unique() %>%
          pull(value) %>%
          unlist()
        
        # predict quantile for each year and site
        quant_reg_pars = this_est_of_qnt_mods %>%
          arrange(tau)
        
        res = c()
        for(q in seq_along(quantiles_to_estimate_bulk)){
          
          qpars = quant_reg_pars[q,]
          res = rbind(res,
                      tibble(quantile =  quantiles_to_estimate_bulk[q],
                             year = temporal_covariates$year,
                             loess_glob_temp_anom = temporal_covariates$loess_glob_temp_anom,
                             loess_temp_anom = temporal_covariates$loess_temp_anom,
                             quant_value = qpars$beta_0 + qpars$beta_1*clim_vals[q] + qpars$beta_2*temporal_covariates$loess_temp_anom))
        }
        
        # interpolate quantiles over tau for each year
        res %>%
          group_by(year) %>%
          group_map(~{
            tibble(year = .x$year[1],
                   quant_spline = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
            # tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
            # temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
          }, .keep = T) %>%
          plyr::rbind.fill() %>%
          as_tibble() %>%
          mutate(Station = .x$Station[1])
        
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    
    # save this 
    obs_smoothed_quantiles %>% saveRDS(paste0("output/bulk_model_fits/", marg_mod, "_num_quantiles_",num_quantiles, "_bts_", bts_num))

    # look at threshold exceedance for this model
    obs_data %>%
      group_by(Station) %>%
      group_map(~{
        
        thresh_exceedance = obs_smoothed_quantiles%>%
          filter(Station == .x$Station[1]) %>%
          pull(quant_spline) %>%
          sapply(function(x) sapply(.x$threshold_9[1], x))

        tibble(bts_num = bts_num,
               Station = .x$Station[1],
               year = temporal_covariates$year,
               thresh_exceedance = 1-thresh_exceedance)  %>%
          write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",bts_num,".csv"), append = T)
      }, .keep = T)
  }
}


# marg_mod = 'mod_1'
# num_quantiles = 15
# 
# for(bts_num in seq(1, 500)){
#   print(bts_num)
#   read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",bts_num,".csv"),
#            col_names = F) %>% tail(n = 16744) %>% 
#     write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",bts_num,".csv"))
# }
# 
# 
# marg_mod = 'mod_2'
# num_quantiles = 15
# 
# for(bts_num in seq(1, 500)){
#   print(bts_num)
#   read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",bts_num,".csv"),
#            col_names = F) %>% tail(n = 16744) %>% 
#     write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",bts_num,".csv"))
# }



  



# job::job({calc_lambda_bts(seq(1,100),  'mod_1', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(101,200),  'mod_1', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(201,300),  'mod_1', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(301,400),  'mod_1', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(401,500),  'mod_1', 15)}) # 47 mins

  # job::job({calc_lambda_bts(seq(1,100),  'mod_1', 15)}) # 47 mins
  # job::job({calc_lambda_bts(seq(21,40), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(41,60), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(61,80), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(81,100),'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(101,120), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(121,140), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(141,160), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(161,180), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(181,200), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(201,220), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(221,240), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(241,260), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(261,280), 'mod_1', 15)})
  # job::job({calc_lambda_bts(seq(281,300), 'mod_1', 15)})
  
  
  
#   
# job::job({calc_lambda_bts(seq(1,100),  'mod_2', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(101,200),  'mod_2', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(201,300),  'mod_2', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(301,400),  'mod_2', 15)}) # 47 mins
# job::job({calc_lambda_bts(seq(401,500),  'mod_2', 15)}) # 47 mins
# 
#   # job::job({calc_lambda_bts(seq(1,20),  'mod_2', 15)}) 
#   job::job({calc_lambda_bts(seq(21,40), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(41,60), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(61,80), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(81,100),'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(101,120), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(121,140), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(141,160), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(161,180), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(181,200), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(201,220), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(221,240), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(241,260), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(261,280), 'mod_2', 15)})
#   job::job({calc_lambda_bts(seq(281,300), 'mod_2', 15)})
#   
  
  

job::job({calc_lambda_bts(seq(1,60),    'mod_4', 25)})
job::job({calc_lambda_bts(seq(61, 120),   'mod_4', 25)})
job::job({calc_lambda_bts(seq(121,180),   'mod_4', 25)})
job::job({calc_lambda_bts(seq(181,240),   'mod_4', 25)})
job::job({calc_lambda_bts(seq(241,300),  'mod_4', 25)})
job::job({calc_lambda_bts(seq(301,360),    'mod_4', 25)})
job::job({calc_lambda_bts(seq(361, 420),   'mod_4', 25)})
job::job({calc_lambda_bts(seq(421,480),   'mod_4', 25)})
job::job({calc_lambda_bts(seq(481,500),   'mod_4', 25)})


# job::job({calc_lambda_bts(seq(201,220), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(221,240), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(241,260), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(261,280), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(281,300), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(301,320), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(321,340), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(341,360), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(361,380), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(381,400), 'mod_4', 25)})
# job::job({calc_lambda_bts(seq(401,420), 'mod_4', 25)}) 
# job::job({calc_lambda_bts(seq(421,450), 'mod_4', 25)}) 
# job::job({calc_lambda_bts(seq(451,500), 'mod_4', 25)}) 
job::job({calc_lambda_bts(seq(493,500), 'mod_4', 25)}) 


job::job({calc_lambda_bts(seq(454,480), 'mod_4', 25)}) 
job::job({calc_lambda_bts(seq(492,500), 'mod_4', 25)}) 












# --- on clim grid

get_bulk_bts_on_clim_grid = function(bts_range, marg_mod, num_quantiles){
  for(this_bts in bts_range){
    
    obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))
    quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
    
    temporal_covariates = obs_data %>%
      dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
      unique() %>%
      arrange(year)
    
    
    
    quant_reg_model_pars =  read_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"),
                                     col_names = c("bts", "tau", "beta_0", "beta_1", "beta_2")) %>%
      filter(bts == this_bts)
    
    
    
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
      saveRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
    
    # thresold at these points 
    
    qunatile_model_fit_9 = readRDS("output/threshold_model_9")
    clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
    
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
      write_csv(paste0("data/processed/climate_thresh_exceedance_lambda_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
    
    clim_date_w_quantile_mod %>% 
      left_join(lambda_thresh_ex) %>%
      saveRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",this_bts,".csv"))
  }
}








# job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_4', 25)}) 
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_4', 25)})

# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_4', 25)})

# job::job({get_bulk_bts_on_clim_grid(seq(201,220), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(221,240), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(241,260), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(261,280), 'mod_4', 25)})
#job::job({get_bulk_bts_on_clim_grid(seq(281,300), 'mod_4', 25)})

# job::job({get_bulk_bts_on_clim_grid(seq(301,340),  'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(341,380), 'mod_4', 25)})
# job::job({get_bulk_bts_on_clim_grid(seq(381,420), 'mod_4', 25)})


job::job({get_bulk_bts_on_clim_grid(seq(421,440), 'mod_4', 25)})
job::job({get_bulk_bts_on_clim_grid(seq(441,480), 'mod_4', 25)})
job::job({get_bulk_bts_on_clim_grid(seq(481,520), 'mod_4', 25)})





# job::job({get_bulk_bts_on_clim_grid(seq(1,100),  'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,200), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,300), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(301,400), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(401,500),'mod_1', 15)})


# 
# job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_1', 15)}) 
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_1', 15)})
# 
# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_1', 15)})
# 
# job::job({get_bulk_bts_on_clim_grid(seq(201,220), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(221,240), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(241,260), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(261,280), 'mod_1', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(281,300), 'mod_1', 15)})



# job::job({get_bulk_bts_on_clim_grid(seq(1,100),  'mod_2', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,200), 'mod_2', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,300), 'mod_2', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(301,400), 'mod_2', 15)})
# job::job({get_bulk_bts_on_clim_grid(seq(401,500),'mod_2', 15)})




# ------------------- PLOTS 

marg_mod = 'mod_4'
num_quantiles = 25

true = read_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                col_names = c("tau", "beta_0", "beta_1", "beta_2"))

files = list.files("output/bts_thresh_ex_lambda/")
files = files[grepl(marg_mod, files) & grepl(paste0("num_quantiles_",num_quantiles), files)]

all_bts_dat = c()
for(f in files){
  print(f)
  all_bts_dat = rbind(all_bts_dat,read_csv(paste0("output/bts_thresh_ex_lambda/", f),
                                           col_names = c("bts", "Station", "year", "thresh_exceedance")))
}




all_bts_dat = rbind(all_bts_dat,read_csv(paste0("output/bts_thresh_ex_lambda/", f),
                                         col_names = c("bts", "Station", "year", "thresh_exceedance")))

true_res = read_csv("data/processed/thresh_exceedance_lambda_num_quantiles_15.csv") %>%
  group_by(year) %>%
  summarise(thresh_exceedance = mean(thresh_exceedance_9))

res = read_csv("output/bts_quant_reg_mod_4_num_quantiles_25.csv",
               col_names = c('bts', 'tau', 'b0', 'b1', 'b2'))



# --- read in fitted quantile regression coefficients
true_pars = read_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",25,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2'))


plts=gridExtra::grid.arrange(res %>%
                          group_by(tau) %>%
                          summarise(upper = quantile(b2, 0.975),
                                    lower = quantile(b2, 0.025)) %>%
                          ggplot()+
                          geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3)+
                          geom_line(data = true_pars, aes(tau, beta_2))+
                          xlim(c(0.08, 0.9))+
                          ylim(c(-1.5, 3))+
                          theme_minimal(12)+
                          labs(y = expression(beta[2]),
                               x = expression(tau))+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                        all_bts_dat %>% 
                          group_by(bts, year) %>%
                          summarise(thresh_exceedance = mean(thresh_exceedance))  %>%
                          group_by(year) %>%
                          summarise(upper = quantile(thresh_exceedance, 0.975),
                                    lower = quantile(thresh_exceedance, 0.025)) %>%
                          ungroup() %>%
                          ggplot()+
                          geom_ribbon(aes(year, ymin = lower, ymax = upper), alpha = 0.3) +
                          geom_line(data = true_res, aes(year, thresh_exceedance))+
                          theme_minimal(12)+
                          labs(y = expression(lambda[t]),
                               x = "Year")+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
                          xlim(c(1931, 2020))+
                          ylim(c(0.05, 0.19)), nrow = 1)

ggsave(plts,filename = "output/figs/bulk_quantile_regression.pdf", height = 3, width = 7.5)


# true_lambda$Station%>% unique
# all_sites = true_lambda %>%
#   ggplot()+
#   geom_ribbon(data = bts_lam, aes(x = year, ymin = thresh_exceedance_c_l, ymax = thresh_exceedance_c_u), alpha = 0.25)+
#   geom_line(data = true_lambda %>% group_by(year) %>% summarise(thresh_exceedance = mean(thresh_exceedance)), aes(year, thresh_exceedance))+
#   geom_line(data = true_lambda %>% filter(Station == 'claremorris'), aes(year, thresh_exceedance, col = 'claremorris', linetype = 'claremorris'))+
#   geom_line(data = true_lambda %>% filter(Station == 'mullingar'), aes(year, thresh_exceedance, col = 'mullingar', linetype = 'mullingar'))+
#   geom_line(data = true_lambda %>% filter(Station == 'phoenixpark'), aes(year, thresh_exceedance, col = 'phoenixpark', linetype = 'phoenixpark'))+
#   geom_line(data = true_lambda %>% filter(Station == 'malin_head'), aes(year, thresh_exceedance, col = 'malin_head', linetype = 'malin_head'))+
#   geom_line(data = true_lambda %>% filter(Station == 'roches_point'), aes(year, thresh_exceedance, col = 'roches_point', linetype = 'roches_point'))+
#   theme_minimal()+
#   labs(colour = expression(s),
#        linetype =  expression(s),
#        y = expression(lambda),
#        x = "Year")+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
# 
# 
# ggsave(all_sites, filename = "lambda_all_sites.pdf", width =5, height = 2.75)
# 
# 
# read_csv("corrections/marginal_models/bulk/bootstrapped_lam.csv",
#          col_names = c("bts_num", "Station", "year", "thresh_exceedance_c")) %>%
#   group_by(bts_num, Station, year) %>%
#   summarise(thresh_exceedance_c = mean(thresh_exceedance_c)) %>%
#   ungroup() %>%
#   group_by(year, Station)%>%
#   summarise(thresh_exceedance_c_u = quantile(thresh_exceedance_c, 0.975),
#             thresh_exceedance_c_l = quantile(thresh_exceedance_c, 0.025),
#             thresh_exceedance_c = mean(thresh_exceedance_c)) %>%
#   ggplot()+
#   geom_ribbon(aes(x = year, ymin = thresh_exceedance_c_l, ymax = thresh_exceedance_c_u), alpha = 0.25)+
#   geom_line(aes(year, thresh_exceedance_c))+
#   theme_minimal()+
#   labs(colour = "Model",
#        fill = "Model",
#        y = expression(lambda),
#        x = "Year")+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
# 
# 

# 
# read_csv("corrections/marginal_models/bulk/bootstrapped_lam.csv",
#          col_names = c("bts_num", "Station", "year", "thresh_exceedance_c")) %>%
#   group_by(bts_num, year)%>% 
#   summarise(thresh_exceedance_c = mean(thresh_exceedance_c)) %>%
#   ungroup() %>%
#   group_by(year)%>%
#   summarise(thresh_exceedance_c_u = quantile(thresh_exceedance_c, 0.975),
#             thresh_exceedance_c_l = quantile(thresh_exceedance_c, 0.025),
#             thresh_exceedance_c = mean(thresh_exceedance_c)) %>%
#   ggplot()+
#   geom_ribbon(aes(x = year, ymin = thresh_exceedance_c_l, ymax = thresh_exceedance_c_u), alpha = 0.25)+
#   geom_line(aes(year, thresh_exceedance_c))+
#   theme_minimal()+
#   labs(colour = "Model",
#        fill = "Model",
#        y = expression(lambda),
#        x = "Year")+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
# 
# 
#   group_by(year) %>%
#   #filter(Station == "claremorris") %>%
#   summarise(thresh_exceedance_c_u = quantile(thresh_exceedance_c, 0.975),
#             thresh_exceedance_c_l = quantile(thresh_exceedance_c, 0.025),
#             thresh_exceedance_c = mean(thresh_exceedance_c)
#             # thresh_exceedance_a_u = quantile(thresh_exceedance_a, 0.975),
#             # thresh_exceedance_a_l = quantile(thresh_exceedance_a, 0.025),
#             # thresh_exceedance_a = mean(thresh_exceedance_a),
#             # thresh_exceedance_b_u = quantile(thresh_exceedance_b, 0.975),
#             # thresh_exceedance_b_l = quantile(thresh_exceedance_b, 0.025),
#             # thresh_exceedance_b = mean(thresh_exceedance_b)
#             ) %>%
#   ggplot()+
#   geom_ribbon(aes(x = year, ymin = thresh_exceedance_c_l, ymax = thresh_exceedance_c_u), alpha = 0.25)+
#   geom_line(aes(year, thresh_exceedance_c))+
#   # geom_ribbon(aes(x = year, ymin = thresh_exceedance_b_l, ymax = thresh_exceedance_b_u, fill = "M^G"), alpha = 0.25)+
#   # geom_line(aes(year, thresh_exceedance_b, col = "M^G"))+
#   # geom_ribbon(aes(x = year, ymin = thresh_exceedance_a_l, ymax = thresh_exceedance_a_u, fill = "year"), alpha = 0.25)+
#   # geom_point(aes(year, thresh_exceedance_a, col = "year"))+
#   theme_minimal()+
#   labs(colour = "Model",
#        fill = "Model",
#        y = expression(lambda),
#        x = "Year")+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
# # # 
# # # 
# # # 
# # # 
# # 
# bts_quant_reg = bts_quant_reg %>%
#   group_by(tau) %>%
#   summarise(beta_a_2_u = quantile(beta_a_2, 0.975),
#             beta_a_2_l = quantile(beta_a_2, 0.025)) 
# 
# 
# # 
bts_quant_reg = read_csv("corrections/marginal_models/bulk/bts_quant_reg.csv",
                         col_names = c("bts",
                                       "tau",
                                       "beta_a_0",
                                       "beta_a_1",
                                       "beta_a_2"))

# 
# bts_quant_reg %>% 
#   ggplot()+
#   geom_ribbon(aes(x = tau, ymin = beta_a_2_l, ymax = beta_a_2_u), alpha = 0.25)+
#   geom_line(data = quant_reg_pars, aes(tau,beta_a_2))+
#   xlim(0.1, 0.9)
# 
# beta_2_plot=bts_quant_reg %>%
#   filter(bts != 96) %>%
#   group_by(tau) %>%
#   summarise(beta_a_2_u = quantile(beta_a_2, 0.975),
#             beta_a_2_l = quantile(beta_a_2, 0.025)) %>%
#   ggplot()+
#   #geom_segment(aes(x = tau, xend = tau,  y = beta_a_2_l, yend = beta_a_2_u), alpha = 0.25)+
#   geom_ribbon(aes(x = tau, ymin = beta_a_2_l, ymax = beta_a_2_u), alpha = 0.25)+
#   geom_line(data = quant_reg_pars %>% filter(tau<0.9, tau> 0.05), aes(tau,beta_a_2))+
#   ylim(-2, 2.7)+
#   theme_minimal(12)+
#   labs(y = expression(beta[2]),
#        x = expression(tau))+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
#   scale_x_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 0.95))
# 
# 
# lambda_plot = true_lambda %>%
#   ggplot()+
#   geom_ribbon(data = bts_lam, aes(x = year, ymin = thresh_exceedance_c_l, ymax = thresh_exceedance_c_u), alpha = 0.25)+
#   geom_line(data = true_lambda %>% group_by(year) %>% summarise(thresh_exceedance = mean(thresh_exceedance)), aes(year, thresh_exceedance))+
#   theme_minimal(12)+
#   labs(colour = expression(s),
#        linetype =  expression(s),
#        y = expression(lambda),
#        x = "Year")+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
# 
# 
# 
# plot_fig_3 = gridExtra::grid.arrange(beta_2_plot, lambda_plot, nrow =1)
# # 
# # 
# ggsave(plot_fig_3,filename =  'corrections/FIG3.pdf', width = 6, height = 2.5)
# # # 
# 
# tst = read_csv("corrections/marginal_models/bulk/bts_quant_reg.csv",
#          col_names = c(
#            "bts",
#            "tau",
#            "beta_a_0",
#            "beta_a_1",
#            "beta_a_2",
#            "beta_b_0",
#            "beta_b_1",
#            "beta_b_2",
#            "beta_c_0",
#            "beta_c_1",
#            "beta_c_2"))
# 
# 
# # # # ----- read in bootstraps and get CI
# qunat_reg_ci = read_csv("corrections/marginal_models/bulk/bts_quant_reg.csv",
#          col_names = c(
#            "bts",
#            "tau",
#           "beta_a_0",
#           "beta_a_1",
#           "beta_a_2",
#           "beta_b_0",
#           "beta_b_1",
#           "beta_b_2",
#           "beta_c_0",
#           "beta_c_1",
#           "beta_c_2")) %>%
#   group_by(tau) %>%
#   summarise(beta_a_0_l = quantile(beta_a_0, 0.025),
#             beta_a_0_u = quantile(beta_a_0, 0.975),
#             beta_a_1_l = quantile(beta_a_1, 0.025),
#             beta_a_1_u = quantile(beta_a_1, 0.975),
#             beta_a_2_l = quantile(beta_a_2, 0.025),
#             beta_a_2_u = quantile(beta_a_2, 0.975),
#             beta_b_0_l = quantile(beta_b_0, 0.025),
#             beta_b_0_u = quantile(beta_b_0, 0.975),
#             beta_b_1_l = quantile(beta_b_1, 0.025),
#             beta_b_1_u = quantile(beta_b_1, 0.975),
#             beta_b_2_l = quantile(beta_b_2, 0.025),
#             beta_b_2_u = quantile(beta_b_2, 0.975),
#             beta_c_0_l = quantile(beta_c_0, 0.025),
#             beta_c_0_u = quantile(beta_c_0, 0.975),
#             beta_c_1_l = quantile(beta_c_1, 0.025),
#             beta_c_1_u = quantile(beta_c_1, 0.975),
#             beta_c_2_l = quantile(beta_c_2, 0.025),
#             beta_c_2_u = quantile(beta_c_2, 0.975))
# 
# 
# 
# 
# # create bootstraps within years?????
# 
# 
# quant_reg_pars = quant_reg_pars %>%
#   left_join(qunat_reg_ci)
# 
# 
# 
# gridExtra::grid.arrange(quant_reg_pars %>%
#                           ggplot()+
#                           geom_ribbon(aes(x = tau, ymin = beta_c_0_l, ymax = beta_c_0_u, fill = "beta_b_2"), alpha = 0.25)+
#                           geom_line(aes(tau, beta_c_0, col = "beta_0"))+
#                           theme_minimal()+
#                           theme(legend.position = 'none'),
#                         
#                         
#                         
#                         quant_reg_pars %>%
#                           ggplot()+
#                           geom_ribbon(aes(x = tau, ymin = beta_c_1_l, ymax = beta_c_1_u, fill = "beta_b_2"), alpha = 0.25)+
#                           geom_line(aes(tau, beta_c_1, col = "beta_1"))+
#                           theme_minimal()+
#                           theme(legend.position = 'none'),
#                         
#               
#                         quant_reg_pars %>%
#                           ggplot()+
#                           geom_ribbon(aes(x = tau, ymin = beta_c_2_l, ymax = beta_c_2_u, fill = "beta_b_2"), alpha = 0.25)+
#                           geom_line(aes(tau, beta_c_2, col = "beta_2"))+
#                           theme_minimal()+
#                           theme(legend.position = 'none'), nrow = 1
# )
# 









# #   ggplot()+
# #   geom_ribbon(aes(x = tau, ymin = beta_c_2_l, ymax = beta_c_2_u, fill = "beta_b_2"), alpha = 0.25)+
# #   geom_line(aes(tau, beta_c_2, col = "beta_b_2"))+
# #   # geom_ribbon(aes(x = tau, ymin = beta_c_2_l, ymax = beta_c_2_u, fill = "beta_c_2"), alpha = 0.25)+
# #   # geom_line(aes(tau, beta_c_2, col = "beta_c_2"))+
# #   # geom_ribbon(aes(x = tau, ymin = beta_c_0_l, ymax = beta_c_0_u, fill = "beta_c_2"), alpha = 0.25)+
# #   # geom_line(aes(tau, beta_c_0, col = "beta_c_2"))+
# #   # geom_ribbon(aes(x = tau, ymin = beta_a_1_l, ymax = beta_a_1_u, fill = "beta_c_2"), alpha = 0.25)+
# #   # geom_line(aes(tau, beta_a_1, col = "beta_c_2"))+
# #   theme_minimal()
# 



# 
# 
# # # --- creates a tibble with each station and its quantule model
# obs_smoothed_quantiles = obs_data %>%
#   group_by(Station) %>%
#   group_map(~{
#     
# 
#     # --- get the climate qunatile estimates closest to current station
#     clim_vals = obs_data %>%
#       filter(Station == .x$Station[1]) %>%
#       dplyr::select(quantile, value) %>%
#       unique() %>%
#       pull(value) %>%
#       unlist()
# 
#     # predict quantile for each year and site
#     quant_reg_pars = quant_reg_pars %>%
#       arrange(tau)
#     
#     res = c()
#     for(q in seq_along(quantiles_to_estimate_bulk)){
#       
#       qpars = quant_reg_pars[q,]
#       res = rbind(res,
#                   tibble(quantile =  quantiles_to_estimate_bulk[q],
#                          year = temporal_covariates$year,
#                          loess_glob_temp_anom = temporal_covariates$loess_glob_temp_anom,loess_temp_anom = temporal_covariates$loess_temp_anom,
#                          quant_value_a = qpars$beta_a_0 + qpars$beta_a_1*clim_vals[q] + qpars$beta_a_2*temporal_covariates$year,
#                          quant_value_b = qpars$beta_b_0 + qpars$beta_b_1*clim_vals[q] + qpars$beta_b_2*temporal_covariates$loess_glob_temp_anom,
#                          quant_value_c = qpars$beta_c_0 + qpars$beta_c_1*clim_vals[q] + qpars$beta_c_2*temporal_covariates$loess_temp_anom,
#                          quant_value_a_u = qpars$beta_a_0_u + qpars$beta_a_1_u*clim_vals[q] + qpars$beta_a_2_u*temporal_covariates$year,
#                          quant_value_b_u = qpars$beta_b_0_u + qpars$beta_b_1_u*clim_vals[q] + qpars$beta_b_2_u*temporal_covariates$loess_glob_temp_anom,
#                          quant_value_c_u = qpars$beta_c_0_u + qpars$beta_c_1_u*clim_vals[q] + qpars$beta_c_2_u*temporal_covariates$loess_temp_anom,
#                          quant_value_a_l = qpars$beta_a_0_l + qpars$beta_a_1_l*clim_vals[q] + qpars$beta_a_2_l*temporal_covariates$year,
#                          quant_value_b_l = qpars$beta_b_0_l + qpars$beta_b_1_l*clim_vals[q] + qpars$beta_b_2_l*temporal_covariates$loess_glob_temp_anom,
#                          quant_value_c_l = qpars$beta_c_0_l + qpars$beta_c_1_l*clim_vals[q] + qpars$beta_c_2_l*temporal_covariates$loess_temp_anom))
#       
#       # res = rbind(res,tibble(year = temporal_covariates$year,
#       #                        loess_glob_temp_anom = temporal_covariates$loess_glob_temp_anom,
#       #                        loess_temp_anom = temporal_covariates$loess_temp_anom)
#       #                        quantile =  quantiles_to_estimate_bulk[q],
#       #                        value_a = obs_q_mod[q,]$quant_mod_a %>% lapply(predict, tibble(value = clim_vals[q], year = seq(1942, 2020))) %>% .[[1]] %>% .$location,
#       #                        #value_a_se = obs_q_mod[q,]$quant_mod %>% lapply(predict, tibble(value = clim_vals[q], year = seq(1942, 2020)), type = 'response', se.fit = T) %>% .[[1]] %>% .$se.fit %>% .$location,
#       #                        value_b = obs_q_mod[q,]$quant_mod_b %>% lapply(predict, tibble(value = clim_vals[q], loess_glob_temp_anom = loess_temp_glob_anom_vals)) %>% .[[1]] %>% .$location,
#       #                        #value_b_se = obs_q_mod[q,]$quant_mod %>% lapply(predict, tibble(value = clim_vals[q], loess_glob_temp_anom = loess_temp_glob_anom_vals), type = 'response', se.fit = T) %>% .[[1]] %>% .$se.fit %>% .$location,
#       #                        value_c = obs_q_mod[q,]$quant_mod_c %>% lapply(predict, tibble(value = clim_vals[q], loess_temp_anom = loess_temp_anom_vals)) %>% .[[1]] %>% .$location
#       #                        #value_c_se = obs_q_mod[q,]$quant_mod %>% lapply(predict, tibble(value = 




#       #))
#     }
# 
#     print(paste0("Interpolating quantile estimates for ", .x$Station[1]))
# 
#     # interpolate quantiles over tau for each year
#     res %>%
#       group_by(year) %>%
#       group_map(~{
#         tibble(year = .x$year[1],
#                quant_spline_a = list(splinefun(.x$quant_value_a,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_b = list(splinefun(.x$quant_value_b,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_c = list(splinefun(.x$quant_value_c,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_a_u = list(splinefun(.x$quant_value_a_u,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_a_l = list(splinefun(.x$quant_value_a_l,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_b_u = list(splinefun(.x$quant_value_b_u,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_b_l = list(splinefun(.x$quant_value_b_l,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_c_u = list(splinefun(.x$quant_value_c_u,.x$quantile,  method = 'monoH.FC')),
#                quant_spline_c_l = list(splinefun(.x$quant_value_c_l,.x$quantile,  method = 'monoH.FC'))
#         )
#       }, .keep = T) %>%
#       plyr::rbind.fill() %>%
#       as_tibble() %>%
#       mutate(Station = .x$Station[1])
# 
#   }, .keep = T) %>%
#   plyr::rbind.fill() %>%
#   as_tibble()
# 
# 
# 
# 
# obs_smoothed_quantiles %>%
#   filter(Station %in% c('claremorris')) 
# 
# 
# 
# obs_sites_quantile_estimtes %>%
#   filter(Station %in% c('claremorris', "mullingar", "malin_head", "armagh")) %>%
#   ggplot()+
#   geom_line(aes(quantile, quantile_value, col = year, group = year))+
#   theme_minimal()+
#   facet_wrap(~Station, nrow = 1)+
#   ggplot2::scale_colour_gradientn(colours = my_pal)
# 
# 
# 
# 
# # fc_0 = obs_smoothed_quantiles%>%
# #   filter(Station == "claremorris") %>%
# #   filter(year == 1970) %>%
# #   pull(quant_spline_c) %>% .[[1]]
# 
# 
# 
# fc_1 = obs_smoothed_quantiles%>%
#   filter(Station == "claremorris") %>%
#   filter(year == 2010) %>%
#   pull(quant_spline_c) %>% .[[1]]
# 
# 
# 
# fc_u = obs_smoothed_quantiles%>%
#   filter(Station == "claremorris") %>%
#   filter(year == 2010) %>%
#   pull(quant_spline_c_u) %>% .[[1]]
# 
# 
# fc_l = obs_smoothed_quantiles%>%
#   filter(Station == "claremorris") %>%
#   filter(year == 2010) %>%
#   pull(quant_spline_c_l) %>% .[[1]]
# 
# tibble(temp = seq(14, 25), 
#        tau0 = fc_1(seq(14, 25)),
#        tau_u = fc_u(seq(14, 25)),
#        tau_l = fc_l(seq(14, 25))) %>%
#   ggplot()+
#   geom_line(aes(tau0, temp , col = "1970"))+
#   geom_line(aes(tau_u,temp,  col = "u"))+
#   geom_line(aes(tau_l,temp,  col = "l"))
# 
# 
# 
# 
# 
# 
# # look at threshold exceedance for this model
# thresh_exceedance_prob = obs_data %>%
#   group_by(Station) %>%
#   group_map(~{
# 
#     thresh_exceedance_a = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_a) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#     
#     thresh_exceedance_a_u = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_a_u) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#     
#     thresh_exceedance_a_l = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_a_l) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
# 
#     thresh_exceedance_b = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_b) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#     
#     thresh_exceedance_b_u = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_b_u) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#     
#     thresh_exceedance_b_l = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_b_l) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#   
#     thresh_exceedance_c = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_c) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#     
#     thresh_exceedance_c_u = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_c_u) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#     
#     thresh_exceedance_c_l = obs_smoothed_quantiles%>%
#       filter(Station == .x$Station[1]) %>%
#       pull(quant_spline_c_l) %>%
#       sapply(function(x) sapply(.x$threshold[1], x))
#     
#     tibble(Station = .x$Station[1],year = seq(1942, 2020),
#            thresh_exceedance_a = 1-thresh_exceedance_a,
#            thresh_exceedance_b = 1-thresh_exceedance_b,
#            thresh_exceedance_c = 1-thresh_exceedance_c,
#            thresh_exceedance_a_u = 1-thresh_exceedance_a_u,
#            thresh_exceedance_b_u = 1-thresh_exceedance_b_u,
#            thresh_exceedance_c_u = 1-thresh_exceedance_c_u,
#            thresh_exceedance_a_l = 1-thresh_exceedance_a_l,
#            thresh_exceedance_b_l = 1-thresh_exceedance_b_l,
#            thresh_exceedance_c_l = 1-thresh_exceedance_c_l)
#   }, .keep = T) %>%
#   plyr::rbind.fill() %>%
#   as_tibble()
# 
# thresh_exceedance_prob %>%
#   group_by(year) %>%
#   summarise(thresh_exceedance_a = mean(thresh_exceedance_a),
#             thresh_exceedance_b = mean(thresh_exceedance_b),
#             thresh_exceedance_c = mean(thresh_exceedance_c),
#             thresh_exceedance_a_u = mean(thresh_exceedance_a_u),
#             thresh_exceedance_a_l = mean(thresh_exceedance_a_l),
#             thresh_exceedance_c_u = mean(thresh_exceedance_c_u),
#             thresh_exceedance_c_l = mean(thresh_exceedance_c_l)) %>%
#   ggplot()+
#   # geom_line(aes(year, thresh_exceedance_a, col = "year(t)")) +
#   # geom_line(aes(year, thresh_exceedance_b, col = "M^G")) +
#   #geom_line(aes(year, thresh_exceedance_a, col = "M^I")) +
#   geom_line(aes(year, thresh_exceedance_a, col = "M^I")) +
#   # geom_line(aes(year, thresh_exceedance_a_u, col = "M^I")) +
#   # geom_line(aes(year, thresh_exceedance_a_l, col = "M^I")) +
#   theme_minimal(12)+
#   labs(x = "Year",
#        y = expression(lambda[t]))
# # 
# # 
# # 
# # 
# # 
# # quant_mod_terms = gridExtra::grid.arrange(tibble(quantiles_to_estimate_bulk, beta_1, beta_2, beta_2_se, beta_1_se)%>%
# #                                             mutate(beta_2_u = beta_2 + 1.96*beta_2_se,
# #                                                    beta_2_l = beta_2 - 1.96*beta_2_se,
# #                                                    beta_1_u = beta_1 + 1.96*beta_1_se,
# #                                                    beta_1_l = beta_1 - 1.96*beta_1_se) %>%
# #                                             ggplot()+
# #                                             geom_line(aes(quantiles_to_estimate_bulk, beta_1))+
# #                                             geom_line(aes(quantiles_to_estimate_bulk, beta_1_u))+
# #                                             geom_line(aes(quantiles_to_estimate_bulk, beta_1_l))+
# #                                             theme_minimal(12)+
# #                                             labs(y = expression(beta[1]),
# #                                                  x = expression(tau)),
# #                                           tibble(quantiles_to_estimate_bulk, beta_1, beta_2, beta_2_se, beta_1_se)%>%
# #                                             mutate(beta_2_u = beta_2 + 1.96*beta_2_se,
# #                                                    beta_2_l = beta_2 - 1.96*beta_2_se,
# #                                                    beta_1_u = beta_1 + 1.96*beta_1_se,
# #                                                    beta_1_l = beta_1 - 1.96*beta_1_se) %>%
# #                                             ggplot()+
# #                                             geom_line(aes(quantiles_to_estimate_bulk, beta_2))+
# #                                             geom_line(aes(quantiles_to_estimate_bulk, beta_2_u))+
# #                                             geom_line(aes(quantiles_to_estimate_bulk, beta_2_l))+
# #                                             theme_minimal(12)+
# #                                             labs(y = expression(beta[2]),
# #                                                  x = expression(tau)),
# #                                           thresh_exceedance_prob %>%
# #                                             group_by(year) %>%
# #                                             summarise(thresh_exceedance = mean(thresh_exceedance),
# #                                                       thresh_exceedance_l = mean(thresh_exceedance_l),
# #                                                       thresh_exceedance_u = mean(thresh_exceedance_u)) %>%
# #                                             ggplot()+
# #                                             geom_line(aes(year, thresh_exceedance)) +
# #                                             geom_line(aes(year, thresh_exceedance_u)) +
# #                                             geom_line(aes(year, thresh_exceedance_l)) +
# #                                             theme_minimal(12)+
# #                                             labs(x = "Year",
# #                                                  y = expression(lambda[t])),nrow = 1)
# # 
# # ggsave(quant_mod_terms, filename = "corrections/quant_mod_terms.pdf", height = 3, width = 10)


# 
# 
# 5/107
# 29/80
# 36/71
# 
# 
# read_csv("corrections/data/processed_data/data_for_r_pareto/bootstrapped/non_stationary_rpareto_fit.csv", 
#          col_names = c("bts",'model',
#                        "alpha",
#                        'beta',
#                        "t")) %>%
#   filter(t<100) %>%
#   filter(model == 'model.4b') %>%
#   filter(t>=0)
# pull(t) %>%
#   quantile(c(0.025, 0.975))
# ggplot()+
#   geom_histogram(aes(t))
