gc()
rm(list = ls())
library(tidyverse)
library(spatialsample)
library(evgam)
setwd("~/Extreme-Irish-Summer-Temperatures/")

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
    if(marg_mod == "mod_0"){
      dat <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name)) %>%
        rename(maxtp = maxtp_0)
    }
    
    if(marg_mod == "mod_1"){
      dat <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name)) %>%
        rename(maxtp = maxtp_1)
    }
    
    if(marg_mod == "mod_2"){
      dat <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name)) %>%
        rename(maxtp = maxtp_2)
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
             beta_0 = qunatile_model_fit$location$coefficients[2],
             beta_1 = qunatile_model_fit$location$coefficients[3])%>%
        write_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"), append = T)
    }
  }
}





# job::job({fit_quant_regression(seq(1,20),  'mod_0', 30)}) 
# job::job({fit_quant_regression(seq(21,40), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(41,60), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(61,80), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(81,100),'mod_0', 30)})
# job::job({fit_quant_regression(seq(101,120), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(121,140), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(141,160), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(161,180), 'mod_0', 30)})
# job::job({fit_quant_regression(seq(181,200), 'mod_0', 30)})

# job::job({fit_quant_regression(seq(1,20),  'mod_1', 30)}) 
# job::job({fit_quant_regression(seq(21,40), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(41,60), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(61,80), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(81,100),'mod_1', 30)})
# job::job({fit_quant_regression(seq(101,120), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(121,140), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(141,160), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(161,180), 'mod_1', 30)})
# job::job({fit_quant_regression(seq(181,200), 'mod_1', 30)})
  
# job::job({fit_quant_regression(seq(1,20),    'mod_2', 30)}) 
# job::job({fit_quant_regression(seq(21,40),   'mod_2', 30)})
# job::job({fit_quant_regression(seq(41,60),   'mod_2', 30)})
# job::job({fit_quant_regression(seq(61,80),   'mod_2', 30)})
# job::job({fit_quant_regression(seq(81,100),  'mod_2', 30)})
# job::job({fit_quant_regression(seq(101,120), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(121,140), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(141,160), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(161,180), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(181,200), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(201,220), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(221,240), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(241,260), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(261,280), 'mod_2', 30)})
# job::job({fit_quant_regression(seq(281,300), 'mod_2', 30)})



calc_lambda_bts = function(bts_range, marg_mod, num_quantiles){

  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))
  quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)
  
  # # # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year, loess_temp_anom, loess_glob_temp_anom) %>%
    unique() %>%
    arrange(year)
  
  bootstrapped_models = read_csv(paste0("output/bts_quant_reg_", marg_mod, "_num_quantiles_", num_quantiles, ".csv"),
                                 col_names = c("bts", "tau", "beta_0", "beta_0", "beta_1"))
  
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
                             quant_value = qpars$beta_0 + qpars$beta_0*clim_vals[q] + qpars$beta_1*temporal_covariates$loess_temp_anom))
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


# job::job({calc_lambda_bts(seq(1,100),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(101,200),  'mod_0', 30)})
# job::job({calc_lambda_bts(seq(201,300),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(301,400),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(401,500),  'mod_0', 30)})
# job::job({calc_lambda_bts(seq(1,100),  'mod_0', 30)}) 
# job::job({calc_lambda_bts(seq(21,40), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(41,60), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(61,80), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(81,100),'mod_0', 30)})
# job::job({calc_lambda_bts(seq(101,120), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(121,140), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(141,160), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(161,180), 'mod_0', 30)})
# job::job({calc_lambda_bts(seq(181,200), 'mod_0', 30)})

# job::job({calc_lambda_bts(seq(1,100),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(101,200),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(201,300),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(301,400),  'mod_1', 30)})
# job::job({calc_lambda_bts(seq(401,500),  'mod_1', 30)}) 
# job::job({calc_lambda_bts(seq(1,20),  'mod_1', 30)})
# job::job({calc_lambda_bts(seq(21,40), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(41,60), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(61,80), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(81,100),'mod_1', 30)})
# job::job({calc_lambda_bts(seq(101,120), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(121,140), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(141,160), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(161,180), 'mod_1', 30)})
# job::job({calc_lambda_bts(seq(181,200), 'mod_1', 30)})

# job::job({calc_lambda_bts(seq(1,60),    'mod_2', 30)})
# job::job({calc_lambda_bts(seq(61, 120),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(121,180),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(181,240),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(241,300),  'mod_2', 30)})
# job::job({calc_lambda_bts(seq(301,360),    'mod_2', 30)})
# job::job({calc_lambda_bts(seq(361, 420),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(421,480),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(481,500),   'mod_2', 30)})
# job::job({calc_lambda_bts(seq(201,220), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(221,240), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(241,260), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(261,280), 'mod_2', 30)})
# job::job({calc_lambda_bts(seq(281,300), 'mod_2', 30)})


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
                                     col_names = c("bts", "tau", "beta_0", "beta_0", "beta_1")) %>%
      filter(bts == this_bts)
    
    
    
    clim_grid = read_csv("data/processed/clim_scale_grid.csv")
    clim_dat_full = read_csv("~/Extreme-Irish-Summer-Temperatures/data/processed_data/full_clim_data.csv")
    
    
    
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
                           quant_value = qpars$beta_0 + qpars$beta_0*clim_vals[q] + (qpars$beta_1)*(this_data_for_pred$loess_temp_anom)))
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








# job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_2', 30)}) 
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,220), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(221,240), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(241,260), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(261,280), 'mod_2', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(281,300), 'mod_2', 30)})


# job::job({get_bulk_bts_on_clim_grid(seq(1,100),  'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,200), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,300), 'mod_0', 30)}
# job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_0', 30)}) 
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_0', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_0', 30)})


# job::job({get_bulk_bts_on_clim_grid(seq(1,100),  'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,200), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(201,300), 'mod_1', 30)}
# job::job({get_bulk_bts_on_clim_grid(seq(1,20),  'mod_1', 30)}) 
# job::job({get_bulk_bts_on_clim_grid(seq(21,40), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(41,60), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(61,80), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(81,100),'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(101,120), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(121,140), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(141,160), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(161,180), 'mod_1', 30)})
# job::job({get_bulk_bts_on_clim_grid(seq(181,200), 'mod_1', 30)})



# ------------------- PLOTS 
marg_mod = 'mod_2'
num_quantiles = 30

true = read_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",num_quantiles,".csv"),
                col_names = c("tau", "beta_0", "beta_0", "beta_1"))

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

true_res = read_csv("data/processed/thresh_exceedance_lambda_num_quantiles_05.csv") %>%
  group_by(year) %>%
  summarise(thresh_exceedance = mean(thresh_exceedance_9))

res = read_csv("output/bts_quant_reg_mod_2_num_quantiles_30.csv",
               col_names = c('bts', 'tau', 'b0', 'b1', 'b2'))

# --- read in fitted quantile regression coefficients
true_pars = read_csv(paste0("data/processed/qunatile_model_fit_pars_num_quantiles_",30,".csv"),
                                col_names = c('tau', 'beta_0', 'beta_0', 'beta_1'))


plts=gridExtra::grid.arrange(res %>%
                          group_by(tau) %>%
                          summarise(upper = quantile(b2, 0.975),
                                    lower = quantile(b2, 0.030)) %>%
                          ggplot()+
                          geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3)+
                          geom_line(data = true_pars, aes(tau, beta_1))+
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
                                    lower = quantile(thresh_exceedance, 0.030)) %>%
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