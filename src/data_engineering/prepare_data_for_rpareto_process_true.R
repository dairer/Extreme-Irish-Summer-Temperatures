rm(list = ls())
library(tidyverse)

setwd("~/Extreme-Irish-Summer-Temperatures/")

prep_rpareto_true = function(marg_mod){
  
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
  
  # get cost threshold
  threshold = quantile(extreme_dates$cost, 0.8) %>% as.numeric
  
  # get extreme dates
  extreme_dates = extreme_dates %>%
    filter(cost > threshold) %>%
    arrange(desc(cost))
  
  # temporally decluster events
  my_data_tmp = extreme_dates %>% arrange(date)
  
  min_dif = my_data_tmp %>% arrange(date) %>% pull(date) %>% diff() %>% min
  
  print(my_data_tmp)
  print("declustering events")
  
  while(min_dif<7){
    i<-1
    while(i < nrow(my_data_tmp)){
      if(diff(c(my_data_tmp[i,]$date, my_data_tmp[i+1,]$date))<7){
        # pick the largest
        if(my_data_tmp[i,]$cost > my_data_tmp[i+1,]$cost){
          my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i+1,]$date,]
        }else{
          my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i,]$date,]
        }
      }
      i<-i+1
    }
    min_dif <- my_data_tmp$date %>% diff() %>% min %>% abs()
  }
  extreme_dates = my_data_tmp
  
  # tibble of extreme events
  obs_data_standardised = obs_data_standardised %>%
    filter(date %in% extreme_dates$date) %>%
    arrange(date)
  
  print("getting data in correct format")
  
  # get observations in list
  exceedances = obs_data_standardised %>%
    group_by(date) %>%
    group_map(~{
      
      # our "conditional" site at the top of the list
      if("aldergrove" %in% .x$Station){
        c(.x %>% filter(Station == "aldergrove") %>% pull(frechet_marg),
          .x %>% filter(Station != "aldergrove") %>% pull(frechet_marg))
      }
    })
  
  # remove observatiosn that dont have valentia_observatory
  to_remove = exceedances %>% map(is.null) %>% unlist %>% which()
  
  exceedances_locs = obs_data_standardised %>%
    group_by(date) %>%
    group_map(~{
      if("aldergrove" %in% .x$Station){
        rbind(.x %>% filter(Station == "aldergrove") %>% dplyr::select(Long.projected, Lat.projected) %>% as.matrix(),
              .x %>% filter(Station != "aldergrove") %>% dplyr::select(Long.projected, Lat.projected ) %>% as.matrix())
      }
    })
  
  if(!is.na(to_remove[1])){
    exceedances = exceedances[-to_remove]
    exceedances_locs = exceedances_locs[-to_remove]
  }
  
  dates_to_rem = extreme_dates %>% arrange(date) %>% .[to_remove,] %>% pull(date)

  extreme_dates = extreme_dates %>% arrange(date)

  list(exceedances = exceedances,
       exceedances_locs = exceedances_locs,
       thresh = threshold,
       extreme_dates = extreme_dates) %>%
    saveRDS(paste0("data/processed/data_for_rpareto/true/data_for_rpareto_", marg_mod))
  
}


job::job({prep_rpareto_true("mod_1")})
job::job({prep_rpareto_true("mod_2")})
job::job({prep_rpareto_true("mod_4")})

