# Description: get bootstrapped data sets under each marginal model
# --- use true gpd and bts models to standardise, sample and back transform 
# --- to create bootstrap sampels

gc() 
rm(list = ls(all.names = TRUE))
library(tidyverse)
setwd("~/Extreme-Irish-Summer-Temperatures/")
source('models/marginal_models/gpd_models.R')

standardise_data = F
if(standardise_data){
  num_quantiles = 30
  obs_data = read_csv("data/processed/obs_data.csv") %>%
    left_join(read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv")))  %>%
    left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea), by = 'Station') %>%
    dplyr::select(-c(loess_glob_temp_anom,id ))
  
  quantile_models = readRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))

  extreme_data = obs_data %>%
    group_by(Station) %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)

  this_fit_mod_0 = fit_mod_0(extreme_data$excess, extreme_data$scale_9)
  this_fit_mod_1 = fit_mod_1(extreme_data$excess, extreme_data$scale_9, extreme_data$loess_temp_anom)
  this_fit_mod_2 = fit_mod_2(extreme_data$excess, extreme_data$scale_9, extreme_data$loess_temp_anom, extreme_data$dist_sea)
  rm(extreme_data)
  
  pred_0 = my_predict_0(this_fit_mod_0, obs_data$scale_9)
  pred_1 = my_predict_1(this_fit_mod_1, obs_data$scale_9, obs_data$loess_temp_anom)
  pred_2 =my_predict_2(this_fit_mod_2, obs_data$scale_9, obs_data$loess_temp_anom, obs_data$dist_sea)
  
  obs_data = obs_data %>%
    mutate(scale_0 = pred_0$scale,
           scale_1 = pred_1$scale,
           scale_2 = pred_2$scale,
           shape_0 = pred_0$shape,
           shape_1 = pred_1$shape,
           shape_2 = pred_2$shape)
  
  rm('this_fit_mod_0', 'this_fit_mod_1',  'this_fit_mod_2',  'pred_0', 'pred_1', 'pred_2')

  print("Joining bulk models")
  obs_data = obs_data %>%
    left_join(quantile_models)
  
  obs_data = obs_data %>%
    group_by(Station) %>%
    group_map(~{
      print(.x$Station[1])
      data = .x$maxtp
      threshold = .x$threshold_9

      scle_0 = .x$scale_0
      shpe_0 = .x$shape_0
      
      scle_1 = .x$scale_1
      shpe_1 = .x$shape_1
      
      scle_2 = .x$scale_2
      shpe_2 = .x$shape_2
      
      yr = .x$year
      my_lambda = .x$thresh_exceedance_9
      
      res_0 = rep(NA, length(data))
      res_1 = rep(NA, length(data))
      res_2 = rep(NA, length(data))

      for(i in seq(length(data))){
        if(data[i] > threshold[i]){ # transform tail
          
          res_0[i] = 1 - (my_lambda[i]) *(1+shpe_0[i]* ((data[i] - threshold[i])/scle_0[i]))^(-1/(shpe_0[i]))
          res_1[i] = 1 - (my_lambda[i]) *(1+shpe_1[i]* ((data[i] - threshold[i])/scle_1[i]))^(-1/(shpe_1[i]))
          res_2[i] = 1 - (my_lambda[i]) *(1+shpe_2[i]* ((data[i] - threshold[i])/scle_2[i]))^(-1/(shpe_2[i]))
          
        }else{ 
          res_0[i] = .x$temp_to_tau[i][[1]](data[i])
          res_1[i] = res_0[i]
          res_2[i] = res_0[i]
        }
      }
      .x$unif_0 = res_0
      .x$unif_1 = res_1
      .x$unif_2 = res_2
      .x
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble() 
  
  obs_data$unif_0[obs_data$unif_0 < 0] = 0
  obs_data$unif_1[obs_data$unif_2 < 0] = 0
  obs_data$unif_2[obs_data$unif_2 < 0] = 0
  
  
  obs_data %>% dplyr::select(-c(tau_to_temp, temp_to_tau)) %>% 
    saveRDS(paste0("data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles))
}



generate_bts = function(bts_rng){

  setwd("~/Extreme-Irish-Summer-Temperatures/")
  source('models/marginal_models/gpd_models.R')
  
  num_quantiles = 30
  obs_data = readRDS(paste0("data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles)) %>% arrange(date)
  quantile_models = readRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))
  obs_data = obs_data %>% left_join(quantile_models) 

  sites_by_date = readRDS("data/processed/bootstrap_data/sites_by_date")
  bootstrap_dates = readRDS("data/processed/bootstrap_data/bootstrapped_dates")

  for(b in bts_rng){
    print(b)
    bts_data <- rlang::duplicate(obs_data, shallow = FALSE) # make a deep copy
    sites_by_date$bts.dat = bootstrap_dates[[b]]

    dates_and_their_data = bts_data %>%
      group_by(date) %>%
      group_map(~{
        tibble(date = .x$date[1],
               Station = list(.x$Station),
               unif_0 = list(.x$unif_0),
               unif_1 = list(.x$unif_1),
               unif_2 = list(.x$unif_2))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()

    dates_and_their_data_copy = duplicate(dates_and_their_data, shallow = FALSE)
    for(i in seq(nrow(sites_by_date))){
      data_to_replace = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]
      sampled_data = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$bts.dat,]
      sites_to_replace = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]$Station[[1]]
      sites_in_bts_date = sampled_data$Station[[1]]
      ind_of_sites_to_rep = which(sites_in_bts_date %in% sites_to_replace)
      dates_and_their_data_copy[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_0 = list(sampled_data$unif_0[[1]][ind_of_sites_to_rep])
      dates_and_their_data_copy[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_1 = list(sampled_data$unif_1[[1]][ind_of_sites_to_rep])
      dates_and_their_data_copy[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_2 = list(sampled_data$unif_2[[1]][ind_of_sites_to_rep])
    }

    expanded_data = dates_and_their_data_copy %>%
      group_by(date) %>%
      group_map(~{
        tibble(date = .x$date,
               Station = .x$Station[[1]],
               unif_0 = .x$unif_0[[1]],
               unif_1 = .x$unif_1[[1]],
               unif_2 = .x$unif_2[[1]])
      }, .keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()

    covars = obs_data %>%
      dplyr::select(Station, year, loess_temp_anom, dist_sea,
                    threshold_9, scale_9,thresh_exceedance_9,
                    "tau_to_temp", 
                    "scale_0", 
                    "scale_1",
                    "scale_2",
                    "shape_0",
                    "shape_1",
                    "shape_2") %>% unique()
    
    expanded_data = expanded_data %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(covars)

    # ---- back transform
    expanded_data%>%
      group_by(Station) %>%
      group_map(~{

        scle_0 = .x$scale_0
        shpe_0 = .x$shape_0

        scle_1 = .x$scale_1
        shpe_1 = .x$shape_1

        scle_2 = .x$scale_2
        shpe_2 = .x$shape_2

        threshold = .x$threshold_9
        my_lambda = .x$thresh_exceedance_9

        .x$maxtp_0 = NA
        .x$maxtp_1 = NA
        .x$maxtp_2 = NA

        for(i in seq(nrow(.x))){
          # --- model 0
          if(.x$unif_0[i] > (1-my_lambda[i])){
            .x$maxtp_0[i] =  evd::qgpd(p =(1 + (.x$unif_0[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_0[i],shape=shpe_0[i])
          }else{
            .x$maxtp_0[i] = .x$tau_to_temp[i][[1]](.x$unif_0[i])
          }

          # --- model 1
          if(.x$unif_1[i] > (1-my_lambda[i])){
            .x$maxtp_1[i] =  evd::qgpd( p =(1 + (.x$unif_1[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_1[i],shape=shpe_1[i])
          }else{
            .x$maxtp_1[i] = .x$tau_to_temp[i][[1]](.x$unif_1[i])
          }

          # --- model 2
          if(.x$unif_2[i] > (1-my_lambda[i])){
            .x$maxtp_2[i] =  evd::qgpd( p =(1 + (.x$unif_2[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_2[i],shape=shpe_2[i])
          }else{
            .x$maxtp_2[i] = .x$tau_to_temp[i][[1]](.x$unif_2[i])
          }
        }
        .x
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      dplyr::select(Station, date, scale_9, threshold_9, maxtp_0, maxtp_1, maxtp_2) %>%
      saveRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_",b))
  }
}


batch_size = 20
for(i in seq(1, 300, by = batch_size)){
  job::job({generate_bts(seq(i,(i+batch_size-1)))})
}
