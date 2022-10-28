# Description: get bootstrapped data sets under each marginal model
# --- use true gpd and bts models to standardise, sample and back transform 
# --- to create bootstrap sampels

gc() 
rm(list = ls(all.names = TRUE))
library(tidyverse)
source('~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/src/models/marginal_models/gpd_models.R')
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")

standardise_data = F
if(standardise_data){
  num_quantiles = 100
  obs_data = read_csv("data/processed/obs_data.csv") %>%
    left_join(read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv")))  %>%
    left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea), by = 'Station') %>%
    dplyr::select(-c(loess_glob_temp_anom,id ))
  
  quantile_models = readRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))
  # obs_data = obs_data %>% left_join(quantile_models) 

  extreme_data = obs_data %>%
    group_by(Station) %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)

  this_fit_mod_1 = fit_mod_1(extreme_data$excess, extreme_data$scale_9)
  this_fit_mod_2 = fit_mod_2(extreme_data$excess, extreme_data$scale_9, extreme_data$loess_temp_anom)
  this_fit_mod_4 = fit_mod_4b(extreme_data$excess, extreme_data$scale_9, extreme_data$loess_temp_anom, extreme_data$dist_sea)
  rm(extreme_data)
  
  pred_1 = my_predict_1(this_fit_mod_1, obs_data$scale_9)
  pred_2 = my_predict_2(this_fit_mod_2, obs_data$scale_9, obs_data$loess_temp_anom)
  pred_4 =my_predict_4b(this_fit_mod_4, obs_data$scale_9, obs_data$loess_temp_anom, obs_data$dist_sea)
  
  obs_data = obs_data %>%
    mutate(scale_1 = pred_1$scale,
           scale_2 = pred_2$scale,
           scale_4 = pred_4$scale,
           shape_1 = pred_1$shape,
           shape_2 = pred_2$shape,
           shape_4 = pred_4$shape)
  
  rm('this_fit_mod_1', 'this_fit_mod_2',  'this_fit_mod_4',  'pred_1', 'pred_2', 'pred_4')

  print("Joining bulk models")
  obs_data = obs_data %>%
    left_join(quantile_models)
  
  obs_data = obs_data %>%
    group_by(Station) %>%
    group_map(~{
      print(.x$Station[1])
      data = .x$maxtp
      threshold = .x$threshold_9
      #threshold = rep(quantile(.x$maxtp, 0.9), length(data))
      scle_1 = .x$scale_1
      shpe_1 = .x$shape_1
      
      scle_2 = .x$scale_2
      shpe_2 = .x$shape_2
      
      scle_4 = .x$scale_4
      shpe_4 = .x$shape_4
      
      yr = .x$year
      my_lambda = .x$thresh_exceedance_9
      
      # qnt_mods = quantile_models %>%
      #   filter(Station == .x$Station[1])

      res_1 = rep(NA, length(data))
      res_2 = rep(NA, length(data))
      res_4 = rep(NA, length(data))

      for(i in seq(length(data))){
        if(data[i] > threshold[i]){ # transform tail
          
          res_1[i] = 1 - (my_lambda[i]) *(1+shpe_1[i]* ((data[i] - threshold[i])/scle_1[i]))^(-1/(shpe_1[i]))
          res_2[i] = 1 - (my_lambda[i]) *(1+shpe_2[i]* ((data[i] - threshold[i])/scle_2[i]))^(-1/(shpe_2[i]))
          res_4[i] = 1 - (my_lambda[i]) *(1+shpe_4[i]* ((data[i] - threshold[i])/scle_4[i]))^(-1/(shpe_4[i]))
          
        }else{ # transform bulk
          # 
          # this_temp_to_tau = qnt_mods %>%
          #   filter(year == yr[i]) %>%
          #   pull(temp_to_tau) %>% .[[1]]

          res_1[i] = .x$temp_to_tau[i][[1]](data[i])
          res_2[i] = res_1[i]
          res_4[i] = res_1[i]
        }
      }
      .x$unif_1 = res_1
      .x$unif_2 = res_2
      .x$unif_4 = res_4
      .x
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble() 
  
  obs_data$unif_1[obs_data$unif_1 < 0] = 0
  obs_data$unif_2[obs_data$unif_4 < 0] = 0
  obs_data$unif_4[obs_data$unif_4 < 0] = 0
  
  
  obs_data %>% dplyr::select(-c(tau_to_temp, temp_to_tau)) %>% 
    saveRDS(paste0("data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles))
}






generate_bts = function(bts_rng){

  setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
  source('~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/src/models/marginal_models/gpd_models.R')
  
  num_quantiles = 50
  obs_data = readRDS(paste0("data/processed/standardised_data_for_bootstrapping_num_quantiles_", num_quantiles)) %>% arrange(date)
  quantile_models = readRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))
  obs_data = obs_data %>% left_join(quantile_models) 

  #
  # obs_data = read_csv("data/processed/obs_data.csv") %>%
  #   left_join(read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_",num_quantiles,".csv")))  %>%
  #   left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea), by = 'Station') %>%
  #   dplyr::select(-c(loess_glob_temp_anom,id ))
  #
  #
  #
  #
  # read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea)
  #
  #
  #
  #
  #
  #
  # quantile_models  = readRDS(paste0("output/quant_models_num_quantiles_",num_quantiles,".csv"))
  # obs_data = obs_data %>% left_join(quantile_models)
  #
  # # fit each marginal model
  # extreme_data = obs_data %>%
  #   group_by(Station) %>%
  #   mutate(excess = maxtp - threshold_9) %>%
  #   filter(excess > 0)
  # extreme_data %>%
  #   filter(is.na(dist_sea))
  #
  # this_fit_mod_1 = fit_mod_1(extreme_data$excess, extreme_data$scale_9)
  # this_fit_mod_2 = fit_mod_2(extreme_data$excess, extreme_data$scale_9, extreme_data$loess_temp_anom)
  # this_fit_mod_4 = fit_mod_4b(extreme_data$excess, extreme_data$scale_9, extreme_data$loess_temp_anom, extreme_data$dist_sea)
  # rm(extreme_data)
  #
  #
  # pred_1 = my_predict_1(this_fit_mod_1, obs_data$scale_9)
  # pred_2 = my_predict_2(this_fit_mod_2, obs_data$scale_9, obs_data$loess_temp_anom)
  # pred_4 =my_predict_4b(this_fit_mod_4, obs_data$scale_9, obs_data$loess_temp_anom, obs_data$dist_sea)
  #
  # obs_data = obs_data %>%
  #   mutate(scale_1 = pred_1$scale,
  #          scale_2 = pred_2$scale,
  #          scale_4 = pred_4$scale,
  #          shape_1 = pred_1$shape,
  #          shape_2 = pred_2$shape,
  #          shape_4 = pred_4$shape)
  #
  #
  # rm('this_fit_mod_1', 'this_fit_mod_2',  'this_fit_mod_4',  'pred_1', 'pred_2', 'pred_4')
  #
  #
  #
  # obs_data = obs_data %>%
  #   group_by(Station) %>%
  #   group_map(~{
  #
  #     data = .x$maxtp
  #     threshold = .x$threshold_9
  #     #threshold = rep(quantile(.x$maxtp, 0.9), length(data))
  #     scle_1 = .x$scale_1
  #     shpe_1 = .x$shape_1
  #
  #     scle_2 = .x$scale_2
  #     shpe_2 = .x$shape_2
  #
  #     scle_4 = .x$scale_4
  #     shpe_4 = .x$shape_4
  #
  #     my_lambda = .x$thresh_exceedance_9
  #     res_1 = rep(NA, length(data))
  #     res_2 = rep(NA, length(data))
  #     res_4 = rep(NA, length(data))
  #
  #
  #     for(i in seq(length(data))){
  #       if(data[i] > threshold[i]){ # transform tail
  #
  #         res_1[i] = 1 - (my_lambda[i]) *(1+shpe_1[i]* ((data[i] - threshold[i])/scle_1[i]))^(-1/(shpe_1[i]))
  #         res_2[i] = 1 - (my_lambda[i]) *(1+shpe_2[i]* ((data[i] - threshold[i])/scle_2[i]))^(-1/(shpe_2[i]))
  #         res_4[i] = 1 - (my_lambda[i]) *(1+shpe_4[i]* ((data[i] - threshold[i])/scle_4[i]))^(-1/(shpe_4[i]))
  #
  #       }else{ # transform bulk
  #         res_1[i] = .x$temp_to_tau[i][[1]](data[i])
  #         res_2[i] = res_1[i]
  #         res_4[i] = res_1[i]
  #       }
  #     }
  #     .x$unif_1 = res_1
  #     .x$unif_2 = res_2
  #     .x$unif_4 = res_4
  #     .x
  #   }, .keep = T) %>%
  #   plyr::rbind.fill() %>%
  #   as_tibble()


  sites_by_date = readRDS("data/processed/bootstrap_data/sites_by_date")
  bootstrap_dates = readRDS("data/processed/bootstrap_data/bootstrapped_dates")
  #obs_data = readRDS("data/processed/standardised_data_for_bootstrapping") %>% arrange(date)

  for(b in bts_rng){
    print(b)
    bts_data <- rlang::duplicate(obs_data, shallow = FALSE) # make a deep copy
    sites_by_date$bts.dat = bootstrap_dates[[b]]

    dates_and_their_data = bts_data %>%
      group_by(date) %>%
      group_map(~{
        tibble(date = .x$date[1],
               Station = list(.x$Station),
               unif_1 = list(.x$unif_1),
               unif_2 = list(.x$unif_2),
               unif_4 = list(.x$unif_4))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()

    for(i in seq(nrow(sites_by_date))){
      data_to_replace = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]
      sampled_data = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$bts.dat,]
      sites_to_replace = dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]$Station[[1]]
      sites_in_bts_date = sampled_data$Station[[1]]
      ind_of_sites_to_rep = which(sites_in_bts_date %in% sites_to_replace)
      #bootstrapped_temp = sampled_data$temp[[1]][ind_of_sites_to_rep]
      #dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]$unif = list(sampled_data$unif[[1]][ind_of_sites_to_rep])

      dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_1 = list(sampled_data$unif_1[[1]][ind_of_sites_to_rep])
      dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_2 = list(sampled_data$unif_2[[1]][ind_of_sites_to_rep])
      dates_and_their_data[dates_and_their_data$date == sites_by_date[i,]$date,]$unif_4 = list(sampled_data$unif_4[[1]][ind_of_sites_to_rep])
    }


    expanded_data = dates_and_their_data %>%
      group_by(date) %>%
      group_map(~{
        tibble(date = .x$date,
               Station = .x$Station[[1]],
               unif_1 = .x$unif_1[[1]],
               unif_2 = .x$unif_2[[1]],
               unif_4 = .x$unif_4[[1]])
      }, .keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()

    covars = obs_data %>%
      dplyr::select(Station, year, loess_temp_anom, dist_sea,
                    threshold_9, scale_9,thresh_exceedance_9,
                    "tau_to_temp", "scale_1",
                    "scale_2","scale_4","shape_1","shape_2","shape_4") %>% unique()
    
    expanded_data = expanded_data %>%
      mutate(year = lubridate::year(date)) %>%
      left_join(covars)

    # ---- back transform
    expanded_data%>%
      group_by(Station) %>%
      group_map(~{

        scle_1 = .x$scale_1
        shpe_1 = .x$shape_1

        scle_2 = .x$scale_2
        shpe_2 = .x$shape_2

        scle_4 = .x$scale_4
        shpe_4 = .x$shape_4

        threshold = .x$threshold_9
        my_lambda = .x$thresh_exceedance_9

        .x$maxtp_1 = NA
        .x$maxtp_2 = NA
        .x$maxtp_4 = NA

        for(i in seq(nrow(.x))){
          # --- model 1
          if(.x$unif_1[i] > (1-my_lambda[i])){
            .x$maxtp_1[i] =  evd::qgpd(p =(1 + (.x$unif_1[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_1[i],shape=shpe_1[i])
          }else{
            .x$maxtp_1[i] = .x$tau_to_temp[i][[1]](.x$unif_1[i])
          }

          # --- model 2
          if(.x$unif_2[i] > (1-my_lambda[i])){
            .x$maxtp_2[i] =  evd::qgpd( p =(1 + (.x$unif_2[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_2[i],shape=shpe_2[i])
          }else{
            .x$maxtp_2[i] = .x$tau_to_temp[i][[1]](.x$unif_2[i])
          }

          # --- model 4
          if(.x$unif_4[i] > (1-my_lambda[i])){
            .x$maxtp_4[i] =  evd::qgpd( p =(1 + (.x$unif_4[i]-1)/(my_lambda[i])), loc=threshold[i],scale=scle_4[i],shape=shpe_4[i])
          }else{
            .x$maxtp_4[i] = .x$tau_to_temp[i][[1]](.x$unif_4[i])
          }
        }
        .x
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      dplyr::select(Station, date, scale_9, threshold_9, maxtp_1, maxtp_2, maxtp_4) %>%
      saveRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_",num_quantiles,"_bts_",b))
  }
}


batch_size = 5
for(i in seq(1, 50, by = batch_size)){
  job::job({generate_bts(seq(i,(i+batch_size-1)))})
}

batch_size = 5
for(i in seq(51, 100, by = batch_size)){
  job::job({generate_bts(seq(i,(i+batch_size-1)))})
}


batch_size = 10
for(i in seq(101, 200, by = batch_size)){
  job::job({generate_bts(seq(i,(i+batch_size-1)))})
}

batch_size = 10
for(i in seq(201, 300, by = batch_size)){
  job::job({generate_bts(seq(i,(i+batch_size-1)))})
}

batch_size = 10
for(i in seq(301, 400, by = batch_size)){
  job::job({generate_bts(seq(i,(i+batch_size-1)))})
}
# ----- here!

batch_size = 10
for(i in seq(401, 500, by = batch_size)){
  job::job({generate_bts(seq(i,(i+batch_size-1)))})
}


job::job({generate_bts(seq(488,494))})
job::job({generate_bts(seq(495,500))})

# # job::job({generate_bts(seq(1,10))}) # ---- 48  to 1hr 40
# # job::job({generate_bts(seq(11,20))})
# # job::job({generate_bts(seq(21,30))})
# # job::job({generate_bts(seq(31,40))})
# # job::job({generate_bts(seq(41,50))})
# # job::job({generate_bts(seq(51,60))})
# # job::job({generate_bts(seq(61,70))})
# # job::job({generate_bts(seq(71,80))})
# # job::job({generate_bts(seq(81,90))})
# # job::job({generate_bts(seq(91,100))})

# # job::job({generate_bts(seq(101,110))})
# # job::job({generate_bts(seq(111,120))})
# # job::job({generate_bts(seq(121,130))})
# # job::job({generate_bts(seq(131,140))})
# # job::job({generate_bts(seq(141,150))})
# # job::job({generate_bts(seq(151,160))})
# # job::job({generate_bts(seq(161,170))})
# # job::job({generate_bts(seq(171,180))})
# # job::job({generate_bts(seq(181,190))})
# # job::job({generate_bts(seq(191,200))})

# # job::job({generate_bts(seq(201,220))})# --- 1hr 45
# # job::job({generate_bts(seq(221,240))})
# # job::job({generate_bts(seq(241,260))})
# # job::job({generate_bts(seq(261,280))})

# # job::job({generate_bts(seq(281,300))})