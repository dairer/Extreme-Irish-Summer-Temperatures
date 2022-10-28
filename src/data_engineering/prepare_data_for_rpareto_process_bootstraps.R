rm(list = ls())
gc()
library(tidyverse)

setwd("~/Extreme-Irish-Summer-Temperatures/")

prep_rpareto_bts = function(marg_mod, bts_range, num_quantiles){

  source('src/models/marginal_models/gpd_models.R')
  
  # # bulk models at each site, year for back transformation
  # obs_smoothed_quantiles = readRDS("output/quant_models")
  for(file_name in bts_range){
    print(file_name)

    
    
    # file_name = 1
    # marg_mod = "mod_1"
    # num_quantiles = 15
    
    bts_lambda = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_",marg_mod,"_num_quantiles_",num_quantiles,"_bts_", file_name,".csv")) %>%
      rename(bts = X1, Station = X2, year = X3, thresh_exceedance_9 = X4) %>%
      mutate(year = as.numeric(year),
             thresh_exceedance_9 = as.numeric(thresh_exceedance_9)) %>%
      dplyr::select(-bts)
    
    obs_data = read_csv("data/processed/obs_data.csv") %>%
      left_join(bts_lambda) %>%
      left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
                  dplyr::select(Station, dist_sea))
    
    covariates = obs_data %>%
      dplyr::select(Station, year, scale_9, threshold_9, thresh_exceedance_9, loess_temp_anom, dist_sea, Long, Lat, Long.projected, Lat.projected) %>%
      unique() 
    
    obs_smoothed_quantiles = readRDS(paste0("output/bulk_model_fits/",marg_mod,"_num_quantiles_",num_quantiles,"_bts_", file_name))
    
    this_data = c()
    pred = c()
    if(marg_mod == 'mod_1'){

      this_data <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name))%>%
        mutate(year = lubridate::year(date)) %>%
        rename(maxtp = maxtp_1) %>%
        dplyr::select(-c(maxtp_2, maxtp_4)) %>%
        dplyr::left_join(covariates %>% dplyr::select(-c(scale_9, threshold_9)), by = c('year', 'Station')) %>%
        filter(maxtp>0)
      
      dates_to_keep = this_data %>%
        group_by(date) %>%
        summarise(check_obs_for_s = 'aldergrove' %in% Station) %>%
        filter(check_obs_for_s) %>%
        pull(date)
      
      this_data <- this_data %>% filter(date %in% dates_to_keep)
      this_data$order = seq(nrow(this_data)) # remember the original order of the data
      
      this_fit_mod = read_csv("output/gpd_model_fits/model_1_pars_corrected.csv",
                              col_names = c("bts", "b0", "b1", "xi")) %>% 
        mutate(bts =  gsub("num_quantiles_50_bts_","", bts) %>% as.numeric()) %>%
        filter(bts == file_name) %>%
        unlist() %>% as.numeric %>% .[2:length(.)]
      
      pred <- my_predict_1(this_fit_mod, this_data$scale_9)
      
    }else if(marg_mod == 'mod_2'){
      
      
      this_data <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name))%>%
        mutate(year = lubridate::year(date)) %>%
        rename(maxtp = maxtp_2) %>%
        dplyr::select(-c(maxtp_1, maxtp_4)) %>%
        left_join(covariates %>% dplyr::select(-c(scale_9, threshold_9)), by = c('year', 'Station'))%>%
        filter(maxtp>0)
      
      dates_to_keep = this_data %>%
        group_by(date) %>%
        summarise(check_obs_for_s = 'aldergrove' %in% Station) %>%
        filter(check_obs_for_s) %>%
        pull(date)
      
      this_data <- this_data %>% filter(date %in% dates_to_keep)
      this_data$order = seq(nrow(this_data)) # remember the original order of the data
      
      this_fit_mod = read_csv("output/gpd_model_fits/model_2_pars_corrected.csv",
                              col_names = c("bts", "b0", "b1", "b2","xi")) %>% 
        mutate(bts =  gsub("num_quantiles_50_bts_","", bts) %>% as.numeric()) %>%
        filter(bts == file_name) %>%
        unlist() %>% as.numeric %>% .[2:length(.)]
      
      pred <- my_predict_2(this_fit_mod, this_data$scale_9, this_data$loess_temp_anom)
      
    }else if(marg_mod == 'mod_4'){

      this_data <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_50_bts_", file_name))%>%
        mutate(year = lubridate::year(date)) %>%
        rename(maxtp = maxtp_4) %>%
        dplyr::select(-c(maxtp_2, maxtp_1)) %>%
        left_join(covariates %>% dplyr::select(-c(scale_9, threshold_9)), by = c('year', 'Station'))%>%
        filter(maxtp>0)
      
      dates_to_keep = this_data %>%
        group_by(date) %>%
        summarise(check_obs_for_s = 'aldergrove' %in% Station) %>%
        filter(check_obs_for_s) %>%
        pull(date)
      
      this_data <<- this_data %>% filter(date %in% dates_to_keep)
      this_data$order = seq(nrow(this_data)) # remember the original order of the data
      
      
      this_fit_mod = read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
                              col_names = c("bts", "b0", "b1", "b2","b3","b4","xi")) %>% 
        mutate(bts =  gsub("num_quantiles_50_bts_","", bts) %>% as.numeric()) %>%
        filter(bts == file_name) %>%
        unlist() %>% as.numeric %>% .[2:length(.)]
      
      pred <- my_predict_4b(this_fit_mod, this_data$scale_9, this_data$loess_temp_anom, this_data$dist_sea)
    }

    this_data$scale = pred$scale
    this_data$shape = pred$shape
    
  obs_data_standardised = this_data %>%
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
          pull(quant_spline) %>%
          .[[1]]
        
        res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])       # alternative ... myecdf(data[data <= threshold])
        res[res<0]=0
        
        .x$unif = res
        .x$frechet_marg = -1/log(.x$unif)
        .x
        
      },.keep = T) %>%
      plyr::rbind.fill()%>%
      as_tibble() 
    
    
    print(paste0("mean: ",mean(obs_data_standardised$unif)))
    print(paste0("min: ",min(obs_data_standardised$unif)))
    print(paste0("max: ",max(obs_data_standardised$unif)))
    
    
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
      saveRDS(paste0("data/processed/data_for_rpareto/bootstraps/data_for_rpareto_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_", file_name))
  }
}


job::job({prep_rpareto_bts('mod_1', seq(1,100), 15)})
job::job({prep_rpareto_bts('mod_1', seq(101,200), 15)})
job::job({prep_rpareto_bts('mod_1', seq(201,300), 15)})
job::job({prep_rpareto_bts('mod_1', seq(301,400), 15)})
job::job({prep_rpareto_bts('mod_1', seq(401,500), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(41,80), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(81,100), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(101,140), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(141,180), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(181,200), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(201,240), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(241,280), 15)})
# job::job({prep_rpareto_bts('mod_1', seq(281,300), 15)})






# job::job({prep_rpareto_bts('mod_2', seq(1,40), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(41,80), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(81,100), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(101,140), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(141,180), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(181,200), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(201,240), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(241,280), 15)})
# job::job({prep_rpareto_bts('mod_2', seq(281,300), 15)})



job::job({prep_rpareto_bts('mod_2', seq(1,100), 15)})
job::job({prep_rpareto_bts('mod_2', seq(101,200), 15)})
job::job({prep_rpareto_bts('mod_2', seq(201,300), 15)})
job::job({prep_rpareto_bts('mod_2', seq(301,400), 15)})
job::job({prep_rpareto_bts('mod_2', seq(401,500), 15)})




job::job({prep_rpareto_bts('mod_4', seq(1,40), 25)})
job::job({prep_rpareto_bts('mod_4', seq(41,80), 25)})
job::job({prep_rpareto_bts('mod_4', seq(81,100), 25)})
job::job({prep_rpareto_bts('mod_4', seq(101,140), 25)})
job::job({prep_rpareto_bts('mod_4', seq(141,180), 25)})
job::job({prep_rpareto_bts('mod_4', seq(181,200), 25)})
job::job({prep_rpareto_bts('mod_4', seq(201,240), 25)})
job::job({prep_rpareto_bts('mod_4', seq(241,280), 25)})
job::job({prep_rpareto_bts('mod_4', seq(281,300), 25)})

job::job({prep_rpareto_bts('mod_4', seq(301,340), 25)})
job::job({prep_rpareto_bts('mod_4', seq(341,380), 25)})
job::job({prep_rpareto_bts('mod_4', seq(381,400), 25)})
job::job({prep_rpareto_bts('mod_4', seq(401,440), 25)})
job::job({prep_rpareto_bts('mod_4', seq(441,480), 25)})
job::job({prep_rpareto_bts('mod_4', seq(481,500), 25)})


job::job({prep_rpareto_bts('mod_4', seq(454,480), 25)})
job::job({prep_rpareto_bts('mod_4', seq(492,500), 25)})












res = readRDS("data/processed/bootstrap_data/sites_by_date") %>%
  left_join(readRDS("data/processed/bootstrap_data/sites_by_date") %>% 
              group_by(Stations) %>%
              group_map(~{
                tibble(Stations = .x$Stations[1], dates = list(.x$date))
              }, .keep = T) %>%
              plyr::rbind.fill() %>%
              as_tibble()) 
res$num_samples_available = res$dates %>% map(length) %>% unlist() 



2000-06-07
res %>%
  filter(lubridate::year(date) == '2000')

res %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(num_samples_available = mean(num_samples_available)) %>%
  ggplot()+
  geom_line(aes(year, (num_samples_available)))

readRDS("data/processed/bootstrap_data/sites_by_date") %>%
  mutate(len = lengths(Stations[[1]])) %>% pull(len) %>% max

res %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(num_samples_available = mean(num_samples_available)) %>%
  pull(num_samples_available)

mutate(size = length(Stations[[1]])) %>%
  tail
  group_by(size) %>%
  summarise(n())

# readRDS("data/processed/data_for_rpareto/bootstraps/data_for_rpareto_mod_1_bts_1")$extreme_dates$cost %>% mean
# readRDS("data/processed/data_for_rpareto/bootstraps/data_for_rpareto_mod_4_bts_2")$extreme_dates$cost %>% mean




