rm(list = ls())
library(tidyverse)
setwd("~/Extreme-Irish-Summer-Temperatures/")

site_pairs = read_csv("data/processed/obs_pairs_with_dist.csv")
site_pairs$V2 %>% unique %>% length
read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station) %>% unique

get_obs_chi = function(my_qnt, num_quantiles, bts_range, marg_mod){
  
  source("src/models/marginal_models/gpd_models.R")
  site_pairs = read_csv("data/processed/obs_pairs_with_dist.csv")
  
  for(file_name in bts_range){
    
    print(file_name)
    
    bts_lambda = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_",marg_mod,"_num_quantiles_",num_quantiles,"_bts_", file_name,".csv"),
                          col_names = c('bts', 'Station', 'year', 'thresh_exceedance_9')) %>%
      filter(bts == file_name) %>%
      dplyr::select(-bts)
    
    
    obs_data = read_csv("data/processed/obs_data.csv") %>%
      left_join(bts_lambda) %>%
      left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea))
    
    covariates = obs_data %>%
      dplyr::select(Station, year, scale_9, threshold_9, thresh_exceedance_9, loess_temp_anom, dist_sea, Long, Lat, Long.projected, Lat.projected) %>%
      unique() 

    this_data <- readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/num_quantiles_30_bts_", file_name))%>%
      mutate(year = lubridate::year(date)) %>%
      rename(maxtp = maxtp_2) %>%
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
    
    
    this_fit_mod = read_csv("output/gpd_model_fits/model_2_pars_corrected.csv",
                            col_names = c("bts", "b0", "b1", "b2","b3","b4","xi")) %>% 
      mutate(bts =  gsub("num_quantiles_30_bts_","", bts) %>% as.numeric()) %>%
      filter(bts == file_name) %>%
      unlist() %>% as.numeric %>% .[2:length(.)]
    
    pred <- my_predict_2b(this_fit_mod, this_data$scale_9, this_data$loess_temp_anom, this_data$dist_sea)
    
    this_data$scale = pred$scale
    this_data$shape = pred$shape
    
    obs_data_standardised = this_data %>%
      group_by(Station, year) %>%
      group_map(~{
        
        data = .x$maxtp
        threshold = .x$threshold_9
        res = rep(NA, length(data))
        num_extremes = sum(data>threshold)
        myecdf = ecdf(data)
        
        
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

        res[data <= threshold] = this_quant_mod(.x$maxtp[data <= threshold])
        res[res<0]=0
        
        .x$unif = res
        .x
        
      },.keep = T) %>%
      plyr::rbind.fill()%>%
      as_tibble()   
    
    chi.emp=function(U,V,q ){
      sum((U>=q)&(V>=q))/sum((U>=q))
    }
    
    
    site_pairs%>%
      mutate(bin = ntile(dist, 25)) %>%
      group_by(bin)%>%
      group_map(~{
        
        print(.x$bin[1])
        x1_weighted = c()
        x2_weighted = c()
        
        for(s in seq(nrow(.x))){
          s1 = obs_data_standardised %>% filter(Station == .x[s,]$V1)
          s2 = obs_data_standardised %>% filter(Station == .x[s,]$V2)
          
          x1_weighted = c(x1_weighted,  s1 %>% filter(date %in% s2$date) %>% arrange(date) %>% pull(unif))
          x2_weighted = c(x2_weighted,  s2 %>% filter(date %in% s1$date) %>% arrange(date) %>% pull(unif))
        }
        
        tibble(bts = file_name, bin = .x$bin[1], dist = mean(.x$dist), 
               chi = chi.emp(x1_weighted, x2_weighted, my_qnt)) %>%
          write_csv(paste0("data/processed/empirical_chi/chi_obs_bts_25_bins", my_qnt,"_bts_number_",file_name,".csv"), append = T)
        
        tibble()
      }, .keep = T)
  }
}




job::job({get_obs_chi(0.9, 25, seq(1, 50), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(51, 100), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(101, 150), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(151, 200), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(201, 250), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(251, 300), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(301, 350), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(351, 400), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(401, 450), 'mod_2')})
job::job({get_obs_chi(0.9, 25, seq(451, 500), 'mod_2')})


job::job({get_obs_chi(0.95, 25, seq(1, 50), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(51, 100), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(101, 150), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(151, 200), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(201, 250), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(251, 300), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(301, 350), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(351, 400), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(401, 450), 'mod_2')})
job::job({get_obs_chi(0.95, 25, seq(451, 500), 'mod_2')})



all_dat = c()
list.files("data/processed/empirical_chi/")[list.files("data/processed/empirical_chi/") %>% map(grepl, pattern='0.95_') %>% unlist] %>%
  map(~{
    all_dat <<- rbind(all_dat,read_csv(paste0("data/processed/empirical_chi/",.x), 
                                     col_names = c('bts', 'bin', 'dist', 'chi')))
    c()
  }) 


all_dat %>%
  group_by(bin) %>%
  summarise(dist = mean(dist),
            upper = quantile(chi, 0.975),
            lower = quantile(chi, 0.025)) %>%
  ggplot()+
  geom_segment(aes(x = dist, xend = dist, y = lower, yend = upper))
