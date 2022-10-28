# Author: Dáire Healy
# Date: 18 Feb 2022

# Description: This script samples blocks of time ans stores them in a list

rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland")


obs_data = read_csv("data/processed/obs_data.csv")

is.subset = function(A, B) all(A %in% B)

sites_by_date = obs_data %>% 
  dplyr::select(Station, date) %>% 
  group_by(date)%>%
  group_map(~{
    tibble(date = .x$date[1], Stations = list(.x$Station))
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

sites_by_date %>% 
  saveRDS(paste0("data/processed/bootstrap_data/sites_by_date"))





run_bts_lookup = function(bts_id){
  # i. ========  ========  Global parameters ========
  # num_bootstraps = 2500
  obs_data = read_csv("data/processed/obs_data.csv")
  
  is.subset = function(A, B) all(A %in% B)
  
  sites_by_date = readRDS("data/processed/bootstrap_data/sites_by_date")
  
  
  # groups of observed sites and dates they were observed
  site_groups_and_dates = sites_by_date %>% 
    group_by(Stations) %>%
    group_map(~{
      tibble(Stations = .x$Stations[1], dates = list(.x$date))
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble()
  
  
  site_groups_and_dates_subsets = c()
  for(i in seq(nrow(site_groups_and_dates))){
    
    group_has_subset = map(.x = site_groups_and_dates$Stations, 
                           .f = ~is.subset(A = site_groups_and_dates[i,]$Stations[[1]], 
                                           .x)) %>% unlist()
    
    all_dates_which_sites_observed = site_groups_and_dates[group_has_subset,]$dates %>% unlist %>% unique() %>% lubridate::as_date()
    site_groups_and_dates_subsets = rbind(site_groups_and_dates_subsets,
                                          tibble(Stations = site_groups_and_dates[i,]$Stations,
                                                 dates = list(all_dates_which_sites_observed)))
  }
  
  
  dates_i_can_sample = function(sites_i_need){
    dts = site_groups_and_dates_subsets[(map(.x = site_groups_and_dates_subsets$Stations, 
                                             .f = ~is.subset(A = sites_i_need, 
                                                             .x)) %>% unlist()),] %>%
      pull(dates) %>% unlist() %>% unique()
    
    
    if(is.null(dts)) NULL
    else dts %>% lubridate::as_date()
  }
  
  
  
  all_bootstraps = c()
  for(nnn in bts_id){
    print(paste0("on bootstrap number ", nnn))
    i=1
    
    dates_samples = c()
    while(length(dates_samples) <= nrow(sites_by_date)){
      
      # sample event size
      block_size = rgeom(n=1, prob = 0.2)
      sites_i_need = sites_by_date[c(i:(i+block_size)),]$Stations %>% unlist() %>% unique()
      dts_to_samp = dates_i_can_sample(sites_i_need) 
      
      if(!(is.null(dts_to_samp) | (length(dts_to_samp) < block_size))){
        
        if(length(dts_to_samp) == block_size){
          dates_samples = c(dates_samples,dts_to_samp)
          
        }else{
          position_of_start_samp = sample(seq(length(dts_to_samp) - block_size), size = 1)
          dates_samples = c(dates_samples,dts_to_samp[position_of_start_samp:(position_of_start_samp+block_size)])
        }
        
        i=length(dates_samples)+1
      }
    }
    all_bootstraps = c(all_bootstraps, list(lubridate::as_date(dates_samples)[1:nrow(sites_by_date)]))
  }
  
  saveRDS(all_bootstraps, paste0("data/processed/bootstrap_data/temp/bootstrapped_dates_",floor(runif(1)*100000)))
}



# --- takes 15 - 25 mins each
job::job({run_bts_lookup(seq(1,10))})
job::job({run_bts_lookup(seq(11,20))})
job::job({run_bts_lookup(seq(21,30))})
job::job({run_bts_lookup(seq(31,40))})
job::job({run_bts_lookup(seq(41,50))})
job::job({run_bts_lookup(seq(51,60))})
job::job({run_bts_lookup(seq(61,70))})
job::job({run_bts_lookup(seq(71,80))})
job::job({run_bts_lookup(seq(81,90))})
job::job({run_bts_lookup(seq(91,100))})

job::job({run_bts_lookup(seq(101,110))})
job::job({run_bts_lookup(seq(111,120))})
job::job({run_bts_lookup(seq(121,130))})
job::job({run_bts_lookup(seq(131,140))})
job::job({run_bts_lookup(seq(141,150))})
job::job({run_bts_lookup(seq(151,160))})
job::job({run_bts_lookup(seq(161,170))})
job::job({run_bts_lookup(seq(171,180))})
job::job({run_bts_lookup(seq(181,190))})
job::job({run_bts_lookup(seq(191,200))})

job::job({run_bts_lookup(seq(201,210))})
job::job({run_bts_lookup(seq(211,220))})
job::job({run_bts_lookup(seq(221,230))})
job::job({run_bts_lookup(seq(231,240))})
job::job({run_bts_lookup(seq(241,250))})
job::job({run_bts_lookup(seq(251,260))})
job::job({run_bts_lookup(seq(261,270))})
job::job({run_bts_lookup(seq(271,280))})
job::job({run_bts_lookup(seq(281,290))})
job::job({run_bts_lookup(seq(291,300))})

job::job({run_bts_lookup(seq(301,310))})
job::job({run_bts_lookup(seq(311,320))})
job::job({run_bts_lookup(seq(321,330))})
job::job({run_bts_lookup(seq(331,340))})
job::job({run_bts_lookup(seq(341,350))})
job::job({run_bts_lookup(seq(351,360))})
job::job({run_bts_lookup(seq(361,370))})
job::job({run_bts_lookup(seq(371,380))})
job::job({run_bts_lookup(seq(381,390))})
job::job({run_bts_lookup(seq(391,400))})


job::job({run_bts_lookup(seq(401,410))})
job::job({run_bts_lookup(seq(411,420))})
job::job({run_bts_lookup(seq(421,430))})
job::job({run_bts_lookup(seq(431,440))})
job::job({run_bts_lookup(seq(441,450))})
job::job({run_bts_lookup(seq(451,460))})
job::job({run_bts_lookup(seq(461,470))})
job::job({run_bts_lookup(seq(471,480))})
job::job({run_bts_lookup(seq(481,490))})
job::job({run_bts_lookup(seq(491,500))})


list.files("data/processed/bootstrap_data/temp/") %>%
  map(~{
    paste0('data/processed/bootstrap_data/temp/',.x) %>% print()
  })



all_bts = c()
files_to_read = list.files("data/processed/bootstrap_data/temp/") 

for(fff in files_to_read){
  all_bts <- c(all_bts, readRDS(paste0("data/processed/bootstrap_data/temp/", fff)))
}
saveRDS(all_bts, "data/processed/bootstrap_data/bootstrapped_dates")


