# Author: Dáire Healy
# Date: 17 Feb 2022

# Description: This script prepares observational temperature data for marginal 
# modelling. This involves incorporating all covariates as descriped in the 
# paper associated to this work.
# Note: bulk modelling in preformed in this script. 
gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(mgcv)
library(raster) # package for netcdf manipulation
# 
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
# 
# # i. ========  ========  Global parameters ========  ========  ========
# marginal threshold 
num_quantiles = 15 
quantiles_to_estimate_bulk = seq(0.001,0.99,length.out = num_quantiles)

obs_data = read_csv("data/processed/obs_data.csv")
clim_data = read_csv("data/processed/clim_data_full.csv")






# estimate empiracle quantules for climate data
clim_quantiles = clim_data %>%
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

clim_quantiles_subset %>%
  saveRDS(paste0("data/processed/clim_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))


obs_data = obs_data %>%
  left_join(clim_quantiles_subset)

# add threshold to obs data
obs_data %>%
  #dplyr::select(-clim_thresh_value) %>%
  saveRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,".csv"))
