
# Description: This script provides comparison between potential marginal GPD
# models using bootstrapped AIC,BIC and log likelihood. Bootstraps are 
# constructed to preserve spatial and temporal dependence.
gc()
rm(list = ls())
setwd("~/Extreme-Irish-Summer-Temperatures/data/processed/")
source('marginal_models/tail/gpd_models.R')

library(tidyverse)
library(spatialsample)
library(scoringRules)
library(job)

# colour palatte 
my_pal = c( 
  '#062c30', # extra dark 
  '#003f5c',
  '#2f4b7c',
  '#665191',
  '#a05195',
  '#d45087',
  '#f95d6a',
  '#ff7c43',
  '#ffa600')


# map of ireland sf
ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]


# ========  ========  Global parameters ========

# ========  ========  Read in data  ========  ========
# Observational data with covariates
obs_data = read_csv("obs_data.csv") %>%
  left_join(read_csv("obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea), by = 'Station')


obs_data %>%
  filter(is.na(dist_sea)) %>%
  pull(Station) %>%
  unique
obs_data = obs_data %>%
  mutate(week = lubridate::week(date))

# onlykeek full years for cv
obs_data = obs_data %>%
  group_by(year, Station) %>%
  summarise(n = n()) %>%
  arrange(n) %>%
  filter(n > 90) %>% 
  left_join(obs_data) %>%
  ungroup()

obs_data$log_dist_sea = log(obs_data$dist_sea)

obs_data$new_covar = obs_data$log_dist_sea*obs_data$loess_temp_anom

# -------- CREATE CV Folds
set.seed(51964)

num_spatial_folds  = 30
clutered = spatial_clustering_cv(read_csv("obs_data_dist_to_sea.csv"), coords = c(Long, Lat), v = num_spatial_folds)

spatial_folds = c()
for(i in seq(num_spatial_folds)){
  spatial_folds = rbind(spatial_folds, 
                        assessment(clutered$splits[i][[1]]) %>% 
                          dplyr::select(Station, Long, Lat) %>% 
                          unique %>%
                          mutate(spatial_fold = i))
}

week_chunks = list(c(22,25,28,31,34),
                   c(23,26,29,32,35),
                   c(24,27,30,33))

obs_data$temporal_fold = 123456789
obs_data[obs_data$week %in% week_chunks[[1]],]$temporal_fold=1
obs_data[obs_data$week %in% week_chunks[[2]],]$temporal_fold=2
obs_data[obs_data$week %in% week_chunks[[3]],]$temporal_fold=3

get_metrics = function(orig_dat, excess, quant, threshold,  pred_scale, pred_shape){
  likelihood_vals = evd::dgpd(x = excess, loc = 0, scale = pred_scale, shape = pred_shape, log = T)
  ll = sum(likelihood_vals[!is.infinite(likelihood_vals)])
  ll_standardised = sum(likelihood_vals[!is.infinite(likelihood_vals)])/length(likelihood_vals[!is.infinite(likelihood_vals)])
  
  # - RMSE
  x = orig_dat
  y = evd::qgpd(p = quant, loc = threshold, scale = pred_scale, shape = pred_shape[1])
  
  my_rmse = sqrt(mean((x-y)^2))
  scoring=0
  scoring = crps(y = excess, family = "gpd", location = 0, scale = pred_scale, shape = pred_shape[1], mass = 0) %>% mean
  return(paste0(ll, ",", ll_standardised, ",",my_rmse, ",", scoring))
}



run_cv = function(cv_method, thresh_qnt, obs_data, spatial_folds, get_metrics, num_spatial_folds = "NA", week_chunks="NA"){
  
  setwd("~/Extreme-Irish-Summer-Temperatures/data/processed/")
  source('marginal_models/tail/gpd_models.R')
  
  get_metrics = function(orig_dat, excess, quant, threshold,  pred_scale, pred_shape){
    likelihood_vals = evd::dgpd(x = excess, loc = 0, scale = pred_scale, shape = pred_shape, log = T)
    ll = sum(likelihood_vals[!is.infinite(likelihood_vals)])
    ll_standardised = sum(likelihood_vals[!is.infinite(likelihood_vals)])/length(likelihood_vals[!is.infinite(likelihood_vals)])
    
    # - RMSE
    x = orig_dat
    y = evd::qgpd(p = quant, loc = threshold, scale = pred_scale, shape = pred_shape[1])
    
    my_rmse = sqrt(mean((x-y)^2))
    scoring = crps(y = excess, family = "gpd", location = 0, scale = pred_scale, shape = pred_shape[1], mass = 0) %>% mean
    return(paste0(ll, ",", ll_standardised, ",",my_rmse, ",", scoring))
  }
  

  if(cv_method == "spatial-temporal"){
    
    extreme_data = obs_data %>%
      mutate(excess = maxtp - threshold_9) %>%
      filter(excess > 0) %>%
      left_join(spatial_folds) %>%
      group_by(year, Station) %>%
      group_map(~{
        
        .x %>%
          mutate(quant = rank(excess)/(length(excess)+1))
        
      },.keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    
    for(s in seq(num_spatial_folds)){
      for(t in week_chunks){

          test = extreme_data %>% 
          filter(spatial_fold == s) %>%
          filter(week %in% t)
        
        train = extreme_data %>%
          anti_join(test)
        
        this_fit_mod_0 = fit_mod_0(train$excess, train$scale_9)
        this_fit_mod_1 = fit_mod_1(train$excess, train$scale_9, train$loess_temp_anom)
        this_fit_mod_2 = fit_mod_2(train$excess, train$scale_9, train$loess_temp_anom, train$dist_sea)

        # calculate scale parameter on climate grid
        pred_0 = my_predict_0(this_fit_mod_0, test$scale_9)
        pred_1 = my_predict_1(this_fit_mod_1, test$scale_9, test$loess_temp_anom)
        pred_2 =my_predict_2(this_fit_mod_2, test$scale_9, test$loess_temp_anom, test$dist_sea)
        
        write.table(paste0("model.0", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_0$scale, pred_0$shape[1]),",",mean(t)),
                    file = paste0('~/Extreme-Irish-Summer-Temperatures/output/cv_res/spatio_temporal_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
        
        write.table(paste0("model.1", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_1$scale, pred_1$shape[1]),",",mean(t)),
                    file = paste0('~/Extreme-Irish-Summer-Temperatures/output/cv_res/spatio_temporal_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)

        write.table(paste0("model.2", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_2$scale, pred_2$shape[1]),",",mean(t)),
                    file = paste0('~/Extreme-Irish-Summer-Temperatures/output/cv_res/spatio_temporal_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
      }
    }
  }
  
  if(cv_method == "10fold"){
    
    extreme_data = obs_data %>%
      mutate(excess = maxtp - threshold_9) %>%
      filter(excess > 0) %>%
      group_by(year, Station) %>%
      group_map(~{
        
        .x %>%
          mutate(quant = rank(excess)/(length(excess)+1))
        
      },.keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    
    #cv_method = '10fold'
    # ---------- define random folds
    num_random_folds = 90
    set.seed(1234567)
    extreme_data$random_fold = sample(seq(num_random_folds), size = nrow(extreme_data), replace = T)
    
    
    for(i in seq(num_random_folds)){
      test = extreme_data %>% filter(random_fold == i)
      train = extreme_data %>% filter(random_fold != i)
      
      this_fit_mod_0 = fit_mod_0(train$excess, train$scale_9)
      this_fit_mod_1 = fit_mod_1(train$excess, train$scale_9, train$loess_temp_anom)
      this_fit_mod_2= fit_mod_2(train$excess, train$scale_9, train$loess_temp_anom, train$dist_sea)
      
      # calculate scale parameter on climate grid
      pred_0 = my_predict_0(this_fit_mod_0, test$scale_9)
      pred_1 = my_predict_1(this_fit_mod_1, test$scale_9, test$loess_temp_anom)
      pred_2 =my_predict_2(this_fit_mod_2, test$scale_9, test$loess_temp_anom, test$dist_sea)
      
      write.table(paste0("model.0", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_0$scale, pred_0$shape[1])),
                  file = paste0('~/Extreme-Irish-Summer-Temperatures/output/cv_res/ten_fold_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                  sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
      
      write.table(paste0("model.1", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_1$scale, pred_1$shape[1])),
                  file = paste0('~/Extreme-Irish-Summer-Temperatures/output/cv_res/ten_fold_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                  sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
      
      write.table(paste0("model.2", ",", get_metrics(test$maxtp, test$excess, test$quant, test$threshold_9,  pred_2$scale, pred_2$shape[1])),
                  file = paste0('~/Extreme-Irish-Summer-Temperatures/output/cv_res/ten_fold_cv_test_mod_n', str_remove(thresh_qnt, "0."), '.csv'),
                  sep = ",", append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
  }
}


job::job({run_cv('10fold', 0.9, obs_data, get_metrics)}, import = c("obs_data", 'run_cv', "get_metrics"))
job::job({run_cv('spatial-temporal', 0.9, obs_data, spatial_folds, get_metrics, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', "spatial_folds",  "get_metrics", "num_spatial_folds", "week_chunks"))
