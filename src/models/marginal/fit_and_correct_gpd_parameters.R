# Author --- Dáire Healy
# Date ---  15/06/2022
# Description --- This script fits all potential marginal gpd models to each 
# bootstrap and predicts scale and shape parameter from the GPD on the climate 
# model grid

rm(list=ls()) 
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
library(tidyverse)
fit_uncorrected_models =T
fit_true_models = T
 
source('src/models/marginal_models/gpd_models.R')

obs_data = vroom::vroom("data/processed/obs_data.csv") %>%
  left_join(vroom::vroom("data/processed/thresh_exceedance_lambda.csv"))  %>%
  left_join(vroom::vroom("data/processed/obs_data_dist_to_sea.csv") %>% dplyr::select(Station, dist_sea), by = 'Station')

covars = obs_data %>%
  dplyr::select(Station, year, scale_9, dist_sea, loess_temp_anom, thresh_exceedance_9, threshold_9) %>%
  unique()

# ----- Fit true model
if(fit_true_models){
  extreme_dat_true = obs_data %>%
    mutate(excess = maxtp - threshold_9) %>%
    filter(excess > 0)
  
  fit_mod_1(extreme_dat_true$excess, extreme_dat_true$scale_9)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_1_true.csv")

  fit_mod_2(extreme_dat_true$excess, extreme_dat_true$scale_9,
            extreme_dat_true$loess_temp_anom)  %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_2_true.csv")

  # fit_mod_4b(extreme_dat_true$excess, extreme_dat_true$scale_9,
  #            extreme_dat_true$loess_temp_anom, extreme_dat_true$dist_sea)%>%
  #   matrix() %>% t() %>% as.data.frame() %>%
  #   write_csv("output/gpd_model_fits/model_4_true.csv")
}


model_1_true = readr::read_csv("output/gpd_model_fits/model_1_true.csv") %>%
  rename(b0 = V1, b1 = V2, xi = V3)

model_2_true = readr::read_csv("output/gpd_model_fits/model_2_true.csv") %>%
  rename(b0 = V1, b1 = V2, b2 = V3, xi = V4)
# 
# model_4_true = readr::read_csv("output/gpd_model_fits/model_4_true.csv") %>%
#   rename(b0 = V1, b1 = V2, b2 = V3, b3 = V4, b4 = V5, xi = V6)



bts_files = list.files("data/processed/bootstrap_data/bts_under_gpd_models/")
# bts_files = bts_files[!(bts_files %in% (fits$bts %>% unique()))]

if(fit_uncorrected_models){ # takes 30 mins
  # file.remove("output/gpd_model_fits/model_1_uncorrected.csv")
  # file.remove("output/gpd_model_fits/model_2_uncorrected.csv")
  # file.remove("output/gpd_model_fits/model_4_uncorrected.csv")
  
  # ----- fit all models
  for(file_name in bts_files){

    print(paste0("fitting to bootstrap ", file_name))
    dat = readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/", file_name))%>%
      mutate(year = lubridate::year(date)) %>%
      left_join(covars %>% dplyr::select(-c(scale_9, threshold_9)), by = c('year', 'Station'))
    
    # --- remove observations poorly interpolated by quantile model (i.e. <0)
    dat = dat %>% filter(maxtp_1 > 0,
                         maxtp_2 > 0,
                         maxtp_4 > 0)
    
    # ----- Model 1
    obs_data_to_pred = dat %>%
      mutate(excess = maxtp_1 - threshold_9) %>%
      filter(excess > 0)

    this_fit_mod_1 = fit_mod_1(obs_data_to_pred$excess, obs_data_to_pred$scale_9,
                              initial_pars = as.numeric(unlist(model_1_true)))

    c(file_name, this_fit_mod_1) %>% matrix() %>% t() %>% as.data.frame() %>%
      write_csv("output/gpd_model_fits/model_1_uncorrected.csv", append = T)

    
    # ----- Model 2
    obs_data_to_pred = dat %>%
      mutate(excess = maxtp_2 - threshold_9) %>%
      filter(excess > 0)
    this_fit_mod_2 = fit_mod_2(obs_data_to_pred$excess, obs_data_to_pred$scale_9,
                               obs_data_to_pred$loess_temp_anom,
                               initial_pars = as.numeric(unlist(model_2_true)))
    c(file_name, this_fit_mod_2)  %>% matrix() %>% t() %>% as.data.frame() %>%
      write_csv("output/gpd_model_fits/model_2_uncorrected.csv", append = T)

    
    # ----- Model 4
    # obs_data_to_pred = dat %>%
    #   mutate(excess = maxtp_4 - threshold_9) %>%
    #   filter(excess > 0)
    # this_fit_mod_4 = fit_mod_4b(obs_data_to_pred$excess, obs_data_to_pred$scale_9,
    #                             obs_data_to_pred$loess_temp_anom,
    #                             obs_data_to_pred$dist_sea,
    #                             initial_pars = as.numeric(unlist(model_4_true)))
    # c(file_name, this_fit_mod_4)  %>% matrix() %>% t() %>% as.data.frame() %>%
    #   write_csv("output/gpd_model_fits/model_4_uncorrected.csv", append = T)
  }
}


read_csv("output/gpd_model_fits/model_1_uncorrected.csv")



# ----  calculate correction terms
 model_1_xi = model_1_true %>% pull(xi)
 model_2_xi = model_2_true %>% pull(xi)
#model_4_xi = model_4_true %>% pull(xi)


model_1_correction = model_1_xi - (read_csv("output/gpd_model_fits/model_1_uncorrected.csv",
                                            col_names = c('bts', 'b0', 'b1', 'xi')) %>%
                                     pull(xi) %>% mean)
#
model_2_correction = model_2_xi - (read_csv("output/gpd_model_fits/model_2_uncorrected.csv",
                                            col_names = c('bts', 'b0', 'b1',  'b2', 'xi')) %>%
                                     pull(xi) %>%
                                     mean)

model_4_correction = model_4_xi - (read_csv("output/gpd_model_fits/model_4_uncorrected.csv",
                                            col_names = c('bts', 'b0', 'b1', 'b2','b3','b4', 'xi')) %>%
                                     pull(xi) %>% mean)



# ------- fit corrected
for(file_name in bts_files){

  print(paste0("fitting to bootstrap ", file_name))
  dat = readRDS(paste0("data/processed/bootstrap_data/bts_under_gpd_models/", file_name))%>%
    mutate(year = lubridate::year(date)) %>%
    left_join(covars %>% dplyr::select(-c(scale_9, threshold_9)), by = c('year', 'Station'))
  
  # --- remove observations poorly interpolated by quantile model (i.e. <0)
  dat = dat %>% filter(maxtp_1 > 0,
                       maxtp_2 > 0,
                       maxtp_4 > 0)
  
  # # ----- Model 1
  obs_data_to_pred = dat %>%
    mutate(excess = maxtp_1 - threshold_9) %>%
    filter(excess > 0)

  uncorrected_1_xi = read_csv("output/gpd_model_fits/model_1_uncorrected.csv",
                              col_names = c('bts', 'b0', 'b1', 'xi'))  %>%
    filter(bts == file_name) %>%
    pull(xi)

  this_fit_mod_1_corrected = fit_mod_1_fix_shape(obs_data_to_pred$excess,
                                                 obs_data_to_pred$scale_9,
                                                 this_shape_est = (uncorrected_1_xi + model_1_correction),
                                                 initial_pars = as.numeric(unlist(model_1_true))[-ncol(model_1_true)] + rnorm(1, 0, 0.01))
  c(file_name, this_fit_mod_1_corrected, (uncorrected_1_xi + model_1_correction)) %>%
    matrix() %>% t() %>% as.data.frame() %>%
    write_csv("output/gpd_model_fits/model_1_pars_corrected.csv", append = T)
  
  
  # # ----- Model 2
  obs_data_to_pred = dat %>%
    mutate(excess = maxtp_2 - threshold_9) %>%
    filter(excess > 0)

  uncorrected_2_xi = read_csv("output/gpd_model_fits/model_2_uncorrected.csv",
                              col_names = c('bts', 'b0', 'b1', 'b2', 'xi'))  %>%
    filter(bts == file_name) %>%
    pull(xi)

  this_fit_mod_2_corrected = fit_mod_2_fix_shape(obs_data_to_pred$excess,
                                                 obs_data_to_pred$scale_9,
                                                 obs_data_to_pred$loess_temp_anom,
                                                 this_shape_est = (uncorrected_2_xi + model_2_correction),
                                                 initial_pars =as.numeric(unlist(model_2_true))[-ncol(model_2_true)])

  c(file_name, this_fit_mod_2_corrected, (uncorrected_2_xi +  model_2_correction)) %>%
    matrix() %>% t() %>% as.data.frame() %>% write_csv("output/gpd_model_fits/model_2_pars_corrected.csv", append = T)
  # 
  # ----- Model 4
  # obs_data_to_pred = dat %>%
  #   mutate(excess = maxtp_4 - threshold_9) %>%
  #   filter(excess > 0)
  # 
  # uncorrected_4_xi = read_csv("output/gpd_model_fits/model_4_uncorrected.csv",
  #                             col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'xi'))  %>%
  #   filter(bts == file_name) %>% pull(xi)
  # 
  # this_fit_mod_4_corrected = fit_mod_4b_fix_shape(obs_data_to_pred$excess,
  #                                                 obs_data_to_pred$scale_9,
  #                                                 obs_data_to_pred$loess_temp_anom,
  #                                                 obs_data_to_pred$dist_sea,
  #                                                 this_shape_est = (uncorrected_4_xi + model_4_correction),
  #                                                 initial_pars =as.numeric(unlist(model_4_true))[-ncol(model_4_true)])
  # 
  # c(file_name, this_fit_mod_4_corrected, (uncorrected_4_xi +  model_4_correction)) %>%
  #   matrix() %>% t() %>% as.data.frame() %>% write_csv("output/gpd_model_fits/model_4_pars_corrected.csv", append = T)
}


# 
# # # # ----- Visualise results




bias_correction_mod_2 = gridExtra::grid.arrange(read_csv("output/gpd_model_fits/model_4_uncorrected.csv",
         col_names = c('bts', 'β0', 'β1', 'β2', 'β3', 'β4', 'ξ')) %>%
  pivot_longer(-bts) %>%
  ggplot()+
  geom_density(aes(value))+
  facet_wrap(~name, scales = 'free', nrow = 1)+
  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_4_true.csv") %>%
               rename(`β0` = V1, `β1` = V2, `β2` = V3, `β3` = V4, `β4` = V5, `ξ` = V6) %>%
               pivot_longer(everything()),
             aes(xintercept = value))+
  facet_wrap(~name, scales = 'free', nrow = 1)+
    labs(y= "Density", x = "Parameter value")+
    theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1),
                            axis.title.x=element_blank()),
read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
         col_names = c('bts', 'β0', 'β1', 'β2', 'β3', 'β4', 'ξ')) %>%
  pivot_longer(-bts) %>%
  ggplot()+
  geom_density(aes(value))+
  facet_wrap(~name, scales = 'free', nrow = 1)+
  geom_vline(data = readr::read_csv("output/gpd_model_fits/model_4_true.csv") %>%
               rename(`β0` = V1, `β1` = V2, `β2` = V3, `β3` = V4, `β4` = V5, `ξ` = V6) %>%
               pivot_longer(everything()),
             aes(xintercept = value))+
  facet_wrap(~name, scales = 'free', nrow = 1)+
  labs(y= "Density", x = "Parameter value")+
  theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)), heights = c(0.87,1))




bias_correction_mod_0 = gridExtra::grid.arrange(read_csv("output/gpd_model_fits/model_1_uncorrected.csv",
                                 col_names = c('bts', 'β0', 'β1', 'ξ')) %>%
                          pivot_longer(-bts) %>%
                          ggplot()+
                          geom_density(aes(value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          geom_vline(data = readr::read_csv("output/gpd_model_fits/model_1_true.csv") %>%
                                       rename(`β0` = V1, `β1` = V2, `ξ` = V3) %>%
                                       pivot_longer(everything()),
                                     aes(xintercept = value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          labs(y= "Density", x = "")+
                          theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1),
                                                  axis.title.x=element_blank()),
                        read_csv("output/gpd_model_fits/model_1_pars_corrected.csv",
                                 col_names = c('bts', 'β0', 'β1', 'ξ')) %>%
                          pivot_longer(-bts) %>%
                          ggplot()+
                          geom_density(aes(value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          geom_vline(data = readr::read_csv("output/gpd_model_fits/model_1_true.csv") %>%
                                       rename(`β0` = V1, `β1` = V2,  `ξ` = V3) %>%
                                       pivot_longer(everything()),
                                     aes(xintercept = value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          labs(y= "Density", x = "Parameter value")+
                          theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)), heights = c(0.87,1))




bias_correction_mod_1 = gridExtra::grid.arrange(read_csv("output/gpd_model_fits/model_2_uncorrected.csv",
                                 col_names = c('bts', 'β0', 'β1', 'β2', 'ξ')) %>%
                          pivot_longer(-bts) %>%
                          ggplot()+
                          geom_density(aes(value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          geom_vline(data = readr::read_csv("output/gpd_model_fits/model_2_true.csv") %>%
                                       rename(`β0` = V1, `β1` = V2, `β2` = V3,  `ξ` = V4) %>%
                                       pivot_longer(everything()),
                                     aes(xintercept = value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          labs(y= "Density", x = "")+
                          theme_minimal(10)+
                            theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1),
                                  axis.title.x=element_blank()),
                        read_csv("output/gpd_model_fits/model_2_pars_corrected.csv",
                                 col_names = c('bts', 'β0', 'β1', 'β2', 'ξ')) %>%
                          pivot_longer(-bts) %>%
                          ggplot()+
                          geom_density(aes(value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          geom_vline(data = readr::read_csv("output/gpd_model_fits/model_2_true.csv") %>%
                                       rename(`β0` = V1, `β1` = V2, `β2` = V3,  `ξ` = V4) %>%
                                       pivot_longer(everything()),
                                     aes(xintercept = value))+
                          facet_wrap(~name, scales = 'free', nrow = 1)+
                          labs(y= "Density", x = "Parameter value")+
                          theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)), heights = c(0.87,1))





ggsave(plot = bias_correction_mod_2, filename = "output/figs/bias_correction_mod_2.png", height = 3.75, width = 8)
ggsave(plot = bias_correction_mod_1, filename = "output/figs/bias_correction_mod_1.png", height = 3.75, width = 5.33333)
ggsave(plot = bias_correction_mod_0, filename = "output/figs/bias_correction_mod_0.png", height = 3.75, width = 4)





