rm(list=ls())
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
source('src/models/marginal_models/gpd_models.R')
library(tidyverse)

# install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")
# install.packages('rnaturalearth')
library(rnaturalearth)

# ----- need 
# 100 year return level 2022
# difference between the 100-year return level estimated for the years 2022 and 1931
# Lower 95% CI for change in 100-year return level from 1931 to 2022
# Upper 95% CI for change in 100-year return level from 1931 to 2022

# map for plotting
ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object, change this to "sp" if needed, this will break the plot though
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]


# colours for plotting
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


# 
# clim_date_w_quantile_mod = readRDS("output/quant_models_clim_num_quantiles_100.csv") %>%
#   dplyr::select(-c(loess_glob_temp_anom, quantile, value)) %>% filter(year %in% c(1931, 1942, 2020, 2022))

# --- read this in seperately for each bootstrap

# # ----- plot thresh exceedance on a map?
# clim_date_w_quantile_mod %>%
#   dplyr::select(id, Long, Lat, year, threshold_9, thresh_exceedance_9) %>%
#   filter(year == 1931) %>% 
#   ggplot()+
#   geom_point(aes(Long, Lat, col = thresh_exceedance_9))+
#   ggplot2::scale_color_gradientn(colours = my_pal)+
#   geom_sf(data = ireland_sf, alpha = 0, col = 'black')
# 
# clim_date_w_quantile_mod %>%
#   dplyr::select(id, Long, Lat, year, threshold_9, thresh_exceedance_9) %>%
#   filter(year == 2022) %>%
#   ggplot()+
#   geom_point(aes(Long, Lat, col = thresh_exceedance_9))+
#   ggplot2::scale_color_gradientn(colours = my_pal)+
#   geom_sf(data = ireland_sf, alpha = 0, col = 'black')
# 


# read in corrected fits

bootstrap_model_fits = read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
         col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'xi'))

# --- get scale prediction for each fit
bts_res = c()
marg_mod = "mod_4"
num_quantiles  =25
for(i in seq(300)){
  clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",i,".csv")) %>%
    dplyr::select(-c(loess_glob_temp_anom, quantile, value)) %>% filter(year %in% c(1931, 1942, 2020, 2022))
  
  
  grid_pred = clim_date_w_quantile_mod %>%
    dplyr::select(id, Long, Lat, year, threshold_9, thresh_exceedance_9, loess_temp_anom) %>%
    left_join(read_csv('data/processed/clim_scale_grid.csv')) %>%
    #left_join(covars %>% dplyr::select(year, loess_temp_anom) %>% unique)%>%
    left_join(read_csv("data/processed/sites_clim_sea_dist.csv"))
  
  grid_pred$rl_qnt = 1 - (1/grid_pred$thresh_exceedance_9)/(100*92)
  
  this_fit_pars = bootstrap_model_fits[i,] %>% unlist %>% .[-1] %>% as.numeric()
  
  this_scale_fit = my_predict_4b(this_fit_pars, grid_pred$scale_9, grid_pred$loess_temp_anom, grid_pred$dist_sea)
  this_rl_fit = rl_mod_4b(this_fit_pars, grid_pred$rl_qnt, grid_pred$threshold_9, grid_pred$scale_9, grid_pred$loess_temp_anom, grid_pred$dist_sea)

  rbind(bts_res, grid_pred %>%
    dplyr::select(-c(threshold_9, thresh_exceedance_9, loess_temp_anom, Long.projected, Lat.projected, scale_8, scale_9, dist_sea, rl_qnt)) %>%
    mutate(bts = i, scale = this_scale_fit$scale, rl = this_rl_fit)) %>%
    write_csv("output/grid_thresh_exceedence_for_plot.csv", append = T)
}



bts_res = read_csv("output/grid_thresh_exceedence_for_plot.csv",
         col_names = c('id', 'Long', 'Lat', 'year', 'bts', 'scale', 'rl'))

read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
         col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'xi')) %>%
  filter(bts == 3)

read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
         col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'xi')) %>%
  pull(xi) %>% hist



bts_res %>%
  filter(year %in% c(1931, 2022)) %>%
  filter(id == 123) %>%
  ggplot()+
 # geom_line(aes(year, rl, group = bts), alpha = 0.25)+
  geom_line(data = true_change %>% filter(id == 123) %>% filter(year %in% c(1931, 2022)), aes(year, rl_100, col = 'true'))+
  geom_line(data = rbind(bts_res %>%
                    filter(year %in% c(2022)) %>%
                    filter(id == 123) %>%
                    summarise(year,
                              upper = quantile(rl, 0.975)) %>%
                    unique(),
                  bts_res %>%
                    filter(year %in% c(1931)) %>%
                    filter(id == 123) %>%
                    summarise(year,
                              upper = quantile(rl, 0.975)) %>%
                    unique()), aes(year, upper), linetype = 'dashed')+
  geom_line(data = rbind(bts_res %>%
                           filter(year %in% c(2022)) %>%
                           filter(id == 123) %>%
                           summarise(year,
                                     lower = quantile(rl, 0.025)) %>%
                           unique(),
                         bts_res %>%
                           filter(year %in% c(1931)) %>%
                           filter(id == 123) %>%
                           summarise(year,
                                     lower = quantile(rl, 0.025)) %>%
                           unique()), aes(year, lower), linetype = 'dashed')+
  theme_minimal()+
  labs("year")+
  scale_x_continuous(breaks = c(1931, 2020))
  







sbts_to_rem = bts_res %>%
  filter(is.na(rl)) %>%
  pull(bts) %>% unique

bts_res = bts_res %>% 
  filter(!(bts %in% bts_to_rem))

# # --- point estimate
clim_date_w_quantile_mod = readRDS("output/quant_models_clim_num_quantiles_100.csv") %>%
  dplyr::select(-c(loess_glob_temp_anom, quantile, value)) %>% filter(year %in% c(1931, 1942, 2020, 2022))
true_change = clim_date_w_quantile_mod %>%
  dplyr::select(id, Long, Lat, year, threshold_9, thresh_exceedance_9, loess_temp_anom) %>%
  left_join(read_csv('data/processed/clim_scale_grid.csv')) %>%
  #left_join(covars %>% dplyr::select(year, loess_temp_anom) %>% unique)%>%
  left_join(read_csv("data/processed/sites_clim_sea_dist.csv"))
true_change$rl_qnt = 1 - (1/true_change$thresh_exceedance_9)/(100*92)


true_change$scale = my_predict_4b((read_csv("output/gpd_model_fits/model_4_true.csv") %>% unlist %>% as.numeric),
                                  true_change$scale_9, true_change$loess_temp_anom, true_change$dist_sea)$scale
true_change$shape = my_predict_4b((read_csv("output/gpd_model_fits/model_4_true.csv") %>% unlist %>% as.numeric),
                                  true_change$scale_9, true_change$loess_temp_anom, true_change$dist_sea)$shape
true_change$rl_100 = rl_mod_4b((read_csv("output/gpd_model_fits/model_4_true.csv") %>% unlist %>% as.numeric),
                               true_change$rl_qnt, true_change$threshold_9, true_change$scale_9, true_change$loess_temp_anom, true_change$dist_sea)
# 
# 
# 




gridExtra::grid.arrange(true_change %>%
                          filter(year == 2022) %>%
                          dplyr::select(Long, Lat, id, rl_100) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = rl_100))+
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(paste('°C')),
                               x = "", y = "")+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank()),
                        bts_res %>%
                          filter(year  == 1931) %>%
                          group_by(id) %>%
                          summarise(Long, Lat, id, lower_1940  = quantile(rl, 0.025)) %>% unique() %>%
                          left_join(bts_res %>%
                                      filter(year  == 2020) %>%
                                      group_by(id) %>%
                                      summarise(Long, Lat, id,
                                                lower_2020  = quantile(rl, 0.025)) %>% unique()) %>%
                          mutate(lower_change = lower_2020 - lower_1940) %>%
                          unique() %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = lower_change))+
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(paste(nabla, '°C')),
                               x = "", y = "")+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank()),
                        true_change %>%
                          filter(year == 1931) %>%
                          rename(rl_100_1942 = rl_100) %>%
                          dplyr::select(Long, Lat, id, rl_100_1942) %>%
                          left_join(true_change %>%
                                      filter(year == 2020) %>%
                                      rename(rl_100_2020 = rl_100) %>%
                                      dplyr::select(Long, Lat, id, rl_100_2020)) %>%
                          mutate(rl_100_diff = rl_100_2020 - rl_100_1942) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = rl_100_diff))+
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(paste(nabla, '°C')),
                               x = "", y = "")+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank()),
                        bts_res %>%
                          filter(year %in% c(1931)) %>%
                          group_by(id) %>%
                          summarise(Long, Lat, id,lower_1940  = quantile(rl, 0.975)) %>% unique() %>%
                          left_join(bts_res %>%
                                      filter(year %in% c(2022)) %>%
                                      group_by(id) %>%
                                      summarise(Long, Lat, id,
                                                lower_2020  = quantile(rl, 0.975)) %>% unique()) %>%
                          mutate(lower_change = lower_2020 - lower_1940) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = lower_change))+
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(paste(nabla, '°C')),
                               x = "", y = "")+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank()), nrow = 1)



# ((bts_res %>%
#   filter(year == 2020) %>%
#   pull(rl)) - (bts_res %>%
#                  filter(year == 1942) %>%
#                  pull(rl))) %>% quantile(c(0.025, 0.975))
# 
# 
# 
# bts_res %>%
#   filter(year == 1942) %>%
#   pull(rl) %>% hist
# 
# bts_res %>%
#   filter(year == 2020) %>%
#   pull(rl) %>% hist
# 
# bts_res %>%
#   filter(year == 1942) %>%
#   rename(scale_1931 = scale,
#          rl_1931 = rl) %>%
#   dplyr::select(-year) %>%
#   left_join(bts_res %>%
#               filter(year == 2020)%>%
#               rename(scale_2022 = scale,
#                      rl_2022 = rl)%>%
#               dplyr::select(-year)) %>%
#   mutate(scale_diff = scale_2022 - scale_1931,
#          rl_diff = rl_2022 - rl_1931) %>%
#   group_by(id) %>%
#   summarise(upper_scale = quantile(scale_diff, 0.975),
#             upper_rl = quantile(rl_diff, 0.975),
#             lower_scale = quantile(scale_diff, 0.025),
#             lower_rl = quantile(rl_diff, 0.025)) %>% pull(upper_rl) %>% hist
# 
# 
# 
