rm(list=ls())
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
source('src/models/marginal_models/gpd_models.R')
library(tidyverse)
library(rnaturalearth)

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


# read in corrected fits
bootstrap_model_fits = read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
                                col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'xi'))





# --- get scale prediction for each fit
bts_res = c()
marg_mod = "mod_4"
num_quantiles = 25

# clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv")) %>%
#   dplyr::select(-c(loess_glob_temp_anom, quantile, value)) %>% filter(year %in% c(1931, 1942, 2020, 2022))
# 
# grid_pred = clim_date_w_quantile_mod %>%
#   dplyr::select(id, Long, Lat, year, threshold_9, thresh_exceedance_9, loess_temp_anom) %>%
#   left_join(read_csv('data/processed/clim_scale_grid.csv')) %>%
#   #left_join(covars %>% dplyr::select(year, loess_temp_anom) %>% unique)%>%
#   left_join(read_csv("data/processed/sites_clim_sea_dist.csv"))
# 
# for(i in seq(500)){
#   this_fit_pars = bootstrap_model_fits[i,] %>% unlist %>% .[-1] %>% as.numeric()
#   this_scale_fit = my_predict_4b(this_fit_pars, grid_pred$scale_9, grid_pred$loess_temp_anom, grid_pred$dist_sea)
# 
#   grid_pred %>%
#     dplyr::select(-c(threshold_9, thresh_exceedance_9, loess_temp_anom, Long.projected, Lat.projected, scale_9, dist_sea)) %>%
#     mutate(bts = i, scale = this_scale_fit$scale) %>%
#     write_csv("output/grid_for_sigma_for_plot.csv", append = T)
# }
# 
# true_change = grid_pred
# 
# true_change$scale = my_predict_4b((read_csv("output/gpd_model_fits/model_4_true.csv") %>% unlist %>% as.numeric),
#                                   true_change$scale_9, true_change$loess_temp_anom, true_change$dist_sea)$scale
# true_change$shape = my_predict_4b((read_csv("output/gpd_model_fits/model_4_true.csv") %>% unlist %>% as.numeric),
#                                   true_change$scale_9, true_change$loess_temp_anom, true_change$dist_sea)$shape
# 
# 
# 
# 
# 
# 
# bts_res = read_csv("output/grid_for_sigma_for_plot.csv",
#                    col_names = c("id", "Long", "Lat", "year", "bts", "scale"))
# bts_res %>% filter(year %in% c(1942, 2020)) %>%
#   filter(id == 123) %>%
#   ggplot()+
#   geom_line(data = true_change %>% filter(id == 123) %>% filter(year %in% c(1942, 2020)), aes(year, scale, col = 'true'))+
#   geom_line(data = rbind(bts_res %>%
#                            filter(year %in% c(2020)) %>%
#                            filter(id == 123) %>%
#                            summarise(year,
#                                      upper = quantile(scale, 0.975)) %>%
#                            unique(),
#                          bts_res %>%
#                            filter(year %in% c(1942)) %>%
#                            filter(id == 123) %>%
#                            summarise(year,
#                                      upper = quantile(scale, 0.975)) %>%
#                            unique()), aes(year, upper), linetype = 'dashed')+
#   geom_line(data = rbind(bts_res %>%
#                            filter(year %in% c(2020)) %>%
#                            filter(id == 123) %>%
#                            summarise(year,
#                                      lower = quantile(scale, 0.025)) %>%
#                            unique(),
#                          bts_res %>%
#                            filter(year %in% c(1942)) %>%
#                            filter(id == 123) %>%
#                            summarise(year,
#                                      lower = quantile(scale, 0.025)) %>%
#                            unique()), aes(year, lower), linetype = 'dashed')+
#   theme_minimal()+
#   labs("year")+
#   scale_x_continuous(breaks = c(1942, 2022))
# 





# # --- get scale prediction for each fit
bts_res = c()
marg_mod = "mod_4"
num_quantiles = 25

for(i in seq(500)){
  print(i)
  clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_", marg_mod, "_num_quantiles_",num_quantiles,"_bts_",i,".csv")) %>%
    dplyr::select(-c(loess_glob_temp_anom)) %>% filter(year %in% c(1931, 1942, 2020, 2022))

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
          dplyr::select(-c(threshold_9, thresh_exceedance_9, loess_temp_anom, Long.projected, Lat.projected, scale_9, dist_sea, rl_qnt)) %>%
          mutate(bts = i, rl = this_rl_fit)) %>%
    write_csv("output/scale_and_rl_for_plotting_mod_4.csv", append = T)
}



num_quantiles = 25
clim_date_w_quantile_mod = readRDS(paste0("output/quant_models_clim_num_quantiles_",num_quantiles,".csv")) %>%
  dplyr::select(-c(loess_glob_temp_anom, quantile, value)) %>% filter(year %in% c(1931, 1942, 2020, 2022))

grid_pred = clim_date_w_quantile_mod %>%
  dplyr::select(id, Long, Lat, year, threshold_9, thresh_exceedance_9, loess_temp_anom) %>%
  left_join(read_csv('data/processed/clim_scale_grid.csv')) %>%
  left_join(read_csv("data/processed/sites_clim_sea_dist.csv"))


true_change = grid_pred
true_change$rl_qnt = 1 - (1/grid_pred$thresh_exceedance_9)/(100*92)
true_change$rl = rl_mod_4b((read_csv("output/gpd_model_fits/model_4_true.csv") %>% unlist %>% as.numeric),
                               true_change$rl_qnt, true_change$threshold_9, true_change$scale_9, true_change$loess_temp_anom, true_change$dist_sea)

bts_res = read_csv("output/scale_and_rl_for_plotting_mod_4.csv",
                   col_names = c("id", "Long", "Lat", "year", "bts", "rl"))

rl_plot = gridExtra::grid.arrange(true_change %>%
                          filter(year == 2020) %>%
                          dplyr::select(Long, Lat, id, rl) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = rl))+
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(paste('°C')))+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                        bts_res %>%
                          filter(year  == 1942) %>%
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
                          labs(col = expression(paste(nabla, '°C')))+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                        true_change %>%
                          filter(year == 1942) %>%
                          rename(rl_100_1942 = rl) %>%
                          dplyr::select(Long, Lat, id, rl_100_1942) %>%
                          left_join(true_change %>%
                                      filter(year == 2020) %>%
                                      rename(rl_100_2020 = rl) %>%
                                      dplyr::select(Long, Lat, id, rl_100_2020)) %>%
                          mutate(rl_100_diff = rl_100_2020 - rl_100_1942) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = rl_100_diff))+
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(paste(nabla, '°C')))+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                plot.margin = unit(c(0,-0.1,0,-0.1),"cm")),
                        bts_res %>%
                          filter(year %in% c(1942)) %>%
                          group_by(id) %>%
                          summarise(Long, Lat, id,lower_1940  = quantile(rl, 0.975)) %>% unique() %>%
                          left_join(bts_res %>%
                                      filter(year %in% c(2020)) %>%
                                      group_by(id) %>%
                                      summarise(Long, Lat, id,
                                                lower_2020  = quantile(rl, 0.975)) %>% unique()) %>%
                          mutate(lower_change = lower_2020 - lower_1940) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = lower_change))+
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal, breaks = c(2.5,3.5))+
                          #scale_color_gradientn( colours=my_pal, breaks = c(2.5,3.5,4.5))+
                          labs(col = expression(paste(nabla, '°C')))+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                plot.margin = unit(c(0,-0.1,0,-0.1),"cm")), nrow = 1, widths = c(0.01,0.01,0.01016,0.01))

ggsave(rl_plot, filename ="output/figs/marginal_rl.pdf", height = 2, width = 9)


# ---- plot for presentation

rl_plot = gridExtra::grid.arrange(true_change %>%
                                    filter(year == 2020) %>%
                                    dplyr::select(Long, Lat, id, rl) %>%
                                    ggplot()+
                                    geom_point(aes(Long, Lat, col = rl))+
                                    geom_sf(data = ireland_sf, col = 'white', alpha = 0)+
                                    scale_color_gradientn( colours=my_pal)+
                                    labs(col = expression(paste('°C')),
                                         x = "", y = "")+
                                    scale_x_continuous(breaks = c(-10, -8, -6))+
                                    theme_minimal()+
                                    theme(axis.text.x  = element_blank(),
                                          axis.text.y  = element_blank(),
                                          legend.text = element_text(color = 'white'),
                                          legend.title  = element_text(color = 'white'),
                                          plot.background = element_rect(fill = "transparent", color = NA), 
                                          legend.background = element_rect(fill = "transparent"),
                                          legend.box.background = element_rect(fill = "transparent"),
                                          axis.text = element_text(color = 'white'),
                                          axis.title = element_text(color = 'white'),
                                          panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), 
                                          panel.border = element_blank(), 
                                          axis.line = element_blank()),
                                  bts_res %>%
                                    filter(year  == 1942) %>%
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
                                    geom_sf(data = ireland_sf, col = 'white', alpha = 0)+
                                    scale_color_gradientn( colours=my_pal)+
                                    labs(col = expression(paste(nabla, '°C')),
                                         x = "", y = "")+
                                    scale_x_continuous(breaks = c(-10, -8, -6))+
                                    theme_minimal()+
                                    theme(axis.text.x  = element_blank(),
                                          axis.text.y  = element_blank(),
                                          legend.text = element_text(color = 'white'),
                                          legend.title  = element_text(color = 'white'),
                                          plot.background = element_rect(fill = "transparent", color = NA), 
                                          legend.background = element_rect(fill = "transparent"),
                                          legend.box.background = element_rect(fill = "transparent"),
                                          axis.text = element_text(color = 'white'),
                                          axis.title = element_text(color = 'white'),
                                          panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), 
                                          panel.border = element_blank(), 
                                          axis.line = element_blank()),
                                  true_change %>%
                                    filter(year == 1942) %>%
                                    rename(rl_100_1942 = rl) %>%
                                    dplyr::select(Long, Lat, id, rl_100_1942) %>%
                                    left_join(true_change %>%
                                                filter(year == 2020) %>%
                                                rename(rl_100_2020 = rl) %>%
                                                dplyr::select(Long, Lat, id, rl_100_2020)) %>%
                                    mutate(rl_100_diff = rl_100_2020 - rl_100_1942) %>%
                                    ggplot()+
                                    geom_point(aes(Long, Lat, col = rl_100_diff))+
                                    geom_sf(data = ireland_sf, col = 'white', alpha = 0)+
                                    scale_color_gradientn( colours=my_pal)+
                                    labs(col = expression(paste(nabla, '°C')),
                                         x = "", y = "")+
                                    scale_x_continuous(breaks = c(-10, -8, -6))+
                                    theme_minimal()+
                                    theme(axis.text.x  = element_blank(),
                                          axis.text.y  = element_blank(),
                                          legend.text = element_text(color = 'white'),
                                          legend.title  = element_text(color = 'white'),
                                          plot.background = element_rect(fill = "transparent", color = NA), 
                                          legend.background = element_rect(fill = "transparent"),
                                          legend.box.background = element_rect(fill = "transparent"),
                                          axis.text = element_text(color = 'white'),
                                          axis.title = element_text(color = 'white'),
                                          panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), 
                                          panel.border = element_blank(), 
                                          axis.line = element_blank()),
                                  bts_res %>%
                                    filter(year %in% c(1942)) %>%
                                    group_by(id) %>%
                                    summarise(Long, Lat, id,lower_1940  = quantile(rl, 0.975)) %>% unique() %>%
                                    left_join(bts_res %>%
                                                filter(year %in% c(2020)) %>%
                                                group_by(id) %>%
                                                summarise(Long, Lat, id,
                                                          lower_2020  = quantile(rl, 0.975)) %>% unique()) %>%
                                    mutate(lower_change = lower_2020 - lower_1940) %>%
                                    ggplot()+
                                    geom_point(aes(Long, Lat, col = lower_change))+
                                    geom_sf(data = ireland_sf, col = 'white', alpha = 0)+
                                    scale_color_gradientn( colours=my_pal, breaks = c(2.5,3.5))+
                                    #scale_color_gradientn( colours=my_pal, breaks = c(2.5,3.5,4.5))+
                                    labs(col = expression(paste(nabla, '°C')),
                                         x = "", y = "")+
                                    scale_x_continuous(breaks = c(-10, -8, -6))+
                                    theme_minimal()+
                                    theme(axis.text.x  = element_blank(),
                                          axis.text.y  = element_blank(),
                                          legend.text = element_text(color = 'white'),
                                          legend.title  = element_text(color = 'white'),
                                          plot.background = element_rect(fill = "transparent", color = NA), 
                                          legend.background = element_rect(fill = "transparent"),
                                          legend.box.background = element_rect(fill = "transparent"),
                                          axis.text = element_text(color = 'white'),
                                          axis.title = element_text(color = 'white'),
                                          panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), 
                                          panel.border = element_blank(), 
                                          axis.line = element_blank()), nrow = 1, widths = c(0.01,0.01,0.01,0.01))


ggsave( rl_plot, filename ="marginal_rl.png", bg = 'transparent', height = 2.5, width = 10)

