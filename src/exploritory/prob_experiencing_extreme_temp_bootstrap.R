gc()
rm(list = ls())
setwd("~/Extreme-Irish-Summer-Temperatures/")
source('src/models/marginal_models/gpd_models.R')
library(tidyverse)


# ------ MODEL 4
# for(bts in seq(1, 500)){
# 
#   if(file.exists(paste0("output/simulations/simulations_on_obs_grid/bootstraps/mod_4_bootstrap_",bts,"_run_1"))){
#     obs_sites = read_csv("data/processed/obs_data.csv") %>%
#       dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
#       unique() %>%
#       left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
#                   dplyr::select(Station, dist_sea))
#     
#     obs_grid = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_mod_4_num_quantiles_25_bts_", bts, ".csv"),
#                         col_names = c('bts', 'Station', 'year', 'thresh_exceedance_9'))  %>%
#       left_join(obs_sites) %>%
#       left_join(read_csv("data/processed/obs_data.csv") %>%
#                   dplyr::select(year, loess_temp_anom) %>% unique)%>%
#       filter(year %in% c(1942, 2020))
#     
#     
#     
#     gpd_pars = read_csv("output/gpd_model_fits/model_4_pars_corrected.csv",
#                         col_names = c('bts', 'b0', 'b1', 'b2', 'b3', 'b4', 'b5')) %>% 
#       filter(bts == paste0("num_quantiles_50_bts_", {{bts}})) %>%
#       unlist() %>% .[2:length(.)] %>% as.numeric
#     
#     
#     obs_grid$scale = my_predict_4b(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$scale
#     obs_grid$shape = my_predict_4b(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$shape
#     
#     
#     my_simulations_extremes = c()
#     for(i in seq(1, 100)){
#       my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/bootstraps/mod_4_bootstrap_",bts,"_run_",i)))
#     }
#     
#     sim_costs = map(my_simulations_extremes, mean) %>% unlist
#     
#     
#     
#     q = 0.5
#     res_2020 = c()
#     res_1942 = c()
#     
#     for(temp_i_want in seq(20, 36, by = 1)){
#       
#       obs_grid$frechet = -1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(temp_i_want - obs_grid$threshold_9,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
#       frechet_2020 = obs_grid %>% filter(year == 2020) %>% pull(frechet)
#       frechet_1942 = obs_grid %>% filter(year == 1942) %>% pull(frechet)
#       frechet_1942[frechet_1942 == -Inf] = Inf
#       frechet_2020[frechet_2020 == -Inf] = Inf
#       
#       above_1942 = my_simulations_extremes %>%
#         map(.f = ~{
#           sum(.x >frechet_1942)>=1
#         }) %>% unlist
#       
#       above_2020 = my_simulations_extremes %>%
#         map(.f = ~{
#           sum(.x > frechet_2020)>=1
#         }) %>% unlist
#       
#       thresh = quantile(sim_costs, q)
#       
#       prob_observing_conditional_1942 = sum((sim_costs > thresh) & above_1942)/sum((sim_costs > thresh))
#       prob_observing_conditional_2020 = sum((sim_costs > thresh) & above_2020)/sum((sim_costs > thresh))
#       
#       res_1942 = c(res_1942, prob_observing_conditional_1942)
#       res_2020 = c(res_2020, prob_observing_conditional_2020)
#     }
#     
#     
#     
#     res_2020 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*res_2020
#     res_1942 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*res_1942
#     
#     tibble(bts,
#            temps = seq(20, 36, by = 1), 
#            prob_1942 = res_1942,
#            prob_2020 = res_2020) %>%
#       write_csv("output/bootstrap_prob_observing_T_anywhere.csv", append = T)
#     
#     
#   }
# }
# 
# 
# rl_plot = read_csv("output/bootstrap_prob_observing_T_anywhere.csv",
#          col_names = c('bts', 'temp', 'p1940', 'p2020')) %>%
#   group_by(temp) %>%
#   summarise(upper_2020 = quantile(p2020, 0.9),
#             lower_2020 = quantile(p2020, 0.1),
#             upper_1940 = quantile(p1940, 0.9),
#             lower_1940 = quantile(p1940, 0.1))%>%
#   ggplot()+
#   geom_ribbon(aes(temp,ymin = 1/(92*lower_1940), ymax = 1/(92*upper_1940),  fill = '1940'), alpha = 0.25)+
#   geom_ribbon(aes(temp,ymin = 1/(92*lower_2020), ymax = 1/(92*upper_2020),  fill = '2020'), alpha = 0.25)+
#   geom_line(data = read_csv("output/prob_observing_T_anywhere_clim.csv"), aes(temps, 1/(92*prob_1942), col = '1940' ))+
#   geom_line(data = read_csv("output/prob_observing_T_anywhere_clim.csv"), aes(temps, 1/(92*prob_2020), col = '2020' ), linetype = 'dashed')+
#   theme_minimal()+
#     labs(x = "Temperature",
#          y = "Return period\n(Years)",
#          fill = 'Year', col = "Year")+
#     theme_minimal(12)+
#     theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
#           legend.position = 'none')+
#     scale_x_continuous(limits = c(27, 34),
#                        breaks = c(25, 28,  30,  32, 34),
#                        label = paste0(c(26, 28,  30,  32,34),"°C"))+
#   scale_y_log10(limits = c(0.075, 1000),breaks = c(0.1, 1, 10,  100,  1000),labels = c(0.1, 1, 10,  100,  1000))+ coord_flip()
#  
# ggsave(rl_plot, filename = 'output/figs/spatial_rl.pdf', width = 5, height = 3)








# -------- MODEL 0

# 
# for(bts in seq(301, 500)){
#   
#   
#   
#   if(file.exists(paste0("output/simulations/simulations_on_obs_grid/bootstraps/mod_1_bootstrap_",bts,"_run_1"))){
#     
#     
#     obs_sites = read_csv("data/processed/obs_data.csv") %>%
#       dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
#       unique() %>%
#       left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
#                   dplyr::select(Station, dist_sea))
#     
#     obs_grid = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_mod_1_num_quantiles_15_bts_", bts, ".csv")) %>%
#       rename(bts = X1, Station = X2, year = X3, thresh_exceedance_9 = X4) %>%
#       left_join(obs_sites) %>%
#       left_join(read_csv("data/processed/obs_data.csv") %>%
#                   dplyr::select(year, loess_temp_anom) %>% unique)%>%
#       filter(year %in% c(1942, 2020))
#     
#     
#     gpd_pars = read_csv("output/gpd_model_fits/model_1_pars_corrected.csv",
#                         col_names = c('bts', 'b0', 'b1', 'b2')) %>% 
#       filter(bts == paste0("num_quantiles_50_bts_", {{bts}})) %>%
#       unlist() %>% .[2:length(.)] %>% as.numeric
#     
#     
#     obs_grid$scale = my_predict_1(unlist(gpd_pars), obs_grid$scale_9)$scale
#     obs_grid$shape = my_predict_1(unlist(gpd_pars), obs_grid$scale_9)$shape
#     
#     
#     my_simulations_extremes = c()
#     for(i in seq(1, 100)){
#       my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/bootstraps/mod_1_bootstrap_",bts,"_run_",i)))
#     }
#     
#     sim_costs = map(my_simulations_extremes, mean) %>% unlist
# 
#     q = 0.5
#     res_2020 = c()
#     res_1942 = c()
#     
#     for(temp_i_want in seq(20, 36, by = 1)){
#       
#       obs_grid$frechet = -1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(temp_i_want - obs_grid$threshold_9,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
#       frechet_2020 = obs_grid %>% filter(year == 2020) %>% pull(frechet)
#       frechet_1942 = obs_grid %>% filter(year == 1942) %>% pull(frechet)
#       frechet_1942[frechet_1942 == -Inf] = Inf
#       frechet_2020[frechet_2020 == -Inf] = Inf
#       
#       above_1942 = my_simulations_extremes %>%
#         map(.f = ~{
#           sum(.x >frechet_1942)>=1
#         }) %>% unlist
#       
#       above_2020 = my_simulations_extremes %>%
#         map(.f = ~{
#           sum(.x > frechet_2020)>=1
#         }) %>% unlist
#       
#       thresh = quantile(sim_costs, q)
#       
#       prob_observing_conditional_1942 = sum((sim_costs > thresh) & above_1942)/sum((sim_costs > thresh))
#       prob_observing_conditional_2020 = sum((sim_costs > thresh) & above_2020)/sum((sim_costs > thresh))
#       
#       res_1942 = c(res_1942, prob_observing_conditional_1942)
#       res_2020 = c(res_2020, prob_observing_conditional_2020)
#     }
# 
#     res_2020 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*res_2020
#     res_1942 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*res_1942
#     
#     tibble(bts,
#            temps = seq(20, 36, by = 1), 
#            prob_1942 = res_1942,
#            prob_2020 = res_2020) %>%
#       write_csv("output/bootstrap_prob_observing_T_anywhere_mod_1.csv", append = T)
#   }
# }
# 
# 
# rl_plot = read_csv("output/bootstrap_prob_observing_T_anywhere.csv",
#                    col_names = c('bts', 'temp', 'p1940', 'p2020')) %>%
#   group_by(temp) %>%
#   summarise(upper_2020 = quantile(p2020, 0.9),
#             lower_2020 = quantile(p2020, 0.1),
#             upper_1940 = quantile(p1940, 0.9),
#             lower_1940 = quantile(p1940, 0.1))%>%
#   ggplot()+
#   geom_ribbon(aes(temp,ymin = 1/(92*lower_1940), ymax = 1/(92*upper_1940),  fill = '1940'), alpha = 0.25)+
#   geom_ribbon(aes(temp,ymin = 1/(92*lower_2020), ymax = 1/(92*upper_2020),  fill = '2020'), alpha = 0.25)+
#   geom_line(data = read_csv("output/prob_observing_T_anywhere_clim.csv"), aes(temps, 1/(92*prob_1942), col = '1940' ))+
#   geom_line(data = read_csv("output/prob_observing_T_anywhere_clim.csv"), aes(temps, 1/(92*prob_2020), col = '2020' ), linetype = 'dashed')+
#   theme_minimal()+
#   labs(x = "Temperature",
#        y = "Return period\n(Years)",
#        fill = 'Year', col = "Year")+
#   theme_minimal(12)+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
#         legend.position = 'none')+
#   scale_x_continuous(limits = c(27, 34),
#                      breaks = c(25, 28,  30,  32, 34),
#                      label = paste0(c(26, 28,  30,  32,34),"°C"))+
#   scale_y_log10(limits = c(0.075, 1000),breaks = c(0.1, 1, 10,  100,  1000),labels = c(0.1, 1, 10,  100,  1000))+ coord_flip()
# 
# ggsave(rl_plot, filename = 'output/figs/spatial_rl.pdf', width = 5, height = 3)










# -------- MODEL 1


for(bts in seq(301, 500)){
  
  
  
  if(file.exists(paste0("output/simulations/simulations_on_obs_grid/bootstraps/mod_2_bootstrap_",bts,"_run_1"))){
    
    
    obs_sites = read_csv("data/processed/obs_data.csv") %>%
      dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
      unique() %>%
      left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
                  dplyr::select(Station, dist_sea))
    
    obs_grid = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_mod_2_num_quantiles_15_bts_", bts, ".csv")) %>%
      rename(bts = X1, Station = X2, year = X3, thresh_exceedance_9 = X4) %>%
      left_join(obs_sites) %>%
      left_join(read_csv("data/processed/obs_data.csv") %>%
                  dplyr::select(year, loess_temp_anom) %>% unique)%>%
      filter(year %in% c(1942, 2020))
    
    
    gpd_pars = read_csv("output/gpd_model_fits/model_2_pars_corrected.csv",
                        col_names = c('bts', 'b0', 'b1', 'b2', 'b3')) %>% 
      filter(bts == paste0("num_quantiles_50_bts_", {{bts}})) %>%
      unlist() %>% .[2:length(.)] %>% as.numeric
    
    
    obs_grid$scale = my_predict_2(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom)$scale
    obs_grid$shape = my_predict_2(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom)$shape
    
    
    my_simulations_extremes = c()
    for(i in seq(1, 100)){
      my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/bootstraps/mod_2_bootstrap_",bts,"_run_",i)))
    }
    
    sim_costs = map(my_simulations_extremes, mean) %>% unlist
    
    q = 0.5
    res_2020 = c()
    res_1942 = c()
    
    for(temp_i_want in seq(20, 36, by = 1)){
      
      obs_grid$frechet = -1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(temp_i_want - obs_grid$threshold_9,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
      frechet_2020 = obs_grid %>% filter(year == 2020) %>% pull(frechet)
      frechet_1942 = obs_grid %>% filter(year == 1942) %>% pull(frechet)
      frechet_1942[frechet_1942 == -Inf] = Inf
      frechet_2020[frechet_2020 == -Inf] = Inf
      
      above_1942 = my_simulations_extremes %>%
        map(.f = ~{
          sum(.x >frechet_1942)>=1
        }) %>% unlist
      
      above_2020 = my_simulations_extremes %>%
        map(.f = ~{
          sum(.x > frechet_2020)>=1
        }) %>% unlist
      
      thresh = quantile(sim_costs, q)
      
      prob_observing_conditional_1942 = sum((sim_costs > thresh) & above_1942)/sum((sim_costs > thresh))
      prob_observing_conditional_2020 = sum((sim_costs > thresh) & above_2020)/sum((sim_costs > thresh))
      
      res_1942 = c(res_1942, prob_observing_conditional_1942)
      res_2020 = c(res_2020, prob_observing_conditional_2020)
    }
    
    res_2020 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*res_2020
    res_1942 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*res_1942
    
    tibble(bts,
           temps = seq(20, 36, by = 1), 
           prob_1942 = res_1942,
           prob_2020 = res_2020) %>%
      write_csv("output/bootstrap_prob_observing_T_anywhere_mod_2.csv", append = T)
  }
}




rl_plot_1 = read_csv("output/bootstrap_prob_observing_T_anywhere_mod_1.csv",
                   col_names = c('bts', 'temp', 'p1940', 'p2020')) %>%
  group_by(temp) %>%
  summarise(upper_2020 = quantile(p2020, 0.9),
            lower_2020 = quantile(p2020, 0.1),
            upper_1940 = quantile(p1940, 0.9),
            lower_1940 = quantile(p1940, 0.1))%>%
  ggplot()+
  geom_ribbon(aes(temp,ymin = 1/(92*lower_1940), ymax = 1/(92*upper_1940),  fill = '1940'), alpha = 0.25)+
  geom_ribbon(aes(temp,ymin = 1/(92*lower_2020), ymax = 1/(92*upper_2020),  fill = '2020'), alpha = 0.25)+
  geom_line(data = read_csv("output/prob_observing_T_anywhere_clim_mod_1.csv"), aes(temps, 1/(92*prob_1942), col = '1940' ))+
  geom_line(data = read_csv("output/prob_observing_T_anywhere_clim_mod_1.csv"), aes(temps, 1/(92*prob_2020), col = '2020' ), linetype = 'dashed')+
  theme_minimal()+
  labs(x = "Temperature",
       y = "Return period (Years)",
       fill = 'Year', col = "Year")+
  theme_minimal(12)+
  theme(legend.position = 'none')+
  scale_x_continuous(limits = c(27, 33),
                     breaks = c(25, 28,  30,  32, 34),
                     label = paste0(c(26, 28,  30,  32,34),"°C"))+
  scale_y_log10(limits = c(0.075, 1000),breaks = c(0.1, 1, 10,  100,  1000),labels = c(0.1, 1, 10,  100,  1000))+ coord_flip()



rl_plot_2 = read_csv("output/bootstrap_prob_observing_T_anywhere_mod_2.csv",
                     col_names = c('bts', 'temp', 'p1940', 'p2020')) %>%
  group_by(temp) %>%
  summarise(upper_2020 = quantile(p2020, 0.9),
            lower_2020 = quantile(p2020, 0.1),
            upper_1940 = quantile(p1940, 0.9),
            lower_1940 = quantile(p1940, 0.1))%>%
  ggplot()+
  geom_ribbon(aes(temp,ymin = 1/(92*lower_1940), ymax = 1/(92*upper_1940),  fill = '1940'), alpha = 0.25)+
  geom_ribbon(aes(temp,ymin = 1/(92*lower_2020), ymax = 1/(92*upper_2020),  fill = '2020'), alpha = 0.25)+
  geom_line(data = read_csv("output/prob_observing_T_anywhere_clim_mod_2.csv"), aes(temps, 1/(92*prob_1942), col = '1940' ))+
  geom_line(data = read_csv("output/prob_observing_T_anywhere_clim_mod_2.csv"), aes(temps, 1/(92*prob_2020), col = '2020' ), linetype = 'dashed')+
  theme_minimal()+
  labs(x = "Temperature",
       y = "Return period (Years)",
       fill = 'Year', col = "Year")+
  theme_minimal(12)+
  theme(legend.position = 'none')+
  scale_x_continuous(limits = c(27, 33),
                     breaks = c(25, 28,  30,  32, 34),
                     label = paste0(c(26, 28,  30,  32,34),"°C"))+
  scale_y_log10(limits = c(0.075, 1000),breaks = c(0.1, 1, 10,  100,  1000),labels = c(0.1, 1, 10,  100,  1000))+ coord_flip()



rl_plot = gridExtra::grid.arrange(rl_plot_1, rl_plot_2, nrow = 1)



ggsave(rl_plot, filename = 'output/figs/spatial_rl_mod_1_and_2.pdf', width = 8, height = 3)







rl_plot_4 = read_csv("output/bootstrap_prob_observing_T_anywhere.csv",
                     col_names = c('bts', 'temp', 'p1940', 'p2020')) %>%
  group_by(temp) %>%
  summarise(upper_2020 = quantile(p2020, 0.9),
            lower_2020 = quantile(p2020, 0.1),
            upper_1940 = quantile(p1940, 0.9),
            lower_1940 = quantile(p1940, 0.1))%>%
  ggplot()+
  geom_ribbon(aes(temp,ymin = 1/(92*lower_1940), ymax = 1/(92*upper_1940),  fill = '1940'), alpha = 0.25)+
  geom_ribbon(aes(temp,ymin = 1/(92*lower_2020), ymax = 1/(92*upper_2020),  fill = '2020'), alpha = 0.25)+
  geom_line(data = read_csv("output/prob_observing_T_anywhere_clim_mod_4.csv"), aes(temps, 1/(92*prob_1942), col = '1940' ))+
  geom_line(data = read_csv("output/prob_observing_T_anywhere_clim_mod_4.csv"), aes(temps, 1/(92*prob_2020), col = '2020' ), linetype = 'dashed')+
  theme_minimal()+
  labs(x = "Temperature",
       y = "Return period (Years)",
       fill = 'Year', col = "Year")+
  theme_minimal(12)+
  theme(legend.position = 'none')+
  scale_x_continuous(limits = c(27, 33),
                     breaks = c(25, 28,  30,  32, 34),
                     label = paste0(c(26, 28,  30,  32,34),"°C"))+
  scale_y_log10(limits = c(0.075, 1000),breaks = c(0.1, 1, 10,  100,  1000),labels = c(0.1, 1, 10,  100,  1000))+ coord_flip()

