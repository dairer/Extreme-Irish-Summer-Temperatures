gc()
rm(list = ls())
setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
source('src/models/marginal_models/gpd_models.R')
library(tidyverse)

# # ------ MODEL 0
# 
# 
# obs_sites = read_csv("data/processed/obs_data.csv") %>%
#   dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
#   unique() %>%
#   left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
#               dplyr::select(Station, dist_sea))
# 
# obs_grid = read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_25.csv"))  %>%
#   left_join(obs_sites) %>%
#   left_join(read_csv("data/processed/obs_data.csv") %>%
#               dplyr::select(year, loess_temp_anom) %>% unique)%>%
#   filter(year %in% c(1942, 2020))
# 
# gpd_pars = read_csv("output/gpd_model_fits/model_1_true.csv")
# obs_grid$scale = my_predict_1(unlist(gpd_pars), obs_grid$scale_9)$scale
# obs_grid$shape = my_predict_1(unlist(gpd_pars), obs_grid$scale_9)$shape
# 
# my_simulations_extremes = c()
# for(i in seq(1, 100)){
#   my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true/mod_1_run_",i)))
# }
# 
# sim_costs = map(my_simulations_extremes, mean) %>% unlist
# 
# q = 0.5
# res_2020 = c()
# res_1942 = c()
# 
# resu_1942 = c()
# resl_1942 = c()
# 
# resu_2020 = c()
# resl_2020 = c()
# 
# for(temp_i_want in seq(20, 36, by = 0.25)){
#   
#   obs_grid$frechet = -1/log(1-obs_grid$thresh_exceedance_9*(1-evd::pgpd(temp_i_want - obs_grid$threshold_9,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
#   frechet_2020 = obs_grid %>% filter(year == 2020) %>% pull(frechet)
#   frechet_1942 = obs_grid %>% filter(year == 1942) %>% pull(frechet)
#   frechet_1942[frechet_1942 == -Inf] = Inf
#   frechet_2020[frechet_2020 == -Inf] = Inf
#   
#   above_1942 = my_simulations_extremes %>%
#     map(.f = ~{
#       sum(.x >frechet_1942)>=1
#     }) %>% unlist
#   
#   above_2020 = my_simulations_extremes %>%
#     map(.f = ~{
#       sum(.x > frechet_2020)>=1
#     }) %>% unlist
#   
#   thresh = quantile(sim_costs, q)
#   print(paste0("Cost threshold = ", signif(thresh,3)))
#   
#   prob_observing_conditional_1942 = sum((sim_costs > thresh) & above_1942)/sum((sim_costs > thresh))
#   print(paste0("p = ", signif(prob_observing_conditional_1942,3)))
#   res_1942 = c(res_1942, prob_observing_conditional_1942)
#   
#   cis = qbinom(c(0.025, 0.975), size = sum(sim_costs > thresh), prob = prob_observing_conditional_1942)/sum(sim_costs > thresh)
#   resu_1942 = c(resu_1942, cis[2])
#   resl_1942 = c(resl_1942, cis[1])
#   
#   prob_observing_conditional_2020 = sum((sim_costs > thresh) & above_2020)/sum((sim_costs > thresh))
#   print(paste0("p = ", signif(prob_observing_conditional_2020,3)))
#   res_2020 = c(res_2020, prob_observing_conditional_2020)
#   
#   cis = qbinom(c(0.025, 0.975), size = sum(sim_costs > thresh), prob = prob_observing_conditional_2020)/sum(sim_costs > thresh)
#   resu_2020 = c(resu_2020, cis[2])
#   resl_2020 = c(resl_2020, cis[1])
# }
# 
# 
# res_2020 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*res_2020
# res_1942 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*res_1942
# 
# resu_1942 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*resu_1942
# resl_1942 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*resl_1942
# 
# resu_2020 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*resu_2020
# resl_2020 = readRDS("data/processed/prob_event_greater_1_mod_1")*q*resl_2020
# 
# 
# tibble(temps = seq(20, 36, by = 0.25), 
#        prob_1942 = res_1942,
#        prob_2020 = res_2020,
#        prob_1942_u = resu_1942,
#        prob_1942_l = resl_1942,
#        prob_2020_u = resu_2020,
#        prob_2020_l = resl_2020) %>%
#   write_csv("output/prob_observing_T_anywhere_clim_mod_1.csv")
# 





# ------ MODEL 1
obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
  unique() %>%
  left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>% 
              dplyr::select(Station, dist_sea))

obs_grid = read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_25.csv"))  %>%
  left_join(obs_sites) %>%
  left_join(read_csv("data/processed/obs_data.csv") %>%
              dplyr::select(year, loess_temp_anom) %>% unique)%>%
  filter(year %in% c(1942, 2020))

gpd_pars = read_csv("output/gpd_model_fits/model_2_true.csv")
obs_grid$scale = my_predict_2(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom)$scale
obs_grid$shape = my_predict_2(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom)$shape

my_simulations_extremes = c()
for(i in seq(1, 100)){
  my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true/mod_2_run_",i)))
}

sim_costs = map(my_simulations_extremes, mean) %>% unlist

q = 0.5
res_2020 = c()
res_1942 = c()

resu_1942 = c()
resl_1942 = c()

resu_2020 = c()
resl_2020 = c()

for(temp_i_want in seq(20, 36, by = 0.25)){
  
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
  print(paste0("Cost threshold = ", signif(thresh,3)))
  
  prob_observing_conditional_1942 = sum((sim_costs > thresh) & above_1942)/sum((sim_costs > thresh))
  print(paste0("p = ", signif(prob_observing_conditional_1942,3)))
  res_1942 = c(res_1942, prob_observing_conditional_1942)
  
  cis = qbinom(c(0.025, 0.975), size = sum(sim_costs > thresh), prob = prob_observing_conditional_1942)/sum(sim_costs > thresh)
  resu_1942 = c(resu_1942, cis[2])
  resl_1942 = c(resl_1942, cis[1])
  
  prob_observing_conditional_2020 = sum((sim_costs > thresh) & above_2020)/sum((sim_costs > thresh))
  print(paste0("p = ", signif(prob_observing_conditional_2020,3)))
  res_2020 = c(res_2020, prob_observing_conditional_2020)
  
  cis = qbinom(c(0.025, 0.975), size = sum(sim_costs > thresh), prob = prob_observing_conditional_2020)/sum(sim_costs > thresh)
  resu_2020 = c(resu_2020, cis[2])
  resl_2020 = c(resl_2020, cis[1])
}


res_2020 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*res_2020
res_1942 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*res_1942

resu_1942 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*resu_1942
resl_1942 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*resl_1942

resu_2020 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*resu_2020
resl_2020 = readRDS("data/processed/prob_event_greater_1_mod_2")*q*resl_2020


tibble(temps = seq(20, 36, by = 0.25), 
       prob_1942 = res_1942,
       prob_2020 = res_2020,
       prob_1942_u = resu_1942,
       prob_1942_l = resl_1942,
       prob_2020_u = resu_2020,
       prob_2020_l = resl_2020) %>%
  write_csv("output/prob_observing_T_anywhere_clim_mod_2.csv")



# ---- MODEL 4
obs_sites = read_csv("data/processed/obs_data.csv") %>%
  dplyr::select(Station, scale_9, threshold_9, Long.projected, Lat.projected) %>%
  unique() %>%
  left_join(read_csv("data/processed/obs_data_dist_to_sea.csv") %>%
              dplyr::select(Station, dist_sea))

obs_grid = read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_25.csv"))  %>%
  left_join(obs_sites) %>%
  left_join(read_csv("data/processed/obs_data.csv") %>%
              dplyr::select(year, loess_temp_anom) %>% unique)%>%
              filter(year %in% c(1942, 2020))

gpd_pars = read_csv("output/gpd_model_fits/model_4_true.csv")
obs_grid$scale = my_predict_4b(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$scale
obs_grid$shape = my_predict_4b(unlist(gpd_pars), obs_grid$scale_9, obs_grid$loess_temp_anom, obs_grid$dist_sea)$shape

my_simulations_extremes = c()
for(i in seq(1, 100)){
  my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true/mod_4_run_",i)))
}

sim_costs = map(my_simulations_extremes, mean) %>% unlist

q = 0.5
res_2020 = c()
res_1942 = c()

resu_1942 = c()
resl_1942 = c()

resu_2020 = c()
resl_2020 = c()

for(temp_i_want in seq(20, 36, by = 0.25)){

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
  print(paste0("Cost threshold = ", signif(thresh,3)))

  prob_observing_conditional_1942 = sum((sim_costs > thresh) & above_1942)/sum((sim_costs > thresh))
  print(paste0("p = ", signif(prob_observing_conditional_1942,3)))
  res_1942 = c(res_1942, prob_observing_conditional_1942)

  cis = qbinom(c(0.025, 0.975), size = sum(sim_costs > thresh), prob = prob_observing_conditional_1942)/sum(sim_costs > thresh)
  resu_1942 = c(resu_1942, cis[2])
  resl_1942 = c(resl_1942, cis[1])

  prob_observing_conditional_2020 = sum((sim_costs > thresh) & above_2020)/sum((sim_costs > thresh))
  print(paste0("p = ", signif(prob_observing_conditional_2020,3)))
  res_2020 = c(res_2020, prob_observing_conditional_2020)

  cis = qbinom(c(0.025, 0.975), size = sum(sim_costs > thresh), prob = prob_observing_conditional_2020)/sum(sim_costs > thresh)
  resu_2020 = c(resu_2020, cis[2])
  resl_2020 = c(resl_2020, cis[1])
}



res_2020 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*res_2020
res_1942 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*res_1942

resu_1942 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*resu_1942
resl_1942 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*resl_1942

resu_2020 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*resu_2020
resl_2020 = readRDS("data/processed/prob_event_greater_1_mod_4")*q*resl_2020


#
# library(scales)
# tibble(temps = seq(20, 36, by = 0.25),
#        rl_1942 = 1/(92*res_1942),
#        rl_2020 = 1/(92*res_2020),
#        rl_1942_u = 1/(92*resu_1942),
#        rl_1942_l = 1/(92*resl_1942),
#        rl_2020_u = 1/(92*resu_2020),
#        rl_2020_l = 1/(92*resl_2020))%>%
#   filter(temps > 22,
#          temps < 34.5) %>%
#   ggplot()+
#   geom_ribbon(aes(x = temps, ymin = (rl_1942_l), ymax = (rl_1942_u), fill = "1942"), alpha = 0.5)+
#   geom_ribbon(aes(x = temps, ymin = (rl_2020_l), ymax = (rl_2020_u), fill = "2020"), alpha = 0.5)+
#   geom_line(aes(temps, (rl_1942), col = "1942"), linetype = 'dashed')+
#   geom_line(aes(temps, (rl_2020), col = "2020"))+
#   scale_y_log10()+
#   labs(x = "Temperature (°C)",
#        y = "Return period\n(Years)",
#        fill = 'Year', col = "Year")+
#   theme_minimal(12)+
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
#         legend.position = 'none')

# 1 year event probability 1/92 = 0.01086957
# 100 year event probability 1/(100*92) = 0.0001086957


tibble(temps = seq(20, 36, by = 0.25),
       prob_1942 = res_1942,
       prob_2020 = res_2020,
       prob_1942_u = resu_1942,
       prob_1942_l = resl_1942,
       prob_2020_u = resu_2020,
       prob_2020_l = resl_2020) %>%
  write_csv("output/prob_observing_T_anywhere_clim_mod_4.csv")
