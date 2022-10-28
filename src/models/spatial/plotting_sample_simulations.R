# simulate under true model
# ---- similate on clim scale grid
rm(list=ls())

setwd("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/")
library(tidyverse)
library(doParallel)
library(foreach)




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




marg_mod = 'mod_4'
obs_smoothed_quantiles = readRDS("output/quant_models_clim_num_quantiles_50.csv")
source('src/models/marginal_models/gpd_models.R')

sset.seed(1234)
# --- marginal model names
marg_mod = 'mod_4'
read_csv("data/processed/clim_data_full.csv") %>% dplyr::select(Long.projected, Lat.projected) %>% unique %>% write_csv("data/processed/clim_grid_simulated_on.csv")
locs_to_pred = as.data.frame(read_csv("data/processed/clim_grid_simulated_on.csv"))
read_csv("data/processed/clim_data_full.csv") %>% dplyr::select(Long, Lat) %>% unique %>% write_csv("data/processed/clim_grid_simulated_on_not_projected.csv")

variogram_model = function(h){
  h = sqrt(norm(h,type = "2")^2)
  nu=0.2
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha_var*(1 - ((((h[h != 0] /beta_var)^nu) * besselK(x = (h[h != 0] / beta_var), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}

fit = read_csv(paste0("output/rpareto_model_fits/true_rpareto_fits_model_",marg_mod, ".csv"),
               col_names = F) %>% 
  .[nrow(.),] %>% # get last row in this file
  unlist %>% as.numeric()

alpha_var <<- fit[1]
beta_var <<- fit[2]
nCores <- 25
cl <- parallel::makeCluster(nCores)
clusterSetRNGStream(cl)

# this returns samples in the unit frechet margin
simulations <- mvPot::simulPareto(n = 500,
                                  loc = locs_to_pred,
                                  vario = variogram_model,
                                  nCores = nCores,
                                  cl = cl)


most_ex = map(simulations, mean) %>% unlist %>% sort() %>% tail(5)
ind_of_ex = which((map(simulations, mean) %>% unlist) %in% most_ex)



sims = read_csv('data/processed/clim_grid_simulated_on_not_projected.csv') %>%
  mutate(sim_1 = simulations[[ind_of_ex[1]]],
         sim_2 = simulations[[ind_of_ex[2]]],
         sim_3 = simulations[[ind_of_ex[3]]],
         sim_4 = simulations[[ind_of_ex[4]]],
         sim_5 = simulations[[ind_of_ex[5]]]) %>%
  pivot_longer(-c(Long, Lat)) %>%
  rename(frechet = value)

    
# sims = read_csv('data/processed/clim_grid_simulated_on_not_projected.csv') %>%
#   mutate(sim_1 = simulations[[1]],
#          sim_2 = simulations[[2]],
#          sim_3 = simulations[[3]],
#          sim_4 = simulations[[4]],
#          sim_5 = simulations[[5]]) %>%
#   pivot_longer(-c(Long, Lat)) %>%
#   rename(frechet = value)

sims = sims %>%
  left_join(read_csv('~/JRSS_organised_code/data/processed_data/clim_data_dist_to_sea.csv') %>%
              dplyr::select(Long, Lat, dist_sea)) %>%
  left_join(read_csv("data/processed/clim_scale_grid.csv") %>% 
              dplyr::select(Long, Lat, scale_9))

sims = rbind(sims %>% mutate(year = 2020), sims %>% mutate(year = 1942)) %>%
  left_join(read_csv("data/processed/obs_data.csv") %>%
              dplyr::select(year, loess_temp_anom) %>%
              filter(year %in% c(2020, 1942)) %>% 
              unique()) %>%
  left_join(obs_smoothed_quantiles %>%
              dplyr::select(Long, Lat, tau_to_temp, threshold_9, thresh_exceedance_9, year))



marg_mod = 'mod_4'
this_fit_mod = read_csv("output/gpd_model_fits/model_4_true.csv") %>%
  unlist() %>% as.numeric
sims$shape = this_fit_mod[length(this_fit_mod)]
sims$scale = my_predict_4b(this_fit_mod, sims$scale_9, sims$loess_temp_anom, sims$dist_sea)$scale



marg_mod = 'mod_1'
this_fit_mod = read_csv("output/gpd_model_fits/model_1_true.csv") %>%
  unlist() %>% as.numeric
sims$shape = this_fit_mod[length(this_fit_mod)]
sims$scale = my_predict_1(this_fit_mod, sims$scale_9)$scale




marg_mod = 'mod_2'
this_fit_mod = read_csv("output/gpd_model_fits/model_2_true.csv") %>%
  unlist() %>% as.numeric
sims$shape = this_fit_mod[length(this_fit_mod)]
sims$scale = my_predict_2(this_fit_mod, sims$scale_9, sims$loess_temp_anom)$scale

# min of 5
sims$unif = evd::pgev(sims$frechet, 1,1,1)

extreme_ind = sims$unif > (1-sims$thresh_exceedance_9)
sims$date_scale = NA
sims$date_scale[extreme_ind] = evd::qgpd(p = 1 + (sims$unif[extreme_ind]-1)/sims$thresh_exceedance_9[extreme_ind],
                           shape = sims$shape[1],
                           loc = sims$threshold_9[extreme_ind],
                           scale = sims$scale[extreme_ind])


for(i in seq(nrow(sims))){
  if(!extreme_ind[i]){
    sims$date_scale[i] = sims[i,]$tau_to_temp[[1]](sims[i,]$unif)
  }
}


sims = sims %>%
  rename(lab = name)

plt = gridExtra::grid.arrange(sims %>%
                          dplyr::select(Long, Lat, year, date_scale, lab) %>%
                          filter(year == 2020) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = date_scale))+
                          facet_grid(year~lab)+
                          geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
                          ggplot2::scale_color_gradientn(colors = my_pal)+
                          theme_minimal()+
                          labs(col = '°C')+
                          theme(axis.text = element_blank(),
                                strip.text = element_blank(),
                                axis.title = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                plot.margin=unit(c(0,0,-0.1,0),"cm")),
                        sims %>%
                          dplyr::select(Long, Lat, year, date_scale, lab) %>%
                          filter(year == 2020) %>% 
                          rename(temp_new = date_scale) %>%
                          dplyr::select(-year) %>%
                          left_join(sims %>%
                                      dplyr::select(Long, Lat, year, date_scale, lab) %>%
                                      filter(year == 1942) %>% 
                                      rename(temp_old = date_scale)) %>%
                          mutate(diff = temp_new - temp_old) %>% 
                          ggplot()+
                          geom_point(aes(Long, Lat, col = diff))+
                          geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
                          facet_grid(year~lab)+
                          ggplot2::scale_color_gradientn(colors = my_pal)+
                          theme_minimal()+
                          labs(col = '°C')+
                          theme(axis.text = element_blank(),
                                strip.text = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.title = element_blank(),
                                plot.margin=unit(c(-0.1,0,0,0),"cm")), nrow = 2)


ggsave(plot = plt, filename = paste0('output/figs/',marg_mod,'_sims.pdf'), width = 9, height = 3.9)








sims %>%
  dplyr::select(Long, Lat, year, date_scale, lab,frechet, unif) %>%
  filter(year == 2020) %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = -log(1 - unif)))+
  facet_wrap(~lab, nrow = 1)+
  geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  theme_minimal()+
  labs(col = '')+
  theme(axis.text = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,-0.1,0),"cm"))
